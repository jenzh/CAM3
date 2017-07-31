#include <misc.h>
#include <preproc.h>

module mksrfdatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksrfdatMod
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI and urban fraction.
!
! !USES:
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PRIVATE MEMBER FUNCTIONS:
  public :: mksrfdat        ! create model surface dataset

! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksrfdat
!
! !INTERFACE:
  subroutine mksrfdat(cam_longxy, cam_latixy, cam_numlon, cam_landfrac, &
                      cam_landmask)
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI and urban fraction.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use nanMod
    use clm_varpar
    use clm_varctl
    use clm_varsur
    use pftvarcon   , only : noveg
    use time_manager, only : is_last_step
    use fileutils   , only : putfil, opnfil, getavu, relavu, get_filename
    use areaMod
    use spmdMod
    use mkgridMod
    use ncdio
    use mklai
    use mkpft
!
! !ARGUMENTS:
    implicit none
    real(r8), optional, intent(in) :: cam_longxy(:,:)   !cam lon values
    real(r8), optional, intent(in) :: cam_latixy(:,:)   !cam lat values
    integer , optional, intent(in) :: cam_numlon(:)     !cam number of longitudes
    real(r8), optional, intent(in) :: cam_landfrac(:,:) !cam fractional land
    integer , optional, intent(in) :: cam_landmask(:,:) !cam land mask
!
! !CALLED FROM:
! routine initialize
!
! !REVISION HISTORY:
! Authors: Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k,m             ! indices
    integer  :: ier                 ! error status
    integer  :: ndiag               ! unit number for surface data summary
    integer  :: ncid                ! netcdf id
    integer  :: omode               ! netcdf output mode
    integer  :: ret                 ! netcdf return status
    integer  :: values(8)           ! temporary
    integer  :: dimid               ! temporary
    integer  :: pftsize             ! size of lsmpft dimension
    real(r8) :: sum                 ! sum for error check
    real(r8) :: pctveg              ! percentage of gridd cell that is vegetated
    real(r8) :: rmax                ! maximum patch cover
    logical  :: lremov = .false.    ! true => remove file after dispose
    character(len=256) :: loc_fn    ! local file name
    character(len=256) :: rem_dir   ! mass store file name
    character(len=256) :: rem_fn    ! mass store full path name
    character(len=  7) :: resol     ! resolution for file name
    character(len=256) :: str       ! global attribute string
    character(len=256) :: name      ! name of attribute
    character(len=256) :: unit      ! units of attribute
    character(len= 18) :: datetime  ! temporary
    character(len=  8) :: date      ! temporary
    character(len= 10) :: time      ! temporary
    character(len=  5) :: zone      ! temporary
    integer , allocatable :: pft_max(:,:,:)    ! PFT data: PFT values
    real(r8), allocatable :: pctpft_max(:,:,:) ! PFT data: % of vegetated area for PFTs
    real(r8), allocatable :: pctpft(:,:,:)     ! PFT data: % of gridcell for PFTs
    real(r8), allocatable :: pctlnd_pft(:,:)   ! PFT data: % land per gridcell
    real(r8), allocatable :: landfrac_pft(:,:) ! PFT data: land fraction per gridcell
    character(len=32) :: subname = 'mksrfdat'  ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to create surface boundary data .....'
       write (6,'(72a1)') ("-",i=1,60)

       ! ----------------------------------------------------------------------
       ! Allocate dynamic memory
       ! ----------------------------------------------------------------------

       allocate(pctlnd_pft(lsmlon,lsmlat), landfrac_pft(lsmlon,lsmlat), &
            pctpft(lsmlon,lsmlat,0:numpft), stat=ier)
       if (ier/=0) then
          write(6,*)'MKSRFDATMOD: allocation error'; call endrun()
       end if
       pctpft(:,:,:) = 1.e36
       if (.not. mksrf_all_pfts) then
          allocate(pft_max(lsmlon,lsmlat,maxpatch_pft), pctpft_max(lsmlon,lsmlat,maxpatch_pft), stat=ier)
          if (ier/=0) then
             write(6,*)'MKSRFDATMOD: allocation error'; call endrun()
          endif
          pft_max(:,:,:) = 0
          pctpft_max(:,:,:) = 1.e36
       end if

       ! ----------------------------------------------------------------------
       ! Open diagnostic output log file
       ! ----------------------------------------------------------------------

       loc_fn = './surface-data.log'
       ndiag = getavu()
       call opnfil (loc_fn, ndiag, 'f')

       ! ----------------------------------------------------------------------
       ! Create netcdf surface dataset and define dimensions and global attributes
       ! ----------------------------------------------------------------------

#if (defined OFFLINE)
       if (mksrf_offline_fgrid /= ' ') then
          write (ndiag,*)'using fractional land data from file= ', &
               trim(mksrf_offline_fgrid),' to create the surface dataset'
       else
          write (ndiag,*)'using fractional land data from file= ', &
               trim(mksrf_offline_fnavyoro),' to create the surface dataset'
       endif
#elif (defined COUP_CAM)
       write (ndiag,*)'using fractional land data from cam', &
            ' model to create the surface dataset'
#elif (defined COUP_CSM)
       write (ndiag,*)'using fractional land data from csm', &
            ' flux coupler to create the surface dataset'
#endif
       write (ndiag,*) 'PFTs from:         ',trim(mksrf_fvegtyp)
       write (ndiag,*) 'glaciers from:     ',trim(mksrf_fglacier)
       write (ndiag,*) 'urban from:        ',trim(mksrf_furban)
       write (ndiag,*) 'inland water from: ',trim(mksrf_flanwat)
       write (ndiag,*) 'soil texture from: ',trim(mksrf_fsoitex)
       write (ndiag,*) 'soil color from:   ',trim(mksrf_fsoicol)

       ! Create netCDF surface dataset.  File will be in define mode
       ! Set fill mode to "no fill" to optimize performance

       write (resol,'(i3.3,"x",i3.3)') lsmlon,lsmlat
       fsurdat = './surface-data.'//trim(resol)//'.nc'

       call check_ret(nf_create(trim(fsurdat), nf_clobber, ncid), subname)
       call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

       ! Define dimensions.

       if (mksrf_all_pfts) then
          pftsize = numpft + 1
       else
          pftsize = maxpatch_pft
       end if

       call check_ret(nf_def_dim (ncid, 'lsmlon' , lsmlon      , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lsmlat' , lsmlat      , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'nlevsoi', nlevsoi     , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lsmpft' , pftsize     , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'time'   , nf_unlimited, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'nchar'  , 128         , dimid), subname)

       ! Create global attributes.

       str = 'NCAR-CSM'
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Conventions', len_trim(str), trim(str)), subname)

       call date_and_time (date, time, zone, values)
       datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
       datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
       str = 'created on: ' // datetime
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'History_Log', len_trim(str), trim(str)), subname)

       call shr_sys_getenv ('LOGNAME', str, ier)
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Logname', len_trim(str), trim(str)), subname)

       call shr_sys_getenv ('HOST', str, ier)
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Host', len_trim(str), trim(str)), subname)

       str = 'Community Land Model: CLM2'
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Source', len_trim(str), trim(str)), subname)

       str = '$Name: cam3_0_brnchT_release01 $'
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Version', len_trim(str), trim(str)), subname)

       str = '$Id: mksrfdatMod.F90,v 1.12.2.10.4.2 2004/06/21 17:57:53 mvr Exp $'
       call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
            'Revision_Id', len_trim(str), trim(str)), subname)

       if (mksrf_offline_fgrid /= ' ') then
          str = mksrf_offline_fgrid
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
               'Input_grid_dataset', len_trim(str), trim(str)), subname)
       else
          str = mksrf_offline_fnavyoro
          call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
               'Input_navy_oro_dataset', len_trim(str), trim(str)), subname)
       endif

       str = get_filename(mksrf_fvegtyp)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Vegetation_type_raw_data_filename', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fsoitex)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_texture_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fsoicol)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_color_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_flanwat)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Inland_water_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fglacier)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Glacier_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_furban)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Urban_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_flai)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)

#if (defined OFFLINE)
       str = 'offline'
#elif (defined COUP_CAM)
       str = 'run through cam'
#elif (defined COUP_CSM)
       str = 'run through flux coupler'
#endif
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Run_mode', len_trim(str), trim(str)), subname)

       ! ----------------------------------------------------------------------
       ! Write fields other than lai, sai, and heights to netcdf surface dataset
       ! ----------------------------------------------------------------------

#ifdef OFFLINE
       if (mksrf_offline_fgrid == ' ') then
          call ncd_defvar(ncid=ncid, varname='EDGEN', xtype=nf_float, &
               long_name='northern edge of surface grid', units='degrees north')

          call ncd_defvar(ncid=ncid, varname='EDGEE', xtype=nf_float, &
               long_name='eastern edge of surface grid', units='degrees east')

          call ncd_defvar(ncid=ncid, varname='EDGES', xtype=nf_float, &
               long_name='southern edge of surface grid', units='degrees north')

          call ncd_defvar(ncid=ncid, varname='EDGEW', xtype=nf_float, &
               long_name='western edge of surface grid', units='degrees east')
       endif
#endif

       call ncd_defvar(ncid=ncid, varname='NUMLON', xtype=nf_int, &
            dim1name='lsmlat', long_name='number of longitudes for each latitude', units='unitless')

       name = 'longitude'
       if (.not. fullgrid) name = 'rlongitude'
       call ncd_defvar(ncid=ncid, varname='LONGXY', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name=trim(name), units='degrees east')

       call ncd_defvar(ncid=ncid, varname='LATIXY', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='latitude', units='degrees north')

       call ncd_defvar(ncid=ncid, varname='LANDMASK', xtype=nf_int, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='land/ocean mask', units='0=ocean and 1=land')

       call ncd_defvar(ncid=ncid, varname='LANDFRAC', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='land fraction', units='unitless')

       call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='land fraction from pft dataset', units='unitless')

       call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='soil color', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent sand', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
            long_name='percent clay', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_WETLAND', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent wetland', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_LAKE', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent lake', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_GLACIER', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent glacier', units='unitless')

       call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=nf_float, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent urban', units='unitless')

       if (mksrf_all_pfts) then
          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=nf_float, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
               long_name='percent plant functional type of gridcell', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PFT', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
               long_name='plant functional type', units='unitless')

          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=nf_float, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
               long_name='percent plant functional type of vegetated gridcell area', units='unitless')
       end if

       call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly leaf area index', units='unitless')

       call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly stem area index', units='unitless')

       call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly height top', units='meters')

       call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=nf_float,  &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
            long_name='monthly height bottom', units='meters')

       ! End of define mode

       call check_ret(nf_enddef(ncid), subname)

       ! ----------------------------------------------------------------------
       ! Initialize surface variables with unusable values
       ! ----------------------------------------------------------------------

       soic2d(:,:)   = -999
       sand3d(:,:,:) = 1.e36
       clay3d(:,:,:) = 1.e36
       pctlak(:,:)   = 1.e36
       pctwet(:,:)   = 1.e36
       pcturb(:,:)   = 1.e36
       pctgla(:,:)   = 1.e36

       ! ----------------------------------------------------------------------
       ! Determine land model grid, fractional land and land mask
       ! ----------------------------------------------------------------------

       ! Initialize grid variables with unusable values

       numlon(:)     = 0
       latixy(:,:)   = 1.e36
       longxy(:,:)   = 1.e36
       landmask(:,:) = -999
       landfrac(:,:) = 1.e36

#if (defined OFFLINE)
       call mkgrid_offline()
#else
       call mkgrid_cam(cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
#endif

       ! ----------------------------------------------------------------------
       ! Make PFTs [pctpft] from dataset [fvegtyp] (1/2 degree PFT data)
       ! ----------------------------------------------------------------------

       if (mksrf_all_pfts) then
          call mkpfts(mksrf_fvegtyp, ndiag, noveg, pctlnd_pft, pctpft)
       else
          call mkpfts(mksrf_fvegtyp, ndiag, noveg, pctlnd_pft, pctpft, pctpft_max, pft_max)
       end if

       ! ----------------------------------------------------------------------
       ! Make inland water [pctlak, pctwet] from Cogley's one degree data [flanwat]
       ! ----------------------------------------------------------------------

       call mklanwat (mksrf_flanwat, ndiag, pctlak, pctwet)

       ! ----------------------------------------------------------------------
       ! Make glacier fraction [pctgla] from [fglacier] dataset
       ! ----------------------------------------------------------------------

       call mkglacier (mksrf_fglacier, ndiag, pctgla)

       ! ----------------------------------------------------------------------
       ! Make soil texture [sand3d, clay3d] from IGBP 5 minute data [fsoitex]
       ! ----------------------------------------------------------------------

       call mksoitex (mksrf_fsoitex, ndiag, pctgla, sand3d, clay3d)

       ! ----------------------------------------------------------------------
       ! Make soil color classes [soic2d] from BATS T42 data [fsoicol]
       ! ----------------------------------------------------------------------

       call mksoicol (mksrf_fsoicol, ndiag, pctgla, soic2d)

       ! ----------------------------------------------------------------------
       ! Make urban fraction [pcturb] from [furban] dataset
       ! ----------------------------------------------------------------------

       call mkurban (mksrf_furban, ndiag, pcturb)

       ! ----------------------------------------------------------------------
       ! Set LAND values on Ross ice shelf to glacier
       ! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (latixy(i,j) < -79. .and. landmask(i,j) == 1) then
                soic2d(i,j) = 0
                pctlak(i,j) = 0.
                pctwet(i,j) = 0.
                pcturb(i,j) = 0.
                pctgla(i,j) = 100.
                if (mksrf_all_pfts) then
                   pctpft(i,j,0) = 100.
                   pctpft(i,j,1:numpft)  = 0.
                else
                   pft_max(i,j,1:maxpatch_pft) = noveg
                   pctpft_max(i,j,1) = 100.
                   pctpft_max(i,j,2:maxpatch_pft) = 0.
                end if
                sand3d(i,j,1:nlevsoi) = 0.
                clay3d(i,j,1:nlevsoi) = 0.
             end if
          end do
       end do

       ! ----------------------------------------------------------------------
       ! Assume wetland and/or lake when input landmask says land and pft
       ! dataset landmask says ocean (assume medium soil color (4) and loamy
       ! texture).
       ! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             if (landmask(i,j)==1 .and. pctlnd_pft(i,j)==0.) then
                soic2d(i,j) = 4
                pctwet(i,j) = 100. - pctlak(i,j)
                pcturb(i,j) = 0.
                pctgla(i,j) = 0.
                if (mksrf_all_pfts) then
                   pctpft(i,j,0) = 100.
                   pctpft(i,j,1:numpft) = 0.
                else
                   pft_max(i,j,1:maxpatch_pft) = noveg
                   pctpft_max(i,j,1) = 100.
                   pctpft_max(i,j,2:maxpatch_pft) = 0.
                end if
                sand3d(i,j,1:nlevsoi) = 43.
                clay3d(i,j,1:nlevsoi) = 18.
             end if
          end do
       end do

       ! ----------------------------------------------------------------------
       ! If have pole points on grid - set south pole to glacier
       ! north pole is as assumed as non-land
       ! ----------------------------------------------------------------------

#if (!defined OFFLINE)
       if (pole_points) then
          do i = 1,numlon(1)
             soic2d(i,1)   = 0
             pctlak(i,1)   = 0.
             pctwet(i,1)   = 0.
             pcturb(i,1)   = 0;
             sand3d(i,1,:) = 0.
             clay3d(i,1,:) = 0.
             pctgla(i,1)   = 100.
             if (mksrf_all_pfts) then
                pctpft(i,1,0) = 100.
                pctpft(i,1,1:numpft) = 0.
             else
                pft_max(i,1,1:maxpatch_pft) = noveg
                pctpft_max(i,1,1) = 100.
                pctpft_max(i,1,2:maxpatch_pft) = 0.
             end if
          end do
       end if
#endif

       ! ----------------------------------------------------------------------
       ! Truncate all percentage fields on output grid. This is needed to
       ! insure that wt is not nonzero (i.e. a very small number such as
       ! 1e-16) where it really should be zero
       ! ----------------------------------------------------------------------

       do j = 1,lsmlat
          do i = 1,numlon(j)
             do k = 1,nlevsoi
                sand3d(i,j,k) = float(nint(sand3d(i,j,k)))
                clay3d(i,j,k) = float(nint(clay3d(i,j,k)))
             end do
             pctlak(i,j) = float(nint(pctlak(i,j)))
             pctwet(i,j) = float(nint(pctwet(i,j)))
             pcturb(i,j) = float(nint(pcturb(i,j)))
             pctgla(i,j) = float(nint(pctgla(i,j)))
             if (.not. mksrf_all_pfts) then
                do m = 1,maxpatch_pft
                   pctpft_max(i,j,m) = float(nint(pctpft_max(i,j,m)))
                end do
             end if
          end do
       end do

       ! ----------------------------------------------------------------------
       ! Make sure sum of land cover types does not exceed 100. If it does,
       ! subtract excess from most dominant land cover.
       ! ----------------------------------------------------------------------

       do j = 1, lsmlat
          do i = 1, numlon(j)

             rmax = -9999.
             k    = -9999
             if (pctlak(i,j) > rmax) then
                k = 1
                rmax = pctlak(i,j)
             end if
             if (pctwet(i,j) > rmax) then
                k = 2
                rmax = pctwet(i,j)
             end if
             if (pcturb(i,j) > rmax) then
                k = 3
                rmax = pcturb(i,j)
             end if
             if (pctgla(i,j) > rmax) then
                k = 4
                rmax = pctgla(i,j)
             end if
             sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
             if (k == -9999) then
                write (6,*) 'MKSRFDAT error: largest patch not found'
                call endrun
             else if (sum > 120.) then
                write (6,*) 'MKSRFDAT error: sum of pctlak, pctwet,', &
                     'pcturb and pctgla is greater than 120%'
                write (6,*)'i,j,pctlak,pctwet,pcturb,pctgla= ', &
                     i,j,pctlak(i,j),pctwet(i,j),pcturb(i,j),pctgla(i,j)
                call endrun
             else if (sum > 100.) then
                if (k==1) pctlak(i,j) = pctlak(i,j) - (sum-100.)
                if (k==2) pctwet(i,j) = pctwet(i,j) - (sum-100.)
                if (k==3) pcturb(i,j) = pcturb(i,j) - (sum-100.)
                if (k==4) pctgla(i,j) = pctgla(i,j) - (sum-100.)
             end if

             if (mksrf_all_pfts) then

                ! Normalize pctpft to be the remainder of [100 - (special landunits)]

                sum = pctlak(i,j)+pctwet(i,j)+pcturb(i,j)+pctgla(i,j)
                do m = 0, numpft
                   pctpft(i,j,m) = 0.01_r8 * pctpft(i,j,m) * (100._r8 - sum)
                end do

             else

                ! If the input landfrac is greater than the pft dataset landfrac, then
                ! convert the discrepancy to wetland. Only convert gridcells where there
                ! is no lake.

!!$             if (landmask(i,j)==1 .and. pctlnd_pft(i,j)>0.) then
!!$                if (((landfrac(i,j) - pctlnd_pft(i,j)/100.) > 0.001) .and. &
!!$                     pctlak(i,j)==0. .and. pctgla(i,j)<100. .and. &
!!$                     pcturb(i,j)<100. .and. pctwet(i,j)<100. ) then
!!$                   pctveg = 100. - (pcturb(i,j)+pctwet(i,j)+pctgla(i,j))
!!$                   pcturb(i,j) = pcturb(i,j) * pctlnd_pft(i,j) / (landfrac(i,j)*100.)
!!$                   pctgla(i,j) = pctgla(i,j) * pctlnd_pft(i,j) / (landfrac(i,j)*100.)
!!$                   pctwet(i,j) = pctwet(i,j) * pctlnd_pft(i,j) / (landfrac(i,j)*100.)
!!$                   pctveg      = pctveg      * pctlnd_pft(i,j) / (landfrac(i,j)*100.)
!!$                   sum = pcturb(i,j)+pctgla(i,j)+pctwet(i,j)+pctveg
!!$                   if (sum < 100.) pctwet(i,j) = pctwet(i,j) + (100.-sum)
!!$                end if
!!$             end if

                ! Make sure sum of PFT cover equals 100 for land points. If it does not,
                ! subtract excess from most dominant PFT.

                rmax = -9999.
                k = -9999
                sum = 0.
                do m = 1, maxpatch_pft
                   sum = sum + pctpft_max(i,j,m)
                   if (pctpft_max(i,j,m) > rmax) then
                      k = m
                      rmax = pctpft_max(i,j,m)
                   end if
                end do
                if (k == -9999) then
                   write (6,*) 'MKSRFDAT error: largest PFT patch not found'
                   call endrun
                else if (landmask(i,j) == 1) then
                   if (sum < 95 .or. sum > 105.) then
                      write (6,*) 'MKSRFDAT error: sum of PFT cover is ',sum
                      call endrun
                   else if (sum /= 100.) then
                      pctpft_max(i,j,k) = pctpft_max(i,j,k) - (sum-100.)
                   endif
                endif

             end if

          end do
       end do

       ! ----------------------------------------------------------------------
       ! Determine fractional land from pft dataset
       ! ----------------------------------------------------------------------

       landfrac_pft(:,:) = pctlnd_pft(:,:)/100.

       ! ----------------------------------------------------------------------
       ! Write fields other than lai, sai, and heights to netcdf surface dataset
       ! ----------------------------------------------------------------------

#ifdef OFFLINE
       if (mksrf_offline_fgrid == ' ') then
          call ncd_ioglobal(varname='EDGEN', data=lsmedge(1), ncid=ncid, flag='write')
          call ncd_ioglobal(varname='EDGEE', data=lsmedge(2), ncid=ncid, flag='write')
          call ncd_ioglobal(varname='EDGES', data=lsmedge(3), ncid=ncid, flag='write')
          call ncd_ioglobal(varname='EDGEW', data=lsmedge(4), ncid=ncid, flag='write')
       endif
#endif
       call ncd_ioglobal(varname='NUMLON'      , data=numlon      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LONGXY'      , data=longxy      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LATIXY'      , data=latixy      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDMASK'    , data=landmask    , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDFRAC'    , data=landfrac    , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDFRAC_PFT', data=landfrac_pft, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='SOIL_COLOR'  , data=soic2d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_SAND'    , data=sand3d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_CLAY'    , data=clay3d      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_WETLAND' , data=pctwet      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_LAKE'    , data=pctlak      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_GLACIER' , data=pctgla      , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='PCT_URBAN'   , data=pcturb      , ncid=ncid, flag='write')
       if (mksrf_all_pfts) then
          call ncd_ioglobal(varname='PCT_PFT', data=pctpft    , ncid=ncid, flag='write')
       else
          call ncd_ioglobal(varname='PFT'    , data=pft_max   , ncid=ncid, flag='write')
          call ncd_ioglobal(varname='PCT_PFT', data=pctpft_max, ncid=ncid, flag='write')
       end if

       ! ----------------------------------------------------------------------
       ! Make LAI and SAI from 1/2 degree data and write to surface dataset
       ! ----------------------------------------------------------------------

       if (mksrf_all_pfts) then
          call mklais(mksrf_flai, ndiag, ncid)
       else
          call mklais(mksrf_flai, ndiag, ncid, pft_max=pft_max)
       end if

       ! ----------------------------------------------------------------------
       ! Close and dispose diagnostic and surface datasets
       ! ----------------------------------------------------------------------

       close (ndiag); call relavu(ndiag)
       call check_ret(nf_close(ncid), subname)

       if (mss_irt > 0) then
          rem_dir = trim(archive_dir) // '/surf/'
          rem_fn  = trim(rem_dir)//'surface-data.'//trim(resol)//'.nc'
          call putfil (fsurdat, rem_fn, mss_wpass, mss_irt, lremov)
          rem_dir = trim(archive_dir) // '/surf/'
          rem_fn  = trim(rem_dir) // 'surface-data.log'
          call putfil (loc_fn, rem_fn, mss_wpass, mss_irt, lremov)
       end if

       write (6,'(72a1)') ("-",i=1,60)
       write (6,'(a46,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
            360./lsmlon,' by ',180./lsmlat,' grid'

       write (6,*)
       write (6,*) 'Surface data output file = ',trim(fsurdat)
       write (6,*) '   This file contains the land model surface data'
       write (6,*) 'Diagnostic log file      = ',trim(loc_fn)
       write (6,*) '   See this file for a summary of the dataset'
       write (6,*)

       ! ----------------------------------------------------------------------
       ! Deallocate dynamic memory
       ! ----------------------------------------------------------------------

       deallocate(pctpft, pctlnd_pft, landfrac_pft)
       if (.not. mksrf_all_pfts) deallocate(pft_max, pctpft_max)

    end if   ! end of if-masterproc block

#if (!defined COUP_CSM)

    ! ----------------------------------------------------------------------
    ! End run if only making surface dataset
    ! ----------------------------------------------------------------------

    ! Note that nestep is determined by the flux coupler and not by
    ! the namelist for a coupled model run

    if (is_last_step()) then
       write (6,*)
       write (6,*)'model stopped because run length is zero'
       write (6,*)'offline case triggers this with nelapse=1'
       call endrun()
    end if
#endif

    ! ----------------------------------------------------------------------
    ! Reset real arrays to 1.e36 and integer arrays to -999 since all
    ! these arrays will be read back in to insure that bit for bit results
    ! are obtained for a run where a surface dataset file is generated and
    ! a run where a surface dataset is read in
    ! ----------------------------------------------------------------------

    lsmedge(:) = inf
    lats(:) = inf
    lonw(:,:) = inf
    numlon(:) =  -999
    latixy(:,:) = 1.e36
    longxy(:,:) = 1.e36
    landmask(:,:) =  -999
    landfrac(:,:) = 1.e36
    soic2d(:,:) = -999
    sand3d(:,:,:) = 1.e36
    clay3d(:,:,:) = 1.e36
    pctwet(:,:) = 1.e36
    pctlak(:,:) = 1.e36
    pctgla(:,:) = 1.e36
    pcturb(:,:) = 1.e36

  end subroutine mksrfdat

end module mksrfdatMod
