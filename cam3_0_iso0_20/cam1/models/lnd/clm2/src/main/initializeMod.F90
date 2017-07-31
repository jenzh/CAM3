#include <misc.h>
#include <preproc.h>

module initializeMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initializeMod
!
! !DESCRIPTION:
! Performs land model initialization
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private header    ! echo version numbers
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize
!
! !INTERFACE:
#if (defined OFFLINE) || (defined COUP_CSM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp   )
#elif (defined COUP_CAM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp   , &
                        cam_caseid  , cam_ctitle  , cam_nsrest, cam_nstep, &
                        cam_irad    , cam_crtinic , cam_nhtfrq, cam_mfilt, &
                        cam_longxy  , cam_latixy  , cam_numlon, &
                        cam_landmask, cam_landfrac, cam_irt   )
#endif
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clmexp] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use spmdMod         , only : masterproc
    use clmtypeInitMod  , only : initClmtype
    use initGridCellsMod, only : initGridCells
    use clm_varpar      , only : lsmlon, lsmlat, maxpatch
    use clm_varsur      , only : varsur_alloc, varsur_dealloc
    use clm_varctl      , only : fsurdat, finidat, nsrest, irad, &
                                 mksrf_offline_fgrid, mksrf_offline_fnavyoro
    use controlMod      , only : control_init, control_print
    use filterMod       , only : initFilters
    use decompMod       , only : initDecomp, get_proc_clumps, get_clump_bounds
    use histFldsMod     , only : initHistFlds
    use restFileMod     , only : restart
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use mksrfdatMod     , only : mksrfdat
    use surfFileMod     , only : surfrd
    use pftvarcon       , only : pftconrd
#if (defined DGVM)
    use DGVMEcosystemDynMod, only : DGVMEcosystemDynini
#else
    use STATICEcosysDynMod , only : EcosystemDynini
#endif
#if (defined DGVM)
    use DGVMMod         , only : resetTimeConstDGVM
#endif
#if (defined RTM)
    use RtmMod          , only : Rtmgridini, Rtmlandini
#endif
#if (defined OFFLINE)
    use atmdrvMod       , only : atm_getgrid
#endif
#if (defined COUP_CSM)
    use clm_csmMod      , only : csm_recvgrid, csm_initialize, csm_sendalb
#endif
#if (defined COUP_CAM)
    use lp_coupling     , only : lp_coupling_init
#endif
#if (defined COUP_CAM)
    use time_manager    , only : get_curr_date, get_nstep
#else
    use time_manager    , only : get_curr_date, get_nstep, advance_timestep, timemgr_init
#endif
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(inout) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(inout) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(inout) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
#if (defined COUP_CAM)
    character(len=*),  intent(in) :: cam_caseid   !cam caseid
    character(len=*),  intent(in) :: cam_ctitle   !cam title
    character(len=*),  intent(in) :: cam_crtinic  !cam initial dataset generation frequency
    integer ,  intent(in) :: cam_irad             !cam radiation frequency
    integer ,  intent(in) :: cam_nsrest           !cam 0=initial run, > 0=continuation run
    integer ,  intent(in) :: cam_nstep            !cam current time index
    integer ,  intent(in) :: cam_nhtfrq           !cam history write freq for tape 1
    integer ,  intent(in) :: cam_mfilt            !cam number of files per tape for tape 1
    integer ,  intent(in) :: cam_irt              !cam mss retention time
    integer ,  intent(in) :: cam_numlon(:)        !cam number of longitudes
    real(r8),  intent(in) :: cam_longxy(:,:)      !cam lon values
    real(r8),  intent(in) :: cam_latixy(:,:)      !cam lat values
    real(r8),  intent(in) :: cam_landfrac(:,:)    !cam fractional land
    integer ,  intent(in) :: cam_landmask(:,:)    !cam land mask
#endif
!
! !CALLED FROM:
! routine program_off if cpp token OFFLINE is defined
! routine program_csm if cpp token COUP_CSM is defined
! routine atmlnd_ini in module atm_lndMod if cpp token COUP_CAM is defined
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k                 !indices
    integer  :: yr                    !current year (0, ...)
    integer  :: mon                   !current month (1 -> 12)
    integer  :: day                   !current day (1 -> 31)
    integer  :: ncsec                 !current time of day [seconds]
    logical  :: readini               !true if read in initial data set
    integer  :: vegxy(lsmlon,lsmlat,maxpatch) !vegetation type
    real(r8) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid weights
#ifdef DGVM
    integer  :: nc            ! clump index
    integer  :: nclumps       ! number of clumps on this processor
    integer  :: begp, endp    ! clump beginning and ending pft indices
    integer  :: begc, endc    ! clump beginning and ending column indices
    integer  :: begl, endl    ! clump beginning and ending landunit indices
    integer  :: begg, endg    ! clump beginning and ending gridcell indices
#endif
#if (defined COUP_CSM)
    integer  :: cam_numlon(lsmlat)           !cam number of longitudes
    real(r8) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8) :: cam_latixy(lsmlon,lsmlat)    !cam lat values
    real(r8) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer  :: cam_landmask(lsmlon,lsmlat)  !cam land mask
#endif
    integer  :: ier
!-----------------------------------------------------------------------

    ! Initialize run control variables, time manager, timestep

    call header()

    if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
#if (defined COUP_CSM) || (defined OFFLINE)
       write (6,*) 'Preset Fortran unit numbers:'
       write (6,*) '   unit  5 = standard input'
       write (6,*) '   unit  6 = standard output'
#endif
       write (6,*)
    endif

#if (defined COUP_CAM)
    call control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                       cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )
#else
    call control_init ()
#endif
    if (masterproc) call control_print()

    ! Allocate surface grid dynamic memory

    call varsur_alloc ()

#if (defined OFFLINE) || (defined COUP_CSM)
    ! Initialize time manager for initial run

    if (nsrest == 0) call timemgr_init()
#endif

#if (defined OFFLINE)
    ! Start at nstep = 1 for an initial offline run

    if (nsrest == 0) call advance_timestep()
#endif

#if (defined RTM)
    ! Initialize RTM river routing grid and mask

     call Rtmgridini()
#endif

#if (defined COUP_CSM)

     ! Get grid and land mask back from flux coupler

     call csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
#endif

     ! Read list of PFTs and their corresponding parameter values
     ! This is independent of the model resolution

     call pftconrd ()

     ! If no surface dataset name is specified then make surface dataset
     ! from original data sources. Always read surface boundary data in.
     ! This insures that bit for bit results are obtained for a run where a
     ! surface dataset file is generated and a run where a surface dataset
     ! is specified and read in. Set up vegetation type [veg] and weight [wt]
     ! arrays for [maxpatch] subgrid patches.

#if (defined OFFLINE)
     if (fsurdat == ' ') then
        call mksrfdat ()
     endif
     call surfrd (vegxy, wtxy)
#else
     if (fsurdat == ' ') then
        call mksrfdat (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
     endif
     call surfrd (vegxy, wtxy, cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
#endif

    ! Initialize clump and processor decomposition

    call initDecomp(wtxy)

#if (defined COUP_CAM)
    ! Initialize mapping between the atmosphere physics chunks and the land clumps

    call lp_coupling_init()
#endif

    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    ! Build hierarchy and topological info for derived typees

    call initGridCells(vegxy, wtxy)

    ! Initialize filters

    call initFilters()

     ! Initialize Ecosystem Dynamics

#if (defined DGVM)
    call DGVMEcosystemDynini()
#else
    call EcosystemDynini()
#endif

     ! Initialize time constant variables

     if (masterproc) write (6,*) 'Attempting to initialize time invariant variables'
     call iniTimeConst()
     if (masterproc) write (6,*) 'Successfully initialized time invariant variables'

#if (defined RTM)
     ! Initialize river routing model

     if (masterproc) write(6,*)'Attempting to initialize RTM'
     call Rtmlandini()
     if (masterproc) write(6,*)'Successfully initialized RTM'
#endif

#if (defined COUP_CSM)
     ! Initialize flux coupler communication

     call csm_initialize(irad,eccen, obliqr, lambm0, mvelpp)
#endif

#if ( !defined SCAM)
    ! Read restart files if continuation run

    if (nsrest > 0) call restart('read')

    ! Initialize master history list. Note, routine initHistFlds will
    ! initialize active history fields if not a restart run.

    call initHistFlds() 
#endif

#if (defined DGVM)
    ! Initialize DGVM and reset DGVM time constant variables

    if (nsrest > 0) then
       nclumps = get_proc_clumps()
!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
       do nc = 1,nclumps
          call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
          call resetTimeConstDGVM(begp, endp)
       end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO
    end if
#endif

#if (defined OFFLINE)
    ! Read atmospheric forcing dataset one time to obtain the longitudes
    ! and latitudes of the atmospheric dataset, as well as the edges. When
    ! coupled to atm model, these are input variables. If no
    ! atmospheric data files are provided, model uses dummy atmospheric
    ! forcing and sets atmospheric grid to land grid.

    if (masterproc) write (6,*) 'Attempting to set up atmospheric grid '
    call atm_getgrid()
    if (masterproc) write (6,*) 'Successfully set up atmospheric grid '
#endif

    ! Initialize accumulator fields to be time accumulated for various purposes.

    if (nsrest == 0) call initAccFlds()

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0 for cam and csm mode
    ! and at nstep=1 for offline mode. This routine is also always called for a
    ! restart run and must therefore be called after the restart file is read in

    call initAccClmtype()

    ! If initial run: initialize time-varying data
    ! If continuation run: end of initialization because time varying data is read
    ! from the restart file

    if (nsrest == 0) then
       if (masterproc) write (6,*) 'Attempting to initialize time variant variables '
       readini = .true.
       if (finidat == ' ') readini = .false.
       call iniTimeVar(readini, eccen, obliqr, lambm0 , mvelpp)
       if (masterproc) then
          write (6,*) 'Successfully initialized time variant variables'
          write (6,*)
       endif
    endif

#if (defined COUP_CSM)
    ! Send first land model data to flux coupler.

    call csm_sendalb()
#endif

    ! Deallocate surface grid dynamic memory

    call varsur_dealloc()

    ! End initialization

    if (masterproc) then
       write (6,*) 'Successfully initialized the land model'
       if (nsrest == 0) then
          write (6,*) 'begin initial run at: '
       else
          write (6,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write (6,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write (6,*)
       write (6,'(72a1)') ("*",i=1,60)
       write (6,*)
    endif

  end subroutine initialize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: header
!
! !INTERFACE:
  subroutine header()
!
! !DESCRIPTION:
! Echo and save model version number
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl  , only : version
    use spmdMod     , only : masterproc
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

    version = 'CLM MODEL version 3.0'
    if ( masterproc )then
      write (6,*) trim(version)
      write (6,*)
    end if

  end subroutine header

end module initializeMod
