#include <misc.h>
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute sea ice to atmosphere surface fluxes; then compute
! sea ice temperature change.
!
! Method: 
! Temperatures over sea-ice surfaces are specified in 'plevmx' layers of
! fixed thickness and thermal properties.  The forecast temperatures are
! determined from a backward/implicit diffusion calculation using
! linearized sensible/latent heat fluxes. The bottom ocean temperature
! is fixed at -2C, allowing heat flux exchange with underlying ocean.
! Temperature over sea ice is not allowed to exceed melting temperature.
! 
! The spectral components of shortwave and albedos must be sent to 
! this sea ice model because the sea ice extinction coefficient depends
! on the wavelength.
!
! Author: C.M. Bitz
! 
!-----------------------------------------------------------------------
  subroutine seaice (c, ncol, dtime, aice, Tsice,          &
                     hi, snowh, ubot, vbot, tbot,          &
                     qbot, thbot, zbot, pbot, flwds,       &
                     swvdr, swidr, swvdf, swidf, alvdr,    &
                     alidr, alvdf, alidf, snowfall,        &
                     sst, frzmlt, Focn,                    &
                     tssub, qflx, taux, tauy, ts, shflx,   &
                     lhflx, lwup, tref, fsns, evap,        &
                     aiceinit, wthi, wtsnowh, wtsnowfall)

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols
  use comsrf,        only: plevmx
  use constituents,  only: pcnst, pnats
  use ice_constants, only: puny, Tfrez, TfrezK, qqqice, TTTice, emissivity_ice, Tffresh, &
                           tmelz, rhofresh, rhos, rhoi, c0, c1, snowpatch, saltz, rLfs, &
                           Lvap, ni, hi_min
  use ice_sfc_flux,  only: ice_atm_flux
  use ice_tstm,      only: tstm
  use ice_dh,        only: dh, prognostic_icesnow, fixice, energ
  use history,       only: outfld
#if ( defined COUP_SOM )
!JR Changed lateral_growth_melt to lateral_melt
  use ice_ocn_flux,  only: init_frzmlt, lateral_melt, init_Tprofile
#endif

  use water_tracers, only: trace_water, ixh2oq, wtrc_qchk1


  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: c               ! chunk index
  integer , intent(in) :: ncol           ! number of columns this chunk

  real(r8), intent(in) :: dtime              ! Land/ocean/seaice flag
  real(r8), intent(inout) :: aice(pcols)    ! Land/ocean/seaice flag
  real(r8), intent(inout) :: tsice(pcols)   ! ice/snow surface temperature (K)
  real(r8), intent(inout) :: snowh(pcols)   ! Snow depth (liquid water equivalent)
  real(r8), intent(inout) :: hi(pcols)      ! Ice thickness
  real(r8), intent(in) :: ubot(pcols)          ! Bottom level u wind
  real(r8), intent(in) :: vbot(pcols)          ! Bottom level v wind
  real(r8), intent(in) :: tbot(pcols)          ! Bottom level temperature
  real(r8), intent(in) :: qbot(pcols,pcnst+pnats)  ! Bottom level specific humidity
  real(r8), intent(in) :: thbot(pcols)         ! Bottom level potential temperature
  real(r8), intent(in) :: zbot(pcols)          ! Bottom level height above surface
  real(r8), intent(in) :: pbot(pcols)        ! Bottom level pressure
  real(r8), intent(in) :: flwds(pcols)   ! net down longwave radiation at surface
  real(r8), intent(in) :: swvdr(pcols)   ! direct beam solar radiation onto srf (sw)
  real(r8), intent(in) :: swidr(pcols)   ! direct beam solar radiation onto srf (lw)
  real(r8), intent(in) :: swvdf(pcols)   ! diffuse solar radiation onto srf (sw)
  real(r8), intent(in) :: swidf(pcols)   ! diffuse solar radiation onto srf (lw)
  real(r8), intent(in) :: alvdr(pcols)   ! ocean + ice albedo: shortwave, direct
  real(r8), intent(in) :: alvdf(pcols)   ! ocean + ice albedo: shortwave, diffuse
  real(r8), intent(in) :: alidr(pcols)   ! ocean + ice albedo: longwave, direct
  real(r8), intent(in) :: alidf(pcols)   ! ocean + ice albedo: longwave, diffuse
  real(r8), intent(in) :: snowfall(pcols)  ! total snow rate (m h2o/s) 
  real(r8), intent(in) :: sst(pcols)       ! sea surface temperature (C)
  real(r8), intent(in) :: frzmlt(pcols)    ! freeze/melt potential  (W/m**2) (if >0 then freeze)

  real(r8), intent(inout) :: Focn(pcols)   ! ocean-ice heat flux for basal and lateral melt (<0)
  real(r8), intent(inout):: tssub(pcols,plevmx)  ! Surface/sub-surface temperatures

! fluxes/quantities summed over surface types
  real(r8), intent(inout):: qflx(pcols,pcnst+pnats)    ! Constituent flux (kg/m2/s)
  real(r8), intent(inout):: taux(pcols)          ! X surface stress (N/m2)
  real(r8), intent(inout):: tauy(pcols)          ! Y surface stress (N/m2)
  real(r8), intent(inout):: ts(pcols)            ! surface temperature (K)
  real(r8), intent(inout):: shflx(pcols)         ! Surface sensible heat flux (J/m2/s)
  real(r8), intent(inout):: lhflx(pcols)         ! Surface latent   heat flux (J/m2/s)
  real(r8), intent(inout):: lwup(pcols)          ! surface longwave up flux (W/m2)
  real(r8), intent(inout):: tref(pcols)          ! 2m reference temperature
  real(r8), intent(inout):: fsns(pcols)          ! SW absorbed in ice
  real(r8), intent(inout):: evap(pcols)
  real(r8), intent(inout):: aiceinit(pcols)

! water tracer variables
  real(r8), intent(in)    :: wtsnowfall(pcols,pcnst+pnats)! total snow rate (m h2o/s) 
  real(r8), intent(inout) :: wtsnowh(pcols,pcnst+pnats)   ! Snow depth (liquid water equivalent)
  real(r8), intent(inout) :: wthi(pcols,pcnst+pnats)      ! Ice thickness


!---------------------------Local variables-----------------------------
! fluxes/quantities over sea ice only
  real(r8) :: tauxice(pcols)       ! X surface stress (N/m2)
  real(r8) :: tauyice(pcols)       ! Y surface stress (N/m2)
  real(r8) :: shflxice(pcols)      ! Surface sensible heat flux (J/m2/s)
  real(r8) :: lhflxice(pcols)      ! Surface latent   heat flux (J/m2/s)
  real(r8) :: trefice(pcols)       ! 2m reference temperature
  real(r8) :: flwup(pcols)

  real(r8) :: Fnet
  real(r8) :: condb         ! conductive flux at the bottom of the ice
  real(r8) :: swbot         ! shortwave passing through the bottom of the ice
  real(r8), dimension(pcols) :: Flwdabs
  real(r8), dimension(pcols) :: Fswabs
  real(r8), dimension(pcols) :: Fswabsv
  real(r8), dimension(pcols) :: Fswabsi
  real(r8), dimension(pcols) :: dflhdT
  real(r8), dimension(pcols) :: dfshdT
  real(r8), dimension(pcols) :: dflwdT
  real(r8) :: Tiz(0:plevmx) ! local 1D ice temperature profile (C)
  real(r8), dimension(pcols) :: Tsfc          ! local ice surface temperature (C)
  real(r8) :: Tbasal        ! ice bottom temp (C)
  real(r8) :: asnow ! snow fractional coverage
  real(r8) :: hs   ! snow thickness 
  real(r8) :: dhs  ! change in snow thickness from melting
  real(r8) :: subs ! change in snow thickness from sublimation/condensation
  real(r8) :: qqq
  real(r8) :: TTT
  real(r8) :: emissivity
  real(r8) :: hi0

  integer :: npts ! number of gridcells with sea ice
  integer :: linpts ! counter for number of ice points reset to linear profile
  integer :: indx(pcols)
  integer :: layer

  real(r8) :: dhi
#if ( defined COUP_SOM )
  real(r8) :: Fbot   ! heat flx to ice bottom only (W/m**2) (<0)
  real(r8) :: Rside
  real(r8) :: vice
! diagnostics for history file 
  real(r8) :: meltb(pcols)
  real(r8) :: meltt(pcols)
  real(r8) :: meltl(pcols)
  real(r8) :: growb(pcols)
  real(r8) :: frazil(pcols)
  real(r8) :: flood(pcols)
#else
! Sea ice thickness is fixed, so
! bottom melt/growth is not computed and the ice-ocean flux 
! is ignored. A flag is set in ice_dh.F accordingly
  real(r8), parameter :: Fbot = 0. ! heat flx to ice bottom  (W/m**2)
#endif

! The following variables are dummies without SOM
  real(r8)  :: dhib, dhit, subi, dhif, dhsf ! changes in ice/snow thickness 
  real(r8)  :: qi(plevmx) ! ice enthalpy

  real(r8)  :: initnrg(pcols)
  real(r8)  :: finlnrg(pcols)
  real(r8)  :: fice2ocn(pcols)
  real(r8)  :: fatm2ice(pcols)
  real(r8)  :: qisum,vsno
  real(r8)  :: nrgerror(pcols)

  real(r8)  :: hssave

  integer  :: i,k,m,ii

!
! Additional storage for water tracers
!
  real(r8) dhsfm
  real(r8) fhsadj(pcols)
  real(r8) fhsred(pcols)		! hs reductions due to melt of ice
  real(r8) wtsubi(pcols)
  real(r8) wtsubs(pcols)
  real(r8) wtdhs(pcols)
  real(r8) wtdhsf(pcols)
  real(r8) wtdhsm(pcols)

  real(r8) wtshsv(pcols)                ! saved total snow height

!-----------------------------------------------------------------------
  call t_startf ('seaice')

!
! Check inputs
!
  if (trace_water) then
!!         write (*,*) '----- SEA ICE -------'
!!         write(*,*) 'WTRCICEini:',wtsnowh(3,4),wtsnowh(3,1),snowh(3)
         call wtrc_qchk1('SEAICE_inp','qbot ',ncol, qbot(:,ixh2oq) , qbot(:,1))
         call wtrc_qchk1('SEAICE_inp','snof1',ncol, wtsnowfall(:,ixh2oq), wtsnowfall(:,1))
         call wtrc_qchk1('SEAICE_inp','snofl',ncol, wtsnowfall(:,ixh2oq), snowfall(:))
         call wtrc_qchk1('SEAICE_inp','sict1',ncol, wthi(:,ixh2oq) , wthi(:,1))
         call wtrc_qchk1('SEAICE_inp','sicth',ncol, wthi(:,ixh2oq) , hi(:))
         call wtrc_qchk1('SEAICE_inp','snoh1',ncol, wtsnowh(:,ixh2oq) , wtsnowh(:,1))
         call wtrc_qchk1('SEAICE_inp','snohw',ncol, wtsnowh(:,ixh2oq) , snowh(:))
  end if


#if ( defined COUP_SOM )

!JR Zero out history diagnostics so non-ice points are zero.

  meltb(:) = 0.
  meltt(:) = 0.
  meltl(:) = 0.
  growb(:) = 0.
  frazil(:) = 0.
  flood(:) = 0.
  nrgerror(:) = 0.
  aiceinit(:) = 0.
#endif

  npts = 0
  do i=1,ncol
     if (aice(i) > puny .or. frzmlt(i) > 0.) then
!        write(6,*) '(ice_srf)',c,i,aice(i),hi(i)
!        write(6,*) '(ice_srf) tssub:', i,c,(tssub(i,k),k=1,plevmx)
        npts = npts + 1
        indx(npts) = i
     end if

     if (aice(i) <= puny) then
        aice(i) = 0.
        Tsice(i) = TfrezK
        hi(i) = 0.
!!!!!!! Make sure the following does not 
!!!!!!! modify variables over land.
!!!!!!! Be sure these variable here are for snow 
!!!!!!! and temperature on sea ice. 
!JR Yup they are.
        snowh(i) = 0.  
        wtsnowh(i,:) = 0.               ! also for water tracers

	do k=1,plevmx
          tssub(i,k) = TfrezK 
        end do
     end if

! Save the input snow height for tracer calultaions
     wtshsv(i) = snowh(i)

  end do

  flwup(:)=0.
  shflxice(:)=0.
  lhflxice(:)=0.
  tauxice(:)=0.
  tauyice(:)=0.
  initnrg(:) = 0.
  fice2ocn(:) = 0.

  Focn(:)=0.
  fsns(:) = 0.

  qqq=qqqice
  TTT=TTTice
  emissivity=emissivity_ice

  if ( npts > 0 ) then

     do ii=1,npts
        i = indx(ii)

! Save initial totals for paralleltracer calculations
        wthi   (i,1) = hi(i)
        wtsnowh(i,1) = snowh(i)

        ! Convert temperatures to C, use 1D array
        Tsfc(i) = Tsice(i) - Tffresh
     end do

     call ice_atm_flux (Tsfc, ubot, vbot, tbot,                       &
                     qbot(:,1),thbot   ,zbot    ,pbot  , flwds ,    &
                     swvdr   ,swidr   ,swvdf   ,swidf ,            &
                     alvdr   ,alidr   ,alvdf   ,alidf ,            &
                     qqq, TTT, emissivity,                         & 
                     tauxice ,tauyice ,flwup,                      &
                     shflxice, lhflxice, trefice,                  &
                     Flwdabs, Fswabs, Fswabsv, Fswabsi,            &
                     dflhdT, dfshdT, dflwdT,                       &
                     npts, indx)
     linpts=0
     do ii=1,npts
        i = indx(ii)
        ! Convert temperatures to C, use 1D array
        do k=1,plevmx
           Tiz(k) = tssub(i,k) - Tffresh
           Tiz(k) = min(Tmelz(k),Tiz(k))
        end do
        ! snow temperature is diagnostic so init to surf
        Tiz(0) = Tsfc(i)
     
        ! snow lands on ice no matter what its temperature
        if (prognostic_icesnow) then
           hs =  (snowh(i) + snowfall(i)*dtime)*rhofresh/rhos 
           fhsadj(i) = hs

           if (fixice) hs = min(hs,0.5_r8)  ! do not let the snow get out of hand

           if (fhsadj(i) > 1.e-22) then		! SOMETHING
               fhsadj(i) = hs / fhsadj(i) ! save for tracers (dcn)
           else
               fhsadj(i)= 0.
           end if

        else
           hs =  snowh(i)*rhofresh/rhos 
           if (trace_water) write(*,*) 'ICE_FLX **WARNING** - must have prognostic snoe for isotopes.'
        endif
     
        ! 1 - snow covered area fraction
        asnow = c1-hs/(hs + snowpatch)
     
        !-----------------------------------------------------------------
        ! compute air to ice heat, momentum, radiative and water fluxes 
        !-----------------------------------------------------------------
!        write(6,*) '(ice_srf) T in C',i,c,Tsfc(i),(Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) b4 ice_sfc_flux',Tsfc(i), ubot(i), vbot(i), tbot(i), &
!             qbot(i)    ,thbot(i)   ,zbot(i)  ,pbot(i)  , flwds(i) ,&
!             swvdr(i)   ,swidr(i)   ,swvdf(i) ,swidf(i) , &
!             alvdr(i)   ,alidr(i)   ,alvdf(i) ,alidf(i)
!        write(6,*) '(ice_srf) ',c,i,alvdr(i),alidr(i)   ,alvdf(i)   ,alidf(i)

      dhsfm = 0.
#if ( defined COUP_SOM )
      !-----------------------------------------------------------------
      ! distribute variable frzmlt into basal melt (Fbot<0), 
      ! new ice growth and lateral melt (0<Rside<1). 
      ! Note that Tiz is the state variable (not qi) at this time
      !-----------------------------------------------------------------

        qisum = 0.
        do k=1,plevmx
           qisum = qisum + energ (tiz(k), saltz(k))
        end do
        qisum = qisum/plevmx
        initnrg(i) = (-rLfs*hs + qisum*hi(i))*aice(i) ! negative
!    write(6,*)'i,c,snow1a=',i,c,hs
!    write(6,*)'i,c,snowterm1a=',i,c,-rLfs*hs*aice(i)
!    write(6,*)'i,c,iceterm1a=',i,c,qisum*hi(i)*aice(i)

! if frzmlt >0 grow new ice
! This is saved for isotope tracers.
        dhsfm = hs
        call init_frzmlt (sst(i), tbot(i), frzmlt(i), aice(i), hi(i), hs, &
                       Tsfc(i), Tiz, saltz, Fbot, frazil(i), Rside, dtime)
        dhsfm = hs - dhsfm
!!!  Rside=0.   ! constraining rside=0 neglects lateral melt 
                   ! but still conserves energy 
#endif

!        write(6,*) '(ice_srf) after ice_sfc_flux', &
!             tauxice(i) ,tauyice(i) ,flwup(i), &
!             shflxice(i),lhflxice(i),trefice(i), &
!             Flwdabs(i),Fswabs(i),Fswabsv(i),Fswabsi(i), &
!             dflhdT(i),dfshdT(i),dflwdT(i)

      !-----------------------------------------------------------------
      ! solve heat equation
      !-----------------------------------------------------------------
!        write(6,*) '(ice_srf) Tiz', (Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) sensible heat flux',shflxice(i)
!        write(6,*) '(ice_srf) latent heat flux',lhflxice(i)

        call tstm( dtime, Tmelz, saltz, Tfrez &
                     , aice(i), hi(i), hs &
                     , fswabs(i), Fswabsv(i), Fswabsi(i) &
                     , Flwdabs(i), dflwdT(i), dflhdT(i), dfshdT(i) &
                     , asnow,  Tbasal &
                     , swbot, Fnet, condb, Tsfc(i), Tiz &
                     , flwup(i), lhflxice(i), shflxice(i),linpts)

!        write(6,*) '(ice_srf) Tiz', (Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) sensible heat flux',shflxice(i)
!        write(6,*) '(ice_srf) latent heat flux',lhflxice(i)
      !-----------------------------------------------------------------
      ! compute snow melt and sublimation
      !-----------------------------------------------------------------

        call dh  ( dtime, saltz, Tiz     &
                     , Tbasal, hi(i), hs, Fbot &
                     , Fnet, condb, lhflxice(i) &
                     , dhib, dhit, dhs, subi &
                     , subs, dhif, dhsf, qi, Focn(i), i,c)


!
! If we are not fixing ice thickness then adjust height 
!
        if (.not. fixice) then
           dhi   = subi + dhit + dhib + dhif
           hi(i) = max (0._r8, hi(i)+dhi)
        end if
!
! Save new value of hi for tracers: This is only needed with SOM.
! need to do this AFTER prognostic equation, but ALSO after
! adjustemnts. The assumption here is that any way of growing
! sea ice must have increment at isotope ratio of frozen sea water.
! (one could imagine this is not the case for flooded snow [dhsf])
!
        wthi(i,1) = hi(i)
!
! If we are prognosing snow then adjust the height due to snow melt and
! sublimation
!

        if (prognostic_icesnow) then
           if (.not. fixice) then
              hs = hs + subs + dhs + dhsf
           else 
              hs = hs  + subs + dhs 
           endif
        endif


! For water tracers, save this to compute "post budget" corrections
       hssave = hs

#if ( defined COUP_SOM )
      !-----------------------------------------------------------------
      ! add shortwave passing through ice to Focn and
      ! multiply Focn by aice to convert from flux per unit ice area 
      ! to flux per unit area
      !-----------------------------------------------------------------
!!!!!! EVAP WAS A BUG
        evap(i) = (rhoi*subi + rhos*subs)/dtime 

! fatm2ice does not include the latent heat flux intentionally,
! instead it has the evap*Lvap (yes Lvap not Lvap+Lfus).
! note that mulitplying the terms infatm2ice by aice
! here is a little deceptive
! the atmosphere model is going to multiply this flux by the ice fraction 
! that is spat out after lateral melt and other adjustments
! Hence the full model energy conservations is a bit worse 
        fatm2ice(i) = fswabs(i) + Flwdabs(i) + flwup(i) + shflxice(i)
!      write(6,*)'i,c,aice1=',i,c,aice(i)
!      write(6,*)'i,c,fswabs1=',i,c,fswabs(i)
!      write(6,*)'i,c,lwsum1=',i,c,Flwdabs(i)+flwup(i)
!      write(6,*)'i,c,shflx1=',i,c,shflxice(i)
!      write(6,*)'i,c,evaplvap1=',i,c,evap(i)*Lvap

        fatm2ice(i) = (fatm2ice(i) + evap(i)*Lvap)*aice(i)

        Focn(i) = (Focn(i) + swbot)*aice(i)
        aiceinit(i) = aice(i)
!JRBPB diagnostic for sw absorbed in ice
        fsns(i) = fswabs(i)

        !-----------------------------------------------------------------
        ! PLEASE send the following variables to the history file 
        !-----------------------------------------------------------------
        growb(i) = aice(i)*max(0._r8,dhib)/dtime
        meltb(i) = aice(i)*min(0._r8,dhib)/dtime
        meltt(i) = aice(i)*dhit/dtime
        flood(i) = aice(i)*dhif/dtime

        !-----------------------------------------------------------------
        ! Reduce ice concentration according to an assumed thickness 
        ! distribution. Mimics melting away thin ice (sort of) by
        ! assuming ice is linearly distributed in thickness between 0 and 2h.
        ! Then melting by dhi results in reducing the area by aice*dhi/hi0.
        ! (Would be okay to skip this.)
        !-----------------------------------------------------------------
        if( dhi .lt. c0 ) then  
           hi0 = hi(i)-dhi
           if (hi0.gt.c0) then 
              vsno = aice(i)*hs  ! ADDED
              aice(i) = aice(i) * sqrt( hi(i)/hi0)
              hi(i) = sqrt(hi0*hi(i))

              !JR Modified 11/20 per C. Bitz to fix pathological hs = inf

              if (aice(i).gt.c0) then
                 hs = vsno/aice(i) 
              end if
           endif
        endif

        !-----------------------------------------------------------------
        ! This routine will alter Focn if there is lateral melt.
        ! It is okay to eliminate lateral melt by setting Rside=0 above
        ! Note qi is the current state variable and not Tiz 
        !-----------------------------------------------------------------
!JR Added subscript to Focn
        call lateral_melt (Rside, aice(i), hi(i), hs, qi, &
                         Focn(i), meltl(i), dtime)

        !-----------------------------------------------------------------
        ! PLEASE send meltl, frazil, Focn, & frzmlt to the history file 
        !-----------------------------------------------------------------

        !-----------------------------------------------------------------
        ! Nuke small ice/open water fractions.
        ! Make sure ice thickness is at least hi_min
        ! (has side effect of reducing ice fraction sometimes as
        ! necessary for numerical reasons, so do not remove!)
        ! Finally, recompute Tiz from qi.
        !-----------------------------------------------------------------
        if (aice(i) <= puny) then
           aice(i) = 0.
           hi(i) = 0.
           hs = 0.       ! Added 11/20 per C. Bitz to fix hs = inf
           Tsfc(i) = Tfrez
           do layer=1,ni
              Tiz(layer) = Tfrez
           end do
        else
           vice       = hi(i)*aice(i)
           vsno       = aice(i)*hs  ! ADDED
           hi(i)      = max(hi(i),hi_min)
           aice(i) = vice/hi(i)
           aice(i) = min(aice(i),1.0_r8)
           hs = vsno/aice(i) ! ADDED
           call init_Tprofile(qi,tiz) 
        endif
!
! Save adjusted ice thickness for tracer calulation
!
        wthi(i,1) = hi(i)

        qisum = 0.
        do k=1,plevmx
           qisum = qisum + energ (tiz(k), saltz(k))
        end do
        qisum = qisum/plevmx
        finlnrg(i)  = (-rLfs*hs + qisum*hi(i))*aice(i) ! negative
!    write(6,*)'i,c,snow1b=',i,c,hs
!    write(6,*)'i,c,snowterm1b=',i,c,-rLfs*hs*aice(i)
!    write(6,*)'i,c,iceterm1b=',i,c,qisum*hi(i)*aice(i)

        fice2ocn(i) = -focn(i) - max(frzmlt(i),0._r8) 
!     write(6,*)'i,c,focn1=',i,c,Focn(i)
!     write(6,*)'i,c,frzmlt1=',i,c,max(frzmlt(i),0._r8)

        nrgerror(i) = (finlnrg(i)-initnrg(i))/dtime - fatm2ice(i) - fice2ocn(i)

!     write(6,*) c,i,nrgerror(i),(finlnrg(i)-initnrg(i))/dtime,fatm2ice(i), &
!          fice2ocn(i)

#else
!JR evap bugfix from C. Bitz
!JR adds sublimation over sea ice
        evap(i) = (rhoi*subi + rhos*subs)/dtime   
!JR OLD code     evap(i) = rhos*subs/dtime  
#endif

        snowh(i) = hs*rhos/rhofresh

        ! Convert temperatures to K, filling 2D arrays
        Tsice(i) = Tsfc(i) + Tffresh
        do k=1,plevmx
           tssub(i,k) = Tiz(k) + Tffresh
        end do

! Save terms for water tracer code
       wtsubi(i) = subi
       wtsubs(i) = subs
       wtdhs(i)  = dhs
       wtdhsf(i) = dhsf
       wtdhsm(i) = dhsfm
       if (hssave > 1.e-30) then
         fhsred(i) = hs/hssave		! normally < 0, but > 1 if aice change
       else
         fhsred(i) = 0.
       end if


     end do
     if (linpts.gt.0) then
        write(6,*)'WARNING: ice_tstm ::profile reset ',linpts,&
           ' points at chunck ',c,' see NOTE in ice_tstm.F for more info'
     end if
     ! Update fluxes, sum ice percentages into total flux values 
     ! ALL FLUXES IN ICE MODEL ARE POSITIVE DOWN
     ! CCM DOES NOT USE THIS STANDARD 
     do ii=1,npts
        i = indx(ii)
        Ts(i)     = Tsice(i)
        tref(i)   = trefice(i)
        lwup(i)   = -1.*(flwup(i)-(1-emissivity)*flwds(i))
        shflx(i)  = -1.*shflxice(i)
        lhflx(i)  = -1.*lhflxice(i)
        taux(i)   = -1.*tauxice(i)
        tauy(i)   = -1.*tauyice(i)
        qflx(i,1) = -1.*evap(i)
     end do

  end if  ! if ( npts > 0 ) then

! Set non-water constituent fluxes to zero
! (these will be overwritten below for water tracers)

  do m=2,pcnst+pnats
     do ii=1,npts
        i = indx(ii)
        qflx(i,m) = 0.
     end do
  end do

! Water tracers (overwrite appropriate water tracers):
!   compute new ice thickness, qflx, and update snow depth
!
  if (trace_water) then
         call wtrc_qchk1('SEAICE_wti','snohi',ncol, wtsnowh(:,ixh2oq) , wtsnowh(:,1) )
         call wtrc_qchk1('SEAICE_wti','snof1',ncol, wtsnowfall(:,ixh2oq), wtsnowfall(:,1))
         call wtrc_qchk1('SEAICE_wti','snofl',ncol, wtsnowfall(:,ixh2oq), snowfall(:))
         call wtrc_qchk1('SEAICE_wti','qbot ',ncol, qbot(:,ixh2oq) , qbot(:,1))


     call wtrc_ice_srf( npts    , indx   , dtime  , &
               wtsubi , wtsubs  , wtdhs  , wtdhsf , wtdhsm  , wtsnowfall, &
               fhsadj , fhsred  , qbot   , wtsnowh, wthi   , qflx  )

         call wtrc_qchk1('SEAICE_wto','qbot ',ncol, qbot(:,ixh2oq) , qbot(:,1))
         call wtrc_qchk1('SEAICE_wto','snof1',ncol, wtsnowfall(:,ixh2oq), wtsnowfall(:,1))
         call wtrc_qchk1('SEAICE_wto','sict1',ncol, wthi(:,ixh2oq) , wthi(:,1))
         call wtrc_qchk1('SEAICE_wto','snoh1',ncol, wtsnowh(:,ixh2oq) , wtsnowh(:,1), 1.e-14)

! Rescale the prognostcs, just to be sure  (this is only for numerical precision)
     call wtrc_ice_checkscale(npts,indx,snowh,wtsnowh)

         call wtrc_qchk1('SEAICE_wto','snohw',ncol, wtsnowh(:,ixh2oq) , snowh(:))
         call wtrc_qchk1('SEAICE_wto','sicth',ncol, wthi(:,ixh2oq) , hi(:))
         call wtrc_qchk1('SEAICE_wto','snofl',ncol, wtsnowfall(:,ixh2oq), snowfall(:))

!!     write(*,*) 'WTRCICEout:',wtsnowh(3,4),wtsnowh(3,1),snowh(3)
  end if


#if ( defined COUP_SOM )
  call outfld ('MELTB   ', meltb,    pcols, c)
  call outfld ('MELTT   ', meltt,    pcols, c)
  call outfld ('MELTL   ', meltl,    pcols, c)
  call outfld ('GROWB   ', growb,    pcols, c)
  call outfld ('FRAZIL  ', frazil,   pcols, c)
  call outfld ('FLOOD   ', flood,    pcols, c)
  call outfld ('FRZMLT  ', frzmlt,   pcols, c)
  call outfld ('NRGERROR', nrgerror, pcols, c)
#endif
  call outfld ('TSICE   ', tsice,  pcols, c)
  call t_stopf ('seaice')

  return

end subroutine seaice
