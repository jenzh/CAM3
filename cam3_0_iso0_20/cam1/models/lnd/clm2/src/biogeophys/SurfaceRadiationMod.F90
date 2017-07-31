#include <misc.h>
#include <preproc.h>

module SurfaceRadiationMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceRadiationMod
!
! !DESCRIPTION:
! Calculate solar fluxes absorbed by vegetation and ground surface
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: SurfaceRadiation ! Solar fluxes absorbed by veg and ground surface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRadiation
!
! !INTERFACE:
   subroutine SurfaceRadiation(lbp, ubp)
!
! !DESCRIPTION:
! Solar fluxes absorbed by vegetation and ground surface.  Note possible
! problem when land is on different grid than atmosphere.  Land may have
! sun above the horizon (coszen $< 0$) but atmosphere may have sun below
! the horizon (forc\_solad $= 0$ and forc\_solai $= 0$). This is okay
! because all fluxes (absorbed, reflected, transmitted) are multiplied
! by the incoming flux and all will equal zero.  Atmosphere may have sun
! above horizon (forc\_solad $> 0$ and forc\_solai $> 0$) but land may
! have sun below horizon. This is okay because fabd, fabi, ftdd, ftid,
! and ftii all equal zero so that sabv = sabg = fsa = 0. Also, albd and
! albi equal one so that fsr = forc\_solad+forc\_solai. In other words,
! all the radiation is reflected. However, the way the code is currently
! implemented, this is only true if (forc\_solad + forc\_solai)|vis =
! (forc\_solad + forc\_solai)|nir.
!
! !USES:
     use clmtype
     use clm_varpar  , only : numrad
     use clm_varcon  , only : spval
     use time_manager, only : get_curr_date, get_step_size
!
! !ARGUMENTS:
     implicit none
     integer, intent(in) :: lbp, ubp      ! pft upper and lower bounds
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/18/02, Peter Thornton: Migrated to new data structures. Added a pft loop.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
     integer , pointer :: pcolumn(:)       ! pft's column index
     integer , pointer :: pgridcell(:)     ! pft's gridcell index
     real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
     real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
     real(r8), pointer :: esai(:)          ! one-sided stem area index with burying by snow
     real(r8), pointer :: londeg(:)        ! longitude (degrees)
     real(r8), pointer :: latdeg(:)        ! latitude (degrees)
     real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (W/m**2)
     real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation (W/m**2)
     real(r8), pointer :: fabd(:,:)        ! flux absorbed by veg per unit direct flux
     real(r8), pointer :: fabi(:,:)        ! flux absorbed by veg per unit diffuse flux
     real(r8), pointer :: ftdd(:,:)        ! down direct flux below veg per unit dir flx
     real(r8), pointer :: ftid(:,:)        ! down diffuse flux below veg per unit dir flx
     real(r8), pointer :: ftii(:,:)        ! down diffuse flux below veg per unit dif flx
     real(r8), pointer :: albgrd(:,:)      ! ground albedo (direct)
     real(r8), pointer :: albgri(:,:)      ! ground albedo (diffuse)
     real(r8), pointer :: albd(:,:)        ! surface albedo (direct)
     real(r8), pointer :: albi(:,:)        ! surface albedo (diffuse)
!
! local pointers to original implicit out arguments
!
     real(r8), pointer :: laisun(:)        ! sunlit leaf area
     real(r8), pointer :: laisha(:)        ! shaded leaf area
     real(r8), pointer :: sabg(:)          ! solar radiation absorbed by ground (W/m**2)
     real(r8), pointer :: sabv(:)          ! solar radiation absorbed by vegetation (W/m**2)
     real(r8), pointer :: fsa(:)           ! solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: parsun(:)        ! average absorbed PAR for sunlit leaves (W/m**2)
     real(r8), pointer :: parsha(:)        ! average absorbed PAR for shaded leaves (W/m**2)
     real(r8), pointer :: fsr(:)           ! solar radiation reflected (W/m**2)
     real(r8), pointer :: fsds_vis_d(:)    ! incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_d(:)    ! incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_i(:)    ! incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i(:)    ! incident diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_d(:)     ! reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_d(:)     ! reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_vis_i(:)     ! reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_i(:)     ! reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_vis_d_ln(:) ! incident direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsds_nir_d_ln(:) ! incident direct beam nir solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_vis_d_ln(:)  ! reflected direct beam vis solar rad at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_ln(:)  ! reflected direct beam nir solar rad at local noon (W/m**2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
     integer , parameter :: nband = numrad ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06   ! prevents overflow for division by zero
     integer  :: p                   ! pft index
     integer  :: c                   ! column index
     integer  :: g                   ! grid cell index
     integer  :: ib                  ! waveband number (1=vis, 2=nir)
     real(r8) :: abs                 ! absorbed solar radiation (W/m**2)
     real(r8) :: rnir                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: laifra              ! leaf area fraction of canopy
     real(r8) :: trd                 ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri                 ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(lbp:ubp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(lbp:ubp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     real(r8) :: fsha                ! shaded fraction of canopy
     real(r8) :: vai                 ! total leaf area index + stem area index, one sided
     integer  :: local_secp1         ! seconds into current date in local time
     real(r8) :: dtime               ! land model time step (sec)
     integer  :: year,month,day      ! temporaries (not used)
     integer  :: secs                ! seconds into current date
!------------------------------------------------------------------------------

     ! Assign local pointers to multi-level derived type members (gridcell level)

     londeg        => clm3%g%londeg
     latdeg        => clm3%g%latdeg
     forc_solad    => clm3%g%a2lf%forc_solad
     forc_solai    => clm3%g%a2lf%forc_solai

     ! Assign local pointers to multi-level derived type members (column level)

     albgrd        => clm3%g%l%c%cps%albgrd
     albgri        => clm3%g%l%c%cps%albgri

     ! Assign local pointers to derived type members (pft-level)

     pcolumn       => clm3%g%l%c%p%column
     pgridcell     => clm3%g%l%c%p%gridcell
     fsun          => clm3%g%l%c%p%pps%fsun
     elai          => clm3%g%l%c%p%pps%elai
     esai          => clm3%g%l%c%p%pps%esai
     laisun        => clm3%g%l%c%p%pps%laisun
     laisha        => clm3%g%l%c%p%pps%laisha
     fabd          => clm3%g%l%c%p%pps%fabd
     fabi          => clm3%g%l%c%p%pps%fabi
     ftdd          => clm3%g%l%c%p%pps%ftdd
     ftid          => clm3%g%l%c%p%pps%ftid
     ftii          => clm3%g%l%c%p%pps%ftii
     albd          => clm3%g%l%c%p%pps%albd
     albi          => clm3%g%l%c%p%pps%albi
     sabg          => clm3%g%l%c%p%pef%sabg
     sabv          => clm3%g%l%c%p%pef%sabv
     fsa           => clm3%g%l%c%p%pef%fsa
     fsr           => clm3%g%l%c%p%pef%fsr
     parsun        => clm3%g%l%c%p%pef%parsun
     parsha        => clm3%g%l%c%p%pef%parsha
     fsds_vis_d    => clm3%g%l%c%p%pef%fsds_vis_d
     fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
     fsds_vis_i    => clm3%g%l%c%p%pef%fsds_vis_i
     fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
     fsr_vis_d     => clm3%g%l%c%p%pef%fsr_vis_d
     fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
     fsr_vis_i     => clm3%g%l%c%p%pef%fsr_vis_i
     fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
     fsds_vis_d_ln => clm3%g%l%c%p%pef%fsds_vis_d_ln
     fsds_nir_d_ln => clm3%g%l%c%p%pef%fsds_nir_d_ln
     fsr_vis_d_ln  => clm3%g%l%c%p%pef%fsr_vis_d_ln
     fsr_nir_d_ln  => clm3%g%l%c%p%pef%fsr_nir_d_ln

     ! Determine seconds off current time step

     dtime = get_step_size()
     call get_curr_date (year, month, day, secs)

     ! Determine fluxes

!dir$ concurrent
!cdir nodep
     do p = lbp,ubp
        sabg(p) = 0._r8
        sabv(p) = 0._r8
        fsa(p)  = 0._r8
     end do

     ! Loop over nband wavebands
     do ib = 1, nband
!dir$ concurrent
!cdir nodep
        do p = lbp,ubp
           c = pcolumn(p)
           g = pgridcell(p)

           ! Absorbed by canopy

           cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
           cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
           sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
           fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)

           ! Transmitted = solar fluxes incident on ground

           trd = forc_solad(g,ib)*ftdd(p,ib)
           tri = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)

           ! Solar radiation absorbed by ground surface

           abs = trd*(1.-albgrd(c,ib)) + tri*(1.-albgri(c,ib))
           sabg(p) = sabg(p) + abs
           fsa(p)  = fsa(p)  + abs

        end do ! end of pft loop
     end do ! end of nband loop

!dir$ concurrent
!cdir nodep
     do p = lbp,ubp
        g = pgridcell(p)

        ! Begin calculations

        fsha = 1.-fsun(p)
        laisun(p) = elai(p)*fsun(p)
        laisha(p) = elai(p)*fsha
        vai = elai(p)+ esai(p)

        ! Partition visible canopy absorption to sunlit and shaded fractions
        ! to get average absorbed par for sunlit and shaded leaves

        laifra = elai(p) / max(vai,mpe)
        if (fsun(p) > 0._r8) then
           parsun(p) = (cad(p,1) + cai(p,1)) * laifra
        else
           parsun(p) = 0._r8
        end if
        parsha(p) = 0._r8

        ! Reflected solar radiation

        rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
        rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
        fsr(p) = rvis + rnir

        fsds_vis_d(p) = forc_solad(g,1)
        fsds_nir_d(p) = forc_solad(g,2)
        fsds_vis_i(p) = forc_solai(g,1)
        fsds_nir_i(p) = forc_solai(g,2)
        fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
        fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
        fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
        fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)

        local_secp1 = secs + nint((londeg(g)/15.*3600.)/dtime)*dtime
        local_secp1 = mod(local_secp1,86400)
        if (local_secp1 == 43200) then
           fsds_vis_d_ln(p) = forc_solad(g,1)
           fsds_nir_d_ln(p) = forc_solad(g,2)
           fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
           fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
        else
           fsds_vis_d_ln(p) = spval
           fsds_nir_d_ln(p) = spval
           fsr_vis_d_ln(p) = spval
           fsr_nir_d_ln(p) = spval
        end if

     end do

   end subroutine SurfaceRadiation

end module SurfaceRadiationMod
