!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModConductance

  use ModIonosphere

  implicit none
  save

  ! Parameters:
  real, parameter :: &
       IONO_Min_EFlux = 0.1e-16, &  ! W/m2
       IONO_Min_Ave_E = 0.5,     &  ! keV
       cFACFloor      = 1.0E-12, &  ! A/m2
       SigmaPar       = 1000.0      ! Parallel conductance, Siemens

  ! Logicals to control what conductance sources are used.
  logical :: DoUseEuvCond=.false., DoUseAurora=.false.

  ! Name of auroral model to use, defaults to 'FAC2FLUX'
  character(len=8) :: NameAuroraMod = 'FAC2FLUX'

  ! Name of empirical conductance model, if used:
  character(len=4) :: NameEmpModel='RLM5' ! Legacy iModel=3,4,5 or CMEE

  ! Background & constant conductance values:
  real :: f107_flux=150, SigmaHalConst=0, SigmaPedConst=0, &
       StarLightCond=0.25, & !replaces starlightpedconductance
       PolarCapPedCond=0.25  !replaces PolarCapPedConductance

  ! Floor values for GM density and pressure, SI units:
  real, parameter :: GmRhoFloor = 1E-21, GmPFloor = 1E-13

  ! Arrays to hold components of conductance - 1 array
  ! for each source of conductance.  Long term storage for output files.
  ! THESE ARE ONLY NEEDED FOR PUTTING OUTPUT IN OLD FILE TYPES.
  ! If we make new output files that are written here, we could avoid this.
  real, dimension(IONO_nTheta,IONO_nPsi) :: &
       SigHalEuv_NORTH=0.0, SigHalEuv_SOUTH=0.0, &
       SigPedEuv_NORTH=0.0, SigPedEuv_SOUTH=0.0
       ! Existing variables that should be moved here
       ! Also, no reason to allocate (we think).
       !IONO_NORTH_DIFF_Ave_E=0,0, IONO_SOUTH_DIFF_Ave_E=0.0  &
       !IONO_NORTH_MONO_Ave_E=0.0, IONO_SOUTH_MONO_Ave_E=0.0, &
       !IONO_NORTH_DIFF_EFlux=0.0, IONO_SOUTH_DIFF_EFlux=0.0, &
       !IONO_NORTH_MONO_EFlux=0.0, IONO_SOUTH_MONO_EFlux=0.0

contains
  !===========================================================================
  subroutine generate_conductance(NameHemiIn)
    ! Calculate conductance from all sources and sum into the IONO_*_Sigma[PH]
    ! arrays, and compute the spatial derivatives.
    ! Arguments:
    !     NameHemiIn -- Select hemisphere (north, south).

    ! Dummy variables:
    character(len=5), intent(in) :: NameHemiIn

    ! Local variables:
    integer :: iBlockNow, nTheta, nPsi, i, j

    ! Local containers for grid in current hemisphere:
    real, dimension(IONO_nTheta, IONO_nPsi) :: psi, theta

    ! Local containers for precipitation:
    real, dimension(IONO_nTheta, IONO_nPsi) :: &
         AvgEDiffe_II=0.0, EfluxDiffe_II=0.0,  &
         AvgEDiffi_II=0.0, EfluxDiffi_II=0.0,  &
         AvgEMono_II=0.0,  EfluxMono_II=0.0,   &
         AvgEBbnd_II=0.0,  EfluxBbnd_II=0.0

    ! Local container for hemispheric grid spacings:
    real, dimension(IONO_nTheta) :: dTheta
    real, dimension(IONO_nPsi)   :: dPsi

    ! Local containers for conductance components:
    real, dimension(IONO_nTheta,IONO_nPsi) ::      &
         SigmaHalEuv_II =0.0, SigmaPedEuv_II =0.0, &
         SigmaHalBbnd_II=0.0, SigmaPedBbnd_II=0.0, &
         SigmaHalAur_II =0.0, SigmaPedAur_II =0.0

    ! Local containers for spatial derivatives:
    real, dimension(IONO_nTheta,IONO_nPsi) :: &
         SigmaH, SigmaP, &
         SigmaThTh, SigmaThPs, SigmaPsPs,     &
         dSigmaThTh_dTheta, dSigmaThTh_dPsi, dSigmaThPs_dTheta, &
         dSigmaThPs_dPsi, dSigmaPsPs_dTheta, dSigmaPsPs_dPsi

    ! Intermediate values:
    real, dimension(IONO_nTheta,IONO_nPsi)  :: &
         sn, cs, sn2, cs2, cs3, cs4, C

    ! Debug variables:
    character(len=*), parameter :: NameSub='generate_conductance'
    logical :: DoTestMe, DoTest
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe) &
         write(*,*)'IE: '//NameSub//' called for hemisphere='//NameHemiIn

    ! Set some convenience variables to prevent very long lines:
    nPsi = IONO_nPsi
    nTheta = IONO_nTheta

    ! Configure approach to match hemisphere:
    if(NameHemiIn .eq. 'north') then
       iBlockNow = 1
       dTheta = dTheta_North
       dPsi   = dPsi_North
    else if (NameHemiIn .eq. 'south') then
       iBlockNow = 2
       dTheta = dTheta_South
       dPsi   = dPsi_South
    else
       call CON_stop(NameSub//': Unrecognized hemisphere -- '//NameHemiIn)
    end if

    ! Create temporary arrays to hold colat/lon for current hemisphere:
    if(NameHemiIn .eq. 'north')then
       theta = IONO_NORTH_Theta
       psi = IONO_NORTH_Psi
    else
       theta = IONO_SOUTH_Theta
       psi = IONO_SOUTH_Psi
    end if

    ! Calculate each source/contribution to conductance individually:
    ! Add EUV conductance
    if (DoUseEuvCond) &
         call calc_euv_cond(NameHemiIn, SigmaHalEuv_II, SigmaPedEuv_II)

    ! Add aurora conductances.  Aurora models obtain/calculate average energy
    ! and energy flux for diffuse, discrete, and broadband (all default to
    ! zero).
    if (DoUseAurora) then
       ! Obtain average energy and energy flux values based on
       ! selected aurora model:
       select case (NameAuroraMod)
       case('FAC2FLUX')
          ! RLM/CMEE fills the diffuse ONLY
          call facs_to_fluxes(NameEmpModel, NameHemiIn, &
               AvgEMono_II, EfluxMono_II)
       case('MAGNIT')
          !call magnit_fluxes(NameHemiIn)
       case default
          call CON_stop(NameSub//': Unrecognized auroral model - ' &
               //NameAuroraMod)
       end select

       ! use robinson here to fill aurora conductance.
       ! SigmaHalAur = ConvertFluxSigma(AvgEDiff, EfluxDiff, 'robinson') + gallrich(AvgEIdif,
    end if

    ! Sum conductance into the correct hemisphere.
    if(NameHemiIn .eq. 'north')then
       IONO_NORTH_Sigma0 = SigmaPar
       IONO_NORTH_SigmaH = sqrt(SigmaHalConst**2 + SigmaHalEuv_II**2 + &
            (2.*StarLightCond)**2 + SigmaHalAur_II**2)
       IONO_NORTH_SigmaP = sqrt(SigmaPedConst**2 + StarLightCond**2 +&
            SigmaPedEuv_II**2 + SigmaPedAur_II**2)
       ! Add broadband conductance:
       IONO_NORTH_SigmaH = IONO_NORTH_SigmaH + SigmaHalBbnd_II
       IONO_NORTH_SigmaP = IONO_NORTH_SigmaP + SigmaPedBbnd_II
       ! Place values into convenience arrays to calculate derivatives:
       SigmaH = IONO_NORTH_SigmaH
       sigmaP = IONO_NORTH_SigmaP
    else
       IONO_SOUTH_Sigma0 = SigmaPar
       IONO_SOUTH_SigmaH = sqrt(SigmaHalConst**2 + SigmaHalEuv_II**2 + &
            (2.*StarLightCond)**2)
       IONO_SOUTH_SigmaP = sqrt(SigmaPedConst**2 + SigmaPedEuv_II**2 + &
            StarLightCond **2)
       ! Add broadband conductance:
       IONO_SOUTH_SigmaH = IONO_SOUTH_SigmaH + SigmaHalBbnd_II
       IONO_SOUTH_SigmaP = IONO_SOUTH_SigmaP + SigmaPedBbnd_II
       ! Place values into convenience arrays to calculate derivatives:
       SigmaH = IONO_SOUTH_SigmaH
       sigmaP = IONO_SOUTH_SigmaP
    end if

    ! Calculate of-diagonal conductance values.
    ! Start with some "helper" values:
    sn = sin(Theta)
    cs = cos(Theta)
    sn2= sn**2
    cs2 = cs**2
    cs3 = 1.00 + 3.00*cs2
    cs4 = sqrt(cs3)
    C = 4.00*SigmaPar*cs2 + SigmaP*sn2

    ! Off-diagonal vals:
    SigmaThTh = SigmaPar*SigmaP*cs3/C
    SigmaThPs = 2.00*SigmaPar*SigmaH*cs*cs4/C
    SigmaPsPs = SigmaP + SigmaH*SigmaH*sn2/C

    ! Intialize derivatives to zero
    dSigmaThTh_dTheta = 0.00
    dSigmaThTh_dPsi   = 0.00
    dSigmaThPs_dTheta = 0.00
    dSigmaThPs_dPsi   = 0.00
    dSigmaPsPs_dTheta = 0.00
    dSigmaPsPs_dPsi   = 0.00

    ! Spatial derivatives of conductance.  Note that dTheta and dPsi
    ! are defined to be the separation over two grid cells, so all
    ! derivatives are central difference formula correct.
    ! Poles and equator are excluded from this calculation.
    ! d/dTheta derivatives:
    do j=1, nPsi
       dSigmaThTh_dTheta(2:nTheta-1,j) = &
            (SigmaThTh(3:nTheta,j) - SigmaThTh(1:nTheta-2,j)) / dTheta(2:nTheta-1)
       dSigmaThPs_dTheta(2:nTheta-1,j) = &
            (SigmaThPs(3:nTheta,j) - SigmaThPs(1:nTheta-2,j)) / dTheta(2:nTheta-1)
       dSigmaPsPs_dTheta(2:nTheta-1,j) = &
            (SigmaPsPs(3:nTheta,j) - SigmaPsPs(1:nTheta-2,j)) / dTheta(2:nTheta-1)
    end do

    ! d/dPsi derivatives:
    do i=2, nTheta-1
       dSigmaThTh_dPsi(i,2:nPsi-1) = &
            (SigmaThTh(i,3:nPsi)-SigmaThTh(i,1:nPsi-2)) &
            / dPsi(2:nPsi-1)
       dSigmaThPs_dPsi(i,2:nPsi-1) = &
            (SigmaThPs(i,3:nPsi)-SigmaThPs(i,1:nPsi-2)) &
            / dPsi(2:nPsi-1)
       dSigmaPsPs_dPsi(i,2:nPsi-1) = &
            (SigmaPsPs(i,3:nPsi)-SigmaPsPs(i,1:nPsi-2)) &
            / dPsi(2:nPsi-1)
    end do
    ! Longitude edge cases (compensating for ghost cell) j=1 and j=nPsi:
    dSigmaThTh_dPsi(2:nTheta-1,1) = &
         (SigmaThTh(2:nTheta-1,2)-SigmaThTh(2:nTheta-1,nPsi-1)) / dPsi(1)
    dSigmaThPs_dPsi(2:nTheta-1,1) = &
         (SigmaThPs(2:nTheta-1,2)-SigmaThPs(2:nTheta-1,nPsi-1)) / dPsi(1)
    dSigmaPsPs_dPsi(2:nTheta-1,1) = &
         (SigmaPsPs(2:nTheta-1,2)-SigmaPsPs(2:nTheta-1,nPsi-1)) / dPsi(1)

    dSigmaThTh_dPsi(2:nTheta-1,nPsi) = &
         (SigmaThTh(2:nTheta-1,2)-SigmaThTh(2:nTheta-1,nPsi-1)) / dPsi(nPsi)
    dSigmaThPs_dPsi(2:nTheta-1,nPsi) = &
         (SigmaThPs(2:nTheta-1,2)-SigmaThPs(2:nTheta-1,nPsi-1)) / dPsi(nPsi)
    dSigmaPsPs_dPsi(2:nTheta-1,nPsi) = &
         (SigmaPsPs(2:nTheta-1,2)-SigmaPsPs(2:nTheta-1,nPsi-1)) / dPsi(nPsi)

    ! Place derivative results into correct hemisphere:
    if(NameHemiIn.eq.'north')then
       ! Off-diagonal terms:
       IONO_NORTH_SigmaThTh = SigmaThTh
       IONO_NORTH_SigmaThPs = SigmaThPs
       IONO_NORTH_SigmaPsPs = SigmaPsPs
       ! Spatial derivative terms:
       IONO_NORTH_dSigmaThTh_dTheta = dSigmaThTh_dTheta
       IONO_NORTH_dSigmaThPs_dTheta = dSigmaThPs_dTheta
       IONO_NORTH_dSigmaPsPs_dTheta = dSigmaPsPs_dTheta
       IONO_NORTH_dSigmaThTh_dPsi = dSigmaThTh_dPsi
       IONO_NORTH_dSigmaThPs_dPsi = dSigmaThPs_dPsi
       IONO_NORTH_dSigmaPsPs_dPsi = dSigmaPsPs_dPsi
    else
       ! Off-diagonal terms:
       IONO_SOUTH_SigmaThTh = SigmaThTh
       IONO_SOUTH_SigmaThPs = SigmaThPs
       IONO_SOUTH_SigmaPsPs = SigmaPsPs
       ! Spatial derivative terms:
       IONO_SOUTH_dSigmaThTh_dTheta = dSigmaThTh_dTheta
       IONO_SOUTH_dSigmaThPs_dTheta = dSigmaThPs_dTheta
       IONO_SOUTH_dSigmaPsPs_dTheta = dSigmaPsPs_dTheta
       IONO_SOUTH_dSigmaThTh_dPsi = dSigmaThTh_dPsi
       IONO_SOUTH_dSigmaThPs_dPsi = dSigmaThPs_dPsi
       IONO_SOUTH_dSigmaPsPs_dPsi = dSigmaPsPs_dPsi

    end if

  end subroutine generate_conductance

  !===========================================================================
  subroutine calc_euv_cond(NameHemiIn, CondHalOut_II, CondPedOut_II)

    ! For a given hemisphere, calculate the EUV-driven conductance using
    ! Moen & Brekke, 1992.  Add scattering across the terminator as from
    ! Rasmussen and Schunk model B (JGR, Vol. 92, pp. 4491-4504, 1987).
    ! The two approaches are blended at a set solar zenith angle (blendAngle).
    ! Following Ridley_serial's parallelization via hemisphere, act only on 1
    ! hemisphere at a time.
    ! Starlight conductance is included in this calculation.

    !use IE_ModSize
    use ModNumConst, ONLY: cDegToRad
    use IE_ModMain,  ONLY: CosThetaTilt, SinThetaTilt

    character(len=*), intent(in) :: NameHemiIn
    real, intent(out), dimension(IONO_nTheta,IONO_nPsi) :: &
         CondHalOut_II, CondPedOut_II

    integer :: i, j

    real, parameter :: cBlendAngle = 70.0*cDegToRad
    real :: f107p53, f107p49, cos_limit, meeting_value_p, meeting_value_h
    real, dimension(IONO_nTheta,IONO_nPsi) :: cos_SZA, X_II, Y_II, Z_II, &
         SigmaH_EUV, SigmaH_SCAT, SigmaP_EUV, SigmaP_SCAT

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='EUV_cond_II'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Grab coordinates/locations based on hemisphere.
    if(NameHemiIn .eq. 'north') then
       X_II = IONO_NORTH_X
       Y_II = IONO_NORTH_Y
       Z_II = IONO_NORTH_Z
    else if (NameHemiIn .eq. 'south') then
       X_II = IONO_SOUTH_X
       Y_II = IONO_SOUTH_Y
       Z_II = IONO_SOUTH_Z
    else
       call CON_stop(NameSub//': Unrecognized hemisphere = '//NameHemiIn)
    endif

    ! Calculate the cosine of the solar zenith angle:
    cos_SZA = (X_II*CosThetaTilt - Z_II*SinThetaTilt) &
         / sqrt(X_II**2 + Y_II**2 + Z_II**2)

    ! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
    ! So, let's just store them straight away:
    f107p53 = f107_flux**0.53
    f107p49 = f107_flux**0.49

    ! Point at which to merge two functions:
    cos_limit = cos(cBlendAngle)

    ! "meeting values" are M&B at SZA=blendAngle degrees.
    meeting_value_p = f107p49*(0.34*cos_limit+0.93*sqrt(cos_limit))
    meeting_value_h = f107p53*(0.81*cos_limit+0.54*sqrt(cos_limit))

    dayside: where (cos_SZA > 0)
       SigmaH_EUV=f107p53*(0.81*cos_SZA+0.54*sqrt(cos_SZA)) ! M & B.
       SigmaP_EUV=f107p49*(0.34*cos_SZA+0.93*sqrt(cos_SZA))
       SigmaH_SCAT = 1.00 ! Scattered directly from R & S
       SigmaP_SCAT = 0.50
       blend: where (cos_SZA < cos_limit)
          SigmaH_EUV = (SigmaH_EUV + &
               meeting_value_h * exp(-((cos_SZA-cos_limit)**2.0)*15.0))/2.0
          SigmaP_EUV = (SigmaP_EUV + &
               meeting_value_p * exp(-((cos_SZA-cos_limit)**2.0)*15.0))/2.0
       endwhere blend
    elsewhere ! nightside:
       SigmaH_EUV = meeting_value_h * exp(-((cos_SZA-cos_limit)**2.0)*15.0)
       SigmaP_EUV = meeting_value_p * exp(-((cos_SZA-cos_limit)**2.0)*15.0)
       SigmaH_SCAT = 1.00*(10.00**cos_SZA) ! Follows R & S
       SigmaP_SCAT = 0.50*(10.00**cos_SZA) ! Follows R & S
    endwhere dayside

    ! Sum the EUV and scattering conductances together:
    CondHalOut_II = sqrt(SigmaH_EUV**2 + SigmaH_SCAT**2)
    CondPedOut_II = sqrt(SigmaP_EUV**2 + SigmaP_SCAT**2)

  end subroutine calc_euv_cond

  !===========================================================================
  subroutine smooth_lagrange_polar(a_II, iLatSize, jLonSize, tolIn)
    ! Use a simple sliding window technique to smooth a 2D array in polar
    ! coordinates such that periodicity and continuity are enforced across
    ! the LON=360/0 boundary. The max and min latitude points are not
    ! considered.  Only points that are more than *tol* percent greater than
    ! the average of their neighbors are changed.
    !
    ! Arguments
    !
    ! a_II(iLatSize,jLonSize): A two-dimentional array to smooth in place.
    ! iLatSize, jLonSize     : The size of a_II in lat/lon space

    integer, intent(in)        :: iLatSize, jLonSize
    real, intent(inout)        :: a_II(iLatSize, jLonSize)
    real, intent(in), optional :: tolIn

    integer                    :: i, j
    real                       :: aSmooth_II(iLatSize, jLonSize), tol
    !------------------------------------------------------------------------

    ! Set default tolerance to 15% of original value:
    if (present(tolIn)) then
       tol = tolIn
    else
       tol=0.15
    endif

    ! Construct smoothed values
    ! No smoothing at poles:
    aSmooth_II(1,:) = a_II(1,:)
    aSmooth_II(iLatSize,:) = a_II(iLatSize,:)
    ! Smooth all non-pole points:
    colat: do i = 2, iLatSize-1
       ! Perform averaging at longitude boundaries (azimuthal periodicity)
       aSmooth_II(i,1) = (                  & ! At lon=0 degrees
            sum(a_II(i-1:i+1, jLonSize)) +  &
            sum(a_II(i-1:i+1, 2))        +  &
            a_II(i-1,1) + a_II(i+1,1)       &
            ) / 8.0

       aSmooth_II(i,jLonSize) = (              & ! At lon=360 degrees
            sum(a_II(i-1:i+1, jLonSize -1)) +  &
            sum(a_II(i-1:i+1, 1))           +  &
            a_II(i-1,jLonSize)              +  &
            a_II(i+1,jLonSize) ) / 8.0

       ! Smooth over all non-longitude boundary points:
       lon: do j = 2, jLonSize-1
          aSmooth_II(i,j) = (             & ! For all other (i,j)
               sum(a_II(i-1:i+1, j-1)) +  & !left of point i,j
               sum(a_II(i-1:i+1, j+1)) +  & !right of point i,j
               a_II(i+1,j) + a_II(i-1,j)  & !above and below
               )/ 8.0
       end do lon

    end do colat

    ! Apply smoothing to discrete_nf only where smoothing changes values
    ! by tolerance percentage *tol* or more from their original value:
    where (a_II<=ABS(aSmooth_II - a_II) /tol) a_II = aSmooth_II

  end subroutine smooth_lagrange_polar

  !===========================================================================
  subroutine imodel_legacy(iModelIn, f10In, StarLightIn, PolarCapIn)
    ! Create mapping between new conductance-related PARAMs and legacy
    ! #IONOSPHERE param with the "imodel" parameter.  This maintains
    ! backwards compatability.

    integer, intent(in) :: iModelIn
    real, intent(in) :: f10In, StarLightIn, PolarCapIn

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='imodel_legacy'
    !------------------------------------------------------------------------
    select case(iModelIn)
       case(0) ! Constant conductance, Pedersen only; no background
          SigmaPedConst = StarLightIn
          SigmaHalConst = 0
          StarLightCond  = 0.0
          PolarCapPedCond= 0.0
       case(1) ! Constant conductance; no background
          SigmaHalConst  = PolarCapIn
          SigmaPedConst  = StarLightIn
          StarLightCond  = 0.0
          PolarCapPedCond= 0.0
       case(2) ! Constant conductance + dayside conductance:
          SigmaHalConst = PolarCapIn
          SigmaPedConst = StarLightIn
          DoUseEuvCond = .true.
          f107_flux = f10In
       case(3)
          DoUseAurora     = .true.
          NameEmpModel    = 'RLM3'
          DoUseEuvCond    = .true.
          StarLightCond   = StarLightIn
          PolarCapPedCond = PolarCapIn
          f107_flux       = f10In
       case(4)
          DoUseAurora     = .true.
          NameEmpModel    = 'RLM4'
          DoUseEuvCond    = .true.
          StarLightCond   = StarLightIn
          PolarCapPedCond = PolarCapIn
          f107_flux       = f10In
       case(5)
          DoUseAurora     = .true.
          NameEmpModel    = 'RLM5'
          DoUseEuvCond    = .true.
          StarLightCond   = StarLightIn
          PolarCapPedCond = PolarCapIn
          f107_flux       = f10In
       case(6) ! Legacy conductance from MHD; DOES NOT FUNCTION.
          call CON_stop(NameSub//': iModel 6 no longer available.')
       case(7) ! Legacy conductance from RCM; DOES NOT FUNCTION.
          call CON_stop(NameSub//': iModel 7 no longer available.')
       case(8)
       case default
          write(*,*) NameSub//': iModel = ', iModelIn
          call CON_stop(NameSub//': Invalid imodel value')
    end select

  end subroutine imodel_legacy
  !===========================================================================


end Module ModConductance
