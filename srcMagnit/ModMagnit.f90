! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

module ModMagnit

  use ModConst, ONLY: cBoltzmann, cProtonMass, cElectronMass, cKToKEV, cKEVToK
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi
  use ModConductance, ONLY: DoUseGmPe

  implicit none
  save

  ! Unit conversion factors and other constants:
  real, parameter :: cPtoKProt = cProtonMass/cBoltzmann * cKToKEV
  real, parameter :: cPtoKElec = cElectronMass/cBoltzmann * cKToKEV

  ! Diffuse auroral parameters
  real :: ratioPe = 1./6. ! Ratio of electron P to proton P.

  MinWidth = 5.0 * cPi / 180.0
  nHalfSmooth = 5
!     MulFac_Dae = 1.0e22
!     MulFac_Def = 5.0e19
!     MulFac_ef = 0.2e7
!     MulFac_ae = 1.0 / 1.0e11
  MulFac_Dae = 1.0e22
  MulFac_Def = 1.0e19
  MulFac_ef = 0.3e6! * 1e03
  MulFac_ae = 3.65e-12

  ! reminder: GmRhoFloor = 1E-21, GmPFloor = 1E-13

  contains
  !============================================================================
  subroutine print_magnit_config(iUnitIn)

    integer, intent(in) :: iUnitIn
    !--------------------------------------------------------------------------

    write(iUnitin,'(a)') "MAGNIT Physics-based Aurora"
    write(iUnitin,'(a)') "(Beta-testing Phase)"
    write(iUnitIn,*) '################## CAUTION ##################'
    write(iUnitIn,*) 'MAGNIT is in BETA TESTING and may not be stable.'
    write(iUnitIn,*) 'Proceed with caution. For issues, contact developers.'
    write(iUnitIn,*) '#############################################'
    if(DoUseGmPe)then
      write(iUnitIn,*) "Electron precip obtained directly from Pe."
    else
      write(iUnitIn,*) "Electron precip obtained from LT-flipped protons"
      write(iUnitIn,'(a, 5.3f)') "Pe/P ratio set to ", ratioPe
    end if

  end subroutine print_magnit_config

  !============================================================================
  subroutine magnit_gen_fluxes(NameHemiIn)
    ! Given magnetospheric density, pressure, and FACs, calculate diffuse and
    ! monoenergetic precipitating fluxes.
    ! Fill the following arrays for use by ModConductance:
    ! SigmaPedAur_II, SigmaHalAur_II, SigmaPedBbnd_II, SigmaHalBbnd_II

    character(len=*), intent(in) :: NameHemiIn

    ! Set arrays to hold magnetospheric values.
    real, dimension(IONO_nTheta, IONO_nPsi) :: &
        MagP_II, MagRhop_II, MagPe_II, MagRhoe_II

    ! Set arrays to hold precip values. Magnetospheric pressure and density
    ! from GM (SI units), Average Energy and Energy Flux (units of KeV and
    ! mW/m2 aka ergs/cm2, respectively) for each precip type.
    real, dimension(IONO_nTheta, IONO_nPsi) :: &
        AveEDiffe_II, AveEDiffp_II, AveEMono_II, AveEBbnd_II, &
        EfluxDiffe_II, EfluxDiffp_II, EfluxMono_II, EfluxBbnd_II


    ! Debug variables:
    logical :: DoTestMe, DoTest
    character(len=*), parameter:: NameSub = 'magnit_gen_fluxes'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe) &
         write(*,*)'IE: '//NameSub//' called for hemisphere='//NameHemiIn

    ! Set values based on hemisphere.
    if(NameHemIn == 'north')then
        MagP_II = iono_north_p
        MagRho_II = iono_north_rho
    else if (NameHemiIn.eq.'south')then
        MagP_II = iono_south_p
        MagRho_II = iono_south_rho
    else
      call CON_stop(NameSub//' : unrecognized hemisphere - '//NameHemiIn)
    end if

    ! Set default/background values.

    ! If we need OCFLB calculation, we would put it here.

    ! If not using explicit Pe from GM, obtain values from proton values.
    ! flip across noon midnight, etc. to fill MagPe, MagRhoe
    ! Scale Pe using ratioPe
    if(.not. DoUseGmPe) then
      do j=1, Iono_nPsi
        MagPe_II(:, j) = ratioPe * MagP_II(:, IONO_nPsi-j+1)
        MagRhoe_II(:, j) = MagRhop_II(:, IONO_nPsi-j+1)
      end do
    end if

    ! Calculate diffuse precipitation: protons.
    AveEDiffp = MagP_II / MagRhop_II * cPtoKProt  ! T = P/nk in
    EfluxDiffp = AveEDiffp * (MagRhop_II/cProtonMass*100^3) * sqrt(kBoltz)

    ! Calculate diffuse precipitation: electrons.
    AveEDiffp = MagPe_II / MagRhoe_II * cPtoKElec
    EfluxDiffp = AveEDiffe * (MagRhop_II/cElectronMass*100^3) * sqrt(kBoltz)

    ! Fill appropriate arrays for IE/Ridley_serial calculations.
    ! SigmaPedAur_II, SigmaHalAur_II, EfluxMono_II, AvgEMono_II

  end subroutine magnit_gen_fluxes
  !============================================================================
