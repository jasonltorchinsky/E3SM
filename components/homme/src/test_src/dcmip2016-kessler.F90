!-----------------------------------------------------------------------
!
!  Version:  2.0
!
!  Date:  January 22nd, 2015
!
!  Change log:
!  v2 - Added sub-cycling of rain sedimentation so as not to violate
!       CFL condition.
!
!  The KESSLER subroutine implements the Kessler (1969) microphysics
!  parameterization as described by Soong and Ogura (1973) and Klemp
!  and Wilhelmson (1978, KW). KESSLER is called at the end of each
!  time step and makes the final adjustments to the potential
!  temperature and moisture variables due to microphysical processes
!  occurring during that time step. KESSLER is called once for each
!  vertical column of grid cells. Increments are computed and added
!  into the respective variables. The Kessler scheme contains three
!  moisture categories: water vapor, cloud water (liquid water that
!  moves with the flow), and rain water (liquid water that falls
!  relative to the surrounding air). There  are no ice categories.
!  Variables in the column are ordered from the surface to the top.
!
!  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
!
!  Input variables:
!     theta  - potential temperature (K)
!     qv     - water vapor mixing ratio (gm/gm)
!     qc     - cloud water mixing ratio (gm/gm)
!     qr     - rain  water mixing ratio (gm/gm)
!     rho    - dry air density (not mean state as in KW) (kg/m^3)
!     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
!     dt     - time step (s)
!     z      - heights of thermodynamic levels in the grid column (m)
!     nz     - number of thermodynamic levels in the column
!     precl  - Precipitation rate (m_water/s)
!
! Output variables:
!     Increments are added into t, qv, qc, qr, and rainnc which are
!     returned to the routine from which KESSLER was called. To obtain
!     the total precip qt, after calling the KESSLER routine, compute:
!
!       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
!       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
!
!
!  Authors: Paul Ullrich
!           University of California, Davis
!           Email: paullrich@ucdavis.edu
!
!           Based on a code by Joseph Klemp
!           (National Center for Atmospheric Research)
!
!  Reference:
!
!    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
!    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
!    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
!    doi:10.1002/2015MS000435
!
!=======================================================================

SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, precl)

use physical_constants,   only:  kappa

  IMPLICIT NONE

  !------------------------------------------------
  !   Input / output parameters
  !------------------------------------------------

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta   ,     & ! Potential temperature (K)
            qv      ,     & ! Water vapor mixing ratio (gm/gm)
            qc      ,     & ! Cloud water mixing ratio (gm/gm)
            qr              ! Rain  water mixing ratio (gm/gm)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho             ! Dry air density (not mean state as in KW) (kg/m^3)

  REAL(8), INTENT(OUT) :: &
            precl          ! Precipitation rate (m_water / s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z       ,     & ! Heights of thermo. levels in the grid column (m)
            pk              ! Exner function (p/p0)**(R/cp)

  REAL(8), INTENT(IN) :: & 
            dt              ! Time step (s)

  INTEGER, INTENT(IN) :: nz ! Number of thermodynamic levels in the column

  !------------------------------------------------
  !   Local variables
  !------------------------------------------------
  REAL, DIMENSION(nz) :: r, rhalf, velqr, sed, pc

  REAL(8) :: f5, f2x, xk, ern, qrprod, prod, qvs, psl, rhoqr, dt_max, dt0

  INTEGER :: k, rainsplit, nt

  !------------------------------------------------
  !   Begin calculation
  !------------------------------------------------
  f2x = 17.27d0
  f5 = 237.3d0 * f2x * 2500000.d0 / 1003.d0
  xk = .2875d0      !  kappa (r/cp)
  !xk = kappa      !  kappa (r/cp)

  psl    = 1000.d0  !  pressure at sea level (mb)
  rhoqr  = 1000.d0  !  density of liquid water (kg/m^3)

  do k=1,nz
    r(k)     = 0.001d0*rho(k)
    rhalf(k) = sqrt(rho(1)/rho(k))
    pc(k)    = 3.8d0/(pk(k)**(1./xk)*psl)

    ! Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)

  end do

  ! Maximum time step size in accordance with CFL condition
  if (dt .le. 0.d0) then
    write(*,*) 'kessler.f90 called with nonpositive dt'
    stop
  end if

  dt_max = dt
  do k=1,nz-1
    if (velqr(k) .ne. 0.d0) then
      dt_max = min(dt_max, 0.8d0*(z(k+1)-z(k))/velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit,8)

  ! Subcycle through rain process
  precl = 0.d0

  do nt=1,rainsplit

    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr

    ! Sedimentation term using upstream differencing
    do k=1,nz-1
      sed(k) = dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end do
    sed(nz)  = -dt0*qr(nz)*velqr(nz)/(.5*(z(nz)-z(nz-1)))

    ! Adjustment terms
    do k=1,nz

      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k)-dt0*max(.001*(qc(k)-.001d0),0.d0))/(1.d0+dt0*2.2d0*qr(k)**.875)
      qc(k) = max(qc(k)-qrprod,0.d0)
      qr(k) = max(qr(k)+qrprod+sed(k),0.d0)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      qvs = pc(k)*exp(f2x*(pk(k)*theta(k)-273.d0)   &
             /(pk(k)*theta(k)- 36.d0))
      prod = (qv(k)-qvs)/(1.d0+qvs*f5/(pk(k)*theta(k)-36.d0)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      ern = min(dt0*(((1.6d0+124.9d0*(r(k)*qr(k))**.2046)  &
            *(r(k)*qr(k))**.525)/(2550000d0*pc(k)            &
            /(3.8d0 *qvs)+540000d0))*(dim(qvs,qv(k))         &
            /(r(k)*qvs)),max(-prod-qc(k),0.d0),qr(k))

      ! Saturation adjustment following KW eq. 3.10
      theta(k)= theta(k) + 2500000d0/(1003.d0*pk(k))*(max( prod,-qc(k))-ern)
      qv(k) = max(qv(k)-max(prod,-qc(k))+ern,0.d0)
      qc(k) = qc(k)+max(prod,-qc(k))
      qr(k) = qr(k)-ern
    end do

    ! Recalculate liquid water terminal velocity
    if (nt .ne. rainsplit) then
      do k=1,nz
        velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)
      end do
    end if
  end do

  precl = precl / dble(rainsplit)

END SUBROUTINE KESSLER

!=======================================================================

!=======================================================================
! Kessler microphysics at constant volume
!=======================================================================

SUBROUTINE KESSLER_Z(theta, qv, qc, qr, rho, pk, dt, z, nz, precl)

use physical_constants,   only:  kappa

  IMPLICIT NONE

  !------------------------------------------------
  !   Input / output parameters
  !------------------------------------------------

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta   ,     & ! Potential temperature (K)
            qv      ,     & ! Water vapor mixing ratio (gm/gm)
            qc      ,     & ! Cloud water mixing ratio (gm/gm)
            qr      ,     & ! Rain  water mixing ratio (gm/gm)
            pk              ! Exner function (p/p0)**(R/cp)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho             ! Dry air density (not mean state as in KW) (kg/m^3)

  REAL(8), INTENT(OUT) :: &
            precl          ! Precipitation rate (m_water / s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z              ! Heights of thermo. levels in the grid column (m)

  REAL(8), INTENT(IN) :: & 
            dt              ! Time step (s)

  INTEGER, INTENT(IN) :: nz ! Number of thermodynamic levels in the column

  !------------------------------------------------
  !   Local variables
  !------------------------------------------------
  REAL, DIMENSION(nz) :: &
            r,           & ! Dry air density (not mean state as in KW) (g/m^3)
            rhalf,       & ! Density ratio coeff in KW ep. 2.15
            velqr,       & ! Rain terminal velocity (m/s)
            sed,         & ! Sedimentation term of M_qr (d/dz[rhobar V qr])
            pc,          & ! Pressure (mb)
            T              ! Temperature (C)

  REAL(8), PARAMETER ::      &
            cp  = 1003.0d0,  & ! Specific heat of dry air at constant pressure (J/(kg K))
            cv  = 716.0d0,   & ! Specific heat of dry air at constant volume (J/(kg K))
            lv  = 2.5D6,     & ! Latent heat of vaporization (J / kg)
            lf  = 3.34D5,    & ! Latent heat of freezing/meting (J / kg)
            xk  = .2875d0,   & ! Gas constant (R / c_p)
            psl = 1000.0d0,  & ! Pressure at sea-level (mb)
            rhoqr = 1000.0d0, &  ! Density of liquid water (kg/m^3)
            k1  = 0.001d0,   & ! k_1 in KW eq. 2.13a (1/s)
            k2  = 2.2d0,     & ! k_2 in KW eq. 2.13b (1/s)
            a   = 0.001d0    ! a in KW eq. 2.13a

  REAL(8) ::             &
            evap,        & ! Amount of rainwater evaporated
            qrprod,      & ! Rain production rate (Ar + Cr)
            prod,        &
            qvs,         &
            C,           & !Ventilation factor
            dt_max,      &
            dt0

  INTEGER :: k, rainsplit, nt

  !------------------------------------------------
  !   Begin calculation
  !------------------------------------------------
  
  ! Calculate pressure, terminal velocity
  do k = 1, nz
    r(k)     = 0.001d0 * rho(k)
    rhalf(k) = sqrt(rho(1) / rho(k))
    pc(k)    = pk(k)**(1./xk)*psl

    ! Rain terminal velocity via KW eq. 2.15 (m/s)
    velqr(k)  = 36.34d0*(qr(k) * r(k))**0.1364 * rhalf(k)
  end do

  ! Maximum time step size in accordance with CFL condition
  if (dt .le. 0.d0) then
    write(*,*) 'kessler.f90 called with nonpositive dt'
    stop
  end if

  dt_max = dt
  do k=1,nz-1
    if (velqr(k) .ne. 0.d0) then
      dt_max = min(dt_max, 0.8d0*(z(k+1)-z(k))/velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit,8)

  ! Subcycle through rain process
  precl = 0.d0

  do nt = 1, rainsplit

    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr

    ! Sedimentation term using upstream differencing (d/dz [rho * V * qr])
    do k = 1, nz-1
      sed(k) = dt0 * (r(k+1) * qr(k+1) * velqr(k+1) - r(k) * qr(k) * velqr(k)) &
                   & /(r(k) * (z(k+1) - z(k)))
    end do
   sed(nz)  = -dt0 * qr(nz) * velqr(nz) / (.5 * (z(nz) - z(nz-1)))

    ! Adjustment terms
    do k = 1, nz

      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k) - dt0 * max(k1 * (qc(k) - a), 0.d0)) &
                            & / (1.d0 + dt0 * k2 * qr(k)**.875)
      qc(k) = max(qc(k) - qrprod, 0.d0)
      qr(k) = max(qr(k) + qrprod + sed(k), 0.d0)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      T(k) = pk(k) * theta(k)
      qvs = (3.8d0 / pc(k)) * exp(17.27d0 * (T(k) - 273.d0) / (T(k)- 36.d0))
      prod = (qv(k) - qvs) / (1.d0 + qvs * (4093.d0 * lv / cp) / (T(k) - 36.d0)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      C = 1.6d0 + 124.9d0 * (r(k) * qr(k))**(0.2046d0)
      evap = min(dt0 * ((C * (r(k) * qr(k))**.525) / (5.4D5 +  2.55D6 / (pc(k) * qvs))) &
                & * (dim(qvs, qv(k)) / (r(k)*qvs)), max(-prod - qc(k), 0.d0), qr(k))

      ! Saturation adjustment following KW eq. 3.10
      theta(k) = theta(k) + lv/(cv * pk(k)) * (max(prod, -qc(k)) - evap)
      pk(k) = pk(k) - (pk(k) / (cv * T(k))) &
                          * (lv * ((1.0D0 / cv) - 1.0D0) * (max(prod, -qc(k)) - evap) &
                                 + lf * sed(k))
      qv(k) = max(qv(k) - max(prod, -qc(k)) + evap, 0.d0)
      qc(k) = qc(k) + max(prod, -qc(k))
      qr(k) = qr(k) - evap
    end do

    ! Recalculate liquid water terminal velocity
    if (nt .ne. rainsplit) then
      do k=1,nz
        velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)
      end do
    end if
  end do

  precl = precl / dble(rainsplit)

END SUBROUTINE KESSLER_Z

!=======================================================================

!=======================================================================
! Microphysics by J.L. Torchinsky, Summer 2022
!=======================================================================

SUBROUTINE TORCHINSKY(theta, qv, qc, qr, rho, pk, dt, z, nz, precl)

use physical_constants,   only:  kappa

  IMPLICIT NONE

  !------------------------------------------------
  !   Input / output parameters
  !------------------------------------------------

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta   ,     & ! Potential temperature (K)
            qv      ,     & ! Water vapor mixing ratio (gm/gm)
            qc      ,     & ! Cloud water mixing ratio (gm/gm)
            qr              ! Rain  water mixing ratio (gm/gm)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho             ! Dry air density (not mean state as in KW) (kg/m^3)

  REAL(8), INTENT(OUT) :: &
            precl          ! Precipitation rate (m_water / s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z       ,     & ! Heights of thermo. levels in the grid column (m)
            pk              ! Exner function (p/p0)**(R/cp)

  REAL(8), INTENT(IN) :: & 
            dt              ! Time step (s)

  INTEGER, INTENT(IN) :: nz ! Number of thermodynamic levels in the column

  !------------------------------------------------
  !   Local variables
  !------------------------------------------------
  REAL, DIMENSION(nz) :: r, rhalf, velqr, sed, pc

  REAL(8) :: f5, f2x, xk, ern, qrprod, prod, qvs, psl, rhoqr, dt_max, dt0

  INTEGER :: k, rainsplit, nt

  !------------------------------------------------
  !   Begin calculation
  !------------------------------------------------
  f2x = 17.27d0
  f5 = 237.3d0 * f2x * 2500000.d0 / 1003.d0
  xk = .2875d0      !  kappa (r/cp)
  !xk = kappa      !  kappa (r/cp)

  psl    = 1000.d0  !  pressure at sea level (mb)
  rhoqr  = 1000.d0  !  density of liquid water (kg/m^3)

  do k=1,nz
    r(k)     = 0.001d0*rho(k)
    rhalf(k) = sqrt(rho(1)/rho(k))
    pc(k)    = 3.8d0/(pk(k)**(1./xk)*psl)

    ! Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)

  end do

  ! Maximum time step size in accordance with CFL condition
  if (dt .le. 0.d0) then
    write(*,*) 'kessler.f90 called with nonpositive dt'
    stop
  end if

  dt_max = dt
  do k=1,nz-1
    if (velqr(k) .ne. 0.d0) then
      dt_max = min(dt_max, 0.8d0*(z(k+1)-z(k))/velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit,8)

  ! Subcycle through rain process
  precl = 0.d0

  do nt=1,rainsplit

    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr

    ! Sedimentation term using upstream differencing
    do k=1,nz-1
      sed(k) = dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end do
    sed(nz)  = -dt0*qr(nz)*velqr(nz)/(.5*(z(nz)-z(nz-1)))

    ! Adjustment terms
    do k=1,nz

      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k)-dt0*max(.001*(qc(k)-.001d0),0.d0))/(1.d0+dt0*2.2d0*qr(k)**.875)
      qc(k) = max(qc(k)-qrprod,0.d0)
      qr(k) = max(qr(k)+qrprod+sed(k),0.d0)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      qvs = pc(k)*exp(f2x*(pk(k)*theta(k)-273.d0)   &
             /(pk(k)*theta(k)- 36.d0))
      prod = (qv(k)-qvs)/(1.d0+qvs*f5/(pk(k)*theta(k)-36.d0)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      ern = min(dt0*(((1.6d0+124.9d0*(r(k)*qr(k))**.2046)  &
            *(r(k)*qr(k))**.525)/(2550000d0*pc(k)            &
            /(3.8d0 *qvs)+540000d0))*(dim(qvs,qv(k))         &
            /(r(k)*qvs)),max(-prod-qc(k),0.d0),qr(k))

      ! Saturation adjustment following KW eq. 3.10
      theta(k)= theta(k) + 2500000d0/(1003.d0*pk(k))*(max( prod,-qc(k))-ern)
      qv(k) = max(qv(k)-max(prod,-qc(k))+ern,0.d0)
      qc(k) = qc(k)+max(prod,-qc(k))
      qr(k) = qr(k)-ern
    end do

    ! Recalculate liquid water terminal velocity
    if (nt .ne. rainsplit) then
      do k=1,nz
        velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)
      end do
    end if
  end do

  precl = precl / dble(rainsplit)

END SUBROUTINE TORCHINSKY
