#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module us_standard_atmosphere_1976_mod
  use kinds,                  only: real_kind
  use physical_constants,     only: g,rearth0
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use common_io_mod,          only: infilenames
#endif

implicit none
private
  
  ! U.S. Standard Atmosphere 1976 parameters
  ! Derived from (NOAA 1976), doi: unknown
  real (kind=real_kind), private, parameter :: M_0 = 2.89644D1 ! Mean molecular weight of dry air at sea-level [kg kmol^{-1}]
  real (kind=real_kind), private, parameter :: R_star = 8.31432D3 ! Universal gas constant [N m K^{-1} kmol^{-1}]
  real (kind=real_kind), private, parameter :: H_b(0:7) &
    = (/ 0.0D3, 1.1D4, 2.0D4, 3.2D4, 4.7D4, 5.1D4, 7.1D4, 8.48520D4 /) ! Geopotential height at reference levels [m']
  real (kind=real_kind), private, parameter :: LM_b(0:6) &
    = (/ -6.5D-3, 0.0D-3, 1.0D-3, 2.8D-3, 0.0D-3, -2.8D-3, -2.0D-3 /) ! Molecular-scale temperature gradient between reference levels [K m'^{-1}]
  real (kind=real_kind), private, parameter :: TM_b(0:7) &
    = (/ 288.15D0, 216.65D0, 216.65D0, 228.65D0, 270.65D0, 270.65D0, 214.65D0, 186.946D0 /) ! Molecular-scale temperature at reference levels [K]
  real (kind=real_kind), private, parameter :: p_b(0:7) &
    = (/ 1.01325D5, 2.2632064D4, 5.47488867D3, 8.68018685D2, 1.10906306D2, &
       6.69388731D1, 3.95642043D0, 3.73383590D-1 /) ! Pressure at reference levels [Pa]
   
  public :: geometric_height_from_pressure
  public :: temperature_from_geometric_height
  public :: temperature_from_pressure

contains

  function geometric_height_from_pressure(p) result(z)
    real (kind=real_kind), intent(in)  :: p
    real (kind=real_kind) :: z

    ! local variables
    real (kind=real_kind) :: h

    if ((p <= p_b(0)) .and. (p > p_b(1))) then
         h = H_b(0) + (TM_b(0) / LM_b(0)) * ((p_b(0) / p)**((R_star * LM_b(0)) / (g * M_0)) - 1.0D0)
    else if ((p <= p_b(1)) .and. (p > p_b(2))) then
         h = H_b(1) - ((R_star * TM_b(1)) / (g * M_0)) * LOG(p / p_b(1))
    else if ((p <= p_b(2)) .and. (p > p_b(3))) then
         h = H_b(2) + (TM_b(2) / LM_b(2)) * ((p_b(2) / p)**((R_star * LM_b(2)) / (g * M_0)) - 1.0D0)
    else if ((p <= p_b(3)) .and. (p > p_b(4))) then
         h = H_b(3) + (TM_b(3) / LM_b(3)) * ((p_b(3) / p)**((R_star * LM_b(3)) / (g * M_0)) - 1.0D0)
    else if ((p <= p_b(4)) .and. (p > p_b(5))) then
         h = H_b(4) - ((R_star * TM_b(4)) / (g * M_0)) * LOG(p / p_b(4))
    else if ((p <= p_b(5)) .and. (p > p_b(6))) then
         h = H_b(5) + (TM_b(5) / LM_b(5)) * ((p_b(5) / p)**((R_star * LM_b(5)) / (g * M_0)) - 1.0D0)
    else if ((p <= p_b(6)) .and. (p >= p_b(7))) then
         h = H_b(6) + (TM_b(6) / LM_b(6)) * ((p_b(6) / p)**((R_star * LM_b(6)) / (g * M_0)) - 1.0D0)
    end if

    z = (rearth0 * h) / (rearth0 - h) ! Geometric height [m]

  end function geometric_height_from_pressure

  function temperature_from_geometric_height(z) result (T)
    real (kind=real_kind), intent(in)  :: z
    real (kind=real_kind) :: T

    ! local variables
    real (kind=real_kind) :: h

    h = (rearth0 * z) / (rearth0 + z) ! Geopotential height [m']
    if ((h >= H_b(0)) .and. (h < H_b(1))) then
         T = TM_b(0) + LM_b(0) * (h - H_b(0))
    else if ((h >= H_b(1)) .and. (h < H_b(2))) then
         T = TM_b(1) + LM_b(1) * (h - H_b(1))
    else if ((h >= H_b(2)) .and. (h < H_b(3))) then
         T = TM_b(2) + LM_b(2) * (h - H_b(2))
    else if ((h >= H_b(3)) .and. (h < H_b(4))) then
         T = TM_b(3) + LM_b(3) * (h - H_b(3))
    else if ((h >= H_b(4)) .and. (h < H_b(5))) then
         T = TM_b(4) + LM_b(4) * (h - H_b(4))
    else if ((h >= H_b(5)) .and. (h < H_b(6))) then
         T = TM_b(5) + LM_b(5) * (h - H_b(5))
    else if ((h >= H_b(6)) .and. (h <= H_b(7))) then
         T = TM_b(6) + LM_b(6) * (h - H_b(6))
    end if

  end function temperature_from_geometric_height

  function temperature_from_pressure(p) result (T)
    real (kind=real_kind), intent(in)  :: p
    real (kind=real_kind) :: T

    T = temperature_from_geometric_height(geometric_height_from_pressure(p))

  end function temperature_from_pressure

end module us_standard_atmosphere_1976_mod