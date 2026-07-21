# Standard Library Imports

# Third-Party Library Imports
import numpy as np

# Local Library Imports
from .dtypes import NP_REAL

# Physical constants
mu_d: NP_REAL = NP_REAL(28.9467e-3) # Mean molar mass of dry air - https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html [kg mol^(-1)]
mu_v: NP_REAL = NP_REAL(18.0153e-3) # Molar mass of water - https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185 [kg mol^(-1)]

R: NP_REAL = NP_REAL(8.314462619) # Molar gas constant - https://physics.nist.gov/cgi-bin/cuu/Value?r [J mol^(-1) K^(-1)]
R_d: NP_REAL = R / mu_d # Gas constant for dry air [J kg^(-1) K^(-1)]
R_v: NP_REAL = R / mu_v # Gas constant for water vapor [J kg^(-1) K^(-1)]

cp_d: NP_REAL = NP_REAL(1.0061e3) # Specific heat of dry air at constant pressure [J kg^{-1} K^{-1}]
cp_v: NP_REAL = NP_REAL(1.884e3) # Specific heat of water vapor at constant pressure [J kg^{-1} K^{-1}]
cp_lw: NP_REAL = NP_REAL(4.184e3) # Specific heat of liquid water at constant pressure [J kg^{-1} K^{-1}]
cp_iw: NP_REAL = NP_REAL(2.093e3) # Specific heat of ice water at constant pressure [J kg^{-1} K^{-1}]

L_v: NP_REAL = NP_REAL(2.5009e6) # Latent heat of vaporization of water at 0.01C - https://www.engineeringtoolbox.com/water-properties-d_1573.html [J kg^(-1)]

# Reference values
g: NP_REAL = NP_REAL(9.80665) # Standard acceleration of gravity [m s^{-2}] - https://physics.nist.gov/cgi-bin/cuu/Value?gn
p_0: NP_REAL = NP_REAL(101325.0) # Standard atmospheric pressure [Pa] - https://physics.nist.gov/cgi-bin/cuu/Value?stdatm

rho_sw: NP_REAL = NP_REAL(1.027e3) # Reference density of sea water [kg^{-1} m^{-3}]
cp_sw: NP_REAL = NP_REAL(3986.) # Reference specific heat of seawater at constant pressure [J kg^{-1} K^{-1}]
h_m: NP_REAL = NP_REAL(19.753) # Approximate mixing layer depth of the GATEIII region in August [m]

# Conversion factors
sec_per_hour: NP_REAL = NP_REAL(3600.)
sec_per_day: NP_REAL = 24. * sec_per_hour # Seconds per day [s d^{-1}]

# Constants from Polvani & Kushner (2002)
pk02_tropopause_pressure: NP_REAL = NP_REAL(1.e2) # Polvani & Kushner (2002) nominal tropopause pressure [hPa]
pk02_sponge_pressure: NP_REAL = NP_REAL(5.e-1) # Polvani & Kushner (2002) sponge layer height [hPa]