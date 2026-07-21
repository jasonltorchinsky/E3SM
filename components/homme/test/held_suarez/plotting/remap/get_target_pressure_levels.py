# Standard Library Imports

# Third-Party Imports
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, NP_REAL, NP_ARRAY

def get_target_pressure_levels() -> NP_ARRAY[NP_REAL]:
    # TO-DO: Base this off of something from the HOMME output

    #---------------------------------------------------------------------------
    # Calculate target pressure levels
    #---------------------------------------------------------------------------
    p_min: NP_REAL = NP_REAL(0.1e2) # Minimum pressure level; [Pa]
    p_max: NP_REAL = NP_REAL(1000.e2) # Maximum pressure levels; [Pa]
    n_lev: NP_INT = NP_INT(256) # Number of target levels
    
    p_lev: NP_ARRAY[NP_REAL] = np.logspace(
        np.log10(p_min),
        np.log10(p_max),
        num = n_lev) # Equally log-spaced pressure levels (min -> max); [Pa]

    return  p_lev