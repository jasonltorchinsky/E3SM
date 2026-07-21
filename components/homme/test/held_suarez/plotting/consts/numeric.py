# Standard Library Imports

# Third-Party Library Imports
import numpy as np

# Local Library Imports
from .dtypes import NP_INT, NP_REAL

# Numeric constants
NP_EPS: NP_REAL = np.finfo(NP_REAL).resolution
NP_INF: NP_REAL = np.finfo(NP_REAL).max

NP_SMALL: NP_REAL = np.sqrt(NP_EPS)
NP_LARGE: NP_REAL = np.sqrt(NP_INF)

# MPI constants
MPI_ROOT: NP_INT = NP_INT(0)