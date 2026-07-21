# Standard Library Imports

# Third-Party Library Imports
import numpy as np

# Local Library Imports

# Numeric precision
INT_PREC: int = 64
REAL_PREC: int = 64

assert(INT_PREC in [32, 64])
if INT_PREC == 32:
    from mpi4py.MPI import INT as MPI_INT
    from numpy import int32 as NP_INT
    NC_INT: str = "i4"
elif INT_PREC == 64:
    from mpi4py.MPI import LONG as MPI_INT
    from numpy import int64 as NP_INT
    NC_INT: str = "i8"

assert(REAL_PREC in [32, 64])
if REAL_PREC == 32:
    from mpi4py.MPI import FLOAT as MPI_REAL
    from numpy import float32 as NP_REAL
    NC_REAL: str = "f4"
elif REAL_PREC == 64:
    from mpi4py.MPI import DOUBLE as MPI_REAL
    from numpy import float64 as NP_REAL
    NC_REAL: str = "f8"

# Class aliases
from matplotlib.axes import Axes as MPL_AXES
from matplotlib.colorbar import Colorbar as MPL_COLORBAR
from matplotlib.figure import Figure as MPL_FIGURE
from matplotlib.contour import QuadContourSet as MPL_CONTOUR
from matplotlib.pyplot import pcolormesh as MPL_PCOLORMESH
from mpi4py.MPI import Intracomm as MPI_COMM
from netCDF4 import Dataset as NC_DATASET
from netCDF4._netCDF4 import Dimension as NC_DIMENSION
from netCDF4._netCDF4 import Variable as NC_VARIABLE
from numpy import ndarray as NP_ARRAY
NP_DATETIME = np.dtype('datetime64[ns]')
from xarray.core.dataset import Dataset as XR_DATASET
from xarray.core.dataarray import DataArray as XR_DATAARRAY
NP_BOOL = np.bool_