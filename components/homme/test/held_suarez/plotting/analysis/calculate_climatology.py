# Standard Library Imports

# Third-Party Imports
import xarray as xr
from mpi4py import MPI
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, NP_REAL, NP_ARRAY, XR_DATASET, XR_DATAARRAY, MPI_COMM
from consts.numeric import MPI_ROOT

def calculate_climatology(homme_output: str,
    spinup_days: NP_INT,
    l_field_tgt: XR_DATAARRAY,
    comm: MPI_COMM) -> XR_DATAARRAY:

    #---------------------------------------------------------------------------
    # Retrieve MPI communicator information
    #---------------------------------------------------------------------------
    l_rank: NP_INT = NP_INT(comm.Get_rank())

    #---------------------------------------------------------------------------
    # Extract relevant fields from HOMME file
    #---------------------------------------------------------------------------
    xr_homme: XR_DATASET
    with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as xr_homme:
        time: XR_DATAARRAY = (xr_homme["time"]
            .isel(time = slice(spinup_days, None))
            .load())
        lon: XR_DATAARRAY = (xr_homme["lon"]
            .load())

    #---------------------------------------------------------------------------
    # Sum local field over time, longitude, divide by n_time, n_lon to get
    # contirbution to global climatology
    #---------------------------------------------------------------------------
    g_n_time: NP_INT = NP_INT(time.size)
    g_n_lon: NP_INT = NP_INT(lon.size)

    l_field_clim: XR_DATAARARAY = (
        (1. / (g_n_time * g_n_lon)) 
        * l_field_tgt.sum(dim = ["time", "lon"], skipna = True, keep_attrs = True)
    )

    l_np_field_clim: NP_ARRAY[NP_REAL] = NP_REAL(l_field_clim.to_numpy())
    g_np_field_clim: NP_ARRAY[NP_REAL] = np.empty_like(l_np_field_clim, dtype = NP_REAL)
    comm.Reduce(l_np_field_clim, g_np_field_clim, op = MPI.SUM, root = MPI_ROOT)

    g_field_clim: XR_DATAARRAY = xr.DataArray(
        data = g_np_field_clim,
        coords = l_field_clim.coords,
        dims = l_field_clim.dims,
        name = l_field_tgt.name,
        attrs = l_field_tgt.attrs
    )

    return g_field_clim