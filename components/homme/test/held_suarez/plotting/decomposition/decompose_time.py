# Standard Library Imports

# Third-Party Imports
import xarray as xr
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, NP_ARRAY, XR_DATASET, XR_DATAARRAY

def decompose_time(homme_output: str,
    spinup_days: NP_INT,
    l_rank: NP_INT,
    comm_size: NP_INT) -> slice:

    #---------------------------------------------------------------------------
    # Extract relevant fields from HOMME file
    #---------------------------------------------------------------------------
    xr_homme: XR_DATASET
    with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as xr_homme:
        time: XR_DATAARRAY = (xr_homme["time"]
            .isel(time = slice(spinup_days, None))
            .load())

    #---------------------------------------------------------------------------
    # Calculate indexers for isel for time dimension
    #---------------------------------------------------------------------------
    g_n_time: NP_INT = NP_INT(time.size)
    comm_n_times: NP_ARRAY[NP_INT] = (
        (g_n_time // comm_size) + NP_INT(np.arange(comm_size) < (g_n_time - comm_size * (g_n_time // comm_size)))
    )

    l_n_time: NP_INT = comm_n_times[l_rank]
    l_time_st_idx: NP_INT = spinup_days + np.sum(comm_n_times[:l_rank])
    l_time_end_idx: NP_INT = l_time_st_idx + l_n_time
    
    return slice(l_time_st_idx, l_time_end_idx)