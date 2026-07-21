# Standard Library Imports

# Third-Party Imports
import xarray as xr
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, NP_REAL, NP_BOOL, NP_ARRAY, XR_DATASET, XR_DATAARRAY

def vertical_remap_field(homme_output: str,
    l_time_slice: slice,
    field_key: str,
    p_tgt: NP_ARRAY[NP_REAL]) -> XR_DATAARRAY:

    supported_field_keys: list[str] = ["T", "u", "v", "w"]
    assert(field_key in supported_field_keys)

    #---------------------------------------------------------------------------
    # Extract relevant fields from HOMME file
    #---------------------------------------------------------------------------
    xr_homme: XR_DATASET
    with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as xr_homme:
        l_p_src: XR_DATAARRAY = (xr_homme["p"]
            .isel(time = l_time_slice)
            .load())
        l_field_src: XR_DATAARRAY = (xr_homme[field_key]
            .isel(time = l_time_slice)
            .load())

    #---------------------------------------------------------------------------
    # Vertically remap field to target pressure levels
    #---------------------------------------------------------------------------
    xr_p_tgt: XR_DATAARRAY = xr.DataArray(
        data = p_tgt,
        coords = {"p" : p_tgt},
        dims = ("p",),
        attrs = l_p_src.attrs
    )

    l_field_tgt: XR_DATAARRAY = xr.apply_ufunc(
        _vertical_remap_field_column,
        l_p_src,
        l_field_src,
        xr_p_tgt,
        input_core_dims = [
            ["lev"],
            ["lev"],
            ["p"]
        ],
        output_core_dims = [
            ["p"]
        ],
        exclude_dims = {
            "p"
        },
        vectorize = True,
        dask = "forbidden",
        output_dtypes = [NP_REAL]
    )

    l_field_tgt = (l_field_tgt
        .assign_coords(coords = {"p" : xr_p_tgt})
        .assign_attrs(l_field_src.attrs)
        .rename(field_key)
    )

    return l_field_tgt

def _vertical_remap_field_column(p_column_src: NP_ARRAY[NP_REAL],
    field_column_src: NP_ARRAY[NP_REAL],
    p_tgt: NP_ARRAY[NP_REAL]) -> NP_ARRAY[NP_REAL]:

    # np.interp requires increasing x-coordinates.
    sort_idx: NP_ARRAY[np.int_] = np.argsort(p_column_src)

    p_column_src_sort: NP_ARRAY[NP_REAL] = p_column_src[sort_idx]
    field_column_src_sort: NP_ARRAY[NP_REAL] = field_column_src[sort_idx]

    field_column_tgt: NP_ARRAY[NP_REAL] = NP_REAL(np.interp(
        p_tgt,
        p_column_src,
        field_column_src,
        left = np.nan,
        right = np.nan
    ))

    return field_column_tgt