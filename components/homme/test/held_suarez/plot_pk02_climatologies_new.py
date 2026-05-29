# Library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime
import xarray as xr
from mpi4py import MPI

pk02_p_T = 100. # Polvani-Kushner 2002 nominal tropopause height [hPa]
pk02_p_sp = 0.5 # Polvani-Kushner 2002 sponge layer height [hPa]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--spinup-days", nargs = "?", default = 0, type = int,
        help = "Spin-up days to skip in calculations.")
    parser.add_argument("--homme-output", nargs = "?", required = True, type = str,
        help = "HOMME output file.")
    parser.add_argument("--working-dir", nargs = "?", default = ".working", type = str,
        help = "Working directory to output calculated values.")
    parser.add_argument("--plotting-dir", nargs = "?", default = ".plotting", type = str,
        help = "Directory to save plots.")
    parser.add_argument("--recalculate", nargs = "?", default = False, type = str2bool,
        help = "Re-calculate climatologies.")
    parser.add_argument("--tag", nargs = "?", default = "", type = str,
        help = "Dataset tag.")
    parser.add_argument("--plot-vars", nargs = "?", default = "", type = str,
        help = "Which variable climatologies to plot: u, T, T_eddy")
    args = parser.parse_args()

    spinup_days = args.spinup_days
    homme_output = args.homme_output
    working_dir = args.working_dir
    plotting_dir = args.plotting_dir
    recalculate = args.recalculate
    tag = args.tag
    plot_vars = [str(plot_var) for plot_var in args.plot_vars.split(",") if plot_var]

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        dirs = [working_dir, plotting_dir]
        for dir in dirs:
            if not os.path.exists(dir):
                os.makedirs(dir)

    comm.barrier()

    #---------------------------------------------------------------------------
    # Set up variables shared across plot variables
    #---------------------------------------------------------------------------
    p_tgt = None
    lat_vals = None
    lat_slice = None

    #---------------------------------------------------------------------------
    # Calculate climatologies and create plots
    #---------------------------------------------------------------------------
    if rank == 0:
        msg = "[{}]: Starting climatology calculation and plotting loop.".format(datetime.now().strftime("%H:%M:%S"))
        print(msg, flush = True)

    for plot_var in plot_vars:
        if rank == 0:
            msg = "[{}]: Starting {}.".format(datetime.now().strftime("%H:%M:%S"), plot_var)
            print(msg, flush = True)

        clim_fileroot = plot_var + "_clim"
        if tag:
            clim_fileroot += "_{}".format(tag)

        clim_filepath = os.path.join(working_dir, clim_fileroot + ".nc")

        if recalculate or (not os.path.exists(clim_filepath)):

            if p_tgt is None or lat_vals is None or lat_slice is None:
                p_tgt, lat_vals, lat_slice = setup_pressure_grid_and_decomp(
                    homme_output, spinup_days, comm)

            local_clim = calc_local_clim(
                plot_var, homme_output, spinup_days, p_tgt, lat_slice)

            gathered = comm.gather(local_clim, root = 0)

            if rank == 0:
                clim = xr.concat(gathered, dim = "lat")
                clim = clim.sortby("lat")

                if plot_var in ["pnh"]:
                    clim = clim / 100. # [Pa] => [hPa]
                    clim.attrs["units"] = "hPa"

                if plot_var in ["T_eddy"]:
                    clim = clim.rename("T_eddy")
                    clim.attrs = {'long_name': 'Temperature Eddy Variance at Midpoints', 'units': 'K^{2}'}

                msg = "[{}]: Saving {} to file.".format(datetime.now().strftime("%H:%M:%S"), plot_var)
                print(msg, flush = True)
                clim.to_netcdf(clim_filepath)

        comm.barrier()

        if rank == 0:
            assert(os.path.exists(clim_filepath))
            with xr.open_dataset(clim_filepath, engine = "netcdf4", decode_timedelta = False) as clim_ds:
                clim = clim_ds[plot_var]

                msg = "[{}]: Plotting {}.".format(datetime.now().strftime("%H:%M:%S"), plot_var)
                print(msg, flush = True)
                plot_clim(plot_var, clim, clim_fileroot, tag, plotting_dir)


def setup_pressure_grid_and_decomp(homme_output, spinup_days, comm):
    rank = comm.Get_rank()
    size = comm.Get_size()

    with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as homme_ds:
        p_ds = homme_ds["p"].sel(time = slice(spinup_days, None))
        np_p_tgt = p_ds.mean(dim = ["time", "lat", "lon"]).to_numpy()
        lat_vals = homme_ds["lat"].to_numpy()

    p_tgt = xr.DataArray(np_p_tgt, dims = ["p"],
        attrs = {"units": "Pa", "long_name": "Hydrostatic Pressure"}
    )

    nlat = len(lat_vals)
    lat_edges = np.linspace(0, nlat, size + 1, dtype = int)
    i0 = lat_edges[rank]
    i1 = lat_edges[rank + 1]
    lat_slice = slice(i0, i1)

    return p_tgt, lat_vals, lat_slice


def calc_local_clim(plot_var, homme_output, spinup_days, p_tgt, lat_slice):

    with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as homme_ds:

        if plot_var in ["u", "T", "pnh"]:
            val_ds = homme_ds[plot_var].sel(time = slice(spinup_days, None)).isel(lat = lat_slice)
        elif plot_var in ["T_eddy"]:
            val_ds = homme_ds[plot_var[0]].sel(time = slice(spinup_days, None)).isel(lat = lat_slice)

        p_ds = homme_ds["p"].sel(time = slice(spinup_days, None)).isel(lat = lat_slice)

        val = interp_to_p(val_ds, p_ds, p_tgt)

        if plot_var in ["T_eddy"]:
            val_zonal_mean = val.mean(dim = "lon", skipna = True)
            val = np.power(val - val_zonal_mean, 2)

        clim = val.mean(dim = ["time", "lon"], skipna = True)

        if plot_var in ["u", "T", "pnh"]:
            clim = clim.rename(plot_var)

        return clim


def ufunc_interp_to_p(val_col, p_col, p_tgt):
    # ASSUME: p_col, p_tgt in ascending order
    return np.interp(p_tgt, p_col, val_col,
        left = np.nan, right = np.nan)


def interp_to_p(val_ds, p_ds, p_tgt):

    val_on_p = xr.apply_ufunc(
        ufunc_interp_to_p, val_ds, p_ds, p_tgt,
        input_core_dims = [["lev"], ["lev"], ["p"]],
        output_core_dims = [["p"]],
        vectorize = True,
        output_dtypes = [val_ds.dtype]
    )
    val_on_p = val_on_p.assign_coords(p = p_tgt)

    return val_on_p.transpose("time", "p", "lat", "lon")


# Source - https://stackoverflow.com/a/43357954
# Posted by Maxim, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-20, License - CC BY-SA 4.0

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def plot_clim(plot_var, clim, clim_fileroot, tag, plotting_dir):
    var_labels = {"u" : r"Zonal Wind $\left[ m\,s^{-1} \right]$",
        "T" : r"Temperature $\left[ K \right]$",
        "pnh" : r"Pressure $\left[ hPa \right]$",
        "T_eddy" : r"Temperature Eddy Variation $\left[ K^{2} \right]$"}
    var_cmaps = {"u" : "RdBu",
        "T" : "plasma",
        "pnh" : "viridis",
        "T_eddy" : "plasma"}

    fig, axs = plt.subplots(sharex = True)

    # Color plot
    if plot_var in ["u"]:
        vmax = np.abs(clim).max()
        vmin = -vmax
    elif plot_var in ["T", "pnh"]:
        vmax = clim.max()
        vmin = clim.min()
    elif plot_var in ["T_eddy"]:
        vmax = clim.max()
        vmin = 0.

    lat = clim["lat"]
    p = clim["p"] / 100. # [Pa] => [hPa]

    cmap = var_cmaps[plot_var]
    clim_plt = axs.pcolormesh(lat, p, clim,
        vmin = vmin, vmax = vmax, cmap = cmap)
    axs.axvline([0], color = "grey")
    axs.axhline([pk02_p_sp], color = "grey", linestyle = "dashed")
    axs.axhline([pk02_p_T], color = "grey")

    # Colorbar
    cb = fig.colorbar(clim_plt, ax = axs)

    # Contours
    line_color = "k"
    pos_linestyle = "solid"
    neg_linestyle = "dashed"
    zero_linewidth = 2.0
    nonzero_linewidth = 1.0

    if plot_var in ["u"]:
        zero_levels = [0]
        pos_levels = [10, 20, 30, 40, 50, 60, 70, 80]
        neg_levels = [-80, -70, -60, -50, -40, -30, -20, -10]
    elif plot_var in ["T"]:
        zero_levels = None
        pos_levels = [180, 220, 260, 300]
        neg_levels = None
    else:
        zero_levels = None
        pos_levels = None
        neg_levels = None

    if zero_levels is not None:
        # Zero contour
        axs.contour(lat, p, clim, levels = zero_levels,
            colors = line_color, linewidths = zero_linewidth)
        cb.ax.axhline(zero_levels, color = line_color, linewidth = zero_linewidth)
    if pos_levels is not None:
        # Positive contours
        axs.contour(lat, p, clim, levels = pos_levels,
            colors = line_color, linewidths = nonzero_linewidth,
            linestyles = pos_linestyle)
        for level in pos_levels:
            cb.ax.axhline(level, color = line_color, linestyle = pos_linestyle,
                linewidth = nonzero_linewidth)
    if neg_levels is not None:
        # Negative contours
        axs.contour(lat, p, clim, levels = neg_levels,
            colors = line_color, linewidths = nonzero_linewidth,
            linestyles = neg_linestyle)
        for level in neg_levels:
            cb.ax.axhline(level, color = line_color, linestyle = neg_linestyle,
                linewidth = nonzero_linewidth)

    # Adjust y-axis
    axs.yaxis.set_inverted(True)
    axs.set_yscale("log")

    # Labels
    cb.set_label(var_labels[plot_var])

    fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
    fig.supylabel(r"Hydrostatic Pressure $\left[ hPa \right]$")
    suptitle = "Climatology"
    if tag:
        suptitle += " - {}".format(tag)
    fig.suptitle(suptitle)

    clim_plt_filepath = os.path.join(plotting_dir, clim_fileroot + ".png")
    plt.savefig(clim_plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()


if __name__ == "__main__":
    main()