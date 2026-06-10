# Library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime
import xarray as xr
from mpi4py import MPI

MPI_ROOT = 0
pk02_p_T = 100. # Polvani-Kushner 2002 nominal tropopause height [hPa]
pk02_p_sp = 0.5 # Polvani-Kushner 2002 sponge layer height [hPa]

xr_var_attrs = {"T" : {"units" : "K",
               "description" : "Temperature"},
    "u" : {"units" : "m s^{-1}",
           "description" : "Zonal Wind Velocity"},
    "v" : {"units" : "m s^{-1}",
           "description" : "Meridional Wind Velocity"},
    "w" : {"units" : "m s^{-1}",
           "description" : "Vertical Wind Velocity"},
}

plt_var_attrs = {"T" : {"label" : r"Temperature $\left[ K \right]$",
                        "cmap" : "plasma"},
    "u" : {"label" : r"Zonal Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
    "v" : {"label" : r"Meridional Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
    "w" : {"label" : r"Vertical Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
}

coord_attrs = { "p" : {"units" : "Pa",
                       "description" : "Hydrostatic pressure"},
    "lat" : {"units" : "degrees",
             "description" : "Latitude",
             "range" : "-90 to 90"}
}


def main():
    #---------------------------------------------------------------------------
    # Set up MPI communicator
    #---------------------------------------------------------------------------
    comm = MPI.COMM_WORLD
    l_rank = comm.Get_rank()
    comm_size = comm.Get_size()

    #---------------------------------------------------------------------------
    # Parse command-line input
    #---------------------------------------------------------------------------
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
    homme_output = os.path.normpath(args.homme_output)
    working_dir = os.path.normpath(args.working_dir)
    plotting_dir = os.path.normpath(args.plotting_dir)
    recalculate = args.recalculate
    tag = args.tag
    plot_vars = [str(plot_var) for plot_var in args.plot_vars.split(",") if plot_var]
    
    #---------------------------------------------------------------------------
    # Create working directories
    #---------------------------------------------------------------------------
    if l_rank == MPI_ROOT:
        dirs = [working_dir, plotting_dir]
        for dir in dirs:
            if not os.path.exists(dir):
                os.makedirs(dir)

    #---------------------------------------------------------------------------
    # Calculate climatologies and create plots
    #---------------------------------------------------------------------------
    if l_rank == MPI_ROOT:
        datetime_now = datetime.now().strftime("%H:%M:%S")
        msg = "[{}]: Starting climatology calculation and plotting loop.".format(datetime_now)
        print(msg, flush = True)

    for plot_var in plot_vars:
        if l_rank == MPI_ROOT:
            datetime_now = datetime.now().strftime("%H:%M:%S")
            msg = "[{}]: Starting {}.".format(datetime_now, plot_var)
            print(msg, flush = True)

        clim_fileroot = plot_var + "_clim"
        if tag:
            clim_fileroot += "_{}".format(tag)

        clim_filepath = os.path.join(working_dir, clim_fileroot + ".nc")

        if recalculate or (not os.path.exists(clim_filepath)):
            #-------------------------------------------------------------------
            # Decompose time grid
            #-------------------------------------------------------------------
            if l_rank == MPI_ROOT:
                datetime_now = datetime.now().strftime("%H:%M:%S")
                msg = "[{}]: Setting up time grid decomposition.".format(datetime_now)
                print(msg, flush = True)

            with (xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)
                .isel(time = slice(spinup_days, None))) as homme_ds:
                ntime = homme_ds["time"].size
                lat = homme_ds["lat"].to_numpy()
            rank_ntimes = (ntime // comm_size) \
                + (np.arange(comm_size) < (ntime - comm_size * (ntime // comm_size))).astype(np.int64)
            l_ntime = rank_ntimes[l_rank]
            l_time_st_idx = spinup_days + np.sum(rank_ntimes[:l_rank])
            l_time_end_idx = l_time_st_idx + l_ntime
            l_time_slice = slice(l_time_st_idx, l_time_end_idx)

            #-------------------------------------------------------------------
            # Set up target pressure grid, interpolate vertically to it,
            # and calculate climatology
            #-------------------------------------------------------------------
            if l_rank == MPI_ROOT:
                datetime_now = datetime.now().strftime("%H:%M:%S")
                msg = "[{}]: Creating target pressure grid and calculating climatology.".format(datetime_now)
                print(msg, flush = True)
            with (xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)
                .isel(time = l_time_slice)) as homme_ds:
                p_tgt = get_p_tgt(homme_ds, comm)
                l_field_vremap = vremap_field(homme_ds, p_tgt, plot_var, comm)
            g_field_clim = calc_clim(l_field_vremap, comm)

            #-------------------------------------------------------------------
            # Save climatology to file
            #-------------------------------------------------------------------
            if l_rank == MPI_ROOT:
                datetime_now = datetime.now().strftime("%H:%M:%S")
                msg = "[{}]: Saving climatology to file.".format(datetime_now)
                print(msg, flush = True)

            if l_rank == MPI_ROOT:
                g_clim_ds = xr.Dataset(
                    data_vars = {
                        plot_var : (["p", "lat"], g_field_clim, xr_var_attrs[plot_var])
                    },
                    coords = {
                        "p" : ("p", p_tgt, coord_attrs["p"]),
                        "lat" : ("lat", lat, coord_attrs["lat"])
                    }
                )

                g_clim_ds.to_netcdf(clim_filepath, mode = "w")

        if l_rank == MPI_ROOT:
            assert(os.path.exists(clim_filepath))
            with xr.open_dataset(clim_filepath, engine = "netcdf4", decode_timedelta = False) as clim_ds:
                clim = clim_ds[plot_var]

                datetime_now = datetime.now().strftime("%H:%M:%S")
                msg = "[{}]: Plotting {}.".format(datetime_now, plot_var)
                print(msg, flush = True)

                plot_clim(plot_var, clim, clim_fileroot, tag, plotting_dir)

def get_p_tgt(homme_ds, comm):
    #---------------------------------------------------------------------------
    # Get MPI communicator parameters
    #---------------------------------------------------------------------------
    l_rank = comm.Get_rank()

    #---------------------------------------------------------------------------
    # Set a constant p_tgt - TO-DO: Caluclate this based on time mean of p_src
    #---------------------------------------------------------------------------
    p_tgt = np.logspace(np.log10(0.1e2), np.log10(1000.e2), num = 256) # [Pa]

    return p_tgt

def vremap_field(homme_ds, p_tgt, plot_var, comm):
    supported_plot_vars = ["T", "u", "v", "w"]
    assert(plot_var in supported_plot_vars)

    #---------------------------------------------------------------------------
    # Get MPI communicator parameters
    #---------------------------------------------------------------------------
    l_rank = comm.Get_rank()

    #---------------------------------------------------------------------------
    # Obtain p_src, field value from file
    #---------------------------------------------------------------------------
    if plot_var in ["T", "u", "v", "w"]:
        field_key = plot_var

    p_src = homme_ds["p"].to_numpy() # Pressure [Pa], [l_nt, nz, nlat, nlon]
    l_field_src = homme_ds[field_key].to_numpy() # [l_nt, nz, nlat, nlon]

    # Reshape for better memeory access
    p_src = np.transpose(p_src, axes = [0, 2, 3, 1]) # [Pa], [l_nt, nlat, nlon, nz]
    l_field_src = np.transpose(l_field_src, axes = [0, 2, 3, 1]) # [l_nt, nlat, nlon, nz]

    #---------------------------------------------------------------------------
    # Perform vertical interpolation
    #---------------------------------------------------------------------------
    [l_nt, nlat, nlon, _] = l_field_src.shape
    [nz_tgt] = p_tgt.shape
    l_field_tgt = np.empty([l_nt, nlat, nlon, nz_tgt], dtype = l_field_src.dtype)

    for tt in range(0, l_nt):
        for jj in range(0, nlat):
            for ii in range(0, nlon):
                l_field_tgt[tt,jj,ii,:] = np.interp(p_tgt, p_src[tt,jj,ii,:],
                    l_field_src[tt,jj,ii,:], left = np.nan, right = np.nan)

    l_field_tgt = np.ascontiguousarray(np.transpose(l_field_tgt, axes = [0, 3, 1, 2])) # [l_nt, nz, nlat, nlon]

    return l_field_tgt

def calc_clim(l_field_vremap, comm):
    #---------------------------------------------------------------------------
    # Get MPI communicator parameters
    #---------------------------------------------------------------------------
    l_rank = comm.Get_rank()

    #---------------------------------------------------------------------------
    # Get sums and counts locally
    #---------------------------------------------------------------------------
    l_not_nan_count = np.sum(~np.isnan(l_field_vremap), axis = (0, 3))
    l_field_sum = np.sum(l_field_vremap, axis = (0, 3), where = ~np.isnan(l_field_vremap)) # [nz, nlat]

    #---------------------------------------------------------------------------
    # Root process gathers and calculates climatology
    #---------------------------------------------------------------------------
    g_not_nan_count = comm.gather(l_not_nan_count, root = MPI_ROOT) # On root is [comm_size * [nz, nlat]]
    g_field_sum = comm.gather(l_field_sum, root = MPI_ROOT) # On root is [comm_size * [nz, nlat]]

    g_field_clim = None
    if l_rank == MPI_ROOT:
        g_not_nan_count = np.sum(np.stack(g_not_nan_count), axis = (0)) # [nz, nlat]
        g_field_sum = np.stack(g_field_sum) # [comm_size, nz, nlat]
        g_field_sum = np.sum(g_field_sum, axis = (0), where = ~np.isnan(g_field_sum)) # [nz, nlat]

        nonzero_mask = (g_not_nan_count != 0) # [nz, nlat]

        g_field_clim = np.full_like(g_field_sum, np.nan)
        g_field_clim[nonzero_mask] = g_field_sum[nonzero_mask] / g_not_nan_count[nonzero_mask] # [nz, nlat]

    return g_field_clim

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


def plot_clim(plot_var, clim, clim_fileroot, tag, plotting_dir, set_ylim = True):
    fig, axs = plt.subplots(sharex = True)

    # Color plot
    if plot_var in ["u", "v", "w"]:
        vmax = np.abs(clim).max()
        vmin = -vmax
    elif plot_var in ["T"]:
        vmax = clim.max()
        vmin = clim.min()

    lat = clim["lat"]
    p = clim["p"] / 100. # [Pa] => [hPa]

    cmap = plt_var_attrs[plot_var]["cmap"]
    clim_plt = axs.pcolormesh(lat, p, clim,
        vmin = vmin, vmax = vmax, cmap = cmap)
    axs.axvline([0], color = "grey")
    axs.axhline([pk02_p_sp], color = "grey", linestyle = "dashed", linewidth = 2.0,
        label = "PK 2002 Sponge Layer Height")
    axs.axhline([pk02_p_T], color = "grey", linewidth = 2.0,
        label = "PK02 Nominal Tropopause Height")

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
        pos_levels = np.arange(10, np.ceil(vmax / 10) * 10, 10)
        neg_levels = np.arange(-np.ceil(vmax / 10) * 10, 0, 10)
    elif plot_var in ["v"]:
        zero_levels = [0]
        pos_levels = np.arange(1, np.ceil(vmax))
        neg_levels = np.arange(-np.ceil(vmax), 0, 1)
    elif plot_var in ["w"]:
        zero_levels = [0]
        pos_levels = np.arange(1.e-3, np.ceil(vmax * 1.e3) * 1e-3, 1.e-3)
        neg_levels = np.arange(-np.ceil(vmax * 1.e3) * 1e-3, 0, 1.e-3)
    elif plot_var in ["T"]:
        zero_levels = [np.ceil(vmin / 20) * 20]
        pos_levels = np.arange((np.ceil(vmin / 20) + 1) * 20, np.ceil(vmax / 20) * 20, 20)
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
    if set_ylim:
        axs.set_ylim([p.max(), 0.2])

    # Legend
    #axs.legend()

    # Labels
    cb.set_label(plt_var_attrs[plot_var]["label"])

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