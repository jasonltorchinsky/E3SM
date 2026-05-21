# Library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

pk02_p_T = 100. # Polvani-Kushner 2002 nominal tropopause height [hPa]

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
    plot_vars = [str(plot_var) for plot_var in args.plot_vars.split(",")]

    dirs = [working_dir, plotting_dir]
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)

    #---------------------------------------------------------------------------
    # Set up variables for plotting
    #---------------------------------------------------------------------------
    p_ds = None
    p_tgt = None

    #---------------------------------------------------------------------------
    # Create plots
    #---------------------------------------------------------------------------
    for plot_var in plot_vars:
        clim_fileroot = plot_var + "_clim"
        if tag:
            clim_fileroot += "_{}".format(tag)

        clim_filepath = os.path.join(working_dir, clim_fileroot + ".nc")

        if recalculate or (not os.path.exists(clim_filepath)):
            if plot_var in ["u", "T", "pnh"]:
                with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as homme_ds:
                    val_ds = homme_ds[plot_var].sel(time = slice(spinup_days, None))

            # Interpolate value to fixed pressure levels
            val_on_p = interp_to_p(val_ds, p_ds, p_tgt,
                spinup_days = spinup_days, homme_output = homme_output)
            clim = val_on_p.mean(dim = ["time", "lon"], skipna = True)
            if plot_var in ["pnh"]:
                clim = clim / 100. # [Pa] => [hPa]
                clim.attrs["units"] = "hPa"
            clim.to_netcdf(clim_filepath)

        assert(os.path.exists(clim_filepath))
        with xr.open_dataset(clim_filepath, engine = "netcdf4", decode_timedelta = False) as clim_ds:
            clim = clim_ds[plot_var]

        plot_clim(plot_var, clim, clim_fileroot, tag, plotting_dir)


    # TO-DO: INCLUDE THESE IN TEMPLATE ABOVE

    #---------------------------------------------------------------------------
    # Temperature eddy variance (T*^2) climatology
    #---------------------------------------------------------------------------
    if "T_eddy" in plot_vars:
        # Get T eddy variance climatology
        T_eddy_clim_fileroot = "T_eddy_clim"
        if tag:
            T_eddy_clim_fileroot += "_{}".format(tag)

        T_eddy_clim_filepath = os.path.join(working_dir, T_eddy_clim_fileroot + ".nc")

        if not recalculate and os.path.exists(T_eddy_clim_filepath):
            T_eddy_clim = xr.open_dataset(T_eddy_clim_filepath, engine = "netcdf4", decode_timedelta = False)["T"]
        else:
            T_ds = xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)["T"].sel(time = slice(spinup_days, None))
            T_zonal_mean = T_ds.mean(dim = "lon")

            T_eddy_clim = np.pow(T_ds - T_zonal_mean, 2).mean(dim = ["time", "lon"])
            T_eddy_clim.to_netcdf(T_eddy_clim_filepath)

        fig, axs = plt.subplots(sharex = True)

        vmax = T_eddy_clim.max()
        vmin = 0.0
        cmap = "plasma"
        T_eddy_clim_plt = axs.pcolormesh(T_eddy_clim["lat"], T_eddy_clim["lev"], T_eddy_clim,
            vmin = vmin, vmax = vmax, cmap = cmap)
        axs.axvline([0], color = "grey")

        axs.yaxis.set_inverted(True)

        cb = fig.colorbar(T_eddy_clim_plt, ax = axs)
        cb.set_label(r"Temperature Eddy Variation $\left[ K^{2} \right]$")

        fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
        fig.supylabel("Level")
        suptitle = "Climatology"
        if tag:
            suptitle += " - {}".format(tag)
        fig.suptitle(suptitle)

        T_eddy_clim_plt_filepath = os.path.join(plotting_dir, T_eddy_clim_fileroot + ".png")
        plt.savefig(T_eddy_clim_plt_filepath, dpi = 256, bbox_inches = "tight")

        plt.close()

def ufunc_interp_to_p(val_col, p_col, p_tgt):
    # ASSUME: p_col, p_tgt in ascending order
    return np.interp(p_tgt, p_col, val_col,
        left = np.nan, right = np.nan)


def interp_to_p(val_ds, p_ds, p_tgt, spinup_days = None, homme_output = None):
    assert((p_ds is not None) and (p_tgt is not None)
           or (spinup_days is not None) and (homme_output is not None))
    
    if (p_ds is None) or (p_tgt is None):
        # Read in HOMME data
        with xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False) as homme_ds:
            p_ds = homme_ds["p"].sel(time = slice(spinup_days, None)) # Hydrostatic pressure [Pa]

        # Get target pressure grid
        np_p_tgt = p_ds.mean(dim = ["time", "lat", "lon"]).values
        p_tgt = xr.DataArray(np_p_tgt, dims = ["p"],
            attrs = {"units": "Pa", "long_name": "Hydrostatic Pressure"}
        )

    val_on_p = xr.apply_ufunc(
        ufunc_interp_to_p, val_ds, p_ds, p_tgt,
        input_core_dims = [["lev"], ["lev"], ["p"]],
        output_core_dims = [["p"]],
        vectorize = True,
        dask = "parallelized",
        output_dtypes = [val_ds.dtype]
    )
    val_on_p = val_on_p.assign_coords(p = p_tgt)

    return val_on_p.transpose()

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
        "pnh" : r"Pressure $\left[ hPa \right]$"}
    var_cmaps = {"u" : "RdBu",
        "T" : "plasma",
        "pnh" : "viridis"}
    
    fig, axs = plt.subplots(sharex = True)

    # Color plot
    if plot_var in ["u"]:
        vmax = np.abs(clim).max()
        vmin = -vmax
    elif plot_var in ["T", "pnh"]:
        vmax = clim.max()
        vmin = clim.min()

    lat = clim["lat"]
    p = clim["p"] / 100. # [Pa] => [hPa]

    cmap = var_cmaps[plot_var]
    clim_plt = axs.pcolormesh(lat, p, clim,
        vmin = vmin, vmax = vmax, cmap = cmap)
    axs.axvline([0], color = "grey")
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
        pos_levels = [10, 30, 50]
        neg_levels = [-50, -30, -10]
    elif plot_var in ["T"]:
        zero_levels = None
        pos_levels = [180, 220, 260, 300]
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