# Library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

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
    parser.add_argument("--recalculate", nargs = "?", default = False, type = bool,
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
    # Zonal wind (u) climatology
    #---------------------------------------------------------------------------
    if "u" in plot_vars:
        u_clim_fileroot = "u_clim"
        if tag:
            u_clim_fileroot += "_{}".format(tag)

        u_clim_filepath = os.path.join(working_dir, u_clim_fileroot + ".nc")

        if not recalculate and os.path.exists(u_clim_filepath):
            u_clim = xr.open_dataset(u_clim_filepath, engine = "netcdf4", decode_timedelta = False)["u"]
        else:
            u_ds = xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)["u"].sel(time = slice(spinup_days, None))
            u_clim = u_ds.mean(dim = ["time", "lon"])

            u_clim.to_netcdf(u_clim_filepath)

        fig, axs = plt.subplots(sharex = True)

        vmax = np.abs(u_clim).max()
        vmin = -vmax
        cmap = "RdBu"
        u_clim_plt = axs.pcolormesh(u_clim["lat"], u_clim["lev"], u_clim,
            vmin = vmin, vmax = vmax, cmap = cmap)
        axs.axvline([0], color = "grey")

        levels = [0]
        colors = "k"
        axs.contour(u_clim["lat"], u_clim["lev"], u_clim, levels = levels, colors = colors)

        axs.yaxis.set_inverted(True)

        cb = fig.colorbar(u_clim_plt, ax = axs)
        cb.ax.axhline(levels, color = colors)
        cb.set_label(r"Zonal Wind $\left[ m\,s^{-1} \right]$")

        fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
        fig.supylabel("Level")
        suptitle = "Climatology"
        if tag:
            suptitle += " - {}".format(tag)
        fig.suptitle(suptitle)

        u_clim_plt_filepath = os.path.join(plotting_dir, u_clim_fileroot + ".png")
        plt.savefig(u_clim_plt_filepath, dpi = 256, bbox_inches = "tight")

        plt.close()

    #---------------------------------------------------------------------------
    # Temperature (T) climatology
    #---------------------------------------------------------------------------
    if "T" in plot_vars:
        T_clim_fileroot = "T_clim"
        if tag:
            T_clim_fileroot += "_{}".format(tag)

        T_clim_filepath = os.path.join(working_dir, T_clim_fileroot + ".nc")

        if not recalculate and os.path.exists(T_clim_filepath):
            T_clim = xr.open_dataset(T_clim_filepath, engine = "netcdf4", decode_timedelta = False)["T"]
        else:
            T_ds = xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)["T"].sel(time = slice(spinup_days, None))
            T_clim = T_ds.mean(dim = ["time", "lon"])

            T_clim.to_netcdf(T_clim_filepath)

        fig, axs = plt.subplots(sharex = True)

        vmax = T_clim.max()
        vmin = T_clim.min()
        cmap = "plasma"
        T_clim_plt = axs.pcolormesh(T_clim["lat"], T_clim["lev"], T_clim,
            vmin = vmin, vmax = vmax, cmap = cmap)
        axs.axvline([0], color = "grey")

        axs.yaxis.set_inverted(True)

        cb = fig.colorbar(T_clim_plt, ax = axs)
        cb.set_label(r"Temperature $\left[ K \right]$")

        fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
        fig.supylabel("Level")
        suptitle = "Climatology"
        if tag:
            suptitle += " - {}".format(tag)
        fig.suptitle(suptitle)

        T_clim_plt_filepath = os.path.join(plotting_dir, T_clim_fileroot + ".png")
        plt.savefig(T_clim_plt_filepath, dpi = 256, bbox_inches = "tight")

        plt.close()

    #---------------------------------------------------------------------------
    # Pressure (pnh) climatology
    #---------------------------------------------------------------------------
    if "pnh" in plot_vars:
        pnh_clim_fileroot = "pnh_clim"
        if tag:
            pnh_clim_fileroot += "_{}".format(tag)

        pnh_clim_filepath = os.path.join(working_dir, pnh_clim_fileroot + ".nc")

        if not recalculate and os.path.exists(pnh_clim_filepath):
            pnh_clim = xr.open_dataset(pnh_clim_filepath, engine = "netcdf4", decode_timedelta = False)["pnh"]
        else:
            pnh_ds = xr.open_dataset(homme_output, engine = "netcdf4", decode_timedelta = False)["pnh"].sel(time = slice(spinup_days, None))
            pnh_clim = pnh_ds.mean(dim = ["time", "lon"]) / 100. # [Pa to hPa]

            pnh_clim.to_netcdf(pnh_clim_filepath)

        fig, axs = plt.subplots(sharex = True)

        vmax = pnh_clim.max()
        vmin = pnh_clim.min()
        cmap = "viridis"
        pnh_clim_plt = axs.pcolormesh(pnh_clim["lat"], pnh_clim["lev"], pnh_clim,
            vmin = vmin, vmax = vmax, cmap = cmap)
        axs.axvline([0], color = "grey")

        axs.yaxis.set_inverted(True)

        cb = fig.colorbar(pnh_clim_plt, ax = axs)
        cb.set_label(r"Pressure $\left[ hPa \right]$")

        fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
        fig.supylabel("Level")
        suptitle = "Climatology"
        if tag:
            suptitle += " - {}".format(tag)
        fig.suptitle(suptitle)

        pnh_clim_plt_filepath = os.path.join(plotting_dir, pnh_clim_fileroot + ".png")
        plt.savefig(pnh_clim_plt_filepath, dpi = 256, bbox_inches = "tight")

        plt.close()

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

if __name__ == "__main__":
    main()