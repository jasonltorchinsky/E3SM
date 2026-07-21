# Standard Library Imports
import os

# Third-Party Imports
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, NP_REAL, NP_ARRAY, XR_DATAARRAY, \
    MPL_FIGURE, MPL_AXES, MPL_PCOLORMESH, MPL_COLORBAR
from consts.physical import pk02_tropopause_pressure as pk02_p_T, pk02_sponge_pressure as pk02_p_sp

# Local "global" variables
plt_var_attrs = {"T" : {"label" : r"Temperature $\left[ K \right]$",
                        "cmap" : "plasma"},
    "u" : {"label" : r"Zonal Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
    "v" : {"label" : r"Meridional Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
    "w" : {"label" : r"Vertical Wind Velocity $\left[ m\,s^{-1} \right]$",
           "cmap" : "RdBu"},
}

def plot_climatology(field_clim: XR_DATAARRAY,
    plot_fileroot: str,
    tag: str,
    plotting_directory: str,
    set_ylim: bool = True):

    #---------------------------------------------------------------------------
    # Extract grids for plotting
    #---------------------------------------------------------------------------
    plot_var: str = field_clim.name
    lat: NP_ARRAY[NP_REAL] = NP_REAL(field_clim["lat"].to_numpy())
    p: NP_ARRAY[NP_REAL] = NP_REAL(field_clim["p"].to_numpy()) * 1.e-2 # [Pa] => [hPa]
    np_field_clim: NP_ARRAY[NP_REAL] = NP_REAL(field_clim.to_numpy())

    #---------------------------------------------------------------------------
    # Set plot parameters
    #---------------------------------------------------------------------------
    # Colorbar limits
    vmax: NP_REAL
    vmin: NP_REAL
    if plot_var in ["u", "v", "w"]:
        vmax = np.abs(field_clim).max()
        vmin = -vmax
    elif plot_var in ["T"]:
        vmax = field_clim.max()
        vmin = field_clim.min()

    # Colormap
    cmap: str = plt_var_attrs[plot_var]["cmap"]

    # Contour plot parameters
    contour_color: str = "k"
    pos_contourstyle: str = "solid"
    neg_contourstyle: str = "dashed"
    zero_contourwidth: NP_REAL = NP_REAL(2.0)
    nonzero_contourwidth: NP_REAL = NP_REAL(1.0)

    # Set contour levels
    zero_levels: Optional[NP_ARRAY[NP_INT]] = None
    if plot_var in ["u", "v", "w"]:
        zero_levels = np.zeros(1, dtype = NP_INT)
    elif plot_var in ["T"]:
        zero_levels = np.array([np.ceil(vmin / 20) * 20], dtype = NP_INT)

    pos_levels: Optional[NP_ARRAY[NP_INT]] = None
    neg_levels: Optional[NP_ARRAY[NP_INT]] = None
    if plot_var in ["u"]:
        pos_levels = np.arange(10, np.ceil(vmax / 10) * 10, 10, dtype = NP_INT)
        neg_levels = np.arange(-np.ceil(vmax / 10) * 10, 0, 10, dtype = NP_INT)
    elif plot_var in ["v"]:
        pos_levels = np.arange(1, np.ceil(vmax), dtype = NP_INT)
        neg_levels = np.arange(-np.ceil(vmax), 0, 1, dtype = NP_INT)
    elif plot_var in ["w"]:
        pos_levels = np.arange(1.e-3, np.ceil(vmax * 1.e3) * 1e-3, 1.e-3, dtype = NP_INT)
        neg_levels = np.arange(-np.ceil(vmax * 1.e3) * 1e-3, 0, 1.e-3, dtype = NP_INT)
    elif plot_var in ["T"]:
        pos_levels = np.arange((np.ceil(vmin / 20) + 1) * 20, np.ceil(vmax / 20) * 20, 20, dtype = NP_INT)

    #---------------------------------------------------------------------------
    # Create plot
    #---------------------------------------------------------------------------
    fig: MPL_FIGURE
    axs: MPL_AXES
    fig, axs = plt.subplots(sharex = True)

    # Climatology
    climatology_plot: MPL_PCOLORMESH = axs.pcolormesh(lat, p, np.transpose(np_field_clim),
        vmin = vmin, vmax = vmax, cmap = cmap,
        shading = "nearest",
        zorder = 0)

    # Colorbar
    climatology_colorbar: MPL_COLORBAR = fig.colorbar(climatology_plot, ax = axs)

    # Contour plot
    if zero_levels is not None:
        # Zero contour
        axs.contour(lat, p, np.transpose(np_field_clim), levels = zero_levels,
            colors = contour_color, 
            linewidths = zero_contourwidth,
            zorder = 2)
        climatology_colorbar.ax.axhline(zero_levels, 
            color = contour_color,
            linewidth = zero_contourwidth)
    
    if pos_levels is not None:
        # Positive contours
        axs.contour(lat, p, np.transpose(np_field_clim), levels = pos_levels,
            colors = contour_color, 
            linewidths = nonzero_contourwidth,
            linestyles = pos_contourstyle,
            zorder = 2)
        for level in pos_levels:
            climatology_colorbar.ax.axhline(level,
                color = contour_color,
                linestyle = pos_contourstyle,
                linewidth = nonzero_contourwidth)

    if neg_levels is not None:
        # Negative contours
        axs.contour(lat, p, np.transpose(np_field_clim), levels = neg_levels,
            colors = contour_color, 
            linewidths = nonzero_contourwidth,
            linestyles = neg_contourstyle,
            zorder = 2)
        for level in neg_levels:
            climatology_colorbar.ax.axhline(level,
                color = contour_color,
                linestyle = neg_contourstyle,
                linewidth = nonzero_contourwidth)

    # Lines to guide the eye
    axs.axvline([0], color = "grey")
    axs.axhline([pk02_p_sp],
        color = "grey",
        linestyle = "dashed",
        linewidth = 2.0,
        label = "PK02 Sponge Layer Height",
        zorder = 1)
    axs.axhline([pk02_p_T],
        color = "grey",
        linewidth = 2.0,
        label = "PK02 Nominal Tropopause Height",
        zorder = 1)

    #---------------------------------------------------------------------------
    # Final plot adjustments
    #---------------------------------------------------------------------------
    # y-axis adjustment
    axs.yaxis.set_inverted(True)
    axs.set_yscale("log")
    if set_ylim:
        axs.set_ylim([p.max(), 2.e-1])

    # Labels
    climatology_colorbar.set_label(plt_var_attrs[plot_var]["label"])

    fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
    fig.supylabel(r"Hydrostatic Pressure $\left[ hPa \right]$")
    suptitle = "Climatology"
    if tag:
        suptitle += " - {}".format(tag)
    fig.suptitle(suptitle)

    #---------------------------------------------------------------------------
    # Save plot to file
    #---------------------------------------------------------------------------
    clim_plt_filepath: str = os.path.join(plotting_directory, plot_fileroot + ".png")
    plt.savefig(clim_plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()