import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import getopt, sys
import xarray as xr

import cmaps
import geocat.viz as gv

def main(argv):

    # Physical constants
    g = 9.801 # Acceleration due to gravity [m s^(-2)]
    # Variables reused for each plot
    file_name = './movies/dcmip2016_test21.nc' # Data file
    uv_max = 15  # Colorbar max [m s^(-1)]
    uv_min = 0  # Colorbar min [m s^(-1)]
    uv_levs = np.linspace(uv_min, uv_max, 40)
    uv_ticks = np.linspace(uv_min, uv_max, 5, dtype = np.uint)
    days = [5, 10, 15, 20]
    cmap = cmaps.WhiteYellowOrangeRed

    projection = ccrs.PlateCarree(central_longitude = -150)

    # Set up figure
    nrows = 2
    ncols = 2
    nplots = nrows * ncols
    fig, axs = plt.subplots(nrows = nrows,
                            ncols = ncols,
                            sharex = True,
                            sharey = True,
                            subplot_kw = {'projection': projection})

    for ax_idx in range(nplots):
        ax = axs.flatten()[ax_idx]
        day = days[ax_idx]

        col_idx = ax_idx % ncols
        row_idx = int(ax_idx / ncols)

        # Open datafile
        ds = xr.open_dataset(file_name)
        
        # Extract zonal meridional wind velocities, shift 180 degrees longitude
        u = ds.u.isel({'time' : day, 'lev' : -1}) # [m s^(-1)]
        v = ds.v.isel({'time' : day, 'lev' : -1}) # [m s^(-1)]
        uv = np.sqrt(u**2 + v**2)

        nlat, nlon = np.shape(u.values)
        uv = uv.roll(lon = int(nlon/2))
        
        # Fix the artifact of not-shown-data around 0 and 360-degree longitudes
        uv = gv.xr_add_cyclic_longitudes(uv, "lon")
        
        ax.coastlines(linewidths=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgray")
        
        # Contourf-plot data
        cplt = uv.plot.contourf(ax = ax,
                                transform = projection,
                                levels = uv_levs,
                                vmin = uv_min,
                                vmax = uv_max,
                                cmap = cmap,
                                add_colorbar = False,
                                add_labels = False)
    
        # Set plot range, tick labels, etc.
        ax.set_xlim(-180, 90)
        ax.set_ylim(-40, 90)
        xticks_major = np.arange(-180, 180, 90)
        xticks_minor = np.setdiff1d(np.arange(-180, 180, 30), xticks_major)
        xticks_labels = [r'$30^{\circ} E$', r'$120^{\circ} E$',
                         r'$150^{\circ} W$', r'$60^{\circ} W$']
        yticks_major = np.arange(-30, 90, 30)
        yticks_minor = np.setdiff1d(np.arange(-90, 90, 10), yticks_major)
        yticks_labels = [r'$30^{\circ} S$', r'$EQ$',
                         r'$30^{\circ} N$', r'$60^{\circ} N$']
        if row_idx == nrows - 1:
            ax.set_xticks(ticks = xticks_major,
                          minor = False)
            ax.set_xticks(ticks = xticks_minor,
                          minor = True)
            ax.set_xticklabels(labels = xticks_labels)
                           
        if col_idx == 0:
            ax.set_yticks(ticks = yticks_major,
                          minor = False)
            ax.set_yticks(ticks = yticks_minor,
                          minor = True)
            ax.set_yticklabels(labels = yticks_labels)
        
        # Label date in plot to save space
        day_str = "Day {}".format(day)
        ax.text(0.05, 0.25, day_str,
                fontsize = 8,
                verticalalignment = 'top',
                transform = ax.transAxes,
                bbox = {'facecolor' : 'w',
                        'alpha' : 0.95})
        
    plt.subplots_adjust(left = 0.075,
                        bottom = 0.05,
                        right = 0.975,
                        top = 1.1,
                        wspace = 0.05,
                        hspace = -0.65)
    
    fig.suptitle("DCMIP2016 Test 2 - Tropical Cyclone",
                 x = .5,
                 y = .98,
                 fontsize = 16)
    
    # Add color bar
    cbar = fig.colorbar(cplt,
                        ax = axs[1,:],
                        orientation = 'horizontal',
                        shrink = 0.8,
                        pad = 0.1,
                        extendrect = True,
                        ticks = uv_ticks,
                        label = r"Surface Wind Speed $[m\,s^{-1}]$")
    
    cbar.ax.tick_params(labelsize = 10)

    fig.set_size_inches(6.5, 3.95)
    fig.savefig('uvs.pdf')
    
if __name__ == "__main__":
    main(sys.argv[1:])
