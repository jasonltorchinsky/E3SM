import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import geocat.viz as gv

def main():

    # Variables reused for each plot
    g = 9.801 # Acceleration due to gravity [m s^(-1)]
    styles = ['k--', 'k-.', 'b-', 'b--', 'r--', 'r-']
    
    uv_max = 35 # y-axis max
    uv_min = 0  # y-axis min

    # Set up figure
    fig, ax = plt.subplots(nrows = 1,
                           ncols = 1)

        
    file_name = './movies/dcmip2016_test21.nc' # Data file
    # Get lowest-level winds
    ds = xr.open_dataset(file_name)
    t = ds.time.to_numpy() / (60 * 60 * 24 * 10**(9)) # [s? -> d]
    u = ds.u.isel({'lev' : -1})
    v = ds.v.isel({'lev' : -1})
        
    ndays, nlats, nlons = np.shape(np.asarray(u))
        
    # Calculate greatest surface winds
    uv = np.sqrt(u**2 + v**2).to_numpy()
    max_uv = np.zeros([ndays])
    for day in range(ndays):
        max_uv[day] = np.amax(uv[day,:,:])
            
    ax.plot(t, max_uv, 
            color = 'black',
            linestyle = 'solid')

    # Plot intensity categories
    full_t = np.arange(t[0] - 1, t[-1] + 1)
    trop_dep = 16.9767 * np.ones_like(full_t)
    trop_storm = 24.1789 * np.ones_like(full_t)
    sev_trop_storm = 32.41 * np.ones_like(full_t)

    ax.plot(full_t, trop_dep,
            color = 'black',
            linestyle = 'dotted',
            linewidth = 0.5)

    ax.plot(full_t, trop_storm,
            color = 'black',
            linestyle = 'dotted',
            linewidth = 0.5)

    ax.plot(full_t, sev_trop_storm,
            color = 'black',
            linestyle = 'dotted',
            linewidth = 0.5)

    ax.text(0.0125, 0.975, 'Typhoon',
                fontsize = 8,
                verticalalignment = 'top',
                transform = ax.transAxes)

    ax.text(0.0125, 0.9, 'Severe Tropical Storm',
                fontsize = 8,
                verticalalignment = 'top',
                transform = ax.transAxes)

    ax.text(0.0125, 0.67, 'Tropical Storm',
                fontsize = 8,
                verticalalignment = 'top',
                transform = ax.transAxes)

    ax.text(0.0125, 0.2, 'Tropical Depression',
                fontsize = 8,
                verticalalignment = 'top',
                transform = ax.transAxes)
    
    # Set plot range, tick labels, etc.
    xticks_major = np.arange(0, ndays, 5)
    xticks_minor = np.setdiff1d(np.arange(0, ndays, 1), xticks_major)
    
    yticks_major = np.arange(uv_min, uv_max, 5)
    yticks_minor = np.setdiff1d(np.arange(uv_min, uv_max, 2.5), yticks_major)
    
    ax.set_xticks(ticks = xticks_major,
                  minor = False)
    ax.set_xticks(ticks = xticks_minor,
                  minor = True)
    
    ax.set_yticks(ticks = yticks_major,
                  minor = False)
    ax.set_yticks(ticks = yticks_minor,
                  minor = True)

    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(uv_min, uv_max)
    
    ax.set_xlabel(r'Time $[d]$')
    ax.set_ylabel(r'Maximum Surface Wind Speed $[m\ s^{-1}]$')


    title = "Tropical Cyclone Intensity"
    fig.suptitle(title,
                 x = .5,
                 y = .98,
                 fontsize = 16)
    
    fig.set_size_inches(6.5, 3.75)
    fig.savefig('intensity.pdf')
    
if __name__ == "__main__":
    main()
