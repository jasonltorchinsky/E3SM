# Library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime

# Local imports
from us_std_atm import z_us_std_atm_p

# Constants
pk02_p_T = 100. # Polvani-Kushner 2002 nominal tropopause height [hPa]
pk02_p_sp = 0.5 # Polvani-Kushner 2002 sponge layer height [hPa]
p0 = 1000. # Base-state surface pressure [hPa]
ps = 1013.25 # Reference standard surface prssure [hPa]

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--interface-file", nargs = "?", required = True, type = str,
        help = "'vcoord' interface coordinates.")
    parser.add_argument("--midpoint-file", nargs = "?", required = True, type = str,
        help = "'vcoord' midpoint coordinates.")
    parser.add_argument("--vcoord-name", nargs = "?", default = "vcoord", type = str,
        help = "Directory to save plots.")
    parser.add_argument("--plotting-dir", nargs = "?", default = ".plotting", type = str,
        help = "Directory to save plots.")
    args = parser.parse_args()

    interface_filepath = os.path.normpath(args.interface_file)
    midpoint_filepath = os.path.normpath(args.midpoint_file)
    vcoord_name = args.vcoord_name
    plotting_dirpath = os.path.normpath(args.plotting_dir)

    dirpaths = [plotting_dirpath]
    for dirpath in dirpaths:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    # Recreate verical coordinate
    assert (vcoord_name in ["scream", "sab", "acme", "turbeville", "cam"])
    with open(interface_filepath) as vfile_int:
        interface_str = vfile_int.read()
        if vcoord_name in ["scream", "sab", "acme", "turbeville"]:
            interface_lines = [line for line in interface_str.splitlines() if line]
        elif vcoord_name in ["cam"]:
            interface_lines = [line for line in interface_str.split() if line[0].isdigit()]
        
        nilev = int(interface_lines[0].split()[0])
        hyai = np.array(interface_lines[1:nilev+1], dtype = np.float32) # A coefficients - Sela 2009, url: https://repository.library.noaa.gov/view/noaa/11401/noaa_11401_DS1.pdf
        hybi = np.array(interface_lines[nilev+2:], dtype = np.float32) # B coefficients - Sela 2009, url: https://repository.library.noaa.gov/view/noaa/11401/noaa_11401_DS1.pdf

    with open(midpoint_filepath) as vfile_mid:
        midpoint_str = vfile_mid.read()
        if vcoord_name in ["scream", "sab", "acme", "turbeville"]:
            midpoint_lines = [line for line in midpoint_str.splitlines() if line]
        elif vcoord_name in ["cam"]:
            midpoint_lines = [line for line in midpoint_str.split() if line[0].isdigit()]

        nlev = int(midpoint_lines[0].split()[0])
        hyam = np.array(midpoint_lines[1:nlev+1], dtype = np.float32) # A coefficients - Sela 2009, url: https://repository.library.noaa.gov/view/noaa/11401/noaa_11401_DS1.pdf
        hybm = np.array(midpoint_lines[nlev+2:], dtype = np.float32) # B coefficients - Sela 2009, url: https://repository.library.noaa.gov/view/noaa/11401/noaa_11401_DS1.pdf

    assert(nlev + 1 == nilev)

    pi = hyai * p0 + hybi * ps # Pressure at interfaces [hPa] (calculated as in components/homme/src/share/prim_driver_base.F90)
    zi = z_us_std_atm_p(pi) # Geometric height at interfaces [m]
    dzi = zi[:-1] - zi[1:] # Geometric thickness of layers NOTE: ToA is k = 0 [m]
    pm = hyam * p0 + hybm * ps # Pressure at midpoints [hPa]

    # Plot vertical coordinate
    fig, axs = plt.subplots(figsize = (9, 6.5))

    axs.plot(dzi, pm, color = "k", marker = ".", markersize = 5.0, zorder = 1)

    axs.axhline([pk02_p_sp], color = "grey", linestyle = "dashed", 
        label = "PK 2002 Sponge Layer Height", zorder = 0)
    axs.axhline([pk02_p_T], color = "grey", label = "PK 2002 Nominal Tropopause Height",
        zorder = 0)

    axs.legend()

    # Label axes
    axs.set_xlabel(r"Layer Thickness $\left[ m \right]$")
    axs.set_ylabel(r"Pressure $\left[ hPa \right]$")

    # Adjust y-axis
    axs.yaxis.set_inverted(True)
    axs.set_yscale("log")

    # Plot title
    if vcoord_name == "scream":
        title = "SCREAM - {} Levels".format(nlev)
    elif vcoord_name == "sab":
        title = r"Pure-$\sigma$ (Ignoring 0 $hPa$) - {} Levels".format(nlev)
    if vcoord_name == "acme":
        title = "ACME - {} Levels".format(nlev)
    if vcoord_name == "cam":
        title = "CAM - {} Levels".format(nlev)
    else:
        title = None
    axs.set_title(title)

    clim_plt_filepath = os.path.join(plotting_dirpath, vcoord_name + "_{}.png".format(nlev))
    plt.savefig(clim_plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

if __name__ == "__main__":
    main()