# library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

# Constants
r0 = 6356766. # Mean Radius of Earth [m]
g0 = 9.80665 # Sea-Level Acceleration of Gravity = Relation of Standard Geopotential to Geometric Meter (g0') [m s^{-2}] = [m^2 s^{-2} (m')^{-1}]
R_star = 8.31432e3 # Gas Constant [N m K^{-1} kmol^{-1}]
M0 = 28.9644 # Mean Molecular Weight of Air at Sea-Level [kg kmol^{-1}]

# Defined U.S. Standard Atmosphere Parameters
H_b = np.array([0, 11, 20, 32, 47, 51, 71, 84.8520]) * 1000 # Geopotential Height at Reference Levels [m']
L_Mb = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]) * 0.001 # Molecular-Scale Temperature Gradient Between Reference Levels [K m'^{-1}]

# Derived U.S. Standard Atmosphere Paremeters
Z_b = (r0 * H_b) / (r0 - H_b) # Geometric Height at Reference Levels [m]
T_Mb = np.zeros(8) # Molecular-Scale Temperature at Reference Levels [K]
T_Mb[0] = 288.15
for ii in range(1, 8):
    T_Mb[ii] = T_Mb[ii-1] + L_Mb[ii-1] * (H_b[ii] - H_b[ii-1])

p_b = np.zeros(8) # Pressure at Reference Levels [hPa]
p_b[0] = 1013.25 # [hPa]
for ii in range(1, 8):
    if L_Mb[ii-1] != 0.0:
        p_b[ii] = p_b[ii-1] * (T_Mb[ii-1] / (T_Mb[ii-1] + L_Mb[ii-1] * (H_b[ii] - H_b[ii-1])))**((g0 * M0) / (R_star * L_Mb[ii-1]))
    else:
        p_b[ii] = p_b[ii-1] * np.exp(-(((g0 * M0) * (H_b[ii] - H_b[ii-1])) / (R_star * T_Mb[ii-1])))

def T_us_std_atm_z(z):
    h = (r0 * z) / (r0 + z) # Convert to Geopotential Height [m]
    T = np.piecewise(
        h,
        [
            (h >= H_b[0])  & (h <  H_b[1]),
            (h >= H_b[1])  & (h <  H_b[2]),
            (h >= H_b[2])  & (h <  H_b[3]),
            (h >= H_b[3])  & (h <  H_b[4]),
            (h >= H_b[4])  & (h <  H_b[5]),
            (h >= H_b[5])  & (h <  H_b[6]),
            (h >= H_b[6])  & (h <= H_b[7]),
        ],
        [
            lambda h: T_Mb[0] + L_Mb[0] * (h - H_b[0]), # 0  - 11 km'
            lambda h: T_Mb[1] + L_Mb[1] * (h - H_b[1]), # 11 - 20 km'
            lambda h: T_Mb[2] + L_Mb[2] * (h - H_b[2]), # 20 - 32 km'
            lambda h: T_Mb[3] + L_Mb[3] * (h - H_b[3]), # 32 - 47 km'
            lambda h: T_Mb[4] + L_Mb[4] * (h - H_b[4]), # 47 - 51 km'
            lambda h: T_Mb[5] + L_Mb[5] * (h - H_b[5]), # 51 - 71 km'
            lambda h: T_Mb[6] + L_Mb[6] * (h - H_b[6]), # 71 - 86 km'
        ]
    )

    return T

def p_us_std_atm_z(z):
    h = (r0 * z) / (r0 + z) # Convert to Geopotential Height [m]
    p = np.piecewise(
        h,
        [
            (h >= H_b[0])  & (h <  H_b[1]),
            (h >= H_b[1])  & (h <  H_b[2]),
            (h >= H_b[2])  & (h <  H_b[3]),
            (h >= H_b[3])  & (h <  H_b[4]),
            (h >= H_b[4])  & (h <  H_b[5]),
            (h >= H_b[5])  & (h <  H_b[6]),
            (h >= H_b[6])  & (h <= H_b[7]),
        ],
        [
            lambda h: p_b[0] * (T_Mb[0] / (T_Mb[0] + L_Mb[0] * (h - H_b[0])))**((g0 * M0) / (R_star * L_Mb[0])), # 0  - 11 km'
            lambda h: p_b[1] * np.exp(-((g0 * M0 * (h - H_b[1])) / (R_star * T_Mb[1]))),                         # 11 - 20 km'
            lambda h: p_b[2] * (T_Mb[2] / (T_Mb[2] + L_Mb[2] * (h - H_b[2])))**((g0 * M0) / (R_star * L_Mb[2])), # 20 - 32 km'
            lambda h: p_b[3] * (T_Mb[3] / (T_Mb[3] + L_Mb[3] * (h - H_b[3])))**((g0 * M0) / (R_star * L_Mb[3])), # 32 - 47 km'
            lambda h: p_b[4] * np.exp(-((g0 * M0 * (h - H_b[4])) / (R_star * T_Mb[4]))),                         # 47 - 51 km'
            lambda h: p_b[5] * (T_Mb[5] / (T_Mb[5] + L_Mb[5] * (h - H_b[5])))**((g0 * M0) / (R_star * L_Mb[5])), # 51 - 71 km'
            lambda h: p_b[6] * (T_Mb[6] / (T_Mb[6] + L_Mb[6] * (h - H_b[6])))**((g0 * M0) / (R_star * L_Mb[6])), # 71 - 84.8520 km'
        ]
    )

    return p

def z_us_std_atm_p(p):
    h = np.piecewise(
        p,
        [
            (p <= p_b[0])  & (p >  p_b[1]),
            (p <= p_b[1])  & (p >  p_b[2]),
            (p <= p_b[2])  & (p >  p_b[3]),
            (p <= p_b[3])  & (p >  p_b[4]),
            (p <= p_b[4])  & (p >  p_b[5]),
            (p <= p_b[5])  & (p >  p_b[6]),
            (p <= p_b[6])  & (p >= p_b[7]),
            (p == 0.)
        ],
        [
            lambda p: H_b[0] + (T_Mb[0] / L_Mb[0]) * (np.pow((p_b[0] / p), (R_star * L_Mb[0]) / (g0 * M0)) - 1.), # 0  - 11 km'
            lambda p: H_b[1] - ((R_star * T_Mb[1]) / (g0 * M0)) * np.log(p / p_b[1]),                      # 11 - 20 km'
            lambda p: H_b[2] + (T_Mb[2] / L_Mb[2]) * (np.pow((p_b[2] / p), (R_star * L_Mb[2]) / (g0 * M0)) - 1.), # 20 - 32 km'
            lambda p: H_b[3] + (T_Mb[3] / L_Mb[3]) * (np.pow((p_b[3] / p), (R_star * L_Mb[3]) / (g0 * M0)) - 1.), # 32 - 47 km'
            lambda p: H_b[4] - ((R_star * T_Mb[4]) / (g0 * M0)) * np.log(p / p_b[4]),                      # 47 - 51 km'
            lambda p: H_b[5] + (T_Mb[5] / L_Mb[5]) * (np.pow((p_b[5] / p), (R_star * L_Mb[5]) / (g0 * M0)) - 1.), # 51 - 71 km'
            lambda p: H_b[6] + (T_Mb[6] / L_Mb[6]) * (np.pow((p_b[6] / p), (R_star * L_Mb[6]) / (g0 * M0)) - 1.), # 71 - 84.8520 km'
            lambda p: np.inf
        ]
    )

    z = (r0 * h) / (r0 - h) # Geometric height [m]

    return z

def T_us_std_atm_p(p):
    z = z_us_std_atm_p(p)
    T = T_us_std_atm_z(z)

    return T

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plotting-dir", nargs = "?", default = ".plotting", type = str,
        help = "Directory to save plots.")
    parser.add_argument("--tag", nargs = "?", default = "", type = str,
        help = "Dataset tag.")
    args = parser.parse_args()

    plotting_dir = args.plotting_dir
    tag = args.tag

    dirs = [plotting_dir]
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)

    #---------------------------------------------------------------------------
    # Temperature U.S. Standard Atmosphere (1976) as a Function of Altitude
    #---------------------------------------------------------------------------
    fileroot = "T_us_std_atm_z"

    n_z = 2048
    z = np.linspace(Z_b[0], Z_b[-1], n_z) # Geometric Height [m]
    T = T_us_std_atm_z(z) # Temperature [K]

    fig, axs = plt.subplots(sharex = True)

    axs.plot(T, z * 0.001, color = "k")
    axs.axhline([Z_b[-1] * 0.001], color = "grey")

    fig.supxlabel(r"Temperature $\left[ K \right]$")
    fig.supylabel(r"Geometric Height $\left[ km \right]$")
    fig.suptitle("U.S. Standard Atmosphere (1976)")

    plt_filepath = os.path.join(plotting_dir, fileroot + ".png")
    plt.savefig(plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

    #---------------------------------------------------------------------------
    # Pressure U.S. Standard Atmosphere (1976) as a Function of Altitude
    #---------------------------------------------------------------------------
    fileroot = "p_us_std_atm_z"

    n_z = 2048
    z = np.linspace(Z_b[0], Z_b[-1], n_z) # Geometric Height [m]
    p = p_us_std_atm_z(z) # Pressure [hPa]

    fig, axs = plt.subplots(sharex = True)

    axs.plot(p, z * 0.001, color = "k")
    axs.axhline([Z_b[-1] * 0.001], color = "grey")

    axs.set_xscale("log")

    fig.supxlabel(r"Pressure $\left[ hPa \right]$")
    fig.supylabel(r"Geometric Height $\left[ km \right]$")
    fig.suptitle("U.S. Standard Atmosphere (1976)")

    plt_filepath = os.path.join(plotting_dir, fileroot + ".png")
    plt.savefig(plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

    #---------------------------------------------------------------------------
    # Geometric Height U.S. Standard Atmosphere (1976) as a Function of Pressure
    #---------------------------------------------------------------------------
    fileroot = "z_us_std_atm_p"

    n_p = 2048
    p = np.linspace(p_b[0], p_b[-1], n_p) # Pressure [hPa]
    z = z_us_std_atm_p(p) # Geometric Height [m]

    fig, axs = plt.subplots(sharex = True)

    axs.plot(p, z * 0.001, color = "k")
    axs.axhline([Z_b[-1] * 0.001], color = "grey")

    axs.set_xscale("log")

    fig.supxlabel(r"Pressure $\left[ hPa \right]$")
    fig.supylabel(r"Geometric Height $\left[ km \right]$")
    fig.suptitle("U.S. Standard Atmosphere (1976)")

    plt_filepath = os.path.join(plotting_dir, fileroot + ".png")
    plt.savefig(plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

    #---------------------------------------------------------------------------
    # Temperature U.S. Standard Atmosphere (1976) as a Function of Pressure
    #---------------------------------------------------------------------------
    fileroot = "T_us_std_atm_p"

    n_p = 2048
    p = np.linspace(p_b[0], p_b[-1], n_p) # Pressure [hPa]
    T = T_us_std_atm_p(p) # Temperature [K]

    fig, axs = plt.subplots(sharex = True)

    axs.plot(T, p, color = "k")
    axs.axhline([p_b[-1]], color = "grey")

    axs.set_yscale("log")
    axs.yaxis.set_inverted(True)

    fig.supxlabel(r"Temperature $\left[ K \right]$")
    fig.supylabel(r"Pressure $\left[ hPa \right]$")
    fig.suptitle("U.S. Standard Atmosphere (1976)")

    plt_filepath = os.path.join(plotting_dir, fileroot + ".png")
    plt.savefig(plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

if __name__ == "__main__":
    main()