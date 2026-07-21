# library imports
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

from us_standard_atmosphere_1976 import temperature_from_pressure

# Constants
g = 9.80665 # Sea-Level Acceleration of Gravity [m s^{-2}]
R_d = 287.047 # Gas Constant for Dry Air [J kg^{-1} K^{-1}]

def main():
    #---------------------------------------------------------------------------
    # Parse command-line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("--plotting-dir", nargs = "?", default = ".plotting", type = str,
        help = "Directory to save plots.")
    parser.add_argument("--tag", nargs = "?", required = True, type = str,
        help = "Dataset tag.")
    args = parser.parse_args()

    plotting_dir = args.plotting_dir
    tag = args.tag

    assert(tag in ["hs94", "pk02", "pk02_north"])

    dir_names = [plotting_dir]
    for dir_name in dir_names:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    #---------------------------------------------------------------------------
    # Equilibrium temperature (T_eq) profile
    #---------------------------------------------------------------------------
    T_eq_fileroot = "T_eq"
    if tag:
        T_eq_fileroot += "_{}".format(tag)

    nx = 360
    ny = 512

    ps = 1013.25e2 # [Pa]

    phi_i = np.linspace(-np.pi / 2., np.pi / 2, nx + 1)
    p_i = np.logspace(np.log10(ps), np.log10(1.e2), ny + 1) # [Pa]

    phi_m = (phi_i[1:] + phi_i[:-1]) / 2.
    p_m = (p_i[1:] + p_i[:-1]) / 2.

    if tag == "hs94":
        T_eq_m = T_eq_hs94(p_m, phi_m)
    elif tag == "pk02":
        T_eq_m = T_eq_pk02(p_m, phi_m)
    elif tag == "pk02_north":
        T_eq_m = T_eq_pk02_north(p_m, phi_m)
    else:
        quit()

    fig, axs = plt.subplots(sharex = True)

    vmax = T_eq_m.max()
    vmin = T_eq_m.min()
    cmap = "plasma"
    T_eq_plt = axs.pcolormesh(np.rad2deg(phi_i), p_i * 1.e-2, T_eq_m,
        vmin = vmin, vmax = vmax, cmap = cmap)
    axs.axvline([0], color = "grey")
    axs.axhline([100], color = "grey")
    axs.axhline([1], color = "grey")

    levels = np.arange(np.ceil(vmin / 20) * 20, np.ceil(vmax / 20) * 20, 20)
    cs = axs.contour(np.rad2deg(phi_m), p_m * 1.e-2, T_eq_m, levels = levels, colors = "black",
        linewidths = 1.0)

    axs.yaxis.set_inverted(True)
    axs.set_yscale("log")

    cb = fig.colorbar(T_eq_plt, ax = axs)

    for level in levels:
        cb.ax.axhline(level, color = "black", linewidth = 1.0)

    plt.clabel(cs, inline = True)

    cb.set_label(r"Temperature $\left[ K \right]$")

    fig.supxlabel(r"Latitude $\left[ ^{\circ} \right]$")
    fig.supylabel(r"Pressure $\left[ hPa \right]$")
    suptitle = "Equilibrum Profile"
    if tag:
        suptitle += " - {}".format(tag)
    fig.suptitle(suptitle)

    T_eq_plt_filepath = os.path.join(plotting_dir, T_eq_fileroot + ".png")
    plt.savefig(T_eq_plt_filepath, dpi = 256, bbox_inches = "tight")

    plt.close()

def T_eq_hs94(p, phi):
    p0 = 1000. # [hPa]
    kappa = 2. / 7.
    dT_y = 60. # [K]
    dtheta_z = 10. # [K]

    phi_l = np.expand_dims(phi, 0)
    p_l = np.expand_dims(p, 1)

    T_eq = (315. - dT_y * np.sin(phi_l)**2 - dtheta_z * np.log(p_l / p0) * np.cos(phi_l)**2) * (p_l / p0)**kappa
    T_eq[T_eq < 200] = 200

    return T_eq

def T_eq_pk02(p, phi, gamma = 4.e-3):
    # gamma [K m^{-1}]
    p0 = 1000.e2 # [Pa]
    kappa = 2. / 7.
    dy = 60. # [K]
    dz = 10. # [K]
    epsilon = 10. # [K]
    pT = 100.e2 # [Pa]
    TT = temperature_from_pressure(pT) # [K]
    T0 = 315 # [K]
    phi_0 = np.deg2rad(-50)
    dphi = np.deg2rad(10)

    phi_l = np.expand_dims(phi, 0)
    p_l = np.expand_dims(p, 1)

    def T_eq_trop(p_arg, phi_arg):
        dT = dy * np.sin(phi_arg)**2 - epsilon * np.sin(phi_arg) + dz * np.log(p_arg / p0) * np.cos(phi_arg)**2
        T = (T0 - dT) * (p_arg / p0)**kappa
        T[T < TT] = TT

        return T

    def W(phi_arg):
        return 0.5 * (1 - np.tanh((phi_arg - phi_0) / dphi))

    def T_eq_strat(p_arg, phi_arg):
        T_PV = TT * (p_arg / pT)**(R_d * gamma / g)
        return (1. - W(phi_arg)) * temperature_from_pressure(p_arg) + W(phi_arg) * T_PV

    T_trop = T_eq_trop(p_l, phi_l)
    T_strat = T_eq_strat(p_l, phi_l)

    mask = p_l >= pT
    T_eq = np.where(mask, T_trop, T_strat)

    return T_eq

def T_eq_pk02_north(p, phi, gamma = 4.e-3):
    # gamma [K m^{-1}]
    p0 = 1000.e2 # [Pa]
    kappa = 2. / 7.
    dy = 60. # [K]
    dz = 10. # [K]
    epsilon = 10. # [K]
    pT = 100.e2 # [Pa]
    TT = temperature_from_pressure(pT) # [K]
    T0 = 315 # [K]
    phi_0 = np.deg2rad(50)
    dphi = np.deg2rad(10)

    phi_l = np.expand_dims(phi, 0)
    p_l = np.expand_dims(p, 1)

    def T_eq_trop(p_arg, phi_arg):
        dT = dy * np.sin(phi_arg)**2 + epsilon * np.sin(phi_arg) + dz * np.log(p_arg / p0) * np.cos(phi_arg)**2
        T = (T0 - dT) * (p_arg / p0)**kappa
        T[T < TT] = TT

        return T

    def W(phi_arg):
        return 0.5 * (1 - np.tanh((phi_0 - phi_arg) / dphi))

    def T_eq_strat(p_arg, phi_arg):
        T_PV = TT * (p_arg / pT)**(R_d * gamma / g)
        return (1. - W(phi_arg)) * temperature_from_pressure(p_arg) + W(phi_arg) * T_PV

    T_trop = T_eq_trop(p_l, phi_l)
    T_strat = T_eq_strat(p_l, phi_l)

    mask = p_l >= pT
    T_eq = np.where(mask, T_trop, T_strat)

    return T_eq

if __name__ == "__main__":
    main()