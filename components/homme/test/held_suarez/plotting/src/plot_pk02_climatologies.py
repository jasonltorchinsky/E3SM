# Append the root directory to the PYTHONPATH for future imports
import os, sys
root_dir: str = os.path.normpath( \
    os.path.join(os.path.dirname(__file__), os.pardir))
if root_dir not in sys.path:
    sys.path.append(root_dir)
    
# Standard Library Imports
import argparse

# Third-Party Imports
import xarray as xr
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np

# Local Library Imports
from consts.dtypes import NP_INT, MPI_COMM
from consts.numeric import MPI_ROOT
from consts.physical import pk02_tropopause_pressure as pk02_p_T, pk02_sponge_pressure as pk02_p_sp

from analysis import calculate_climatology
from decomposition import decompose_time
from messages import print_msg
from remap import get_target_pressure_levels, vertical_remap_field
from visualization import plot_climatology


def main():
    #---------------------------------------------------------------------------
    # Set up MPI communicator
    #---------------------------------------------------------------------------
    comm: MPI_COMM = MPI.COMM_WORLD
    l_rank: NP_INT = NP_INT(comm.Get_rank())
    comm_size: NP_INT = NP_INT(comm.Get_size())

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
        help = "Which variable climatologies to plot: T, u, v, w")
    args = parser.parse_args()

    spinup_days: NP_INT = NP_INT(args.spinup_days)
    homme_output: str = os.path.normpath(args.homme_output)
    working_dir: str = os.path.normpath(args.working_dir)
    plotting_dir: str = os.path.normpath(args.plotting_dir)
    recalculate: bool = args.recalculate
    tag : str = args.tag
    plot_vars: list[str] = [str(plot_var) for plot_var in args.plot_vars.split(",") if plot_var]
    
    #---------------------------------------------------------------------------
    # Create working directories
    #---------------------------------------------------------------------------
    if l_rank == MPI_ROOT:
        directories: list[str] = [working_dir, plotting_dir]
        for directory in directories:
            if not os.path.exists(directory):
                os.makedirs(directory)

    #---------------------------------------------------------------------------
    # Calculate climatologies and create plots
    #---------------------------------------------------------------------------
    msg: str = "Starting climatology calculation and plotting loop."
    print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

    for plot_var in plot_vars:
        msg: str = "Starting {}.".format(plot_var)
        print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

        clim_fileroot: str = plot_var + "_clim"
        if tag:
            clim_fileroot += "_{}".format(tag)

        clim_filepath: str = os.path.join(working_dir, clim_fileroot + ".nc")

        # Calculate climatology only if necessary
        if recalculate or (not os.path.exists(clim_filepath)):
            #-------------------------------------------------------------------
            # Decompose time grid
            #-------------------------------------------------------------------
            msg: str = "Setting up time grid decomposition."
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            l_time_slice: slice = decompose_time(
                homme_output = homme_output,
                spinup_days = spinup_days,
                l_rank = l_rank,
                comm_size = comm_size)

            #-------------------------------------------------------------------
            # Set up target pressure grid
            #-------------------------------------------------------------------
            msg: str = "Creating target pressure grid."
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            p_tgt: NP_ARRAY[NP_REAL] = get_target_pressure_levels()

            #-------------------------------------------------------------------
            # Vertically remap to target pressure grid
            #-------------------------------------------------------------------
            msg: str = "Vertically remapping {}.".format(plot_var)
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            l_field_tgt: XR_DATAARRAY = vertical_remap_field(
                homme_output = homme_output,
                l_time_slice = l_time_slice,
                field_key = plot_var,
                p_tgt = p_tgt)

            #-------------------------------------------------------------------
            # Calculate climatology
            #-------------------------------------------------------------------
            msg: str = "Calculating climatology of {}.".format(plot_var)
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            g_field_clim: Optional[XR_DATAARRAY] = calculate_climatology(
                homme_output = homme_output,
                spinup_days = spinup_days,
                l_field_tgt = l_field_tgt,
                comm = comm
            )

            #-------------------------------------------------------------------
            # Save climatology to file
            #-------------------------------------------------------------------
            msg: str = "Saving climatology of {} to file.".format(plot_var)
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            if l_rank == MPI_ROOT:
                g_field_clim.to_dataset().to_netcdf(clim_filepath, mode = "w")

        if l_rank == MPI_ROOT:
            assert(os.path.exists(clim_filepath))
            with xr.open_dataset(clim_filepath, engine = "netcdf4", decode_timedelta = False) as xr_field_clim:
                field_clim: XR_DATAARRAY = (xr_field_clim[plot_var]
                    .load())

            msg: str = "Plotting {}.".format(plot_var)
            print_msg(msg, l_rank = l_rank, print_rank = MPI_ROOT)

            plot_climatology(
                field_clim = field_clim, 
                plot_fileroot = clim_fileroot,
                tag = tag,
                plotting_directory = plotting_dir)

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

if __name__ == "__main__":
    main()