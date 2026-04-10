#!/bin/bash 
#
#   Jobscript for launching DCMIP 2016 Test 2 - Tropical Cyclone on Flight,
#   an SRN Capacity Cluster
#

#SBATCH --job-name=dcmip2016-2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=normal
#SBATCH --nodes=4

#SBATCH --time=0-00:30:00

### This is currently set up for a specific account - jltorch. MUST be changed per-user.
#SBATCH --account=fy260144
#SBATCH --reservation=flight-cldera
#SBATCH --partition=short

### We assume a source code directory structure of ${SRC_ROOT}/e3sm/${BRANCH_NAME},
### which we will recreate in the "work space".
SRC_ROOT=${HOME}/codes
BRANCH_NAME=jasonltorchinsky/flight

SRC_E3SM=${SRC_ROOT}/e3sm/${BRANCH_NAME}

### Load the necessary modules
echo "-- WARNING: Ensure that necessary modules are lodaded by running"
echo "            'source ${SRC_E3SM}/cime/CIME/Tools/get_case_env'"
echo "            We will run this now, but that is only local to this script."
eval `${SRC_E3SM}/cime/CIME/Tools/get_case_env`
module load climate-utils/2023q2 # For NCL

### Process distribution
OMP_NUM_THREADS=1
CORES_PER_NODE=56  # Cores per node on Flight
NCORES=$(( $SLURM_NNODES * $CORES_PER_NODE )) # Physical Cores
CORES_PER_TASK=2 # Flight has hyperthreading enabled, so use 2
NVCORES=$(( CORES_PER_TASK * NCORES )) # Virtual Cores

### Executable location, module for NCL
WORK_HOMME=../../.. # Relative path
TEST_EXECS=${WORK_HOMME}/test_execs
EXEC=${TEST_EXECS}/theta-l-nlev30/theta-l-nlev30 # Default executable

# Function to run the test
function run {
    local NTASKS=$1
    echo "NTASKS = $NTASKS"
    namelist=namelist-$prefix.nl
    \cp -f $namelist input.nl
    srun -K -c $CORES_PER_TASK -n $NTASKS -N $SLURM_NNODES $EXEC < input.nl
    date

    # Plot initial conditions for: u, T, th, q, p_nh, geo, and ps
    #ncl plot-tropical-cyclone-init.ncl
    # Plot horizontal cross-section an 10 days for: u, v, T, Q, precl, geo, and ps
    #ncl plot-horiz-crossx.ncl
    # Plot surface pressure and lowest level horizontal wind magnitude over time.
    #ncl plot-intensity-trace.ncl
    # Plot surface pressure map at days 0, 8, 9, and 10.
    #ncl plot-horiz-ps.ncl

    # Save output to run-specific files
    #\mv -f movies/dcmip2016_test21.nc   movies/${prefix}_dcmip2016_test21.nc

    #\mv -f init.pdf ${prefix}_init.pdf
    #\mv -f x-sections.pdf ${prefix}_x-sections.pdf
    #\mv -f wind.pdf ${prefix}_wind.pdf
    #\mv -f psmap.pdf ${prefix}_psmap.pdf
}

# Max NTASKS is ne*ne*6, with ne specified in the namelist
MAX_NTASKS=$(( 8 * 8 * 6 ))
prefix=r400  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

MAX_NTASKS=$(( 30 * 30 * 6 ))
prefix=r100  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

MAX_NTASKS=$(( 60 * 60 * 6 ))
prefix=r50   ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))
prefix=r50-h ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

