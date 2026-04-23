#!/bin/bash 
#
#   Jobscript for launching Held-Suarez on Flight, an SRN Capacity Cluster
#

#SBATCH --job-name=held-suarez0
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=normal
#SBATCH --nodes=8

#SBATCH --time=0-00:30:00

### This is currently set up for a specific account - jltorch. MUST be changed per-user.
#SBATCH --account=fy260144
#SBATCH --reservation=flight-cldera
#SBATCH --partition=batch

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

### Process distribution
OMP_NUM_THREADS=1
CORES_PER_NODE=56  # Cores per node on Flight
NCORES=$(( $SLURM_NNODES * $CORES_PER_NODE )) # Physical Cores
CORES_PER_TASK=2 # Flight has hyperthreading enabled, so use 2
NVCORES=$(( CORES_PER_TASK * NCORES )) # Virtual Cores

### Create output directory - SHOULD MATCH `output_dir` OF NAMELIST
OUTPUT_DIR="movies"
mkdir -p ${OUTPUT_DIR}

### Executable location, module for NCL
WORK_HOMME=../.. # Relative path
TEST_EXECS=${WORK_HOMME}/test_execs

### Function to run the test
function run {
    local NTASKS=$1
    echo "NTASKS = ${NTASKS}"
    NAMELIST=namelist-${PREFIX}.nl
    \cp -f ${NAMELIST} input.nl
    srun -K -c ${CORES_PER_TASK} -n ${NTASKS} -N ${SLURM_NNODES} ${EXEC} < input.nl
    date

    # Save output to run-specific files
    \mv -f ${OUTPUT_DIR}/held_suarez01.nc   ${OUTPUT_DIR}/${PREFIX}-held_suarez0.nc
}

### Max NTASKS is ne*ne*6, with ne specified in the namelist
MAX_NTASKS=$(( 8 * 8 * 6 ))
MODE=theta-l
NLEVS=(30 72 128)
for NLEV in ${NLEVS[@]}
do
    EXEC=${TEST_EXECS}/${MODE}-nlev${NLEV}/${MODE}-nlev${NLEV} # Executable
    PREFIX=r400-nlev${NLEV}
    run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))
done

### The follow two won't run within the time limit, so we ignore them by default.
#MAX_NTASKS=$(( 30 * 30 * 6 ))
#PREFIX=r100  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

#MAX_NTASKS=$(( 60 * 60 * 6 ))
#PREFIX=r050  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))
