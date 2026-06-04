#!/bin/bash 
#
#   Jobscript for creating plots on Flight, an SRN Capacity Cluster
#

#SBATCH --time=0-02:30:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --account=fy260144
#SBATCH --reservation=flight-cldera
#SBATCH --job-name=pk02-plotting
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# NOTE: This will run the plotting scripts created for the SCREAM-STRAT project
# It will have to be adjusted for individual purposes.

# This copy of the script is meant for Flight, an SRN Capacity Cluster.

HS_DIR=/projects/scream_strat/jltorch/held_suarez
TAG=r400-nlev30

if false; then
    #---------------------------------------------------------------------------
    CURRENT_TIME=$(date +"%T")
    STAGE_NAME="HELD-SUAREZ CLIMATOLOGY VISUALIZATION"
    echo "[${CURRENT_TIME}]: Beginning ${STAGE_NAME}..."

    python plot_hs94_climatologies.py \
        --spinup-days 365 \
        --homme-output ${HS_DIR}/${TAG}-held_suarez.nc \
        --plot-vars u,T,T_eddy,pnh \
        --tag ${TAG} \
        --recalculate true

    CURRENT_TIME=$(date +"%T")
    echo "[${CURRENT_TIME}]: ${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi

PK_DIR=/pscratch/jltorch/scream-strat/e3sm/jasonltorchinsky/held-suarez/homme/test/held_suarez/movies
TAG=ne30-nlev72

if true; then
    #---------------------------------------------------------------------------
    CURRENT_TIME=$(date +"%T")
    STAGE_NAME="POLVANI-KUSHNER CLIMATOLOGY VISUALIZATION"
    echo "[${CURRENT_TIME}]: Beginning ${STAGE_NAME}..."

    srun --mpi=pmi2 -n ${SLURM_NTASKS} --cpu-bind=cores \
      python plot_pk02_climatologies.py \
        --spinup-days 200 \
        --homme-output ${PK_DIR}/${TAG}-held_suarez.nc \
        --plot-vars u \
        --tag ${TAG} \
        --working-dir .polvani_kushner/${TAG} \
        --plotting-dir polvani_kushner/${TAG} \
        --recalculate true

    CURRENT_TIME=$(date +"%T")
    echo "[${CURRENT_TIME}]: ${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi

VCOORD_DIR=/ascldap/users/jltorch/codes/e3sm/jasonltorchinsky/held-suarez/components/homme/test/vcoord
VCOORD_NAME=acme
I_FILE=${VCOORD_DIR}/acme-72i.ascii
M_FILE=${VCOORD_DIR}/acme-72m.ascii

if false; then
    #---------------------------------------------------------------------------
    CURRENT_TIME=$(date +"%T")
    STAGE_NAME="VCOORD VISUALIZATION"
    echo "[${CURRENT_TIME}]: Beginning ${STAGE_NAME}..."
    
    python plot_vcoord.py \
        --interface-file ${I_FILE} \
        --midpoint-file ${M_FILE} \
        --vcoord-name ${VCOORD_NAME} \
        --plotting-dir vcoord
    
    CURRENT_TIME=$(date +"%T")
    echo "[${CURRENT_TIME}]: ${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi