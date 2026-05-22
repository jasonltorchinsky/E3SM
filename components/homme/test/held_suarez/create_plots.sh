#!/bin/bash 
# NOTE: This will run the plotting scripts created for the SCREAM-STRAT project
# It will have to be adjusted for individual purposes.

HS_DIR=/projects/scream_strat/jltorch/held_suarez
TAG=r400-nlev30

if false; then
    #---------------------------------------------------------------------------
    STAGE_NAME="HELD-SUAREZ CLIMATOLOGY VISUALIZATION"
    echo "Beginning ${STAGE_NAME}..."

    python plot_hs94_climatologies.py \
        --spinup-days 365 \
        --homme-output ${HS_DIR}/${TAG}-held_suarez.nc \
        --plot-vars u,T,T_eddy,pnh \
        --tag ${TAG} \
        --recalculate true

    echo "${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi

PK_DIR=/projects/scream_strat/jltorch/polvani_kushner
TAG=ne8-nlev72

if true; then
    #---------------------------------------------------------------------------
    STAGE_NAME="POLVANI-KUSHNER CLIMATOLOGY VISUALIZATION"
    echo "Beginning ${STAGE_NAME}..."

    python plot_pk02_climatologies.py \
        --spinup-days 200 \
        --homme-output ${PK_DIR}/${TAG}-held_suarez.nc \
        --plot-vars u,T \
        --tag ${TAG} \
        --working-dir .polvani_kushner/${TAG} \
        --plotting-dir polvani_kushner/${TAG} \
        --recalculate false

    echo "${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi

VCOORD_DIR=/ascldap/users/jltorch/codes/e3sm/jasonltorchinsky/held-suarez/components/homme/test/vcoord
VCOORD_NAME=cam
I_FILE=${VCOORD_DIR}/cami-30.ascii
M_FILE=${VCOORD_DIR}/camm-30.ascii

if false; then
    #---------------------------------------------------------------------------
    STAGE_NAME="VCOORD VISUALIZATION"
    echo "Beginning ${STAGE_NAME}..."
    
    python plot_vcoord.py \
        --interface-file ${I_FILE} \
        --midpoint-file ${M_FILE} \
        --vcoord-name ${VCOORD_NAME} \
        --plotting-dir vcoord
    
    echo "${STAGE_NAME} complete!"
    #---------------------------------------------------------------------------
fi