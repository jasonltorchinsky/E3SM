#!/bin/bash
# Builds test Standalone HOMME on the SRN Capacity Cluster 'Flight'.
# Must be adjusted for each user.
# Currently set up for account jltorch on Flight.
# Author: Jason Torchinsky

# Parse command-line input
while getopts t: flag
do
	case "${flag}" in
		t) TEST_ID=${OPTARG};;
	esac
done

### Store the current directory to return to at script conclusion
PWD=$(pwd)

### We assume a source code directory structure of ${SRC_ROOT}/e3sm/${BRANCH_NAME},
### which we will recreate in the "work space".
SRC_ROOT=${HOME}/codes
BRANCH_NAME=jasonltorchinsky/flight

SRC_E3SM=${SRC_ROOT}/e3sm/${BRANCH_NAME}
SRC_HOMME=${SRC_E3SM}/components/homme
SRC_TEST=${SRC_HOMME}/test
SRC_VCOORD=${SRC_TEST}/vcoord

### The "work space" will contain the build files, executables, etc.
WORK_ROOT=/pscratch/jltorch/scream-strat
WORK_E3SM=${WORK_ROOT}/e3sm/${BRANCH_NAME}
WORK_HOMME=${WORK_E3SM}/homme
WORK_TEST=${WORK_HOMME}/test
WORK_VCOORD=${WORK_TEST}/vcoord

### Ensure the work space directories exist
mkdir -p ${WORK_ROOT}
mkdir -p ${WORK_E3SM}
mkdir -p ${WORK_HOMME}
mkdir -p ${WORK_VCOORD}
mkdir -p ${WORK_TEST}

### Copy vertical coordinate files (vcoord) to the work space
cp -r ${SRC_VCOORD}/* ${WORK_VCOORD}

### The machine file to use for compilation - MUST be changed for each machine
MACHINE_FILE=${SRC_HOMME}/cmake/machineFiles/flight-intel.cmake

### Load the necessary modules
echo "-- WARNING: Ensure that necessary modules are loaded by running"
echo "            'source ${SRC_E3SM}/cime/CIME/Tools/get_case_env'"
echo "            We will run this now, but that is only local to this script."
eval `${SRC_E3SM}/cime/CIME/Tools/get_case_env`

### Set some common variables to shorten code later
CMAKE_FLAGS=""

### Select test to compile
case "${TEST_ID}" in
	held_suarez0)
		TEST_NAME="Held-Suarez"
		MODE=theta-l
		SRC_TEST_DIR=${SRC_TEST}/held_suarez0
		WORK_TEST_DIR=${WORK_TEST}/held_suarez0
		;;
	*)
		echo "-- Unable to parse test ID, or test ID is unsupported. Aborting..."
		exit 2
		;;
esac

### Set CMake flags
if [[ "${MODE}" == "sweqx" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_SWEQX=TRUE"
elif [[ "${MODE}" == "preqx" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_PREQX=TRUE"
elif [[ "${MODE}" == "theta-l" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_THETA=TRUE"
elif [[ "${MODE}" == "preqx_acc" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_PREQX_ACC=TRUE"
elif [[ "${MODE}" == "theta-l_kokkos" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_THETA_KOKKOS=TRUE"
elif [[ "${MODE}" == "swim" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_SWIM=TRUE"
fi

echo ${CMAKE_FLAGS}

### Build Standalone HOMME
echo "-- Building HOMME..."

#cd ${WORK_HOMME}
#cmake -C ${MACHINE_FILE} ${SRC_HOMME} ${CMAKE_FLAGS}
#make -j ${MODE}

echo "-- HOMME built!"

### Copy namelists, etc. to the work space
echo "-- Copying ${TEST_NAME} files to work space..."

mkdir -p ${WORK_TEST_DIR}
cp ${SRC_TEST_DIR}/*.nl ${WORK_TEST_DIR}
cp ${SRC_TEST_DIR}/*.sh ${WORK_TEST_DIR}

echo "-- Copied ${TEST_NAME} files to work space!"

### Return to directory this script was called from
cd ${PWD}
