#!/bin/bash
# Builds DCMIP tests for Standalone HOMME on the SRN Capacity Cluster 'Flight'.
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

### The "work space" will contain the build files, executables, etc.
WORK_ROOT=/pscratch/jltorch/scream-strat
WORK_E3SM=${WORK_ROOT}/e3sm/${BRANCH_NAME}
WORK_HOMME=${WORK_E3SM}/homme

### Ensure the work space directories exist
mkdir -p ${WORK_ROOT}
mkdir -p ${WORK_E3SM}
mkdir -p ${WORK_HOMME}

### The machine file to use for compilation - MUST be changed for each machine
MACHINE_FILE=${SRC_HOMME}/cmake/machineFiles/flight-intel.cmake

### Load the necessary modules
echo "-- WARNING: Ensure that necessary modules are loaded by running"
echo "            'source ${SRC_E3SM}/cime/CIME/Tools/get_case_env'"
echo "            We will run this now, but that is only local to this script."
eval `${SRC_E3SM}/cime/CIME/Tools/get_case_env`

### Set some common variables to shorten code later - directories here are
### created by CMake
WORK_DCMIP=${WORK_HOMME}/dcmip_tests
CMAKE_FLAGS=""

### Select test to compile
case "${TEST_ID}" in
	dcmip2012_1_1)
		TEST_NAME='DCMIP 2012 Test 1.1 - 3D Deformational Flow'
		MODE=preqx
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test1.1_3d_deformational_flow/${MODE}
		;;
	dcmip2012_1_2)
		TEST_NAME='DCMIP 2012 Test 1.2 - Hadley-Like Meridional Circulation'
		MODE=preqx
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test1.2_hadley_meridional_circulation/${MODE}
		;;
	dcmip2012_1_3)
		TEST_NAME='DCMIP 2012 Test 1.3 - Thin Clouds Over Orography'
		MODE=preqx
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test1.3_thin_clouds_over_orography/${MODE}
		;;
	dcmip2012_2_0)
		TEST_NAME='DCMIP 2012 Test 2.0 - Steady State with Orography'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test2.0_steady_state_with_orography/${MODE}
		;;
	dcmip2012_2_1)
		TEST_NAME='DCMIP 2012 Test 2.1 - Non-Sheared Background Flow with Orography'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test2.1_nh_mountain_waves_no_shear/${MODE}
		;;
	dcmip2012_2_2)
		TEST_NAME='DCMIP 2012 Test 2.2 - Sheared Background Flow with Orography'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test2.2_nh_mountain_waves_with_shear/${MODE}
		;;
	dcmip2012_3)
		TEST_NAME='DCMIP 2012 Test 3 - Non-Hydrostatic Gravity Waves'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test3.1_nh_gravity_waves/${MODE}
		;;
	dcmip2012_4_1)
		TEST_NAME='DCMIP 2012 Test 4.1 - Baroclinic Instability'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2012_test4.1_baroclinic_instability/${MODE}
		;;
	dcmip2016_1)
		TEST_NAME='DCMIP 2016 Test 1 - Moist Baroclinic Wave'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2016_test1_baroclinic_wave/${MODE}
		;;
	dcmip2016_2)
		TEST_NAME='DCMIP 2016 Test 2 - Tropical Cyclone'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2016_test2_tropical_cyclone/${MODE}
		;;
	dcmip2016_3)
		TEST_NAME='DCMIP 2016 Test 3 - Supercell'
		MODE=theta-l
		WORK_TEST_DIR=${WORK_DCMIP}/dcmip2016_test3_supercell/${MODE}
		;;
	*)
		echo "-- Unable to parse test ID, or test ID is unsupported. Aborting..."
		exit 2
		;;
esac

### Set CMake flags
if [[ "${MODE}" == "sweqx" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_SWEQX=true"
elif [[ "${MODE}" == "preqx" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_PREQX=true"
elif [[ "${MODE}" == "theta-l" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_THETA=true"
elif [[ "${MODE}" == "preqx_acc" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_PREQX_ACC=true"
elif [[ "${MODE}" == "theta-l_kokkos" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_THETA_KOKKOS=true"
elif [[ "${MODE}" == "swim" ]]; then
	CMAKE_FLAGS+=" -DBUILD_HOMME_SWIM=true"
fi

echo ${CMAKE_FLAGS}

### Build Standalone HOMME
echo "-- Building HOMME..."

cd ${WORK_HOMME}
cmake -C ${MACHINE_FILE} ${SRC_HOMME} ${CMAKE_FLAGS}
make -j ${MODE}

echo "-- HOMME built!"

### Compile Test
echo "-- Compiling ${TEST_NAME}..."

cd ${WORK_TEST_DIR}
make install
./build.sh

echo "-- ${TEST_NAME} compiled!"

### Return to directory this script was called from
cd ${PWD}
