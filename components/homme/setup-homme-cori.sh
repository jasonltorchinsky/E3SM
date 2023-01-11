#!/bin/bash
# Compiles standalone HOMME for various (but not all) test cases.
# Must be adjusted for each user, and depends on machine.
# Currently set up for account jlturner on Cori.
# Author: Jason Torchinsky

# Parse command-line input
while getopts t: flag
do
	case "${flag}" in
		t) test=${OPTARG};;
	esac
done

# [Optional] Declare development branch for E3SM
# Useful for developing multiple branches simultaneously
branch=jasonltorchinsky/dcmip2016-physics-p

# Declare directory paths
e3sm=$(pwd)/../../
homme=$(pwd)
scratch=/global/cscratch1/sd/jlturner/E3SM
wdir=$scratch/$branch/homme
mach=$homme/cmake/machineFiles/cori-knl.cmake

case "${test}" in
	12_1_1) 
		testname='DCMIP 2012 Test 1.1 - 3D Deformational Flow'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.1_3d_deformational_flow/preqx
		mode='preqx'
		;;
	12_1_2) 
		testname='DCMIP 2012 Test 1.2 - Hadley-Like Meridional Circulation'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.2_hadley_meridional_circulation/preqx
		mode='preqx'
		;;
	12_1_3) 
		testname='DCMIP 2012 Test 1.3 - Thin Clouds Over Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test1.3_thin_clouds_over_orography/preqx
		mode='preqx'
		;;
	12_2_0) 
		testname='DCMIP 2012 Test 2.0 - Steady State with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.0_steady_state_with_orography/theta-l
		mode='theta-l'
		;;
	12_2_1) 
		testname='DCMIP 2012 Test 2.1 - Non-Sheared Background Flow with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.1_nh_mountain_waves_no_shear/theta-l
		mode='theta-l'
		;;
	12_2_2) 
		testname='DCMIP 2012 Test 2.2 - Sheared Background Flow with Orography'
		testdir=$wdir/dcmip_tests/dcmip2012_test2.2_nh_mountain_waves_with_shear/theta-l
		mode='theta-l'
		;;
	12_3) 
		testname='DCMIP 2012 Test 3 - Non-Hydrostatic Gravity Waves'
		testdir=$wdir/dcmip_tests/dcmip2012_test3.1_nh_gravity_waves/theta-l
		mode='theta-l'
		;;
	12_4_1) 
		testname='DCMIP 2012 Test 4.1 - Baroclinic Instability'
		testdir=$wdir/dcmip_tests/dcmip2012_test4.1_baroclinic_instability/theta-l
		mode='theta-l'
		;;
	16_1) 
		testname='DCMIP 2016 Test 1 - Moist Baroclinic Wave'
		testdir=$wdir/dcmip_tests/dcmip2016_test1_baroclinic_wave/theta-l
		mode='theta-l'
		;;
	16_2) 
		testname='DCMIP 2016 Test 2 - Tropical Cyclone'
		testdir=$wdir/dcmip_tests/dcmip2016_test2_tropical_cyclone/theta-l
		mode='theta-l'
		;;
	16_3) 
		testname='DCMIP 2016 Test 3 - Supercell'
		testdir=$wdir/dcmip_tests/dcmip2016_test3_supercell/theta-l
		mode='theta-l'
		;;
	mrb)	
		testname='Moist Rising Bubble'
		testdir=$wdir/test_execs/theta-l-nlev20-native
		mode='theta-l'
		;;
	*)
		echo "-- Unable to parse test number, or test number is unsupported. Aborting..."
                exit 2
		;;
esac

# Ensure directory paths exists
mkdir -p $wdir
mkdir -p $testdir

echo "-- Setting up test: $testname"

export e3sm
export branch
export homme
export scratch
export wdir
export mach
export testdir
export mode

# Load the necessary modules

echo "-- Loading necessary modules..."

eval `$e3sm/cime/CIME/Tools/get_case_env`
# Create GeoCAT conda environment if it doesn't exist
module load python
find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

if find_in_conda_env ".*geocat.*" ; then
	echo '-- Found geocat conda environment.'
	conda activate geocat
else
	echo '-- Failed to find geocat conda environment,  installing...'
	conda create -n geocat -c conda-forge geocat-comp geocat-viz matplotlib cartopy jupyter xarray cmaps -y
	conda activate geocat
fi

# Compile HOMME
echo "-- Compiling HOMME..."

cd $wdir
cmake -C $mach $homme
make -j4 $mode

# Compile Test

cd $testdir
make install
./build.sh

