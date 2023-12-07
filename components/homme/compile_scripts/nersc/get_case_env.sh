#!/bin/bash
# Gets the case environment for standalone HOMME based on the `get_case_env`
# script of cime/CIME/Tools.
# Author: Jason Torchinsky

# Declare directory paths
# E3SM source code directory
e3sm=$(pwd)/../../../..

# Load the necessary modules
echo "-- Loading necessary modules..."
export e3sm
# Manually added modules to load based on current get_case_env script output on Perlmutter
eval `. /usr/share/lmod/8.3.1/init/sh && module unload cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf PrgEnv-gnu PrgEnv-intel PrgEnv-nvidia PrgEnv-cray PrgEnv-aocc intel intel-oneapi cudatoolkit craype-accel-nvidia80 craype-accel-host perftools-base perftools darshan && module load PrgEnv-intel/8.3.3 intel/2023.1.0 craype-accel-host craype/2.7.20 cray-mpich/8.1.25 cray-hdf5-parallel/1.12.2.3 cray-netcdf-hdf5parallel/4.9.0.3 cray-parallel-netcdf/1.12.3.3 cmake/3.24.3 && export MPICH_ENV_DISPLAY=1 && export MPICH_VERSION_DISPLAY=1 && export OMP_STACKSIZE=128M && export OMP_PROC_BIND=spread && export OMP_PLACES=threads && export HDF5_USE_FILE_LOCKING=FALSE && export PERL5LIB=/global/cfs/cdirs/e3sm/perl/lib/perl5-only-switch && export FI_CXI_RX_MATCH_MODE=software && export MPICH_COLL_SYNC=MPI_Bcast && export Albany_ROOT=/global/common/software/e3sm/mali_tpls/albany-e3sm-serial-release-gcc && export Trilinos_ROOT=/global/common/software/e3sm/mali_tpls/trilinos-e3sm-serial-release-gcc && export NETCDF_PATH=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.3/intel/19.0 && export PNETCDF_PATH=/opt/cray/pe/parallel-netcdf/1.12.3.3/intel/19.0 && export ADIOS2_ROOT=/global/cfs/cdirs/e3sm/3rdparty/adios2/2.9.1/cray-mpich-8.1.25/intel-2023.1.0 && export BLA_VENDOR=Intel10_64_dyn`
echo "-- Loaded necessary modules!"

# Unset declared variables
unset e3sm

