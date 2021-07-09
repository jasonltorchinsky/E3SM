#!/bin/bash 
#
#   Jobscript for launching dcmip2012 test3-1 on Cori KNL
#
#SBATCH --job-name=d31-theta
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=45
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=all

#SBATCH --time=0-00:30:00

# Set variables, load modules
OMP_NUM_THREADS=1
EXEC=../../../test_execs/theta-l-nlev20/theta-l-nlev20
module load ncl

echo '-- Loaded modules. Running hydrostatic version...'

#############################################################################
# theta (hydrostatic)
#############################################################################
#namelist=./namelist-h.nl
#cp -f $namelist input.nl
#srun -K -c 4 -N $SLURM_NNODES  $EXEC < input.nl

#echo '-- Completed hydrostatic simulation. Creating hydrostatic plots...'

#ncl plot_omega.ncl
#ncl plot_theta.ncl

#mv -f test31_omega.pdf                 hydro_test31_omega.pdf
#mv -f dcmip2012_test3_theta_diff.pdf   hydro_test3_theta_diff.pdf
#mv -f dcmip2012_test3_theta_diff_last.pdf   hydro_test3_theta_diff_last.pdf
#mv -f movies/dcmip2012_test31.nc       movies/hydro_dcmip2012_test31.nc

#echo '-- Completed hydrostatic plots. Running non-hydrostatic version...'

#############################################################################
# theta-nh (non-hydrostatic)
#############################################################################
namelist=./namelist-nh.nl
cp -f $namelist input.nl
srun -K -c 4 -N $SLURM_NNODES  $EXEC < input.nl

echo '-- Completed non-hydrostatic simulation. Creating non-hydrostatic plots...'

ncl plot_omega.ncl
ncl plot_theta.ncl

mv -f test31_omega.pdf                 nonhydro_test31_omega.pdf
mv -f dcmip2012_test3_theta_diff.pdf   nonhydro_test3_theta_diff.pdf  
mv -f dcmip2012_test3_theta_diff_last.pdf   nonhydro_test3_theta_diff_last.pdf  
mv -f movies/dcmip2012_test31.nc        movies/nonhydro_dcmip2012_test31.nc 

echo '-- Completed non-hydrostatic plots. Job complete!'
