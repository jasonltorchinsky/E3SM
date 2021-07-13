#!/bin/bash 
#
#   Jobscript for launching dcmip2012 test3-1 on Cori KNL
#
#SBATCH --job-name=d20-theta
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=45
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=all

#SBATCH --time=0-00:10:00

# Set variables, load modules
echo '-- Loding modules...'
OMP_NUM_THREADS=1
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
module load ncl

echo '-- Loaded modules.'

#############################################################################
# theta (hydrostatic)
#############################################################################
#echo '-- Running hydrostatic version...'
#namelist=./namelist-h.nl
#cp -f $namelist input.nl
#srun -K -c 4 -N $SLURM_NNODES  $EXEC < input.nl

#echo '-- Completed hydrostatic simulation. Creating hydrostatic plots...'

#ncl plot_z_lon.ncl
#ncl test200-range.ncl

#mv -f dcmip2012_test2_0_u_t6.00.pdf hydro_test2_0_u_z.pdf
#mv -f movies/dcmip2012_test2_01.nc.pdf  hydro_test2_0_u.pdf
#mv -f movies/dcmip2012_test2_01.nc  movies/hydro_dcmip2012_test2_01.nc 

#echo '-- Completed hydrostatic plots.'

#############################################################################
# theta-nh (non-hydrostatic)
#############################################################################
echo '-- Running non-hydrostatic version...'
namelist=./namelist-nh.nl
cp -f $namelist input.nl
srun -K -c 4 -N $SLURM_NNODES  $EXEC < input.nl

echo '-- Completed non-hydrostatic simulation. Creating non-hydrostatic plots...'

ncl plot_z_lon.ncl
ncl test200-range.ncl

mv -f dcmip2012_test2_0_u_t6.00.pdf nonhydro_test2_0_u_t6.00.pdf
mv -f movies/dcmip2012_test2_01.nc.pdf  nonhydro_test2_0_u.pdf
mv -f movies/dcmip2012_test2_01.nc  movies/nonhydro_dcmip2012_test2_01.nc 

echo '-- Completed non-hydrostatic plots. Job complete!'


