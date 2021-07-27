#!/bin/bash 
#
#   Jobscript for launching dcmip2012 test2-1 on Cori KNL
#
#SBATCH --job-name=d21-theta
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=19
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=all

#SBATCH --time=0-00:20:00

# Set variables, load modules
export OMP_NUM_THREADS=1
# compute number of MPI tasks                                                                                               
if [ -n "$SLURM_NNODES" ]; then
  NNODES=$SLURM_NNODES
else
  NNODES=1
fi

EXEC=../../../test_execs/theta-l-nlev60/theta-l-nlev60
module load ncl


# hydrostatic theta
#namelist=namelist-h.nl
#\cp -f $namelist input.nl
#srun -c 1 -N $SLURM_NNODES  $EXEC < input.nl
#ncl plot_lon_vs_z.ncl
#\mv -f movies/dcmip2012_test2_11.nc  movies/hydro_dcmip2012_test2_11.nc
#\mv -f dcmip2012_test2_1_T_t10.pdf   hydro_test2_1_T_t10.pdf

# nonhydrostatic theta
namelist=namelist-nh.nl
\cp -f $namelist input.nl
srun -c 4 -N $SLURM_NNODES $EXEC < input.nl
ncl plot_lon_vs_z.ncl
\mv -f movies/dcmip2012_test2_11.nc  movies/nonhydro_dcmip2012_test2_11.nc
\mv -f dcmip2012_test2_1_T_t10.pdf   nonhydro_test2_1_T_t10.pdf

date
