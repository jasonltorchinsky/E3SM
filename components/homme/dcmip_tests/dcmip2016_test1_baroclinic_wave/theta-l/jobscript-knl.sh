#!/bin/bash 
#
#   Jobscript for launching dcmip2016 test1 on Cori KNL
#
#SBATCH --job-name=dcmip2016-1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=6
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=all

#SBATCH --time=0-00:30:00

# Set variables, load modules
echo '-- Loading modules...'

export OMP_STACKSIZE=16M # Cori has 96GB/node, had to lower to 8M on 3K nodes.
export OMP_NUM_THREADS=2
PER_NODE=64              # MPI per node
VC_PER_MPI=256 # Num. of virtual cores per MPI task
               # Set to 272 if PER_NODE divides 272 instead of 256
let VC_PER_MPI/=$PER_NODE

export KMP_AFFINITY="granularity=core,scatter"
bind=--cpu_bind=core

# compute number of MPI tasks
if [ -n "$SLURM_NNODES" ]; then
  NNODES=$SLURM_NNODES
else
  NNODES=1
fi


NMPI=$NNODES
let NMPI*=$PER_NODE


# Hydrostatic theta-l
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30

module load ncl

echo '-- Loaded modules.'

function run { 
local NMPI=$1

echo NODES =            $NNODES
echo NMPI_PER_NODE =    $PER_NODE
echo NTHREADS_PER_MPI = $OMP_NUM_THREADS
mpirun="srun -n $NMPI -N $NNODES -c $VC_PER_MPI $bind"
echo mpi commnand:
echo $mpirun

namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
$mpirun  $EXEC < input.nl
date

ncl plot-baroclinicwave-init.ncl
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=1'
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=2'
ncl plot-lat-lon-TPLSPS.ncl 'var_choice=3'
ncl plot-all.ncl
\mv -f plot_baroclinicwave_init.pdf  ${prefix}_init.pdf
\mv -f preqx-test16-1latlonT850.pdf  ${prefix}_T850.pdf
\mv -f preqx-test16-1latlonPS.pdf  ${prefix}_PS.pdf
\mv -f preqx-test16-1latlonPRECL.pdf  ${prefix}_PRECL.pdf
\mv -f plot_all.pdf  ${prefix}_all.pdf

\mv -f movies/dcmip2016_test11.nc    movies/${prefix}_dcmip2016_test11.nc
}

prefix=r400    ; run $(($NMPI>384?384:NMPI))


