#!/bin/bash 
#
#   Jobscript for launching dcmip2016 test2 on Cori KNL
#
#SBATCH --job-name=dcmip-2016-2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=32
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jason.torchinsky@wisc.edu
#SBATCH --mail-type=all

#SBATCH --time=0-00:30:00

# Set variables, load modules
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
module load python
conda activate geocat
module load ncl

function run {
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
srun -K -c 4 -N $SLURM_NNODES  $EXEC < input.nl
date

ncl plot-tropical-cyclone-init.ncl  # u,t,th,q,pnh,geo,ps, time=0
ncl plot-horiz-crossx.ncl     # contour plot, time=10d, U,V,T,ps,precl,Q,geo
#ncl plot-intensity-trace.ncl
python plot-ps.py
python plot-intensity.py
python plot-precip.py
python plot-Ts.py
python plot-uvs.py
#ncl plot-horiz-ps.ncl

# save output
\mv -f movies/dcmip2016_test21.nc   movies/${prefix}_dcmip2016_test21.nc

\mv -f init.pdf ${prefix}_init.pdf
\mv -f x-sections.pdf ${prefix}_x-sections.pdf
#\mv -f wind.pdf ${prefix}_wind.pdf
#\mv -f psmap.pdf ${prefix}_psmap.pdf
\mv -f ps.pdf ${prefix}_ps.pdf
\mv -f intensity.pdf ${prefix}_intensity.pdf
\mv -f precip.pdf ${prefix}_precip.pdf
\mv -f Ts.pdf ${prefix}_Ts.pdf
\mv -f uvs.pdf ${prefix}_uvs.pdf
}

prefix=r100;  run

