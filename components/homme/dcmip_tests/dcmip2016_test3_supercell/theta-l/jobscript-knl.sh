#!/bin/bash 
#
#   Jobscript for launching dcmip2016 test3 on Cori KNL
#
#SBATCH --job-name=dcmip-2016-3
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=36
#SBATCH --constraint=knl

#SBATCH --account=e3sm
#SBATCH --mail-user=jltorch@sandia.gov
#SBATCH --mail-type=all

OMP_NUM_THREADS=1
NCPU=8
if [ -n "$SLURM_NNODES" ]; then
   patt='([[:digit:]]+)'
   if [[ $SLURM_TASKS_PER_NODE =~ $patt ]]; then
     PER_NODE=${BASH_REMATCH[1]}
   else
     PER_NODE=16
   fi
   NCPU=$(( $SLURM_NNODES * $PER_NODE ))
fi
EXEC=../../../test_execs/theta-l-nlev40/theta-l-nlev40



function run { 
local NCPU=$1
echo "NCPU = $NCPU"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
date
srun -K -c 1 -n $NCPU -N $SLURM_NNODES  $EXEC < input.nl
date

ncl plot_supercell_wvel.ncl
ncl plot_supercell_2.5km_wvel_xsec.ncl
ncl plot_supercell_5km_xsec.ncl
ncl plot_supercell_prect.ncl

\mv movies/dcmip2016_test31.nc movies/dcmip2016_test3_${prefix}.nc
\mv HommeTime                  HommeTime_${prefix}
\mv max_w.pdf                  max_w_${prefix}.pdf
\mv max_precip.pdf             max_precip_${prefix}.pdf
\mv 2.5km_wvel_xsec.pdf        2.5km_wvel_xsec_${prefix}.pdf
\mv 5km_xsec.pdf               5km_xsec_${prefix}.pdf
\mv measurement_wmax.txt       measurement_wmax_${prefix}.txt
\mv measurement_time.txt       measurement_time_${prefix}.txt
\mv measurement_prect_rate.txt measurement_prect_rate_${prefix}.txt

}

prefix=r400             ; run $(($NCPU>384?384:NCPU))
#prefix=r200             ; run $NCPU
#prefix=r100             ; run $NCPU
#prefix=r50              ; run $NCPU

#prefix=explicit-r100    ; run $NCPU

