#!/bin/bash

echo '-- Setting file paths...'

e3sm=/global/homes/j/jlturner/practicum-2021/E3SM-heightcoord

export e3sm

# Load the necessary modules

echo '-- Loading necessary modules...'

eval `$e3sm/cime/scripts/Tools/get_case_env`

cwd=`pwd`
cd ../../..
  echo "make -j 4 theta-l-nlev20"
  make -j 4 theta-l-nlev20
cd $cwd

# Remove previous output files
echo '-- Removing previous output files...'
rm *.out *.err *.pdf
