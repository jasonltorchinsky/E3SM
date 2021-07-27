#!/bin/bash

echo '-- Setting file paths...'

e3sm=/global/homes/j/jlturner/practicum-2021/E3SM-height+remap

export e3sm

# Load the necessary modules

echo '-- Loading necessary modules...'

eval `$e3sm/cime/scripts/Tools/get_case_env`

# Build the executable
cwd=`pwd`
cd ../../..
  echo "make -j theta-l-nlev30"
  make -j theta-l-nlev30
cd $cwd

# Remove previous output files
echo '-- Removing previous output files...'
rm *.out *.err *.pdf
