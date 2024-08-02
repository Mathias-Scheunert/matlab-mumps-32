#!/bin/bash
set -e

# Default values
: ${MUMPS_VERSION:=5.6.2}
: ${PREFIX:="${PWD}"}

# Build mex-files
cd src/MUMPS_${MUMPS_VERSION}/MATLAB
cp Makefile Makefile.bak
cp ../../../Makefile_matlab ./Makefile
make

# Install
mkdir -p "${PREFIX}/lib/matlab"
cp *.m *.mex* *.mat "${PREFIX}/lib/matlab"

# Test MATLAB
echo " "
echo "running tests ..."
echo " "
cd /tmp
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.5
MATLABPATH="${PREFIX}/lib/matlab:${MATLABPATH}" matlab -nojvm -batch 'simple_example; zsimple_example; schur_example; diagainv_example; multiplerhs_example; sparserhs_example; polyfit(1:10, sin(1:10), 2)'
