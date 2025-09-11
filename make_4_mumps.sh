#!/bin/bash
set -e

# Default values
: ${MUMPS_VERSION:=5.6.2}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
wget -c -O MUMPS_${MUMPS_VERSION}.tar.gz https://coin-or-tools.github.io/ThirdParty-Mumps/MUMPS_${MUMPS_VERSION}.tar.gz
tar -xzf MUMPS_${MUMPS_VERSION}.tar.gz
cd MUMPS_${MUMPS_VERSION}

# Build MUMPS static lib
cp ../../Makefile_mumps.inc Makefile.inc
make d prefix=${PREFIX} openmp=1
make z prefix=${PREFIX} openmp=1

# Test
cd examples
echo " "
echo "running tests ..."
echo " "
cp ../../../check_test .
./dsimpletest < input_simpletest_real | ./check_test
./zsimpletest < input_simpletest_cmplx | ./check_test
cd ..

