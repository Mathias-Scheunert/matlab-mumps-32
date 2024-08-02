#!/bin/bash
set -e

# Default values
: ${METIS_VERSION:=5.1.0}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
wget -c -O metis-${METIS_VERSION}.tar.gz https://sourceforge.net/projects/openfoam-extend/files/foam-extend-3.0/ThirdParty/metis-${METIS_VERSION}.tar.gz/download
tar -xzf metis-${METIS_VERSION}.tar.gz
cd metis-${METIS_VERSION}

# Build METIS static lib
make config prefix=${PREFIX} openmp=1 i64=0 r64=0
make

# Install i.e. copy content to desired paths
make install
