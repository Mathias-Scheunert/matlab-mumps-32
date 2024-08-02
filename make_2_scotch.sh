#!/bin/bash
set -e

# Default values
: ${SCOTCH_VERSION:=v7.0.3}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
git clone https://gitlab.inria.fr/scotch/scotch.git
cd scotch
git checkout tags/${SCOTCH_VERSION}
cd src

# Build SCOTCH static lib (DEBIAN system!)
# -> the linux packages bison and flex are required! (sudo apt-get install bison flex)
cp ./Make.inc/Makefile.inc.i686_pc_linux2 Makefile.inc
make esmumps

# Install i.e. copy content to desired paths
make install prefix=${PREFIX}
