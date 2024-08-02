#!/bin/bash
set -e

# Define paths
: ${OBLAS_VERSION:=v0.3.26}
: ${PREFIX:="${PWD}"}

# Download
mkdir -p src && cd src
git clone https://github.com/OpenMathLib/OpenBLAS.git
cd OpenBLAS
git checkout tags/${OBLAS_VERSION}

# Build 
# -> NUM_PARALLEL=... default automatically detected (if set to $(nproc) test 105/105 fork:safety_after_fork_in_parent will fail!)
# -> DYNAMIC_ARCH=1 DYNAMIC_OLDER=1 support for multiple (older) processor types on runtime in one ninary
# -> NO_SHARED=1 disable generating shared library
# -> PREFIX=... not required as automatically build in current folder
# -> INTERFACE64=1 enable 64-bit interface (not supported by all compilers!)
# -> BINARY=64 64-bit binary generation (seems to be set automatically if INTERFACE64=1)
make INTERFACE64=0 NO_LAPACKE=1 NO_CBLAS=1 NO_SHARED=1 USE_OPENMP=1 DYNAMIC_ARCH=1 DYNAMIC_OLDER=1

# Copy library files
mkdir -p "${PREFIX}/lib"
cp libopenblas.a "${PREFIX}/lib"
