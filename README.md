Forked from: https://github.com/xmjiao/mumps4m-openmp and https://github.com/blechta/mumps-matlab-recipes \
Original Authors: Xiangmin Jiao, Jan Blechta

Recipe tested with Ubuntu 22.04.1, MATLAB 2021b
- dependencies:\
openblas 0.3.26 \
mumps    5.6.2 \
metis    5.1.0 \
scotch   7.0.3

- package requirements:
bison flex zlib-dev

- compiler requirements: \
gcc-9 gfortran-9 \
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021b-supported-compilers.pdf \
may use update-alternatives under Linux: \
https://linuxconfig.org/how-to-switch-between-multiple-gcc-and-g-compiler-versions-on-ubuntu-20-04-lts-focal-fossa \

- prepend before calling MATLAB: \
```export MATLABPATH="<Path-to-mex-files>":<sth-like-/usr/local/MATLAB>"```

- may prepend before calling MATLAB: \
```export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.5```

Simply execute ```make_...``` in order (and hope that it will also work as simple ;-)

## Literature
P. R. Amestoy, I. S. Duff, J. Koster and J.-Y. L’Excellent,
A fully asynchronous multifrontal solver using distributed dynamic
scheduling, SIAM Journal on Matrix Analysis and Applications,
Vol 23, No 1, pp 15-41 (2001).

P. R. Amestoy, A. Buttari, J.-Y. L’Excellent and T. Mary,
Performance and scalability of the block low-rank multifrontal
factorization on multicore architectures,
ACM Transactions on Mathematical Software,
Vol 45, Issue 1, pp 2:1-2:26 (2019)


