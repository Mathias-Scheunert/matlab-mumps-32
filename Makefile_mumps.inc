# Begin orderings
LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

ISCOTCHDIR = $(topdir)/../../include
LSCOTCHDIR = $(topdir)/../../lib
ISCOTCH   = -I$(ISCOTCHDIR)
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr

IMETISDIR = $(topdir)/../../include
LMETISDIR = $(topdir)/../../lib
IMETIS    = -I$(IMETISDIR)
LMETIS    = -L$(LMETISDIR) -lmetis

# Corresponding variables reused later
ORDERINGSF = -Dmetis -Dpord -Dscotch 
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(IMETIS) $(IPORD) $(ISCOTCH)
IORDERINGSC = $(IORDERINGSF)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
CC = gcc -fPIC
FC = gfortran -fPIC 
# -> -fallow-argument-mismatch -fallow-invalid-boz fixes https://github.com/coin-or-tools/ThirdParty-Mumps/issues/4 (only required for gcc/gfortran >= 10)
# -> -fPIC produces adress position independent code
FL = gfortran -fPIC
AR = ar vr 
RANLIB = ranlib

LAPACK = -L$(topdir)/../lib/ -lopenblas

LIBBLAS = -L$(topdir)/../lib/ -lopenblas
LIBOTHERS = -lpthread

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq/ -lmpiseq

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O3 -fopenmp -DBLR_MT
OPTL    = -O3 -fopenmp -DBLR_MT
OPTC    = -O3 -fopenmp -DBLR_MT
#End Optimized options

INCS = $(INCSEQ)
LIBS = -L$(libdir) -lmumps_common $(LPORD) $(LSCOTCH) $(LMETIS) $(LIBSEQ) $(LAPACK) $(LIBBLAS)
LIBSEQNEEDED = libseqneeded
