

# uncomment if compiler supports MPI, to allow for MPI parallelisation
DMPI=-DMPI

# Compiler
F90 = mpif90
#F90 = gfortran

# main compiler flags
F90FLAGS = -O3 -funroll-loops 
#F90FLAGS = -g -Wall -fbacktrace -fbounds-check

# uncomment to use OpenMP parallisation
OMP = -fopenmp

# Normal blas/lapack. If MKL is aliased by these, the PARDISO (if included) will crush
BLAS = -lblas
LAPACK = -llapack 

# uncomment to allow for PARDISO sparse sover, it is available via linking with MKL
#DPARDI=-DPARDI

# MKL sequential blask/lapack, MKL allows for sparse PARDISO solver
#BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm 

# MKL Parallel blask/lapack. MKL allows for sparse PARDISO solver.
# Note the importanse to have a correct threading option of MKL
#BLAS = -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lm 

