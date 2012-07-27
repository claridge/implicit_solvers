
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS =


# ---------------------------------
# List of sources for this program:
# ---------------------------------


CLAW_SOURCES = \
	apply_linearized_pde_operator.f90 \
	apply_pde_operator.f90 \
	bc1.f90 \
	boundary_conditions.f90 \
	qinit.f90 \
	rp1.f90 \
	setprob.f90 \
	src1.f90 \
	../lib/1d/apply_newton_operator.f90 \
	../lib/1d/driver.f90 \
	../lib/1d/inner_product.f90 \
	../lib/1d/solve_newton_system_bicgstab.f90 \
	../lib/1d/take_backward_euler_step.f90 \
	../lib/1d/take_crank_nicolson_step.f90 \
	../lib/1d/take_forward_euler_step.f90 \
	../lib/LAPACK/dgtsv.f \
  ../lib/LAPACK/dgbsv.f \
  ../lib/LAPACK/dgbtf2.f \
  ../lib/LAPACK/dgbtrf.f \
  ../lib/LAPACK/dgbtrs.f \
  ../lib/LAPACK/dlaswp.f \
  ../lib/LAPACK/ilaenv.f \
  ../lib/LAPACK/ieeeck.f \
  ../lib/LAPACK/iparmq.f \
  ../lib/BLAS/idamax.f \
  ../lib/BLAS/lsame.f \
  ../lib/BLAS/dger.f \
  ../lib/BLAS/xerbla.f \
  ../lib/BLAS/dscal.f \
  ../lib/BLAS/dgemm.f \
  ../lib/BLAS/dswap.f \
  ../lib/BLAS/dtbsv.f \
  ../lib/BLAS/dcopy.f \
  ../lib/BLAS/dgemv.f \
  ../lib/BLAS/dtrsm.f



# Clawpack libraries to be used:
CLAW_LIB = $(CLAW)/clawpack/1d/lib

CLAW_LIBSOURCES = \
	$(CLAW_LIB)/b4step1.f \
  $(CLAW_LIB)/claw1.f \
  $(CLAW_LIB)/claw1ez.f \
  $(CLAW_LIB)/copyq1.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/out1.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/step1fw.f



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
