
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
CLAW_OUTDIR = _output               # Directory for output
CLAW_PLOTDIR = _plots               # Directory for plots
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_setplot_file = setplot.py      # File containing function to set plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS ?= -O3 -fopenmp -fno-align-commons


# ---------------------------------
# List of sources for this program:
# ---------------------------------

IMPLICIT_LIB=$(IMPLICIT_CLAW)/implicit_claw

CLAW_SOURCES = \
  set_implicit_boundary_data.f90 \
  qinit.f90 \
  apply_pde_operator.f90 \
  setprob.f90 \
	surface_tension_utils.f90 \
  $(IMPLICIT_LIB)/1d/bc_coefficients.f90 \
  $(IMPLICIT_LIB)/1d/interface_derivatives.f90 \
  $(IMPLICIT_LIB)/2d/apply_linearized_pde_operator.f90 \
  $(IMPLICIT_LIB)/2d/get_backward_euler_lhs.f90 \
  $(IMPLICIT_LIB)/2d/get_backward_euler_rhs.f90 \
  $(IMPLICIT_LIB)/2d/get_laplacian.f90 \
  $(IMPLICIT_LIB)/2d/builtin_bc_routines.f90 \
  $(IMPLICIT_LIB)/2d/driver.f90 \
  $(IMPLICIT_LIB)/2d/inner_product.f90 \
  $(IMPLICIT_LIB)/2d/setprob_implicit.f90 \
  $(IMPLICIT_LIB)/2d/solve_backward_euler_system-bicgstab.f90 \
  $(IMPLICIT_LIB)/2d/src2.f90 \
  $(IMPLICIT_LIB)/2d/take_backward_euler_step.f90 \
  $(IMPLICIT_LIB)/2d/take_crank_nicolson_step.f90 \
  $(IMPLICIT_LIB)/2d/take_forward_euler_step.f90 \
  $(IMPLICIT_LIB)/2d/null_riemann_solver.f90


# Clawpack library to be used:
CLAW_LIB = $(CLAW)/clawpack/2d/lib

CLAW_LIBSOURCES = \
  $(CLAW_LIB)/claw2ez.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/out2.f \
  $(CLAW_LIB)/restart2.f \
  $(CLAW_LIB)/claw2.f \
  $(CLAW_LIB)/step2.f \
  $(CLAW_LIB)/step2ds.f \
  $(CLAW_LIB)/dimsp2.f \
  $(CLAW_LIB)/flux2fw.f \
  $(CLAW_LIB)/copyq2.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
