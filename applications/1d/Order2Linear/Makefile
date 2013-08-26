
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
FFLAGS ?= -O3 -fno-align-commons


# ---------------------------------
# List of sources for this program:
# ---------------------------------

IMPLICIT_LIB = $(IMPLICIT_CLAW)/lib/1d

CLAW_SOURCES = \
  apply_pde_operator.f90 \
	boundary_conditions.f90 \
	initial_conditions.f90 \
  setprob.f90 \
  true_solution.f90 \
  $(IMPLICIT_LIB)/builtin_bc_routines.f90 \
  $(IMPLICIT_LIB)/apply_linearized_pde_operator.f90 \
  $(IMPLICIT_LIB)/bc_coefficients.f90 \
  $(IMPLICIT_LIB)/driver.f90 \
  $(IMPLICIT_LIB)/get_backward_euler_rhs.f90 \
  $(IMPLICIT_LIB)/get_backward_euler_lhs.f90 \
  $(IMPLICIT_LIB)/get_true_solution_by_index_range.f90 \
  $(IMPLICIT_LIB)/inner_product.f90 \
  $(IMPLICIT_LIB)/interface_derivatives.f90 \
  $(IMPLICIT_LIB)/null_riemann_solver.f90 \
  $(IMPLICIT_LIB)/setprob_implicit.f90 \
  $(IMPLICIT_LIB)/solve_backward_euler_system-bicgstab.f90 \
  $(IMPLICIT_LIB)/src1.f90 \
  $(IMPLICIT_LIB)/take_backward_euler_step.f90 \
  $(IMPLICIT_LIB)/take_crank_nicolson_step.f90 \
  $(IMPLICIT_LIB)/take_forward_euler_step.f90


# Clawpack libraries to be used:
CLAW_LIB = $(CLAW)/clawpack/1d/lib

CLAW_LIBSOURCES = \
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
    
