
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
FFLAGS = -O3 -fbounds-check -llapack -lblas
#-fbounds-check


# ---------------------------------
# List of sources for this program:
# ---------------------------------

IMPLICIT_LIB=../../../lib/2d

CLAW_SOURCES = \
  boundary_conditions.f90 \
  initial_conditions.f90 \
  implicit_operators.f90 \
  setprob.f90 \
  true_solution.f90 \
  $(IMPLICIT_LIB)/apply_newton_operator.f90 \
  $(IMPLICIT_LIB)/driver.f90 \
  $(IMPLICIT_LIB)/extend_to_ghost_cells.f90 \
  $(IMPLICIT_LIB)/inner_product.f90 \
  $(IMPLICIT_LIB)/null_riemann_solver.f90 \
  $(IMPLICIT_LIB)/solve_newton_system_bicgstab.f90 \
  $(IMPLICIT_LIB)/src2.f90 \
  $(IMPLICIT_LIB)/take_backward_euler_step.f90 \
  $(IMPLICIT_LIB)/take_crank_nicolson_step.f90 \
  $(IMPLICIT_LIB)/take_forward_euler_step.f90


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
    