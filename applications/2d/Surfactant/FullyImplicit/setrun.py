""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
from pyclaw import data 
from math import pi

#------------------------------
def setrun(claw_pkg='classic'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "classic" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    assert claw_pkg.lower() == 'classic',  "Expected claw_pkg = 'classic'"

    ndim = 2
    rundata = data.ClawRunData(claw_pkg, ndim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
    probdata.add_param('max_time_step_splits', 2,
                       'Max number of times to halve the lenght of the '
                       'implicit time step, should Newton''s method fail to '
                       'converge.')

    probdata.add_param('newton_max_iter', 100,
                       'Max iterations for Newton''s method before enforcing reduction criterion')
    probdata.add_param('newton_reduction_factor', .5,
                       'Required reduction in residual norm for Newton''s method to continue '
                       'beyond newton_max_iter iterations')
    probdata.add_param('newton_tolerance', 1e-4,
                       'Newton''s method stops when norm(delta(iterate)) is below this.')
    probdata.add_param('newton_verbosity', 2, 'Logging level for Newton''s method')
    
    probdata.add_param('cg_tolerance', 1e-4,
                       'CG or BiCGStab terminate when norm(residual) is below this.')
    probdata.add_param('cg_verbosity', 1, 'Logging level for CG/BiCGStab')

    probdata.add_param('num_threads', 4, 'Number of OpenMP threads.')
    
    probdata.add_param('film_bc_options', ['13', '01', 'p', 'p'])
    probdata.add_param('surfactant_bc_options', ['1', '0', 'p', 'p'])
    
    probdata.add_param('beta', 0., 'gravitational constant')  # Estimated value: beta=.271
    probdata.add_param('kappa', 1e-4, 'capillarity')  # Estimated value: kappa=.013
    probdata.add_param('delta', 1e-4, 'surfactant diffusivity')
    probdata.add_param('mu', 0.,  'surface tension parameter: sigma = (1+mu*Gamma)**(-3)')
    probdata.add_param('right_fllm_height', 0.05)
    
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = ndim
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = 0
    clawdata.xupper = 2*pi
    
    clawdata.ylower = -pi
    clawdata.yupper = pi
        

    # Number of grid cells:
    clawdata.mx = 200
    clawdata.my = 200

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 2

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 4
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = 1

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.nout = 20
        clawdata.tfinal = 100.

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.01, 0.02]   # used if outstyle == 2
        clawdata.tout += map(lambda x: 0.02+5*c*x,range(1,10))
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 5
        clawdata.iout = [iout, ntot]
    


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = .1
    
    # Max time step to be allowed if variable dt used:
    # clawdata.dt_max = clawdata.dt_initial
    clawdata.dt_max = .1  # Don't want the implicit steps getting too large.
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 10000 #int(2*(clawdata.tfinal/clawdata.nout)/clawdata.dt_initial)

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 2
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 2
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [4,4]
    
    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used, 
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity
    
    clawdata.mthbc_xlower = 0
    clawdata.mthbc_xupper = 0
    
    clawdata.mthbc_ylower = 2
    clawdata.mthbc_yupper = 2
    
    clawdata.restart = 0
    clawdata.N_restart = 0
    
    return rundata
    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
	rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
    
