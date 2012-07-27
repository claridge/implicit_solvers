""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
from pyclaw import data 
from numpy import pi

#------------------------------
def setrun(claw_pkg='Classic'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "Classic4" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    assert claw_pkg.lower() == 'classic',  "Expected claw_pkg = 'classic'"

    rundata = data.ClawRunData(pkg=claw_pkg, ndim=1)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    
    probdata.add_param('newton_max_iter', 30, 'Max iterations for Newton''s '
                       'method before enforcing reduction criterion')
    probdata.add_param('newton_reduction_factor', .5, 'Required reduction in '
                       'residual norm for Newton''s method to continue beyond '
                       'newton_max_iter iterations')
    probdata.add_param('newton_tolerance', .1e-8, 'Newton''s method stops when '
                       'norm(delta(iterate)) is below this.')
    probdata.add_param('newton_verbosity', 1, 'Logging level for Newton''s method')
    
    probdata.add_param('cg_tolerance', 1e-12, 'CG or BiCGStab terminate when '
                       'norm(residual) is below this.')
    probdata.add_param('cg_verbosity', 1, 'Logging level for CG/BiCGStab')

    probdata.add_param('gamma', 1., 'whatever')
    
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated
    clawdata.ndim = 1
    clawdata.xlower = -pi
    clawdata.xupper = pi
    dx = .05
    clawdata.mx = int(round((clawdata.xupper - clawdata.xlower) / dx))
    clawdata.meqn = 1
    clawdata.maux = 1
    clawdata.mcapa = 0
    clawdata.t0 = 0.0
    clawdata.outstyle = 1
    clawdata.nout = 50
    clawdata.tfinal = .5
    clawdata.verbosity = 1
    clawdata.dt_variable = 0
    clawdata.dt_initial = 1e-3
    clawdata.dt_max = 1e+99
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    clawdata.max_steps = 10000
    clawdata.order = 2
    clawdata.order_trans = 0
    clawdata.mwaves = 1
    clawdata.mthlim = [0]
    clawdata.src_split = 1
    clawdata.mbc = 2
    clawdata.mthbc_xlower = 0  # Not used
    clawdata.mthbc_xupper = 0  # Not used
    
    return rundata



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
	rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
    
