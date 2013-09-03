""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
from pyclaw import data 
from math import pi

def setrun(claw_pkg='classic'):
    
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

    # Problem-specific parameters to be written to setprob.data:
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
    probdata.add_param('newton_max_iter', 10)
    probdata.add_param('newton_tolerance', 1e-8)
    probdata.add_param('newton_verbosity', 0)
    probdata.add_param('linear_solver_tolerance', 1e-8)
    probdata.add_param('linear_solver_verbosity', 0)
    probdata.add_param('num_threads', 1)
    probdata.add_param('bc_options', ['0', '1', '1', '1'])

    clawdata = rundata.clawdata

    clawdata.ndim = 2
    clawdata.meqn = 1

    clawdata.xlower = -1.
    clawdata.xupper = 1
    clawdata.ylower = -1.
    clawdata.yupper = 1
    clawdata.mx = 50
    clawdata.my = 50

    clawdata.tfinal = 1.0
    clawdata.dt_initial = 0.1

    clawdata.nout = 10
    clawdata.verbosity = 1
    clawdata.dt_variable = 0
    clawdata.src_split = 1
    clawdata.mbc = 2

    return rundata


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
	rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
    
