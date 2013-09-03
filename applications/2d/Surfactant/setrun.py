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

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')

    probdata.add_param('newton_max_iter', 30)
    probdata.add_param('newton_tolerance', 1e-4)
    probdata.add_param('newton_verbosity', 2)
    probdata.add_param('linear_solver_tolerance', 1e-4)
    probdata.add_param('linear_solver_verbosity', 1)
    probdata.add_param('num_threads', 4)
    probdata.add_param('film_bc_options', ['13', '01', 'p', 'p'])
    probdata.add_param('surfactant_bc_options', ['1', '0', 'p', 'p'])
    
    probdata.add_param('beta', 0., 'gravitational constant')  # Estimated value: beta=.271
    probdata.add_param('kappa', 1e-4, 'capillarity')  # Estimated value: kappa=.013
    probdata.add_param('delta', 1e-4, 'surfactant diffusivity')
    probdata.add_param('mu', 0.,  'surface tension parameter: sigma = (1+mu*Gamma)**(-3)')
    probdata.add_param('right_fllm_height', 0.05)
    
    
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    clawdata.ndim = ndim
    clawdata.meqn = 2
    clawdata.maux = 4
    clawdata.mcapa = 0

    clawdata.xlower = 0
    clawdata.xupper = 2*pi    
    clawdata.ylower = -pi
    clawdata.yupper = pi
    clawdata.mx = 200
    clawdata.my = 200
    
    clawdata.t0 = 0.0
    clawdata.dt_initial = .1
    clawdata.dt_variable = 0
    clawdata.tfinal = 100
    clawdata.nout = 20

    clawdata.verbosity = 1
    
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 10000

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
    
