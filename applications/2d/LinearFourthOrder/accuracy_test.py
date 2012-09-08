from numpy import *
import accuracy_testing
import pylab
import pyclaw.data


def BuildRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=2)

  probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
  probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
  probdata.add_param('max_time_step_splits', 0)
  probdata.add_param('newton_max_iter', 30)
  probdata.add_param('newton_reduction_factor', .5)
  probdata.add_param('newton_tolerance', 1e-8)
  probdata.add_param('newton_verbosity', 0)
  probdata.add_param('cg_tolerance', 1e-8)
  probdata.add_param('cg_verbosity', 0)
  probdata.add_param('num_threads', 1)
  
  # probdata.add_param('bc_options', ['0', '0', '0', '0'])
  probdata.add_param('bc_options', ['13', '13', '03', '03'])

  clawdata = rundata.clawdata
  clawdata.ndim = 2
  clawdata.xlower = .1
  clawdata.xupper = 1.1
  clawdata.ylower = .1
  clawdata.yupper = 1.1
  clawdata.meqn = 1
  clawdata.maux = 1
  clawdata.mcapa = 0
  clawdata.t0 = 0.0
  clawdata.outstyle = 1
  clawdata.nout = 1
  clawdata.verbosity = 1
  clawdata.dt_variable = 0
  clawdata.dt_max = 100
  clawdata.cfl_desired = 0.9
  clawdata.cfl_max = 1.0
  clawdata.max_steps = 10000
  clawdata.order = 1
  clawdata.order_trans = 1
  clawdata.mwaves = 1
  clawdata.mthlim = [0]
  clawdata.src_split = 1
  clawdata.mbc = 2
  clawdata.mthbc_xlower = 0
  clawdata.mthbc_xupper = 0
  clawdata.mthbc_ylower = 0
  clawdata.mthbc_yupper = 0
  clawdata.restart = 0
  clawdata.N_restart = 0

  # Specified by AccuracyTest.
  clawdata.mx = None
  clawdata.my = None
  clawdata.tfinal = None
  clawdata.dt_initial = None

  return rundata


def TrueSolution(x, y, t):
  return exp(-4*t) * sin(x) * sin(y)


if __name__ == '__main__':
  run_simulations, show_plots = accuracy_testing.ParseFlags()

  t_final = 0.2
  steps1 = t_final / 1e-1
  steps2 = t_final / 1e-2
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 5)][:-1]
  dt_values = array([t_final / n for n in num_steps])
  mx_min = 10

  test = accuracy_testing.AccuracyTest(BuildRunData, TrueSolution, t_final, dt_values, mx_min)

  if run_simulations:
    test.RunSimulations()

  test.CalculateErrors()
  
  for name in sorted(test.errors.keys()):
    print test.CheckConvergenceOrder(name, 1.9)
  
  if show_plots:
    test.errors['L1'].Plot('ro')
    test.errors['L1'].PlotDtFit('r:')

    test.errors['L2'].Plot('go')
    test.errors['L2'].PlotDtFit('g:')

    test.errors['LInfinity'].Plot('bo')
    test.errors['LInfinity'].PlotDtFit('b:')
    
    pylab.legend(loc='lower right')

  raw_input('Press ENTER to finish')
