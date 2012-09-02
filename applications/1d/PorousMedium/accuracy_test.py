from numpy import *
import accuracy_testing
import pylab
import pyclaw.data


def BuildRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=1)

  probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
  probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
  probdata.add_param('max_time_step_splits', 0)
  probdata.add_param('newton_max_iter', 10)
  probdata.add_param('newton_reduction_factor', .5)
  probdata.add_param('newton_tolerance', 1e-8)
  probdata.add_param('newton_verbosity', 0)
  probdata.add_param('cg_tolerance', 1e-8)
  probdata.add_param('cg_verbosity', 0)

  probdata.add_param('bc_options', ['0', '0'])

  clawdata = rundata.clawdata
  clawdata.ndim = 1
  clawdata.xlower = -1.
  clawdata.xupper = 1.
  clawdata.meqn = 1
  clawdata.maux = 1
  clawdata.mcapa = 0
  clawdata.outstyle = 1
  clawdata.nout = 1
  clawdata.verbosity = 0
  clawdata.dt_variable = 0
  clawdata.dt_max = 1e+99
  clawdata.cfl_desired = 0.9
  clawdata.cfl_max = 1.0
  clawdata.max_steps = 10000
  clawdata.order = 1
  clawdata.order_trans = 0
  clawdata.mwaves = 1
  clawdata.mthlim = [0]
  clawdata.src_split = 1
  clawdata.mbc = 2
  clawdata.mthbc_xlower = 0
  clawdata.mthbc_xupper = 0
  clawdata.t0 = 0.0  # Setting this in clawdata doesn't appear to work.

  # These will be specified separately for each test.
  clawdata.mx = None
  clawdata.tfinal = None
  clawdata.dt_initial = None

  return rundata


_MASS = 1.
_T0 = .3


def TrueSolution(x, t):
    tau = t + _T0
    tmp = (_MASS - 1./12. * x**2 / tau**(2./3.)) / tau**(1./3.)
    return tmp if tmp > 0. else 0.


if __name__ == '__main__':
  run_simulations, show_plots = accuracy_testing.ParseFlags()

  t_final = 0.5
  steps1 = t_final / 1e-3
  steps2 = t_final / 1e-4
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 11)]
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
