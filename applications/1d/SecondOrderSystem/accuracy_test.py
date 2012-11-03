from numpy import *
import accuracy_testing
import pylab
import pyclaw.data
import setrun


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
  probdata.add_param('bc_options_1', ['0', 'n'])
  probdata.add_param('bc_options_2', ['0', '0'])
  
  clawdata = rundata.clawdata
  clawdata.ndim = 1
  clawdata.xlower = 0.
  clawdata.xupper = .5
  clawdata.meqn = 2
  clawdata.mcapa = 0
  clawdata.nout = 1
  clawdata.verbosity = 0
  clawdata.dt_variable = 0
  clawdata.src_split = 1
  clawdata.mbc = 2

  # These will be specified separately for each test.
  clawdata.mx = None
  clawdata.tfinal = None
  clawdata.dt_initial = None

  return rundata


_T0 = .3


def TrueSolution(x, t):
    tau = t + _T0
    rho = x / tau**(1./3)
    return (2. * rho, 1 / (6 * tau**(1./3)) * (1. - rho))


if __name__ == '__main__':
  run_simulations, show_plots = accuracy_testing.ParseFlags()

  t_final = 0.5
  steps1 = t_final / 1e-1
  steps2 = t_final / 1e-2
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 11)]

  test = accuracy_testing.AccuracyTest(
    BuildRunData, TrueSolution, t_final,
    dt_values=array([t_final / n for n in num_steps]),
    mx_min=10)

  if run_simulations:
    test.RunSimulations()

  test.CalculateErrors()
  
  for name in sorted(test.errors.keys()):
    print test.CheckConvergenceOrder(name, 1.9)
  
  if show_plots:
    test.PlotErrors()

  raw_input('Press ENTER to finish')
