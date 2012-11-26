from numpy import *
import pyclaw.data
import unittest
import lib.python.convergence_test as convergence_test


def BuildRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=2)

  probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
  probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
  probdata.add_param('max_time_step_splits', 0)
  probdata.add_param('newton_max_iter', 10)
  probdata.add_param('newton_reduction_factor', .5)
  probdata.add_param('newton_tolerance', 1e-8)
  probdata.add_param('newton_verbosity', 0)
  probdata.add_param('cg_tolerance', 1e-8)
  probdata.add_param('cg_verbosity', 0)
  probdata.add_param('num_threads', 1)
  
  probdata.add_param('bc_options_1', ['0', 'n', '0', 'n'])
  probdata.add_param('bc_options_2', ['0', '0', '0', '0'])

  clawdata = rundata.clawdata
  clawdata.ndim = 2
  clawdata.xlower = .1
  clawdata.xupper = .6
  clawdata.ylower = .1
  clawdata.yupper = .6
  clawdata.meqn = 2
  clawdata.nout = 1
  clawdata.verbosity = 0
  clawdata.dt_variable = 0
  clawdata.src_split = 1
  clawdata.mbc = 2

  # These will be specified separately for each test.
  clawdata.mx = None
  clawdata.my = None
  clawdata.tfinal = None
  clawdata.dt_initial = None
  
  return rundata


_T0 = 1.

def TrueSolution(x, y, t):
    tau = t + _T0
    rho = sqrt(x**2 + y**2) / tau**.25
    return (2. * rho**2, -1./(8. * tau**.5) * log(rho))


class ConvergenceTest(convergence_test.ConvergenceTest):
  
  build_rundata = staticmethod(BuildRunData)
  true_solution = staticmethod(TrueSolution)
  t_final = 0.5
  steps1 = t_final / 1e-1
  steps2 = t_final / 1e-2
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 5)]
  dt_values = array([t_final / n for n in num_steps])
  mx_min = 10
  
  def testL1Convergence(self):
    self.assertGreater(self.GetConvergenceOrder('L1'), 1.9)

  def testL2Convergence(self):
    self.assertGreater(self.GetConvergenceOrder('L2'), 1.9)

  def testLInfinityConvergence(self):
    self.assertGreater(self.GetConvergenceOrder('LInfinity'), 1.9)
    
    
if __name__ == '__main__':
  unittest.main()