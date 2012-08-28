import setrun
import setplot
import pyclaw.data
import pyclaw.runclaw
from pyclaw.plotters import plotclaw

from numpy import *
from numpy import linalg
import pylab

import claw_solution_1d
import unittest

_MASS = 1.
_T0 = .3


def TrueSolution(x, t):
    tau = t + _T0
    tmp = (_MASS - 1./12. * x**2 / tau**(2./3.)) / tau**(1./3.)
    return tmp if tmp > 0. else 0.


def BuildDefaultClawRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=1)

  probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
  probdata.add_param('implicit_integration_scheme', 'Crank-Nicolson')
  probdata.add_param('max_time_step_splits', 0,
                     'Max number of times to halve the lenght of the '
                     'implicit time step, should Newton''s method fail to '
                     'converge.')
  probdata.add_param('newton_max_iter', 30,
                     'Max iterations for Newton''s method before enforcing '
                     'reduction criterion')
  probdata.add_param('newton_reduction_factor', .5,
                     'Required reduction in residual norm for Newton''s method '
                     'to continue beyond newton_max_iter iterations')
  probdata.add_param('newton_tolerance', 1e-8, 'Newton''s method stops when '
                     'norm(delta(iterate)) is below this.')
  probdata.add_param('newton_verbosity', 0, 'Logging level for Newton''s method')

  probdata.add_param('cg_tolerance', 1e-8, 'Conjugate gradient stops when '
                     'norm(residual) is below this.')
  probdata.add_param('cg_verbosity', 0, 'Logging level for CG/BiCGStab')



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
  clawdata.dt_initial = 1e-1
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


  # Specify these manually
  clawdata.mx = None
  clawdata.tfinal = None

  return rundata

    
class RefinementTest(object):

  def __init__(self, t_final, dt_values, min_mx):
    self._t_final = t_final
    self._dt_values = dt_values

    max_dt = max(dt_values)
    mx_values = [int(round(min_mx * max_dt / dt)) for dt in dt_values]
    self._mx_values = array(mx_values)  

  def RunSimulations(self):
    rundata = BuildDefaultClawRunData()
    clawdata = rundata.clawdata
    clawdata.tfinal = self._t_final
    
    for i in xrange(len(dt_values)):
      clawdata.dt_initial = self._dt_values[i]
      clawdata.mx = self._mx_values[i]
      rundata.write()
      pyclaw.runclaw.runclaw(xclawcmd='xclaw > /dev/null', outdir='_output%02d' % i)

  def GetNumericalErrors(self):
    numerical_errors = {'l1': [], 'l2': [], 'linf': []}
    for i in xrange(len(self._dt_values)):
      solution = claw_solution_1d.ClawSolution('_output%02d' % i)
      solution.SetFrame(1)
      true_solution = array([TrueSolution(x, self._t_final)
                             for x in solution.GetCellCenters()])
      cellwise_error = abs(solution.q[:,0] - true_solution)

      numerical_errors['l1'].append(sum(cellwise_error) * solution.dx)
      numerical_errors['l2'].append(sqrt(sum(cellwise_error**2)) * solution.dx)
      numerical_errors['linf'].append(max(cellwise_error))

    return numerical_errors


def PowerFit(error_values, dt_values):
  """Fits error values with a curve c*dt**exponent.
  
  Args:
    error_values: Vector of error values.
    dt_values: Vector of dt values, parallel to error_values.
  """
  A = vstack([ones(len(error_values)), log(dt_values)]).T
  c, exponent = linalg.lstsq(A, log(error_values))[0]
  return exp(c), exponent


t_final = 0.5
num_steps = [round(x) for x in logspace(log10(5), log10(50), 11)]
dt_values = array([t_final / n for n in num_steps])

refinement_test = RefinementTest(t_final, dt_values, 10)
refinement_test.RunSimulations()

numerical_errors = refinement_test.GetNumericalErrors()

l1_factor, l1_exponent = PowerFit(numerical_errors['l1'], dt_values)
l2_factor, l2_exponent = PowerFit(numerical_errors['l2'], dt_values)
linf_factor, linf_exponent = PowerFit(numerical_errors['linf'], dt_values)


def GetTestMessages():
  def _OkOrFail(passing):
    return 'ok' if passing else 'FAILURE'

  messages = []
  passing = l1_exponent > 1.9
  messages.append('L1 error ~ dt**%f ... %s' % 
                  (l1_exponent, _OkOrFail(passing)))

  passing = l2_exponent > 1.9
  messages.append('L2 error ~ dt**%f ... %s' % 
                  (l2_exponent, _OkOrFail(passing)))

  passing = linf_exponent > 1.9
  messages.append('Linf error ~ dt**%f ... %s' % 
                  (linf_exponent, _OkOrFail(passing)))
  return messages
    
if __name__ == '__main__':
  print '\n'.join(GetTestMessages())
# pylab.loglog(dt_values, numerical_errors['l1'], 'ro', label='$L^1$ error')
# pylab.loglog(dt_values, l1_factor*dt_values**l1_exponent, 'r:', label='$\sim \Delta t^{%.3f}$' % l1_exponent)
# 
# pylab.loglog(dt_values, numerical_errors['l2'], 'go', label='$L^2$ error')
# pylab.loglog(dt_values, l2_factor*dt_values**l2_exponent, 'g:', label='$\sim \Delta t^{%.3f}$' % l2_exponent)
# 
# pylab.loglog(dt_values, numerical_errors['linf'], 'bo', label='$L^\infty$ error')
# pylab.loglog(dt_values, linf_factor*dt_values**linf_exponent, 'b:', label='$\sim \Delta t^{%.3f}$' % linf_exponent)
# 
# pylab.legend(loc='lower right')
