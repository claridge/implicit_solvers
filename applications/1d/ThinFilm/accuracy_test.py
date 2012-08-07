import setrun
import setplot
import pyclaw.data
import pyclaw.runclaw
from pyclaw.plotters import plotclaw

from numpy import *
from numpy import linalg
import pylab

import claw_solution_1d

def TrueSolution(x, t):
  ell = 2.
  t0 = .3
  tau = (5. * (t + t0))**.2
  eta = x / tau

  if abs(eta) < ell:
    return 1. / (24. * tau) * (ell**2 - eta**2)**2
  else:
    return 0.


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
  probdata.add_param('newton_tolerance', 1e-12, 'Newton''s method stops when '
                     'norm(delta(iterate)) is below this.')
  probdata.add_param('newton_verbosity', 0, 'Logging level for Newton''s method')

  probdata.add_param('cg_tolerance', 1e-12, 'Conjugate gradient stops when '
                     'norm(residual) is below this.')
  probdata.add_param('cg_verbosity', 0, 'Logging level for CG/BiCGStab')

  probdata.add_param('gamma', 1., 'TODO: remove this')


  clawdata = rundata.clawdata
  clawdata.ndim = 1
  clawdata.xlower = -1.
  clawdata.xupper = 1.
  clawdata.meqn = 1
  clawdata.maux = 1
  clawdata.mcapa = 0
  clawdata.outstyle = 1
  clawdata.nout = 1
  clawdata.verbosity = 1
  clawdata.dt_variable = 0
  clawdata.dt_initial = 1e-3
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

    # Using constant mx for now.
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
      pyclaw.runclaw.runclaw(xclawcmd='xclaw', outdir='_output%02d' % i)

  def CalculateErrors(self):
    l2_errors = []
    linfinity_errors = []
    for i in xrange(len(dt_values)):
      solution = claw_solution_1d.ClawSolution('_output%02d' % i)
      solution.SetFrame(1)
      true_solution = array([TrueSolution(x, self._t_final)
                             for x in solution.GetCellCenters()])
      error = abs(solution.q[:,0] - true_solution)
      l2_errors.append(sqrt(sum(error**2)) * solution.dx)
      linfinity_errors.append(max(error))

    return l2_errors, linfinity_errors


if __name__ == '__main__':
  t_final = 0.2
  steps1 = t_final / 1e-3
  steps2 = t_final / 1e-4
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 11)]
  dt_values = array([t_final / n for n in num_steps])

  refinement_test = RefinementTest(t_final, dt_values, 10)
  refinement_test.RunSimulations()

  l2, linf = refinement_test.CalculateErrors()
  pylab.loglog(dt_values, l2, 'r.', label='$L^2$ error')
  pylab.loglog(dt_values, linf, 'b.', label='$L^\infty$ error')
  pylab.loglog(dt_values, dt_values / 1000, 'k--', label='$\sim\Delta t$')
  pylab.loglog(dt_values, dt_values**2, 'k-.', label='$\sim\Delta t^2$')
  pylab.legend(loc='lower right')

