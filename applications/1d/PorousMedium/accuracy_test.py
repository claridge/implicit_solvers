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
  t0 = .3
  tau = t + t0
  one_third = 1./3.
  a = (4.5)**one_third

  if abs(x) <= a * tau**one_third:
    return 1. / (6 * tau) * ((a * tau**one_third)**2 - x**2)
  else:
    return 0.


def BuildDefaultClawRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=1)

  probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
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
  clawdata.verbosity = 1
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


t_final = 0.5
num_steps = [round(x) for x in logspace(log10(5), log10(50), 11)]
dt_values = array([t_final / n for n in num_steps])

refinement_test = RefinementTest(t_final, dt_values, 20)
refinement_test.RunSimulations()

l2, linf = refinement_test.CalculateErrors()
pylab.loglog(dt_values, l2, 'r.', label='$L^2$ error')
pylab.loglog(dt_values, linf, 'b.', label='$L^\infty$ error')
pylab.loglog(dt_values, dt_values/100, 'k--', label='$\sim\Delta t$')
pylab.loglog(dt_values, dt_values**2/80, 'k-.', label='$\sim\Delta t^2$')
legend(loc='lower right')
