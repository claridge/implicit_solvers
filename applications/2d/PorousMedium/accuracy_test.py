import setrun
import pyclaw.data
import pyclaw.runclaw
from pyclaw.plotters import plotclaw

from numpy import *
from numpy import linalg
import pylab

import claw_solution_2d

_MASS = 1.
_T0 = .3

def TrueSolution(x, y, t):
    tau = t + _T0
    tmp = _MASS - 1./16. * (x**2 + y**2) / sqrt(tau)
    return tmp/sqrt(tau) if tmp > 0. else 0.


# def TrueSolution(x, y, t):
#     tau = t + _T0
#     tmp = (_MASS - 1./12. * y**2 / tau**(2./3.)) / tau**(1./3.)
#     return tmp if tmp > 0. else 0.


def BuildDefaultClawRunData():
  rundata = pyclaw.data.ClawRunData(pkg='Classic', ndim=2)

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
  probdata.add_param('newton_verbosity', 1, 'Logging level for Newton''s method')
  probdata.add_param('cg_tolerance', 1e-8, 'Conjugate gradient stops when '
                     'norm(residual) is below this.')
  probdata.add_param('cg_verbosity', 1, 'Logging level for CG/BiCGStab')
  probdata.add_param('num_threads', 4, 'Number of OpenMP threads.')


  clawdata = rundata.clawdata
  clawdata.ndim = 2
  clawdata.xlower = -1.
  clawdata.xupper = 1
  clawdata.ylower = -1.
  clawdata.yupper = 1
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


  clawdata.my = 4

  clawdata.mx = None
  clawdata.tfinal = None
  clawdata.dt_initial = None

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
      clawdata.my = self._mx_values[i]
      rundata.write()
      pyclaw.runclaw.runclaw(xclawcmd='xclaw', outdir='_output%02d' % i)

  def CalculateErrors(self):
    l2_errors = []
    linfinity_errors = []
    
    for i in xrange(len(dt_values)):
      solution = claw_solution_2d.ClawSolution('_output%02d' % i)
      solution.SetFrame(1)

      x, y = solution.cell_centers()
      # true_solution = TrueSolution(x, y, self._t_final)

      true_solution = zeros((solution.mx, solution.my))
      for i in xrange(solution.mx):
        for j in xrange(solution.my):
          true_solution[i,j] = TrueSolution(x[i,j], y[i,j], self._t_final)

      # array([TrueSolution(x, self._t_final)
      #                        for x in solution.cell_centers()])

      pointwise_error = abs(solution.q[:,:,0] - true_solution)
      l2_errors.append(sqrt(sum(pointwise_error**2)) * solution.dx * solution.dy)
      linfinity_errors.append(pointwise_error.max())

    return l2_errors, linfinity_errors


if __name__ == '__main__':
  t_final = 0.2
  steps1 = t_final / 1e-2
  steps2 = t_final / 1e-3
  num_steps = [round(x) for x in logspace(log10(steps1), log10(steps2), 11)]
  dt_values = array([t_final / n for n in num_steps])
  
  refinement_test = RefinementTest(t_final, dt_values, 10)
  refinement_test.RunSimulations()
  
  l2, linf = refinement_test.CalculateErrors()
  pylab.loglog(dt_values, l2, 'r.', label='$L^2$ error')
  pylab.loglog(dt_values, linf, 'b.', label='$L^\infty$ error')
  pylab.loglog(dt_values, dt_values/10, 'k--', label='$\sim\Delta t$')
  pylab.loglog(dt_values, dt_values**2*10, 'k-.', label='$\sim\Delta t^2$')
  #legend(loc='lower right')
