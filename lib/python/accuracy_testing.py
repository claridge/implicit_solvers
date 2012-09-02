import setrun
import pyclaw.data
import pyclaw.runclaw

from numpy import *
from numpy import linalg
import pylab
import sys
import getopt

import claw_solution_1d
import claw_solution_2d


def _GetOutputDirectory(i):
  return '_output%02d' % i


class NumericalError(object):
  
  def __init__(self, name):
    assert name in ('L1', 'L2', 'LInfinity')
    self.name = name
    self.dt_values = []
    self.dx_values = []
    self.error_values = []
    self.scale_factor = None
    self.exponent = None
  
  def AddDataPoint(self, dt, dx, cellwise_error):
    self.dt_values.append(dt)
    self.dx_values.append(dx)
    self.error_values.append(self._CalculateErrorValues(dx, cellwise_error))
  
  def _CalculateErrorValues(self, dx, cellwise_error):
    dv = reduce(lambda x, y: x*y, dx) if isinstance(dx, tuple) else dx
    
    if self.name == 'L1':
      return sum(abs(cellwise_error)) * dv
    elif self.name == 'L2':
      return sqrt(sum(cellwise_error**2)) * dv
    elif self.name == 'LInfinity':
      return abs(cellwise_error).max()
      
  def PowerFit(self):
    """Fits error values with a curve c*dt**exponent.

    Args:
      dt_values: Vector of dt values, parallel to error_values.
      error_values: Vector of error values.
    """
    A = vstack([ones(len(self.error_values)), log(self.dt_values)]).T
    c, exponent = linalg.lstsq(A, log(self.error_values))[0]
    self.scale_factor = exp(c)
    self.exponent = exponent

  # TODO: Pass kwargs?
  def Plot(self, style='ko'):
    if self.name == 'L1':
      legend_name = 'L^1'
    elif self.name == 'L2':
      legend_name = 'L^2'
    elif self.name == 'LInfinity':
      legend_name = 'L^\infty'

    pylab.loglog(self.dt_values, self.error_values, style,
                 label='$%s$ error' % legend_name)
                 
  def PlotDtFit(self, style='k:'):
    pylab.loglog(self.dt_values, self.scale_factor*self.dt_values**self.exponent,
                 style, label='$\sim \Delta t^{%.3f}$' % self.exponent)
    

class AccuracyTest(object):

  def __init__(self, build_rundata, true_solution,
               t_final, dt_values, mx_min):
    self._build_rundata = build_rundata
    self._true_solution = true_solution
    self._t_final = t_final
    self._dt_values = dt_values

    mx_values = [int(round(mx_min * max(dt_values) / dt)) for dt in dt_values]
    self._mx_values = array(mx_values)
    
    # Not the most efficient way of getting the dimension, but it works.
    self.ndim = build_rundata().clawdata.ndim
    
    self.errors = None

  def RunSimulations(self):
    rundata = self._build_rundata()
    clawdata = rundata.clawdata
    clawdata.tfinal = self._t_final

    for i in xrange(len(self._dt_values)):
      clawdata.dt_initial = self._dt_values[i]
      clawdata.mx = self._mx_values[i]
      if self.ndim == 2:
        clawdata.my = clawdata.mx
      rundata.write()
      pyclaw.runclaw.runclaw(xclawcmd='xclaw',
                             outdir=_GetOutputDirectory(i))

  def CalculateErrors(self):
    self.errors = {'L1': NumericalError('L1'),
                   'L2': NumericalError('L2'),
                   'LInfinity': NumericalError('LInfinity')}
    for i in xrange(len(self._dt_values)):
      if self.ndim == 1:
        solution = claw_solution_1d.ClawSolution(_GetOutputDirectory(i))
        solution.SetFrame(1)
        true_solution = array([self._true_solution(x, self._t_final)
                              for x in solution.GetCellCenters()])
        cellwise_error = abs(solution.q[:,0] - true_solution)
        for e in self.errors.itervalues():
          e.AddDataPoint(self._dt_values[i], solution.dx, cellwise_error)
      elif self.ndim == 2:
        solution = claw_solution_2d.ClawSolution(_GetOutputDirectory(i))
        solution.SetFrame(1)
        x, y = solution.GetCellCenters()
        
        true_solution = zeros((solution.mx, solution.my))
        for ix in xrange(solution.mx):
          for iy in xrange(solution.my):
            true_solution[ix,iy] = self._true_solution(x[ix,iy], y[ix,iy], self._t_final)

        cellwise_error = abs(solution.q[:,:,0] - true_solution)
        for e in self.errors.itervalues():
          e.AddDataPoint(self._dt_values[i], (solution.dx, solution.dy), cellwise_error)
        

    for e in self.errors.itervalues():
      e.PowerFit()

  def CheckConvergenceOrder(self, name, target_order):    
    def _Status(passing):
      if passing:
        return '\033[1;32mOK\033[1;m'  # Print green
      else:
        return '\033[1;31mFAILED\033[1;m'  # Print red

    error = self.errors[name]
    exponent = error.exponent
    passing = exponent > target_order
    return '%s error ~ dt**%f ... %s' % (name, exponent, _Status(passing))


def ParseFlags():
  try:
    options, arguments = getopt.getopt(sys.argv[1:], '', ['show_plots', 'norun_simulations'])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

  run_simulations = True
  show_plots = False
  for option, value in options:
    if option == '--norun_simulations':
      run_simulations = False
    elif option == '--show_plots':
      show_plots = True
  
  return run_simulations, show_plots