import setrun
import pyclaw.data
import pyclaw.runclaw

from scipy import interpolate
from numpy import *
from numpy import linalg
import pylab
import sys
import getopt
import subprocess

import claw_solution_1d
import claw_solution_2d


# Print 'OK' or 'FAILED' in green and red, respectively
_OK = '\033[1;32mOK\033[1;m'
_FAILED = '\033[1;31mFAILED\033[1;m' 


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
    
  def _GetOutputDirectory(self, i):
    return '_output%02d' % i    

  def RunSimulations(self):
    """Makes xclaw and runs a simulation for each dt value."""
    
    # TODO: Natural to break here, as the make process carries its own
    # expectation.
    return_code = subprocess.call(['make', 'xclaw'])
    if return_code:
      print '\'make xclaw\' unsuccessful ... ' + _FAILED
      return
    
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
                             outdir=self._GetOutputDirectory(i))
    print 'Simulations executed ... ' + _OK

  def _GetCellwiseError(self, solution):
    """Gets cellwise error for a ClawSolution.
    
    Args:
      solution: 1d or 2d ClawSolution
    Returns:
      error (magnitude) for each cell and solution component.  Shape is the
        same as solution.values.
    """
    
    if self.ndim == 1:
      x_values = solution.GetCellCenters()
      true_values = zeros((solution.mx, solution.meqn))
      for ix in range(solution.mx):
        true_values[ix, :] = self._true_solution(x_values[ix], self._t_final)
    elif self.ndim == 2:
      x, y = solution.GetCellCenters()
      true_values = zeros((solution.mx, solution.my, solution.meqn))
      for ix in xrange(solution.mx):
        for iy in xrange(solution.my):
          true_values[ix,iy,:] = self._true_solution(x[ix,iy], y[ix,iy], self._t_final)
    return abs(solution.values - true_values)

  def CalculateErrors(self):
    """Fills dictionary of numerical errors: L1, L2, and LInfinity.
    
    Call after RunSimulations.
    """

    i_skip = -1
    if self._true_solution is None and self.ndim == 1:
      i_skip = list(self._dt_values).index(min(self._dt_values))
      most_refined = claw_solution_1d.ClawSolution(self._GetOutputDirectory(i_skip))
      most_refined.SetFrame(1)
        
      f = interpolate.interp1d(
          most_refined.GetCellCenters(), most_refined.values[:,0], kind='cubic')
      self._true_solution = lambda x, t: f(x)
    
    self.errors = {'L1': NumericalError('L1'),
                   'L2': NumericalError('L2'),
                   'LInfinity': NumericalError('LInfinity')}
                   
    for i in xrange(len(self._dt_values)):
      if i == i_skip: continue
      
      if self.ndim == 1:
        solution = claw_solution_1d.ClawSolution(self._GetOutputDirectory(i))
        dx_vector = solution.dx
      elif self.ndim == 2:
        solution = claw_solution_2d.ClawSolution(self._GetOutputDirectory(i))
        dx_vector = (solution.dx, solution.dy)

      solution.SetFrame(1)
      cellwise_error = self._GetCellwiseError(solution)

      for e in self.errors.itervalues():
        e.AddDataPoint(self._dt_values[i], dx_vector, cellwise_error)

    for e in self.errors.itervalues():
      e.PowerFit()

  def CheckConvergenceOrder(self, name, target_order):
    """Checks that an error type matches the specified convergence order.
    
    Args:
      name: Name of the error being specified; one of the valid names for
        NumericalError.
      target_order: The minimum order of convergence expected.
    Returns:
      A string indicating success or failure.
    """
    
    def _Status(passing):
      if passing:
        return _OK
      else:
        return _FAILED

    error = self.errors[name]
    exponent = error.exponent
    passing = exponent > target_order
    return '%s error ~ dt**%f ... %s' % (name, exponent, _Status(passing))
    
  def PlotErrors(self):
    self.errors['L1'].Plot('ro')
    self.errors['L1'].PlotDtFit('r:')

    self.errors['L2'].Plot('go')
    self.errors['L2'].PlotDtFit('g:')

    self.errors['LInfinity'].Plot('bo')
    self.errors['LInfinity'].PlotDtFit('b:')
    
    pylab.legend(loc='best')
    pylab.show()
    

def FakeTrueSolution(x, t, refined_solution):
  return interpolate.interp1d(
      x, refined_solution.GetCellCenters(), refined_solution.values[:,0])


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