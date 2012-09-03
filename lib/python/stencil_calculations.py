from numpy import *

def GetBoundaryCoefficients(derivative_orders,
                            input_coordinates,
                            output_coordinates):
  def _DxSuffix(order):
    if order == 0:
      return ''
    elif order == 1:
      return ' * dx'
    else:
      return ' * dx**%d' % order
      
  num_conditions = len(derivative_orders) + len(input_coordinates)
  rows = []
  orders_to_row_indices = {}
  for order in sorted(derivative_orders):
    row = zeros(num_conditions)
    row[order] = math.factorial(order)
    orders_to_row_indices[order] = len(rows)
    rows.append(row)

  coordinates_to_row_indices = {}
  for coordinate in input_coordinates:
    coordinates_to_row_indices[coordinate] = len(rows)
    rows.append([coordinate**i for i in xrange(num_conditions)])
  
  A = array(rows).T

  lines = []

  for i_output, coordinate in enumerate(output_coordinates):
    coefficients = linalg.solve(A, [coordinate**i for i in xrange(num_conditions)])

    for i, o in enumerate(derivative_orders):
      coefficient = coefficients[orders_to_row_indices[o]]
      coeff_string = ('%.15e' % coefficient).replace('e', 'd')
      lines.append('boundary_coefficients(%d, %d) = %s%s' % 
                   (i_output+1, i+1, coeff_string, _DxSuffix(o)))

    for i, c in enumerate(input_coordinates):
      coefficient = coefficients[coordinates_to_row_indices[c]]
      coeff_string = ('%.15e' % coefficient).replace('e', 'd')
      lines.append('cell_coefficients(%d, %d) = %s' % (i_output+1, i+1, coeff_string)) 

  return lines


def GetBcCoefficients():
  cases = [(0,), (1,), (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
  input_coordinates = (.5, 1.5)

  print 'Lower coefficients'
  print '=================='
  first = True
  for case in cases:
    if len(case) == 1:
      output_coordinates = (-.5,)
    else:
      output_coordinates = (-1.5, -.5)
    
    if first:
      print 'if (orders == \'%s\') then' % ''.join(str(o) for o in case)
      first = False
    else:
      print 'else if (orders == \'%s\') then' % ''.join(str(o) for o in case)
      
    lines = GetBoundaryCoefficients(case, input_coordinates, output_coordinates)
    lines = ['    ' + line for line in lines]
    print '\n'.join(lines)
  print 'end if'
  
  print '\n'
  
  print 'Upper coefficients'
  print '=================='
  cases = [(0,), (1,), (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
  input_coordinates = (-1.5, -.5)
  output_coordinates = (.5, 1.5)
  
  first = True
  for case in cases:
    if len(case) == 1:
      output_coordinates = (.5,)
    else:
      output_coordinates = (.5, 1.5)

    if first:
      print 'if (orders == \'%s\') then' % ''.join(str(o) for o in case)
      first = False
    else:
      print 'else if (orders == \'%s\') then' % ''.join(str(o) for o in case)
      
    lines = GetBoundaryCoefficients(case, input_coordinates, output_coordinates)
    lines = ['    ' + line for line in lines]
    print '\n'.join(lines)
  print 'end if'


def GetCubicExtrapolant():
  input_coordinates = (.5, 1.5, 2.5, 3.5)
  output_coordinates = (-1.5, -.5)
  lines = GetBoundaryCoefficients([], input_coordinates, output_coordinates)
  print '\n'.join(lines)


def GetQuarticExtrapolant():
  input_coordinates = (.5, 1.5, 2.5, 3.5, 4.5)
  output_coordinates = (-.5,)
  lines = GetBoundaryCoefficients([], input_coordinates, output_coordinates)
  print '\n'.join(lines)
  
if __name__ == '__main__':
  GetQuarticExtrapolant()