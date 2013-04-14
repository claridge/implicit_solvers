import sys
from lib.python.claw_solution_2d import ClawSolution
from pylab import contour, show

frame_number = 1 if len(sys.argv) == 1 else int(sys.argv[1])
solution_dir = '_output' if len(sys.argv) < 3 else sys.argv[2]
solution = ClawSolution(solution_dir)
solution.SetFrame(frame_number)
x, y = solution.GetCellCenters()
v = [float(i)/15 for i in xrange(16)]
contour(x, y, solution.q[:,:,0], v)
show()