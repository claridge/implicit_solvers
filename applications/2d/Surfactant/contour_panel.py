import sys
from implicit_claw.python.claw_solution_2d import ClawSolution
import pylab

# Set plot parameters
pylab.rcParams.update({
'lines.linewidth': 1.2,
'font.size': 16,
})

# Read solution from output files.
solution_dir = '_output' if len(sys.argv) < 3 else sys.argv[2]
solution = ClawSolution(solution_dir)
x, y = solution.GetCellCenters()

# Plot contours at intervals of 1/15 on [0, 1].
v = [float(i)/15 for i in xrange(16)]

# Plot for times 5, 25, 50, and 100.
frame_numbers = [1, 5, 10, 20]

for i, frame_number in enumerate(frame_numbers):
  solution.SetFrame(frame_number)
  pylab.subplot(2, 2, i + 1)

  pylab.contour(x, y, solution.q[:,:,0], v, colors='k')
  pylab.title('$t = %s$' % int(solution.t))

  # Only label the y-axis for the leftmost subplots.
  if i + 1 in [1, 3]:
    pylab.ylabel('$y$')
    pylab.yticks(range(-3, 4))
  else:
    pylab.yticks([])

  # Only label x-axis for the bottom subplots.
  if i + 1 in [3, 4]:
    pylab.xlabel('$x$')
    pylab.xticks(range(0, 7, 2))
  else:
    pylab.xticks([])

pylab.show()
