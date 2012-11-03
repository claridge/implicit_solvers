# This is a terribly hacky way of running tests.

import os
import subprocess

TEST_FILES = [
  'applications/1d/PorousMedium/accuracy_test.py',
  'applications/1d/ThinFilm/accuracy_test.py',
  'applications/1d/LinearDiffusion/accuracy_test.py',
  'applications/1d/LinearFourthOrder/accuracy_test.py',
  'applications/1d/SecondOrderSystem/accuracy_test.py',
  'applications/2d/PorousMedium/accuracy_test.py',
  'applications/2d/ThinFilm/accuracy_test.py',
  'applications/2d/LinearDiffusion/accuracy_test.py',
  'applications/2d/LinearFourthOrder/accuracy_test.py',
  'applications/2d/SecondOrderSystem/accuracy_test.py',
  ]

HOME = os.getcwd()

for test_file in TEST_FILES:
  os.chdir(HOME)
  directory = os.path.dirname(test_file)
  file_name = os.path.basename(test_file)
  os.chdir(directory)

  return_code = subprocess.call(['make', 'xclaw'])
  if return_code:
    print 'Make failed'
  else:
    subprocess.call(['python', file_name])


for test_file in TEST_FILES:
  os.chdir(HOME)
  directory = os.path.dirname(test_file)
  os.chdir(directory)
  return_code = subprocess.call(['make', 'clobber'])
