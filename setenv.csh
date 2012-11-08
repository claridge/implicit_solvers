#!/bin/csh
# Sets environment variable IMPLICIT_CLAW to the current directory, and appends
# it to PYTHONPATH.  Unlike setenv.bash, this script isn't smart enough to be
# executed from a separate directory.  (Might be possible, but I don't know
# much csh.)

setenv IMPLICIT_CLAW_LIB "$cwd"
setenv PYTHONPATH "$PYTHONPATH":"$IMPLICIT_CLAW_LIB/lib/python"
