#!/bin/bash
# Sets environment variable IMPLICIT_CLAW to the directory in which this
# script resides, and appends it to PYTHONPATH.  You should be able to
# execute this directly from your .bashrc (or Mac equivalent) via
#     source path-to-this-script/setenv.bash

export IMPLICIT_CLAW="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH="${PYTHONPATH}":${IMPLICIT_CLAW}
