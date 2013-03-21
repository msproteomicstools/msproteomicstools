#!/bin/bash
#put this script in same folder as mProphet.R and then add to PATH to simplify mProphet calls (no path hardcoding required)
export MPROPHET_BINDIR=$(dirname $0)
set -x
R --slave --file=$MPROPHET_BINDIR/mProphet.R --args bin_dir=$MPROPHET_BINDIR $@
