#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh

outcome=0
runfino engine-piston.fin > engine-piston.dat
if [ $? -ne 0 ]; then
  outcome=99
fi

if [ "`cat engine-piston.dat`" != "593.4" ]; then
  outcome=99
fi

# exit
exit $outcome
