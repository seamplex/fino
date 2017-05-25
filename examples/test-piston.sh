#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh

outcome=0
runfino engine-piston.fin > engine-piston.dat
if [ $? -ne 0 ]; then
  outcome=99
fi

if [ "`cat engine-piston.dat | cut -c-3`" != "593" ]; then
  outcome=99
fi

# exit
exit $outcome
