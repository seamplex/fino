#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh

outcome=0
runfino piston.fin > piston.dat

if [ ! -e piston.dat ]; then
  exit 99
fi
temp=`grep -v \# piston.dat`

if (( $temp > 590 )) && (( $temp < 600 )); then
  outcome=0
else
  outcome=99
fi

# exit
exit $outcome
