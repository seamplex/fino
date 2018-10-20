#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh


outcome=0
gmsh -v 0 -3 wire-general.geo
if [ $? != 0 ]; then
  exit 99
fi

runfino wire.fin | tee wire.dat
if [ ! -e wire.dat ]; then
  exit 99
fi
relerror=`grep ^1 wire.dat | cut -f 4 | cut -c1`

if [ $relerror -lt 2 ]; then
  outcome=0
else
  outcome=99
fi

# exit
exit $outcome
