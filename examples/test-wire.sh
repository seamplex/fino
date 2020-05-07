#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh

# if Fino is not linked against SLEPc skip the test
${finobin} -i | grep -i slepc | grep -v none > /dev/null
if [ $? -eq 1 ]; then
  exit 77
fi  

outcome=0
if [ ! -e wire-general.msh ]; then
  gmsh -v 0 -3 wire-general.geo
fi  
if [ $? -ne 0 ]; then
  exit 99
fi

runfino wire.fin | tee wire.dat
if [ $? -eq 1 ]; then
  exit 99
fi
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
