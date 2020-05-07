#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh

outcome=0
if [ ! -e cylinder.msh ]; then
  gmsh -3 cylinder.geo
fi

runfino heat-tran-3d.fin

if [ $? -ne 0 ]; then
  outcome=99
fi

# exit
exit $outcome
