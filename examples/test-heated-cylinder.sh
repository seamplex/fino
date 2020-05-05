#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

outcome=0
runfino heat-tran-3d.fin

if [ $? -ne 0 ]; then
  outcome=99
fi

# exit
exit $outcome
