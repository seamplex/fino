#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

outcome=0
runfino tensile-test.fin
if [ $? -ne 0 ]; then
  outcome=99
fi

if [ "`cat tensile-sigma.dat`" != "100" ]; then
  outcome=99
fi

# exit
exit $outcome
