#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

runfino tensile-test.fin
outcome=$?

# exit
exit $outcome
