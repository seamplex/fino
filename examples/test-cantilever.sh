#!/bin/bash

# locate and include runtest.sh
. locateruntest.sh

checkgmsh
checkm4

outcome=0
for order in 1 2; do
 for struct in 0 1; do
  runfino cantilever.fin $order $struct | tee cantilever-${order}-${struct}.dat
  if [ $? -ne 0 ]; then
    outcome=99
  fi
 done
done

plot "set xrange [0:0.5]; set key bottom left; plot 'cantilever-1-0.dat' u 1:5 w lp, 'cantilever-1-1.dat' u 1:5 w lp, 'cantilever-2-0.dat' u 1:5 w lp, 'cantilever-2-1.dat' u 1:5 w lp"

# exit
exit $outcome
