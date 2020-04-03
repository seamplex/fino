% FINO(1) Fino User Manual
% Jeremy Theler
% April 2, 2020

# NAME

fino - a free finite-element thermo-mechanical solver

# SYNOPSIS

fino [*options*] inputfile [*optional_extra_arguments*]...


# DESCRIPTION

Fino is a free and open source tool released under the terms of the GPLv3+ that uses the finite-element method to solve

 * steady-state thermo-mechanical problems, or
 * steady or transient heat conduction problems, or
 * modal analysis problems.
 
# OPTIONS

include(help.md)

# REFERENCE

## Fino keywords

esyscmd([x../wasora/doc/reference.sh parser kw | sed 's/^#/##/' x])


## Mesh keywords

esyscmd([x../wasora/doc/reference.sh parser kw ../wasora/src/mesh  | sed 's/^#/##/' x])


## Special input distributions

TBD.


## Boundary conditions

TBD.


## Result functions

TBD.


## Wasora keywords

esyscmd([x../wasora/doc/reference.sh parser kw ../wasora/src  | sed 's/^#/##/' x])


## Fino variables

esyscmd([x../wasora/doc/reference.sh init va  | sed 's/^#/##/' x])


# SEE ALSO

The Fino web page contains full source code, updates, examples, V&V cases and full reference:
<https://www.seamplex.com/fino>.

