% FINO(1) Fino User Manual
% Jeremy Theler
% April 2, 2020

# NAME

fino - a free finite-element thermo-mechanical solver

# SYNOPSIS

fino [*options*] input-file [*optional_extra_arguments*]...


# DESCRIPTION

Fino is a free and open source tool released under the terms of the GPLv3+ that uses the finite-element method to solve

 * steady-state thermo-mechanical problems, or
 * steady or transient heat conduction problems, or
 * modal analysis problems.

It can a Gmsh(1) mesh file and can write results for post-processing in either `.msh` or `.vtk` files.

 
# EXAMPLES

## Minimum working example

The following is a MWE input file for Fino that reads a Gmsh-generated `.msh` file, solves a linear elastic problem and wriets the results in a `.vtk` file which can be post-processed by Paraview:

```
include(../examples/tensile-mwe.fin)
```

The `.geo` file that generates the mesh with Gmsh and the CAD file in `.step` format can be found in the `examples` directory.

## Extended annotated example

The example above can be extended to give more information as the following annotated input shows:

```
include(../examples/tensile-test.fin)
```

 
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

`gmsh`(1), `paraview`(1)

The Fino Case files at <https://www.seamplex.com/fino/cases/> contains fully-discussed examples.

The Fino web page contains full source code, updates, examples, V&V cases and full reference:
<https://www.seamplex.com/fino>.

# AUTHOR

Jeremy Theler <https://www.seamplex.com>
