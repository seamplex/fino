---
title: Fino
subtitle: a free finite-element thermo-mechanical solver
desc: a completely free-as-in-freedom finite-element thermo-mechancial solver desinged and implemented following the UNIX principles
author: Jeremy Theler
date: April 16th, 2020
version: v0.7
infoname: fino
lang: en-US
---

# Overview


Fino is a free and open source tool released under the terms of the GPLv3+ that uses the finite-element method to solve

 * steady-state thermo-mechanical problems, or
 * steady or transient heat conduction problems, or
 * modal analysis problems.

![Updates, examples, V&V cases and full reference: <https://www.seamplex.com/fino>](fino-logo)


# Running `fino`

## Invocation

The format for running the `fino` program is:

```
fino [options] inputfile [optional_extra_arguments]...
```

The `fino` executable supports the following options:


include(help.md)


## Example input files

### Minimum working example

The following is a MWE input file for Fino that reads a Gmsh-generated `.msh` file, solves a linear elastic problem and wriets the results in a `.vtk` file which can be post-processed by Paraview:

```
include(../examples/tensile-mwe.fin)
```

### Extended annotated example

The example can be extended to give more information as the following annotated input shows:

```
include(../examples/tensile-test.fin)
```



# Test case

include(000-tensile-test.m4)


# Reference

include(reference-manual.md)
