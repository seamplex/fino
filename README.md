---
title: Fino, a free finite element solver
lang: en-US
...


:::{.text-center}
![](doc/fino.svg)

[Download Fino](https://www.seamplex.com/fino/install.html){.btn .btn-primary}
[Test cases](https://www.seamplex.com/fino/cases){.btn .btn-secondary}
[Reference sheet](https://www.seamplex.com/fino/reference.html){.btn .btn-info}
:::


# What

[Fino](https://www.seamplex.com/fino) is a [free](https://www.gnu.org/philosophy/free-sw.en.html) and [open source](https://opensource.com/resources/what-open-source) tool released under the terms of the [GPLv3+](#licensing) that uses the finite-element method to solve

 * steady or quasistatic thermo-mechanical problems, or
 * steady or transient heat conduction problems, or
 * modal analysis problems.

:::{.alert .alert-warning}
Please note that Fino is a [back-end](https://en.wikipedia.org/wiki/Front_and_back_ends) aimed at advanced users. For an easy-to-use web-based front-end with Fino running in the cloud directly from your browser see [CAEplex](https://www.caeplex.com) at <https://www.caeplex.com>.

::::: {.embed-responsive .embed-responsive-16by9 .mb-3}
 <iframe class="embed-responsive-item" src="https://www.youtube.com/embed/kD3tQdq17ZE" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
:::::
:::

## Why another finite-element program? The world is already full of them!

Because almost every piece of FEA software falls in either one of these two categories:

 a. Powerful, flexible and complex advanced numerical solvers of general non-linear partial differential equations written by academics (for academics) distributed under the form of 
     i. libraries, which the user has to compile and link to their own codes, or
     ii. interpreted languages (i.e. Python) wrappers, which the user has to call from their own scripts, or
     iii. input-file reading binaries, which the user needs to fill in with the weak form of the equation they need to solve.
     
    **Examples:** [MoFEM](http://mofem.eng.gla.ac.uk/mofem/html/), [Sparselizard](http://sparselizard.org/), [GetDP](http://getdp.info/), [FEnICS](https://fenicsproject.org/), [MOOSE](https://mooseframework.org/), [FreeFEM](https://freefem.org/), ...  
     
 b. Commercial, non-free (well some of them are free but coded in FORTRAN 77 so the source is unintelligible) and complex GUI-based programs that are
     i. closed-source, so nobody can know what the actual equations are nor how they are solved, and/or
     ii. complicated, so the only way to use them is through their embedded mouse-based GUI, and/or
     iii. expensive and out of the league of many companies and professionals.
     
    **Examples:** [CalculiX](http://calculix.de/), [CodeAster](https://www.code-aster.org), [NASTRAN](https://github.com/nasa/NASTRAN-95)^[We list just the open-source ones because we [at Seamplex do not want to encourage the usage of non-free software](https://www.seamplex.com/mission.html#principles), but any of the commercial packages out there would also apply.] 
     
Hence, Fino tries to fill in the gap between these two worlds with a different design basis.^[Somewhat like [milonga](https://www.seamplex.com/milonga) although the landscape is slightly different.]
Read the foreword of the [tensile-test example](https://seamplex.com/fino/cases/000-tensile-test/) within the [case files](https://seamplex.com/fino/cases) for a deeper insight into Fino’s design and implementation philosophy.
     
## How does Fino fill in the gap?

Fino...

 * is [free and open source software](#licensing). It is free as in “free speech,” it gives every user [the four essential freedoms](https://www.gnu.org/philosophy/free-sw.en.html) and its [source code is published as a Git repository](https://github.com/seamplex/fino).
 * reads an English-like input file with the problem definition. And here is the main thing: **simple problems ought to have simple inputs**.  See [the first example](#tensile-test) and the first cases in the [Fino case files](https://www.seamplex.com/fino/cases)
 * follows, among [others](https://www.seamplex.com/principles.html), the [UNIX philosophy](https://en.wikipedia.org/wiki/Unix_philosophy).
 * according to the UNIX rule of separation, leaves the pre and post-processing steps to software written by professional programmers in the CAD management, meshing and analysis fields. 
 * is written in plain [ANSI C](https://en.wikipedia.org/wiki/ANSI_C) (neither C++ nor Fortran) and uses state-of-the-art libraries and resources---such as [GNU Autotools](https://www.gnu.org/software/automake/manual/automake.html) and [PETSc](https://www.mcs.anl.gov/petsc/)---written by professional programmers. 
 * delegates all the grid management to [Gmsh](http://gmsh.info/) in a way that the input file Fino reads does not (necessarily) contain any reference to the mesh properties. In particular, this means that the very same Fino input file can be used to solve the same problem with different meshing schemes (shapes, sizes, optimizations, etc.).
 * can write [VTK](https://vtk.org/) files to be post-processed with [Paraview](https://www.paraview.org/) or other compatible tool. It can also write [MSH](http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format) files for post-processing with [Gmsh](http://gmsh.info/).
 * is particularly designed to handle complex dependence of material properties (i.e. temperature-dependent properties) and boundary conditions dependent on space over different geometric entities (i.e. volumes, faces, edges and/or vertices).
 * can perform parametric or optimization runs---see [the parametric cantilever study](#cantilever).





# Features

Fino uses a main input file (see below for [examples](#examples)), that in turn instructs Fino to read one or more mesh files in [Gmsh](http://gmsh.info/) format. Fino works on top of the [wasora framework](https://www.seamplex.com/wasora) so it shares its [design basis](https://www.seamplex.com/docs/SP-FI-17-BM-5460-A.pdf) and inherits all of its features:


  * evaluation of [algebraic expressions](https://www.seamplex.com/wasora/doc/realbook/002-expressions)
  * [one](https://www.seamplex.com/wasora/doc/realbook/007-functions) and [multi-dimensional](https://www.seamplex.com/wasora/doc/realbook/010-2dfunctions) function interpolation
  * [scalar](https://www.seamplex.com/wasora/doc/realbook/000-hello), [vector](https://www.seamplex.com/wasora/doc/realbook/006-fibonacci) and matrix operations
  * numerical [integration](https://www.seamplex.com/wasora/doc/realbook/008-integrals), [differentiation](https://www.seamplex.com/wasora/doc/realbook/012-mechanics) and root finding of functions
  * possibility to solve iterative and/or [time-dependent](https://www.seamplex.com/wasora/doc/realbook/003-lag) problems
  * adaptive [integration of systems of differential-algebraic equations](https://www.seamplex.com/wasora/realbook/real-018-waterwheel.html)
  * I/O from files and shared-memory objects (with optional synchronization using semaphores) for external coupling
  * execution of arbitrary code provided as shared object files
  * parametric runs using quasi-random sequence numbers to efficiently sweep a sub-space of parameter space 
  * solution of systems of non-linear algebraic equations
  * non-linear fit of scattered data to one or multidimensional functions
  * non-linear multidimensional optimization
  * management of unstructured grids, definition and operation of functions and smooth interpolation between meshes

Output is 100% defined in the input file. If no explicit output instruction is provided, Fino remains silent (as required by the [UNIX rule of silence](http://www.linfo.org/rule_of_silence.html)). Besides terminal and plain-text files (e.g. user-defined results in [JSON](https://en.wikipedia.org/wiki/JSON)), post-processing files in [VTK](http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html) o [MSH](http://gmsh.info/doc/texinfo/gmsh.html#File-formats) formats can be generated.

# Examples

This section lists a few relevant examples.
See the [Fino cases list](https://www.seamplex.com/fino/cases) for a full list of annotated examples, verification cases and discussion of results.
See also the directory [examples](https://github.com/seamplex/fino/tree/master/examples) in the source repository for further problems.


## Tensile test 

Let us consider a tensile test specimen like [this one](https://caeplex.com/p/41dd1)


![](examples/tensile-test-cad.png){.img-fluid}\ 


One end of the specimen is fixed and the other one has a tension load of 10 kN.
We would like to obtain the displacements and stresses distribution within the geometry.
Given that the problem is relatively simple, the input file [`examples/tensile-test.fin`](examples/tensile-test.fin) should also be rather simple:

```fino
# tensile test example for Fino, see https://caeplex.com/p/41dd1
MESH FILE_PATH tensile-test.msh  # mesh file in Gmsh format (either version 2.2 or 4)

# uniform properties given as scalar variables
E = 200e3   # [ MPa ] Young modulus = 200 GPa
nu = 0.3    # Poisson’s ratio

# boundary conditions ("left" and "right" come from the names in the mesh)
PHYSICAL_ENTITY left  BC fixed       # fixed end
PHYSICAL_ENTITY right BC Fx=1e4      # [ N ] load in x+

FINO_SOLVER PROGRESS_ASCII  # print ascii progress bars (optional) 
FINO_STEP                   # solve

# compute reaction force at fixed end
FINO_REACTION PHYSICAL_ENTITY left RESULT R

# write results (Von Mises, principal and displacements) in a VTK file
MESH_POST FILE_PATH tensile-test.vtk sigma sigma1 sigma2 sigma3 VECTOR u v w

# print some results (otherwise output will be null)
PRINT "displ_max = " %.3f displ_max "mm"
PRINT "sigma_max = " %.1f sigma_max "MPa"
PRINT "principal1 at center = " %.8f sigma1(0,0,0) "MPa"
PRINT "reaction  = [" %.3e R "] Newtons"
```

If we ran this example from a terminal, we would get something like this:

```bash
$ fino tensile-test.fin
....................................................................................................
----------------------------------------------------------------------------------------------------
====================================================================================================
displ_max =     0.076   mm
sigma_max =     160.5   MPa
principal1 at center =  99.99948213     MPa
reaction  = [   -1.000e+04      1.622e-03       6.226e-03       ] Newtons
$
```

 * The three lines with the dots, dashes and double dashes are ASCII progress bars for the assembly of the stiffness matrix, the solution of the linear system and the computation of stresses, respectively.
 * There is no need to have a node at the origin to know what the stress is at $\vec{x}=(0,0,0)$. Fino (actually [wasora](https://www.seamplex.com/wasora)) can evaluate functions at any arbitrary point.
 * Almost any location in the input file where a numerical value is expected can be replaced by an algebraic expression, including standard functions like `log`, `exp`, `sin`, etc. See [wasora’s reference](https://www.seamplex.com/wasora/reference.html#functions).
 * If the `MESH_POST` and `PRINT` instructions were not included, there would not be any default output of the execution ([UNIX rule of silence](http://www.linfo.org/rule_of_silence.html)). This should be emphasized, as I have recently (i.e. twelve years after the commercial introduction of smartphones) stumbled upon a the output file of a classical FEM program that seems to be written in the 1970s:  paginated ASCII text apparently ready to be printed with all the possible numerical output because the CPU cost of re-running the case of course overwhelms the hourly rate of the engineer that has to understand the results. 
 * The output is 100% controlled by the user, including the precision of the printed results with [`printf` format specifiers](https://en.wikipedia.org/wiki/Printf_format_string). Note the eight decimal positions in the evaluation of $\sigma_1$ at the origin, whilst the expected value was $100~\text{MPa}$ (the load is $F_x=10^4~\text{N}$ and the cross-sectional area is $100~\text{mm}^2$).
 * The generation of the mesh `tensile-test.msh` is not covered in this example. See [`tensile-test.geo`](https://github.com/seamplex/fino/blob/master/examples/tensile-test.geo). Yet, please do consider the comparison of two few-lines [syntactically-sweetened](https://en.wikipedia.org/wiki/Syntactic_sugar) plain-text files which live near the English language with a similar input file for a certain open-source solver (whose input is in turn inspired by another non-free solver) that condenses both the mesh and the problem definition in a single 5Mb file that lives far away from plain English. 
 * The `VTK` output can be post-processed with the free tool [ParaView](http://www.paraview.org/):

![Tensile test results obtained by [Fino](https://www.seamplex.com/fino) and post-processed by [ParaView](http://www.paraview.org/).](examples/tensile-test.png){.img-fluid #fig:tensile-vtk}


## Cantilever beam with first & second order elements {#cantilever}

This example is far more complex as it studies the shear-locking effect of first order elements under bending loads, how displacements compares to second-order elements and how they depend on mesh size.
Hence, the input file is also more complex. 

```fino
DEFAULT_ARGUMENT_VALUE 1 1        # use first (1) or second (2) order elements
DEFAULT_ARGUMENT_VALUE 2 0        # use structured (1) or unstructured (0) tets
DEFAULT_ARGUMENT_VALUE 3 5        # maximum number of elements along h

h = 10   # beam width and height
l = 50   # beam length

PARAMETRIC n MIN 2 MAX $3 STEP 1
n = 2

OUTPUT_FILE geo  cantilever-$1-$2-%g.geo n
M4 {
 INPUT_FILE_PATH  cantilever.geo.m4
 OUTPUT_FILE geo
 MACRO h      h
 MACRO l      l
 MACRO lc     $1*h/n
 MACRO order  $1
 MACRO struct $2
}

SHELL "gmsh -v 0 -3 cantilever-$1-$2-%g.geo" n

INPUT_FILE mesh cantilever-$1-$2-%g.msh n
MESH FILE mesh DIMENSIONS 3

P = 1000   # load in [ N ]
E = 200e3  # Young modulus in [ MPa ]
nu = 0.3   # Poisson’s ratio

PHYSICAL_ENTITY NAME bulk
PHYSICAL_ENTITY NAME left  BC fixed
PHYSICAL_ENTITY NAME right BC Fz=-1000


FINO_STEP

# this is fun but takes a lot of time for a test
# energy_density(x,y,z) := 0.5*{( 
#  sigmax(x,y,z)*dudx(x,y,z) +
#  sigmay(x,y,z)*dvdy(x,y,z) +
#  sigmaz(x,y,z)*dwdz(x,y,z) +
#  tauxy(x,y,z)*(dudy(x,y,z) + dvdx(x,y,z)) +
#  tauyz(x,y,z)*(dvdz(x,y,z) + dwdy(x,y,z)) +
#  tauzx(x,y,z)*(dwdx(x,y,z) + dudz(x,y,z))
# )}
# MESH_INTEGRATE FUNCTION energy_density OVER bulk RESULT integrated_energy

FINO_REACTION PHYSICAL_ENTITY left RESULT R

# reference max deflection according to Euler-Bernoulli
# https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory#Cantilever_beams
wc = P*l^3/(3*E*(h^4)/12)

PRINT %.3e 1/nodes %g n nodes elements \
      %.5g -w(l,0,0) wc sigma_max \
      %.3f time_cpu_build time_cpu_solve time_cpu_stress %.0f memory/1e6 \
      %g R(3) strain_energy #integrated_energy

OUTPUT_FILE vtk cantilever-$1-$2-%g.vtk n
MESH_POST FILE vtk sigma sigma1 sigma2 sigma3 VECTOR u v w sigma
```

![Cantilever beam displacement for different grids and element order.](examples/cantilever.svg){.img-fluid width=100% #fig:cantilever-displ}

![Cantilever beam strain energy for different grids and element order.](examples/cantilever-energy.svg){.img-fluid width=100% #fig:cantilever-energy}


## Thermal conduction in a piston engine

Problem taken from [Simscale’s thermal tutorial](https://www.simscale.com/docs/content/tutorials/tutorial_heat-transfer.html):


```fino
# thermal conductivity in an engine piston as in
# https://www.simscale.com/docs/content/tutorials/tutorial_heat-transfer.html

SHELL "if [ ! -e piston.msh ]; then gmsh -v 0 -3 piston.geo; fi"
MESH FILE_PATH piston.msh        # the mesh is in mm
FINO_PROBLEM HEAT DIMENSIONS 3

f = 1e-3   # factor to convert from m to mm
# thermal conductivity numerically in W/(m*K) converted to W/(mm*K)
k = 160*f

# heat transfer coefficient in W/(m^2*K) converted to W/(mm^2*K)
# note that the names contain spaces so they must be quoted
PHYSICAL_ENTITY "top"                BC   h=450*f^2   Tref=1400
PHYSICAL_ENTITY "ring 1"             BC   h=150*f^2   Tref=450
PHYSICAL_ENTITY "ring 1 groove"      BC   h=1e3*f^2   Tref=450
PHYSICAL_ENTITY "ring 2"             BC   h=150*f^2   Tref=450
PHYSICAL_ENTITY "ring 2 groove"      BC   h=400*f^2   Tref=380
PHYSICAL_ENTITY "ring 3"             BC   h=150*f^2   Tref=380
PHYSICAL_ENTITY "ring 3 groove"      BC   h=400*f^2   Tref=380
PHYSICAL_ENTITY "interior and skirt" BC   h=650*f^2   Tref=380

FINO_STEP

MESH_POST FILE_PATH piston-temp.vtk T
MESH_POST FILE_PATH piston-temp.msh T

PRINT "\# cpu time [sec] = "  %.2f time_cpu_build "(build) "  time_cpu_solve "(solve)"  SEP " "
PRINT "\# memory [Mb]    = "  %.0f memory/1024^2
PRINT %.0f T(0,0,0)
```

![Fino results.](examples/engine-piston.png){.img-fluid #fig:engine-piston}

![Simscale (CalculiX) results.](examples/piston-simscale.png){.img-fluid #fig:piston-simscale}


## Conic valve

Can your solver constrain your model faces to algebraically-defined surfaces such as cones? Ours can (and it is open source):

[![Conic valve](examples/conic_valve.png){.img-fluid}](https://twitter.com/seamplex/status/789440535329181696)\ 

```fino
# can you fem solver constrain your model faces
# to algebraically-defined surfaces such as cones?
# ours can! (and it is open source)
# <https://twitter.com/seamplex/status/789440535329181696>

SHELL "gmsh -v 0 -3 conic_valve.geo"
MESH FILE_PATH conic_valve.msh DIMENSIONS 3

FINO_SOLVER PROGRESS_ASCII

E = 200e3
nu = 0.3

PHYSICAL_ENTITY NAME base  BC u=0 v=1e-2 w=0
PHYSICAL_ENTITY NAME top   BC u=0 v=1e-2 w=0

# the cone equation
x1 = -4
y1 = 2

x2 = -2
y2 = 4

f(x) := (y2-y1)/(x2-x1)*(x-x1) + y1
h = f(0)
r = root(f(x), x, -10, 0) 

PHYSICAL_ENTITY NAME cone  BC 0=((x+u)^2+(z+w)^2)/(r/h)^2-(y+v-h)^2

FINO_STEP
MESH_POST FILE_PATH conic_valve.vtk sigma VECTOR u v w #dudx dvdx dwdx dudy dvdy dwdy dudz dvdz dwdz 
```

See the original tweet at <https://twitter.com/seamplex/status/789440535329181696>

## Thermal expansion of finite cylinders

![Veeder Benchmark problem](examples/veeder.png){.img-fluid #fig:veeder}


See <https://www.seamplex.com/docs/SP-FI-17-BM-5460-A.pdf>.

# Licensing

Fino is distributed under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) version 3 or (at your option) any later version. The following text was borrowed from the [Gmsh documentation](http://gmsh.info/doc/texinfo/gmsh.html#Copying-conditions). Replacing “Gmsh” with “Fino” gives:

::: {.alert .alert-light}
Fino is “free software”; this means that everyone is free to use it and to redistribute it on a free basis. Fino is not in the public domain; it is copyrighted and there are restrictions on its distribution, but these restrictions are designed to permit everything that a good cooperating citizen would want to do. What is not allowed is to try to prevent others from further sharing any version of Fino that they might get from you.

Specifically, we want to make sure that you have the right to give away copies of Fino, that you receive source code or else can get it if you want it, that you can change Fino or use pieces of Fino in new free programs, and that you know you can do these things.

To make sure that everyone has such rights, we have to forbid you to deprive anyone else of these rights. For example, if you distribute copies of Fino, you must give the recipients all the rights that you have. You must make sure that they, too, receive or can get the source code. And you must tell them their rights.

Also, for our own protection, we must make certain that everyone finds out that there is no warranty for Fino. If Fino is modified by someone else and passed on, we want their recipients to know that what they have is not what we distributed, so that any problems introduced by others will not reflect on our reputation.

The precise conditions of the license for Fino are found in the [General Public License](https://github.com/seamplex/fino/blob/master/COPYING) that accompanies the source code. Further information about this license is available from the GNU Project webpage <http://www.gnu.org/copyleft/gpl-faq.html>.
:::

# Further information

Home page: <https://www.seamplex.com/fino>  
Repository: <https://github.com/seamplex/fino.git>  
Mailing list and bug reports: <wasora@seamplex.com>  (you need to subscribe first at <wasora+subscribe@seamplex.com>)  
Web interface for mailing list: <https://www.seamplex.com/lists.html>  
Follow us: [Twitter](https://twitter.com/seamplex/) [YouTube](https://www.youtube.com/channel/UCC6SzVLxO8h6j5rLlfCQPhA) [LinkedIn](https://www.linkedin.com/company/seamplex/) [Github](https://github.com/seamplex)

---------------------------

Fino is copyright ©2014-2020 [Seamplex](https://www.seamplex.com)  
Fino is licensed under [GNU GPL version 3](http://www.gnu.org/copyleft/gpl.html) or (at your option) any later version.  
Fino is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  
See the [copying conditions](COPYING).  
