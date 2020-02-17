% Fino reference sheet

This reference sheet is for [Fino](https://www.seamplex.com/fino) v0.6.88-g90f0ab5
. 
Note that Fino works on top of [wasora](https://www.seamplex.com/wasora), so you should also check the [wasora reference sheet](https://www.seamplex.com/wasora/reference.html) also---not to mention the [wasora RealBook](https://www.seamplex.com/wasora/realbook).
See Fino in action at the [Fino case files](https://www.seamplex.com/fino).

-   [Keywords](#keywords)
    -   [FINO\_DEBUG](#fino_debug)
    -   [FINO\_LINEARIZE](#fino_linearize)
    -   [FINO\_PROBLEM](#fino_problem)
    -   [FINO\_REACTION](#fino_reaction)
    -   [FINO\_SOLVER](#fino_solver)
    -   [FINO\_STEP](#fino_step)
-   [Variables](#variables)
    -   [displ\_max](#displ_max)
    -   [displ\_max\_x](#displ_max_x)
    -   [displ\_max\_y](#displ_max_y)
    -   [displ\_max\_z](#displ_max_z)
    -   [fino\_abstol](#fino_abstol)
    -   [fino\_divtol](#fino_divtol)
    -   [fino\_gamg\_threshold](#fino_gamg_threshold)
    -   [fino\_iterations](#fino_iterations)
    -   [fino\_max\_iterations](#fino_max_iterations)
    -   [fino\_penalty\_weight](#fino_penalty_weight)
    -   [fino\_reltol](#fino_reltol)
    -   [fino\_residual\_norm](#fino_residual_norm)
    -   [lambda](#lambda)
    -   [memory](#memory)
    -   [memory\_available](#memory_available)
    -   [memory\_petsc](#memory_petsc)
    -   [petsc\_flops](#petsc_flops)
    -   [sigma\_max](#sigma_max)
    -   [sigma\_max\_x](#sigma_max_x)
    -   [sigma\_max\_y](#sigma_max_y)
    -   [sigma\_max\_z](#sigma_max_z)
    -   [strain\_energy](#strain_energy)
    -   [time\_cpu\_build](#time_cpu_build)
    -   [time\_cpu\_solve](#time_cpu_solve)
    -   [time\_cpu\_stress](#time_cpu_stress)
    -   [time\_petsc\_build](#time_petsc_build)
    -   [time\_petsc\_solve](#time_petsc_solve)
    -   [time\_petsc\_stress](#time_petsc_stress)
    -   [time\_wall\_build](#time_wall_build)
    -   [time\_wall\_solve](#time_wall_solve)
    -   [time\_wall\_stress](#time_wall_stress)
    -   [time\_wall\_total](#time_wall_total)
    -   [T\_max](#t_max)
    -   [T\_min](#t_min)
    -   [u\_at\_displ\_max](#u_at_displ_max)
    -   [u\_at\_sigma\_max](#u_at_sigma_max)
    -   [v\_at\_displ\_max](#v_at_displ_max)
    -   [v\_at\_sigma\_max](#v_at_sigma_max)
    -   [w\_at\_displ\_max](#w_at_displ_max)
    -   [w\_at\_sigma\_max](#w_at_sigma_max)


# Keywords

##  FINO_DEBUG

> Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.

~~~wasora
FINO_DEBUG [ FILE <file_id> | FILE_PATH <file_path> ]
 [ MATRICES_ASCII ]
 [ MATRICES_ASCII_STRUCTURE ]
 [ MATRICES_PETSC_BINARY ]
 [ MATRICES_PETSC_COMPRESSED_BINARY ]
 [ MATRICES_PETSC_ASCII ]
 [ MATRICES_PETSC_OCTAVE ]
 [ MATRICES_PETSC_DENSE ]
 [ MATRICES_X ]
 [ MATRICES_SNG ]
 [ MATRICES_SNG_STRUCT ]
 [ MATRICES_SIZE <expr> ]
 [ MATRICES_STRIDE <expr> ]
 [ INCLUDE_INPUT ]

~~~



##  FINO_LINEARIZE

> Performs stress linearization according to ASME VII-Sec 5 over a
Stress Classification Line given either as a one-dimensional physical entity in the
mesh or as the (continuous) spatial coordinates of two end-points.

~~~wasora
FINO_LINEARIZE { PHYSICAL_ENTITY <physical_entity_name> | START_POINT <x1> <y1> <z1> END_POINT <x2> <y2> <z2> }
 [ FILE <file_id> | FILE_PATH <file_path> ]
 [ TOTAL { vonmises tresca | tresca | principal1 | principal2 | principal3 }
 [ M <variable> ]
 [ MB <variable> ]
 [ PEAK <variable> ]

~~~


If the SCL is given as a `PHYSICAL_ENTITY`, the entity should be one-dimensional (i.e a line)
independently of the dimension of the problem.
If the SCL is given with `START_POINT` and `END_POINT`, the number of coordinates given should
match the problem dimension (i.e three coordinates for full\ 3D problems and two coordinates for
axisymmetric or plane problems).
Coordinates can be given algebraic expressions that will be evaluated at the time of the linearization.
If either a `FILE` or a `FILE_PATH` is given, the total, membrane and membrane plus bending
stresses are written as a function of a scalar $t \in [0,1]$.
Moreover, the individual elements of the membrane and bending stress tensors are written
within comments (i.e. lines starting with the hash symbol `#`).
By default, the linearization uses the Von\ Mises criterion for the composition of stresses.
The definition of what _total stress_ means can be changed using the `TOTAL` keyword.
The membrane, bending and peak stress tensor elements are combined using the
Von\ Mises criterion and stored as variables.
If no name for any of the variables is given, they are stored in
`M_entity`, `B_entity` and `P_entity` respectively if there is a physical entity.
Otherwise `M_1`, `B_1` and `P_1` for the first instruction, `M_2`... etc.

##  FINO_PROBLEM

> Sets the problem type that Fino has to solve.      

~~~wasora
FINO_PROBLEM [ mechanical | thermal | modal ]
 [ AXISYMMETRIC | PLANE_STRESS | PLANE_STRAIN ] [ SYMMETRY_AXIS { x | y } ] [ LINEAR | NON_LINEAR ]
 [ DIMENSIONS <expr> ] [ MESH <identifier> ] 
 [ N_MODES <expr> ] 

~~~


 * `mechanical` (or `elastic` or `break`) solves the mechanical elastic problem (default).
 * `thermal` (or `heat` or `bake`) solves the heat conduction problem.
 * `modal` (or `shake`) computes the natural frequencies and oscillation modes.

If the `AXISYMMETRIC` keyword is given, the mesh is expected to be two-dimensional in the $x$-$y$ plane
and the problem is assumed to be axi-symmetric around the axis given by `SYMMETRY_AXIS` (default is $y$).
If the problem type is mechanical and the mesh is two-dimensional on the $x$-$y$ plane and no
axisymmetry is given, either `PLANE_STRESS` and `PLAIN_STRAIN` can be provided (default is plane stress).
By default Fino tries to detect wheter the computation should be linear or non-linear.
An explicit mode can be set with either `LINEAR` on `NON_LINEAR`.
The number of spatial dimensions of the problem needs to be given either with the keyword `DIMENSIONS`
or by defining a `MESH` (with an explicit `DIMENSIONS` keyword) before `FINO_PROBLEM`.
If there are more than one `MESH`es define, the one over which the problem is to be solved
can be defined by giving the explicit mesh name with `MESH`. By default, the first mesh to be
defined in the input file is the one over which the problem is solved.
The number of modes to be computed in the modal problem. The default is DEFAULT_NMODES.

##  FINO_REACTION

> Computes the reaction at the selected physical entity.

~~~wasora
FINO_REACTION PHYSICAL_ENTITY <physical_entity_name> [ RESULT { <variable> | <vector> } ]
~~~


The result is stored in the variable or vector provided, depending on the number of degrees of freedoms of the problem.
If the object passed as `RESULT` does not exist, an appropriate object (scalar variable or vector) is created.
For the elastic problem, the components of the total reaction force are stored in the result vector.
For the thermal problem, the total power passing through the entity is computed as an scalar.

##  FINO_SOLVER

> Sets options related to the solver and the computation of gradients.

~~~wasora
FINO_SOLVER [ PROGRESS_ASCII ]
 [ PC_TYPE { gamg | mumps | lu | hypre | sor | bjacobi | cholesky | ... } ]
 [ KSP_TYPE { gmres | mumps | bcgs | bicg | richardson | chebyshev | ... } ]
 [ SNES_TYPE { newtonls | newtontr | nrichardson | ngmres | qn | ngs | ... } ]
 [ GRADIENT { gauss | nodes | none } ]
 [ GRADIENT_HIGHER { average | nodes | } ]
 [ SMOOTH | ROUGH ]
 [ GRADIENT_ELEMENT_WEIGHT { volume | flat | quality | volume_times_quality } ]
 [ GRADIENT_QUALITY_THRESHOLD <expr> ]

~~~


If the keyword `PROGRESS_ASCII` is given, three ASCII lines will show in the terminal the
progress of the ensamble of the stiffness matrix, the solution of the linear system and the
computation of gradients (stresses).
The preconditioner, linear and non-linear solver might be any of those available in PETSc:

 * List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>. 
 * List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>. 
* List of `SNES_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESType.html>. 

If either `PC_TYPE` or `KSP_TYPE` is set to `mumps` (and PETSc is compiled with MUMPS support) then this direct solver is used.

##  FINO_STEP

> Solves the linear eigenvalue problem.

~~~wasora
FINO_STEP [ JUST_BUILD | JUST_SOLVE ]
~~~






# Variables

##  displ_max

> The module of the maximum displacement of the elastic problem.



##  displ_max_x

> The\ $x$ coordinate of the maximum displacement of the elastic problem.



##  displ_max_y

> The\ $y$ coordinate of the maximum displacement of the elastic problem.



##  displ_max_z

> The\ $z$ coordinate of the maximum displacement of the elastic problem.



##  fino_abstol

> Absolute tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html)
Default `1e-50`.



##  fino_divtol

> Divergence tolerance,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e+4`.  



##  fino_gamg_threshold

> Relative threshold to use for dropping edges in aggregation graph for the
[Geometric Algebraic Multigrid Preconditioner](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html)
as passed to PETSc’s
[`PCGAMGSetThreshold`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMGSetThreshold.html).
A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
Default `0.01`.  



##  fino_iterations

> This variable contains the actual number of iterations used
by the solver. It is set after `FINO_STEP`.



##  fino_max_iterations

> Number of maximum iterations before diverging,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `10000`.



##  fino_penalty_weight

> The weight $w$ used when setting multi-freedom boundary conditions.
Higher values mean better precision in the constrain but distort
the matrix condition number. 
Default is `1e8`.



##  fino_reltol

> Relative tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e-6`.



##  fino_residual_norm

> This variable contains the residual obtained
by the solver. It is set after `FINO_STEP`.



##  lambda

> 
Requested eigenvalue. It is equal to 1.0 until
`FINO_STEP` is executed.  



##  memory

> Maximum resident set size (global memory used), in bytes.



##  memory_available

> Total available memory, in bytes.



##  memory_petsc

> Maximum resident set size (memory used by PETSc), in bytes.



##  petsc_flops

> Number of floating point operations performed by PETSc/SLEPc.



##  sigma_max

> The maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_x

> The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_y

> The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_z

> The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  strain_energy

> The strain energy stored in the solid, computed as

[\ \frac{1}{2} \vec{u}^T \cdot K \cdot \vec{u}]

where $\vec{u}$ is the displacements vector and $K$ is the stiffness matrix.



##  time_cpu_build

> CPU time insumed to build the problem matrices, in seconds.



##  time_cpu_solve

> CPU time insumed to solve the problem, in seconds.



##  time_cpu_stress

> CPU time insumed to compute the stresses from the displacements, in seconds.



##  time_petsc_build

> CPU time insumed by PETSc to build the problem matrices, in seconds.



##  time_petsc_solve

> CPU time insumed by PETSc to solve the eigen-problem, in seconds.



##  time_petsc_stress

> CPU time insumed by PETSc to compute the stresses, in seconds.



##  time_wall_build

> Wall time insumed to build the problem matrices, in seconds.



##  time_wall_solve

> Wall time insumed to solve the problem, in seconds.



##  time_wall_stress

> Wall time insumed to compute the stresses, in seconds.



##  time_wall_total

> Wall time insumed to initialize, build and solve, in seconds.
CPU time insumed to initialize, build and solve, in seconds.
CPU time insumed by PETSc to initialize, build and solve, in seconds.



##  T_max

> The maximum temperature\ $T_\text{max}$ of the thermal problem.



##  T_min

> The minimum temperature\ $T_\text{min}$ of the thermal problem.



##  u_at_displ_max

> The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.



##  u_at_sigma_max

> The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  v_at_displ_max

> The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.



##  v_at_sigma_max

> The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  w_at_displ_max

> The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.



##  w_at_sigma_max

> The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.






