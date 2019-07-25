% Fino reference sheet
% Jeremy Theler

This reference sheet is for [Fino](index.html) v0.6.13-gead0afa
. 
Note that Fino works on top of [wasora](https://www.seamplex.com/wasora), so you should also check the [wasora reference sheet](https://www.seamplex.com/wasora/reference.html) also---not to mention the [wasora RealBook](https://www.seamplex.com/wasora/realbook).

# Keywords

##  FINO_DEBUG

Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.

~~~wasora
FINO_DEBUG [ FILE <file_id> | FILE_PATH <file_path> ] [ MATRICES_ASCII ] [ MATRICES_ASCII_STRUCTURE ] [ MATRICES_PETSC_BINARY ] [ MATRICES_PETSC_COMPRESSED_BINARY ] [ MATRICES_PETSC_ASCII ] [ MATRICES_PETSC_OCTAVE ] [ MATRICES_PETSC_DENSE ] [ MATRICES_X ] [ MATRICES_SNG ] [ MATRICES_SNG_STRUCT ] [ MATRICES_SIZE <expr> ] [ MATRICES_STRIDE <expr> ] [ INCLUDE_INPUT ]
~~~



##  FINO_LINEARIZE

Performs stress linearization according to ASME VII-Sec 5 over a
Stress Classification Line given either as a one-dimensional physical entity in the
mesh or as the (continuous) spatial coordinates of two end-points.
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
Von\  Mises criterion and stored as variables.
If no name for any of the variables is given, they are stored in
`M_entity`, `B_entity` and `P_entity` respectively if there is a physical entity.
Otherwise `M_1`, `B_1` and `P_1` for the first instruction, `M_2`... etc.

~~~wasora
FINO_LINEARIZE { PHYSICAL_ENTITY <physical_entity_name> | START_POINT <x1> <y1> <z1> END_POINT <x2> <y2> <z2> } [ FILE <file_id> | FILE_PATH <file_path> ] [ M <variable> ] [ MB <variable> ] [ PEAK <variable> ]
~~~



##  FINO_PROBLEM

Sets the problem type that Fino has to solve.      

~~~wasora
FINO_PROBLEM [ BAKE | SHAKE | BREAK | HEAT_AXISYMMETRIC | PLANE_STRESS | PLANE_STRAIN | ELASTIC_AXISYMMETRIC ] [ DIMENSIONS <expr> ] [ DEGREES <expr> ] [ SYMMETRY_AXIS { x | y } ] [ MESH <identifier> ] [ N_EIGEN <expr> ] [ UNKNOWNS <name1> <name2> ... <name_degrees> ]
~~~



 * `BAKE` (or `HEAT` or `THERMAL`) solves the heat conduction problem.
 * `SHAKE` (or `MODAL`) computes the natural frequencies and modes.
 * `BREAK` (or `ELASTIC` or `MECHANICAL`) solves the elastic problem.
 * `HEAT_AXISYMMETRIC` solves the heat conduction problem in an axysimmetric way.
 * `PLANE_STRESS` solves the plane stress elastic problem.
 * `PLANE_STRAIN` solves the plane strain elastic problem.
 * `ELASTIC_AXISYMMETRIC` solves the elastic problem in an axysimmetric way.

For the heat conduction problem the number of dimensions needs to be given either with the keyword `DIMENSIONS`
or by defining a `MESH` (with an explicit `DIMENSIONS` keyword) before `FINO_PROBLEM`.

##  FINO_SOLVER

Sets options related to the eigen-solver.

~~~wasora
FINO_SOLVER [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ] [ PC_TYPE { lu | gamg | hypre | sor | bjacobi | cholesky | ... } ] [ SET_NEAR_NULLSPACE { rigidbody | fino | none } ] [ DO_NOT_SET_BLOCK_SIZE | SET_BLOCK_SIZE ] [ GRADIENT_EVALUATION { mass_matrix_consistent | mass_matrix_row_sum | mass_matrix_lobatto | mass_matrix_diagonal | node_average_all | node_average_corner | gauss_average | none } ] [ GRADIENT_JACOBIAN_THRESHOLD <expr> ] [ PROGRESS_ASCII ] ///kw+FINO_SOLVER+usage [ PROGRESS_BUILD_SHM <shmobject> ] [ PROGRESS_SOLVE_SHMEM <shmobject> ] [ MEMORY_USAGE_SHMEM <shmobject> ] [ TOTAL {
~~~


List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>
List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>

##  FINO_STEP

Solves the linear eigenvalue problem.

~~~wasora
FINO_STEP [ JUST_BUILD | JUST_SOLVE ] [ DUMP_FILE_PATH <filepath> ]
~~~






# Variables

##  displ_max

The module of the maximum displacement of the elastic problem.



##  displ_max_x

The\ $x$ coordinate of the maximum displacement of the elastic problem.



##  displ_max_y

The\ $y$ coordinate of the maximum displacement of the elastic problem.



##  displ_max_z

The\ $z$ coordinate of the maximum displacement of the elastic problem.



##  fino_abstol

Absolute tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html)
Default `1e-50`.



##  fino_divtol

Divergence tolerance,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e+4`.  



##  fino_gamg_threshold

Relative threshold to use for dropping edges in aggregation graph for the
[Geometric Algebraic Multigrid Preconditioner](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html)
as passed to PETSc’s
[`PCGAMGSetThreshold`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMGSetThreshold.html).
A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
Default `0.01`.  



##  fino_iterations

This variable contains the actual number of iterations used
by the solver. It is set after `FINO_STEP`.



##  fino_max_iterations

Number of maximum iterations before diverging,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `10000`.



##  fino_penalty_weight

The weight $w$ used when setting multi-freedom boundary conditions.
Higher values mean better precision in the constrain but distort
the matrix condition number. 
Default is `1e8`.



##  fino_reltol

Relative tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e-6`.



##  fino_residual_norm

This variable contains the residual obtained
by the solver. It is set after `FINO_STEP`.



##  lambda


Requested eigenvalue. It is equal to 1.0 until
`FINO_STEP` is executed.  



##  memory

Maximum resident set size (global memory used), in bytes.



##  memory_available

Total available memory, in bytes.



##  memory_petsc

Maximum resident set size (memory used by PETSc), in bytes.



##  petsc_flops

Number of floating point operations performed by PETSc/SLEPc.



##  sigma_max

The maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_x

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_y

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  sigma_max_z

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  time_cpu_build

CPU time insumed to build the problem matrices, in seconds.



##  time_cpu_solve

CPU time insumed to solve the problem, in seconds.



##  time_cpu_stress

CPU time insumed to compute the stresses from the displacements, in seconds.



##  time_petsc_build

CPU time insumed by PETSc to build the problem matrices, in seconds.



##  time_petsc_solve

CPU time insumed by PETSc to solve the eigen-problem, in seconds.



##  time_petsc_stress

CPU time insumed by PETSc to compute the stresses, in seconds.



##  time_wall_build

Wall time insumed to build the problem matrices, in seconds.



##  time_wall_solve

Wall time insumed to solve the problem, in seconds.



##  time_wall_stress

Wall time insumed to compute the stresses, in seconds.



##  time_wall_total

Wall time insumed to initialize, build and solve, in seconds.
CPU time insumed to initialize, build and solve, in seconds.
CPU time insumed by PETSc to initialize, build and solve, in seconds.



##  T_max

The maximum temperature\ $T_\text{max}$ of the thermal problem.



##  T_min

The minimum temperature\ $T_\text{min}$ of the thermal problem.



##  u_at_displ_max

The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.



##  u_at_sigma_max

The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  v_at_displ_max

The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.



##  v_at_sigma_max

The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  w_at_displ_max

The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.



##  w_at_sigma_max

The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.






