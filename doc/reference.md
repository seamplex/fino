% Fino reference sheet
% Jeremy Theler

This reference sheet is for [Fino](index.html) v0.5.102-g640bd66
. 

~~~
$ fino
fino v0.5.98-g882bf1e 
a partial differential equation solver based on the finite element method
$
~~~

Note that Fino works on top of [wasora](/wasora), so you should also check the [wasora reference sheet](/wasora/reference.html) also---not to mention the [wasora RealBook](/wasora/realbook).

# Keywords

##  `FINO_DEBUG`

Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.

~~~wasora
FINO_DEBUG [ FILE <file_id> | [ FILE_PATH <file_path> ] [ MATRICES_ASCII ] [ MATRICES_ASCII_STRUCTURE ] [ MATRICES_PETSC_BINARY ] [ MATRICES_PETSC_COMPRESSED_BINARY ] [ MATRICES_PETSC_ASCII ] [ MATRICES_PETSC_OCTAVE ] [ MATRICES_PETSC_DENSE ] [ MATRICES_X ] [ MATRICES_SNG ] [ MATRICES_SNG_STRUCT ] [ MATRICES_SIZE <expr> ] [ MATRICES_STRIDE <expr> ] [ INCLUDE_INPUT ]
~~~



##  `FINO_PROBLEM`


~~~wasora
FINO_PROBLEM [ BAKE | SHAKE | BREAK | HEAT_AXISYMMETRIC | PLANE_STRESS | PLANE_STRAIN | ELASTIC_AXISYMMETRIC | [ DIMENSIONS <expr> ] [ DEGREES <expr> ] [ SYMMETRY_AXIS { x | y } ] [ MESH <identifier> ] [ N_EIGEN <expr> ] [ UNKNOWNS <name1> <name2> ... <name_degrees> ]
~~~



##  `FINO_REACTION`

Asks Fino to compute the reactions at physical entities with Dirichlet boundary conditions.

~~~wasora
FINO_REACTION PHYSICAL_ENTITY <physical_entity> [ NAME_ROOT <name> ]
~~~



##  `FINO_SOLVER`

Sets options related to the eigen-solver.

~~~wasora
FINO_SOLVER [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ] [ PC_TYPE { lu | gamg | hypre | sor | bjacobi | cholesky | ... } ] [ SET_NEAR_NULLSPACE { rigidbody | fino | none } ] [ DO_NOT_SET_BLOCK_SIZE | SET_BLOCK_SIZE ] [ GRADIENT_EVALUATION { mass_matrix_consistent mass_matrix_row_sum mass_matrix_lobatto mass_matrix_diagonal node_average_all node_average_corner gauss_average none } ] [ GRADIENT_JACOBIAN_THRESHOLD <expr> ] [ SHMEM_PROGRESS_BUILD <shmobject> ] [ SHMEM_PROGRESS_SOLVE <shmobject> ] [ MEMORY_SHMEM <shmobject> ]
~~~



List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>
         
List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>
         

##  `FINO_STEP`

Solves the linear eigenvalue problem.

~~~wasora
FINO_STEP [ JUST_BUILD | JUST_SOLVE ] [ DO_NOT_COMPUTE_GRADIENTS | COMPUTE_GRADIENTS ] [ DUMP_FILE_PATH <filepath> ]
~~~






# Variables

##  `displ_max`

The module of the maximum displacement of the elastic problem.



##  `displ_max_x`

The\ $x$ coordinate of the maximum displacement of the elastic problem.



##  `displ_max_y`

The\ $y$ coordinate of the maximum displacement of the elastic problem.



##  `displ_max_z`

The\ $z$ coordinate of the maximum displacement of the elastic problem.



##  `fino_abstol`

Absolute tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html)
Default `1e-50`.



##  `fino_dirichlet_diagonal`

Value that is inserted in the diagonal of the rows
that correspond to Dirichlet boundary conditions.
Default is one, but PETSc internally scales it up
automatically to keep a good condition number.



##  `fino_divtol`

Divergence tolerance,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e+4`.  



##  `fino_gamg_threshold`

Relative threshold to use for dropping edges in aggregation graph for the
[Geometric Algebraic Multigrid Preconditioner](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html)
as passed to PETSc’s
[`PCGAMGSetThreshold`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMGSetThreshold.html).
A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
Default `0.01`.  



##  `fino_iterations`

This variable contains the actual number of iterations used
by the solver. It is set after `FINO_STEP`.



##  `fino_max_iterations`

Number of maximum iterations before diverging,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `10000`.



##  `fino_reltol`

Relative tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
Default `1e-6`.



##  `fino_residual_norm`

This variable contains the residual obtained
by the solver. It is set after `FINO_STEP`.



##  `lambda`


Requested eigenvalue. It is equal to 1.0 until
`FINO_STEP` is executed.  



##  `memory_usage_global`

Maximum resident set size (global memory used), in bytes.



##  `memory_usage_petsc`

Maximum resident set size (memory used by PETSc), in bytes.



##  `memory_use`

Total available memory, in bytes.



##  `petsc_flops`

Number of floating point operations performed by PETSc/SLEPc.



##  `sigma_max`

The maximum von Mises stress\ $\sigma$ of the elastic problem.



##  `sigma_max_x`

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  `sigma_max_y`

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  `sigma_max_z`

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.



##  `time_cpu_build`

CPU time insumed to build the problem matrices, in seconds.



##  `time_cpu_solve`

CPU time insumed to solve the eigen-problem, in seconds.



##  `time_petsc_build`

CPU time insumed by PETSc to build the problem matrices, in seconds.



##  `time_petsc_solve`

CPU time insumed by PETSc to solve the eigen-problem, in seconds.



##  `time_wall_build`

Wall time insumed to build the problem matrices, in seconds.



##  `time_wall_solve`

Wall time insumed to solve the eigen-problem, in seconds.



##  `time_wall_total`

Wall time insumed to initialize, build and solve, in seconds.
CPU time insumed to initialize, build and solve, in seconds.
CPU time insumed by PETSc to initialize, build and solve, in seconds.



##  `T_max`

The maximum temperature\ $T_\text{max}$ of the thermal problem.



##  `T_min`

The minimum temperature\ $T_\text{min}$ of the thermal problem.



##  `u_at_displ_max`

The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.



##  `u_at_sigma_max`

The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  `v_at_displ_max`

The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.



##  `v_at_sigma_max`

The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.



##  `w_at_displ_max`

The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.



##  `w_at_sigma_max`

The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.






