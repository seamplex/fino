% Fino reference sheet
% Jeremy Theler

This reference sheet is for [Fino](index.html) v0.5.19-g3537a8b
. 

~~~
$ fino
fino v0.5.19-g3537a8b 
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
FINO_PROBLEM [ BAKE | SHAKE | BREAK ] [ DIMENSIONS <expr> ] [ DEGREES <expr> ] [ MESH <identifier> ] [ N_EIGEN <expr> ] [ UNKNOWNS <name1> <name2> ... <name_degrees> ] [ SHAPE <name> ] [ MATRIX <name> ] [ VECTOR <name> ] [ EIGEN_MATRIX <name> ]
~~~



##  `FINO_SOLVER`

Sets options related to the eigen-solver.

~~~wasora
FINO_SOLVER [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ] [ PC_TYPE { lu | gamg | hypre | sor | bjacobi | cholesky | ... } ] [ SET_NEAR_NULLSPACE ] [ USE_PCSETCOORDINATES ] [ SET_BLOCK_SIZE ]
~~~



List of `KSP_TYPE`s http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html
         
List of `PC_TYPE`s http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
         

##  `FINO_STEP`

Solves the linear eigenvalue problem.

~~~wasora
FINO_STEP [ JUST_BUILD | JUST_SOLVE ] [ DO_NOT_COMPUTE_GRADIENTS | COMPUTE_GRADIENTS ] [ DUMP_FILE_PATH <filepath> ]
~~~






# Variables

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






