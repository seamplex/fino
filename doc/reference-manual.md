## Fino keywords

###  FINO_LINEARIZE

> Performs stress linearization according to ASME VII-Sec 5 over a
Stress Classification Line

~~~wasora
FINO_LINEARIZE { PHYSICAL_GROUP <physical_group> | START_POINT <x1> <y1> <z1> END_POINT <x2> <y2> <z2> }
 [ FILE <file_id> | FILE_PATH <file_path> ]
 [ TOTAL { vonmises tresca | tresca | principal1 | principal2 | principal3 }
 [ M <variable> ]
 [ MB <variable> ]
 [ PEAK <variable> ]

~~~


The Stress Classification Line (SCL) may be given either as a one-dimensional physical entity
in the mesh or as the (continuous) spatial coordinates of two end-points.
If the SCL is given as a `PHYSICAL_GROUP`, the entity should be one-dimensional (i.e a line)
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

###  FINO_PROBLEM

> Sets the problem type that Fino has to solve.      

~~~wasora
FINO_PROBLEM [ mechanical | thermal | modal ]
 [ AXISYMMETRIC | PLANE_STRESS | PLANE_STRAIN ] [ SYMMETRY_AXIS { x | y } ] [ LINEAR | NON_LINEAR ]
 [ DIMENSIONS <expr> ] [ MESH <identifier> ] 
 [ N_MODES <expr> ] 

~~~


 * `mechanical` (or `elastic` or `break`, default) solves the mechanical elastic problem (default).
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

###  FINO_REACTION

> Computes the reaction at the selected physical entity.

~~~wasora
FINO_REACTION PHYSICAL_GROUP <physical_group> [ RESULT { <variable> | <vector> } ]
~~~


The result is stored in the variable or vector provided, depending on the number of degrees of freedoms of the problem.
If the object passed as `RESULT` does not exist, an appropriate object (scalar variable or vector) is created.
For the elastic problem, the components of the total reaction force are stored in the result vector.
For the thermal problem, the total power passing through the entity is computed as an scalar.

###  FINO_SOLVER

> Sets options related to the solver and the computation of gradients.

~~~wasora
FINO_SOLVER [ PROGRESS_ASCII ]
 [ PC_TYPE { gamg | mumps | lu | hypre | sor | bjacobi | cholesky | ... } ]
 [ KSP_TYPE { gmres | mumps | bcgs | bicg | richardson | chebyshev | ... } ]
 [ SNES_TYPE { newtonls | newtontr | nrichardson | ngmres | qn | ngs | ... } ]
 [ GRADIENT { gauss | nodes | none } ]
 [ GRADIENT_HIGHER { average | nodes } ]
 [ SMOOTH { always | never | material } ]
 [ ELEMENT_WEIGHT { volume_times_quality volume | quality | flat } ]

~~~


If the keyword `PROGRESS_ASCII` is given, three ASCII lines will show in the terminal the
progress of the ensamble of the stiffness matrix (or matrices), the solution of the system of equations
and the computation of gradients (stresses).
The preconditioner, linear and non-linear solver might be any of those available in PETSc:

 * List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>.
 * List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>.
 * List of `SNES_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESType.html>.

If either `PC_TYPE` or `KSP_TYPE` is set to `mumps` (and PETSc is compiled with MUMPS support) then this direct solver is used.
For the mechanical problem, the default is to use GAMG as the preconditioner and PETSc’s default solver (GMRES).
For the thermal problem, the default is to use the default PETSc settings.
For the modal problem, the default is to use the default SLEPc settings.
The `GRADIENT` keyword controls how the derivatives (i.e. strains) at the first-order nodes
are to be computed out of the primary unknowns (i.e. displacements).

 * `gauss` (default) computes the derivatives at the gauss points and the extrapolates the values to the nodes
 * `nodes` computes the derivatives direcetly at the nodes
 * `none` does not compute any derivative at all

The way derivatives are computed at high-order nodes (i.e. those at the middle of edges or faces)
is controlled with `GRADIENT_HIGHER`:

 * `average` (default) assigns the plain average of the first-order nodes that surrond each high-order node
 * `none` computes the derivatives at the location of the high-order nodes

The keyword `SMOOTH` controls how the gradient-based functions (i.e. strains, stresses, etc) are
smoothed---or not---to obtain nodal values out of data which primarily comes from element-wise evaluations at the Gauss points.

 * `always` (default) computes a single value for each node by averaging the contributions of individual elements.
 * `never` keeps the contribution of each individual element separate.
This option implies that the output mesh is different from the input mesh as each element now
has a “copy” of the original shared nodes.
 * `material` averages element contribution only for those elements that belong to the same material (i.e. physical group).
As with `never`, a new output mesh is created where the nodes are duplicated even for those elements which belong to the same physical group.

The way individual contributions of different elements to the same node are averaged is controlled by `ELEMENT_WEIGHT`:

 * `volume_times_quality` (default) weights each element by the product of its volume times its quality
 * `volume` weights each element by the its volume
 * `quality` weights each element by the its quality
 * `flat` performs plain averages (i.e. the same weight for all elements)


###  FINO_STEP

> Ask Fino to solve the problem and advance one step.

~~~wasora
FINO_STEP [ JUST_BUILD | JUST_SOLVE ]
~~~


The location of the `FINO_STEP` keyword within the input file marks the logical location where
the problem is solved and the result functions (displacements, temperatures, stresses, etc.) are available
for output or further computation.





## Mesh keywords

###  MATERIAL

> 
~~~wasora
MATERIAL <name> [ MESH <name> ] [ PHYSICAL_GROUP <name_1> [ PHYSICAL_GROUP <name_2> [ ... ] ] ] [ <property_name_1> <expr_1> [ <property_name_2> <expr_2> [ ... ] ] ]
~~~



###  MESH

> Reads an unstructured mesh from an external file in MSH, VTK or FRD format.

~~~wasora
MESH [ NAME <name> ] { FILE <file_id> | FILE_PATH <file_path> } [ DIMENSIONS <num_expr> ]
 [ SCALE <expr> ] [ OFFSET <expr_x> <expr_y> <expr_z> ]
 [ READ_SCALAR <name_in_mesh> AS <function_name> ] [...]
 [ READ_FUNCTION <function_name> ] [...]
~~~


If there will be only one mesh in the input file, the `NAME` is optional.
Yet it might be needed in cases where there are many meshes and one needs to refer to a particular mesh,
such as in `MESH_POST` or `MESH_INTEGRATE`.
When solving PDEs (such as in Fino or milonga), the first mesh is the problem mesh.
Either a file identifier (defined previously with a `FILE` keyword) or a file path should be given.
The format is read from the extension, which should be either

 * `.msh` [Gmsh ASCII format](http:\//gmsh.info/doc/texinfo/gmsh.html#MSH-file-format), versions 2.2, 4.0 or 4.1
 * `.vtk` [ASCII legacy VTK](https:\//lorensen.github.io/VTKExamples/site/VTKFileFormats/)
 * `.frd` [CalculiX’s FRD ASCII output](https:\//web.mit.edu/calculix_v2.7/CalculiX/cgx_2.7/doc/cgx/node4.html))

Note than only MSH is suitable for defining PDE domains, as it is the only one that provides information about physical groups.
The spatial dimensions should be given with `DIMENSION`. If material properties are uniform and
given with variables, the dimensions are not needed and will be read from the file.
But if spatial functions are needed (either for properties or read from the mesh file), an
explicit value for the mesh dimensions is needed.
If either `SCALE` or `OFFSET` are given, the node position if first shifted and then scaled by the provided amounts.
For each `READ_SCALAR` keyword, a point-wise defined function of space named `<function_name>` is
defined and filled with the scalar data named `<name_in_mesh>` contained in the mesh file.
The `READ_FUNCTION` keyword is a shortcut when the scalar name and the to-be-defined function are the same.
If no `NAME` is given, the first mesh to be defined is called `first`.

###  MESH_FILL_VECTOR

> Fills the elements of a vector with data evaluated at the nodes or the cells of a mesh.

~~~wasora
MESH_FILL_VECTOR VECTOR <vector> { FUNCTION <function> | EXPRESSION <expr> } 
 [ MESH <name> ] [ NODES | CELLS ]
~~~


The vector to be filled needs to be already defined and to have the appropriate size,
either the number of nodes or cells of the mesh depending on `NODES` or `CELLS` (default is nodes).
The elements of the vectors will be either the `FUNCTION` or the `EXPRESSION` of $x$, $y$ and $z$
evaluated at the nodes or cells of the provided mesh.
If there is more than one mesh, the name has to be given.

###  MESH_FIND_MINMAX

> Finds absolute extrema of a function or expression within a mesh-based domain.

~~~wasora
MESH_FIND_MINMAX { FUNCTION <function> | EXPRESSION <expr> }
 [ MESH <name> ] [ NODES | CELLS ]
 [ MIN <variable> ] [ I_MIN <variable> ] [ X_MIN <variable> ] [ Y_MIN <variable> ] [Z_MIN <variable> ]
 [ MAX <variable> ] [ I_MAX <variable> ] [ X_MAX <variable> ] [ Y_MAX <variable> ] [Z_MAX <variable> ]

~~~


Either a `FUNCTION` or an `EXPRESSION` should be given.
In the first case, just the function name is expected (i.e. not its arguments).

###  MESH_INTEGRATE

> Performs a spatial integration of a function or expression over a mesh.

~~~wasora
MESH_INTEGRATE { FUNCTION <function> | EXPRESSION <expr> }
 [ MESH <mesh_identifier> ] [ OVER <physical_group> ] [ NODES | CELLS ]
 RESULT <variable>

~~~


The integrand may be either a `FUNCTION` or an `EXPRESSION`.
In the first case, just the function name is expected (i.e. not its arguments).
In the second case, a full algebraic expression including the arguments is expected.
If the expression is just `1` then the volume (or area or length) of the domain is computed.
Note that arguments ought to be `x`, `y` and/or `z`.
If there are more than one mesh defined, an explicit one has to be given with `MESH`.
By default the integration is performed over the highest-dimensional elements of the mesh.
If the integration is to be carried out over just a physical group, it has to be given in `OVER`.
Either `NODES` or `CELLS` define how the integration is to be performed.
In the first case a the integration is performed using the Gauss points and weights associated to each element type.
In the second case, the integral is computed as the sum of the product of the function evaluated at the center of each cell (element) and the cell’s volume.
The scalar result of the integration is stored in the variable given by `RESULT`.
If the variable does not exist, it is created.
In the second case, a full algebraic expression including the arguments is expected.

###  MESH_MAIN

> 
~~~wasora
MESH_MAIN [ <name> ]
~~~



###  MESH_POST

> 
~~~wasora
MESH_POST [ MESH <mesh_identifier> ] { FILE <name> | FILE_PATH <file_path> } [ NO_MESH ] [ FORMAT { gmsh | vtk } ] [ CELLS | ] NODES ] [ NO_PHYSICAL_NAMES ] [ VECTOR <function1_x> <function1_y> <function1_z> ] [...] [ <scalar_function_1> ] [ <scalar_function_2> ] ...
~~~



###  PHYSICAL_GROUP

> Defines a physical group of elements within a mesh file.

~~~wasora
PHYSICAL_GROUP <name> [ MESH <name> ] [ DIMENSION <expr> ]
 [ MATERIAL <name> ]
 [ BC <bc_1> <bc_2> ... ]

~~~


A name is mandatory for each physical group defined within the input file.
If there is no physical group with the provided name in the mesh, this instruction makes no effect.
If there are many meshes, an explicit mesh can be given with `MESH`.
Otherwise, the physical group is defined on the main mesh.
An explicit dimension of the physical group can be provided with `DIMENSION`.
For volumetric elements, physical groups can be linked to materials using `MATERIAL`.
Note that if a material is created with the same name as a physical group in the mesh,
they will be linked automatically. The `MATERIAL` keyword in `PHYSICAL_GROUP` is used
to link a physical group in a mesh file and a material in the wasora input file with
different names.
For non-volumetric elements, boundary conditions can be assigned by using the `BC` keyword.
This should be the last keyword of the line, and any token afterwards is treated
specially by the underlying solver (i.e. Fino or milonga).

###  PHYSICAL_PROPERTY

> 
~~~wasora
PHYSICAL_PROPERTY <name> [ <material_name1> <expr1> [ <material_name2> <expr2> ] ... ]
~~~







## Special input distributions

TBD.


## Boundary conditions

TBD.


## Result functions

TBD.


## Wasora keywords

###  =

> Assign an expression to a variable, a vector or a matrix.

~~~wasora
<var>[ [<expr_tmin>, <expr_tmax>] | 
<expr_t> ] = <expr> <vector>(<expr_i>)[<expr_i_min, expr_i_max>] [ [<expr_tmin>, <expr_tmax>] | 
<expr_t> ] = <expr> <matrix>(<expr_i>,<expr_j>)[<expr_i_min, expr_i_max; expr_j_min, expr_j_max>] [ [<expr_tmin>, <expr_tmax>] | 
<expr_t> ] = <expr>
~~~



###  _.=

> Add an equation to the DAE system to be solved in the phase space spanned by `PHASE_SPACE`.

~~~wasora
{ 0[(i[,j]][<imin:imax[;jmin:jmax]>] | <expr1> } .= <expr2>
~~~



###  ABORT

> Catastrophically abort the execution and quit wasora.

~~~wasora
ABORT
~~~


Whenever the instruction `ABORT` is executed, wasora quits without closing files
or unlocking shared memory objects. The objective of this instruction is, as
illustrated in the examples, either to debug complex input files and check the
values of certain variables or to conditionally abort the execution using `IF`
clauses.

###  ALIAS

> Define a scalar alias of an already-defined indentifier.

~~~wasora
ALIAS { <new_var_name> IS <existing_object> | <existing_object> AS <new_name> }
~~~


The existing object can be a variable, a vector element or a matrix element.
In the first case, the name of the variable should be given as the existing object.
In the second case, to alias the second element of vector `v` to the new name `new`, `v(2)` should be given as the existing object.
In the third case, to alias second element (2,3) of matrix `M` to the new name `new`, `M(2,3)` should be given as the existing object.

###  CALL

> Call a previously dynamically-loaded user-provided routine.

~~~wasora
CALL <name> [ expr_1 expr_2 ... expr_n ]
~~~



###  CLOSE

> Explicitly close an already-`OPEN`ed file.

~~~wasora
CLOSE
~~~



###  CONST

> Mark a scalar variable, vector or matrix as a constant.

~~~wasora
CONST name_1 [ <name_2> ] ... [ <name_n> ]
~~~



###  DEFAULT_ARGUMENT_VALUE

> Give a default value for an optional commandline argument.

~~~wasora
DEFAULT_ARGUMENT_VALUE <constant> <string>
~~~


If a `$n` construction is found in the input file but the
commandline argument was not given, the default behavior is to
fail complaining that an extra argument has to be given in the
commandline. With this keyword, a default value can be assigned if
no argument is given, thus avoiding the failure and making the argument
optional.

###  DIFFERENTIAL

> Explicitly mark variables, vectors or matrices as “differential” to compute intial conditions of DAE systems.

~~~wasora
DIFFERENTIAL { <var_1> <var_2> ... | <vector_1> <vector_2> ... | <matrix_1> <matrix_2> ... }
~~~



###  DO_NOT_EVALUATE_AT_PARSE_TIME

> Ask wasora not to evaluate assignments at parse time.

~~~wasora
DO_NOT_EVALUATE_AT_PARSE_TIME
~~~



###  FILE

> Define a file, either as input or as output, for further usage.

~~~wasora
< FILE | OUTPUT_FILE | INPUT_FILE > <name> <printf_format> [ expr_1 expr_2 ... expr_n ] [ INPUT | OUTPUT | MODE <fopen_mode> ] [ OPEN | DO_NOT_OPEN ]
~~~



###  FIT

> Fit a function of one or more arguments to a set of data.

~~~wasora
FIT <function_to_be_fitted> TO <function_with_data> VIA <var_1> <var_2> ... <var_n> [ GRADIENT <expr_1> <expr_2> ... <expr_n> ] [ RANGE_MIN <expr_1> <expr_2> ... <expr_n> ] [ RANGE_MAX <expr_1> <expr_2> ... <expr_n> ] [ DELTAEPSREL <expr> ] [ DELTAEPSABS <expr> ] [ MAX_ITER <expr> ] [ VERBOSE ] [ RERUN | DO_NOT_RERUN ]
~~~


The function with the data has to be point-wise defined.
The function to be fitted hast to be parametrized with at least one of the variables provided after the `VIA` keyword.
Only the names of the functions have to be given.
Both functions have to have the same number of arguments.
The initial guess of the solution is given by the initial value of the variables listed in the `VIA` keyword.
Analytical expressions for the gradient of the function to be fitted with respect to the parameters to be fitted can be optionally given with the `GRADIENT` keyword.
If none is provided, the gradient will be computed numerically using finite differences.
A range over which the residuals are to be minimized can be given with `RANGE_MIN` and `RANGE_MAX`.
For multidimensional fits, the range is an hypercube.
If no range is given, all the definition points of the function witht the data are used for the fit.
Convergence can be controlled by given the relative and absolute tolreances with
`DELTAEPSREL` (default 1e-4) and `DELTAEPSABS` (default 1e-6),
and with the maximum number of iterations `MAX_ITER` (default 100).
If the optional keyword `VERBOSE` is given, some data of the intermediate steps is written in the standard output.

###  FUNCTION

> Define a function of one or more variables.

~~~wasora
FUNCTION <name>(<var_1>[,var2,...,var_n]) { [ = <expr> | FILE_PATH <file_path> | ROUTINE <name> | | MESH <name> { DATA <new_vector_name> | VECTOR <existing_vector_name> } { NODES | CELLS } | [ VECTOR_DATA <vector_1> <vector_2> ... <vector_n> <vector_n+1> ] } [COLUMNS <expr_1> <expr_2> ... <expr_n> <expr_n+1> ] [ INTERPOLATION { linear | polynomial | spline | spline_periodic | akima | akima_periodic | steffen | nearest | shepard | shepard_kd | bilinear } ] [ INTERPOLATION_THRESHOLD <expr> ] [ SHEPARD_RADIUS <expr> ] [ SHEPARD_EXPONENT <expr> ] [ SIZES <expr_1> <expr_2> ... <expr_n> ] [ X_INCREASES_FIRST <expr> ] [ DATA <num_1> <num_2> ... <num_N> ]
~~~


The number of variables $n$ is given by the number of arguments given between parenthesis after the function name.
The arguments are defined as new variables if they had not been already defined as variables.
If the function is given as an algebraic expression, the short-hand operator `:=` can be used.
That is to say, `FUNCTION f(x) = x^2` is equivalent to `f(x) := x^2`.
If a `FILE_PATH` is given, an ASCII file containing at least $n+1$ columns is expected.
By default, the first $n$ columns are the values of the arguments and the last column is the value of the function at those points.
The order of the columns can be changed with the keyword `COLUMNS`, which expects $n+1$ expressions corresponding to the column numbers.
A function of type `ROUTINE` calls an already-defined user-provided routine using the `CALL` keyword and passes the values of the variables in each required evaluation as a `double *` argument.
If `MESH` is given, the definition points are the nodes or the cells of the mesh.
The function arguments should be $(x)$, $(x,y)$ or $(x,y,z)$ matching the dimension the mesh.
If the keyword `DATA` is used, a new empty vector of the appropriate size is defined.
The elements of this new vector can be assigned to the values of the function at the $i$-th node or cell.
If the keyword `VECTOR` is used, the values of the dependent variable are taken to be the values of the already-existing vector.
Note that this vector should have the size of the number of nodes or cells the mesh has, depending on whether `NODES` or `CELLS` is given.
If `VECTOR_DATA` is given, a set of $n+1$ vectors of the same size is expected.
The first $n+1$ correspond to the arguments and the last one is the function value.
Interpolation schemes can be given for either one or multi-dimensional functions with `INTERPOLATION`.
Available schemes for $n=1$ are:

 * linear
 * polynomial, the grade is equal to the number of data minus one
 * spline, cubic (needs at least 3 points)
 * spline_periodic
 * akima (needs at least 5 points)
 * akima_periodic (needs at least 5 points)
 * steffen, always-monotonic splines-like (available only with GSL >= 2.0)

Default interpolation scheme for one-dimensional functions is `(*gsl_interp_linear)`.

Available schemes for $n>1$ are:

 * nearest, $f(\vec{x})$ is equal to the value of the closest definition point
 * shepard, [inverse distance weighted average definition points](https://en.wikipedia.org/wiki/Inverse_distance_weighting) (might lead to inefficient evaluation)
 * shepard_kd, [average of definition points within a kd-tree](https://en.wikipedia.org/wiki/Inverse_distance_weighting#Modified_Shepard&#39;s_method) (more efficient evaluation provided `SHEPARD_RADIUS` is set to a proper value)
 * bilinear, only available if the definition points configure an structured hypercube-like grid. If $n>3$, `SIZES` should be given.

For $n>1$, if the euclidean distance between the arguments and the definition points is smaller than `INTERPOLATION_THRESHOLD`, the definition point is returned and no interpolation is performed.
Default value is square root of `9.5367431640625e-07`.
The initial radius of points to take into account in `shepard_kd` is given by `SHEPARD_RADIUS`. If no points are found, the radius is double until at least one definition point is found.
The radius is doubled until at least one point is found.
Default is `1.0`.
The exponent of the `shepard` method is given by `SHEPARD_EXPONENT`.
Default is `2`.
When requesting `bilinear` interpolation for $n>3$, the number of definition points for each argument variable has to be given with `SIZES`,
and wether the definition data is sorted with the first argument changing first (`X_INCREASES_FIRST` evaluating to non-zero) or with the last argument changing first (zero).
The function can be pointwise-defined inline in the input using `DATA`. This should be the last keyword of the line, followed by $N=k\cdot (n+1)$ expresions giving $k$ definition points: $n$ arguments and the value of the function.
Multiline continuation using brackets `{` and `}` can be used for a clean data organization. See the examples.

###  HISTORY

> Record the time history of a variable as a function of time.

~~~wasora
HISTORY <variable> <function>
~~~



###  IF

> Begin a conditional block.

~~~wasora
IF expr
<block_of_instructions_if_expr_is_true>
[ ELSE ]
[block_of_instructions_if_expr_is_false]
ENDIF
~~~



###  IMPLICIT

> Define whether implicit declaration of variables is allowed or not.

~~~wasora
IMPLICIT { NONE | ALLOWED }
~~~


By default, wasora allows variables (but not vectors nor matrices) to be
implicitly declared. To avoid introducing errors due to typos, explicit
declaration of variables can be forced by giving `IMPLICIT NONE`.
Whether implicit declaration is allowed or explicit declaration is required
depends on the last `IMPLICIT` keyword given, which by default is `ALLOWED`.

###  INCLUDE

> Include another wasora input file.

~~~wasora
INCLUDE <file_path> [ FROM <num_expr> ] [ TO <num_expr> ]
~~~


Includes the input file located in the string `file_path` at the current location.
The effect is the same as copying and pasting the contents of the included file
at the location of the `INCLUDE` keyword. The path can be relative or absolute.
Note, however, that when including files inside `IF` blocks that instructions are
conditionally-executed but all definitions (such as function definitions) are processed at
parse-time independently from the evaluation of the conditional.
The optional `FROM` and `TO` keywords can be used to include only portions of a file.

###  INITIAL_CONDITIONS_MODE

> Define how initial conditions of DAE problems are computed.

~~~wasora
INITIAL_CONDITIONS_MODE { AS_PROVIDED | FROM_VARIABLES | FROM_DERIVATIVES }
~~~


In DAE problems, initial conditions may be either:

 * equal to the provided expressions (`AS_PROVIDED`)
 * the derivatives computed from the provided phase-space variables (`FROM_VARIABLES`)
 * the phase-space variables computed from the provided derivatives (`FROM_DERIVATIVES`)

In the first case, it is up to the user to fulfill the DAE system at\ $t = 0$.
If the residuals are not small enough, a convergence error will occur.
The `FROM_VARIABLES` option means calling IDA’s `IDACalcIC` routine with the parameter `IDA_YA_YDP_INIT`.
The `FROM_DERIVATIVES` option means calling IDA’s `IDACalcIC` routine with the parameter IDA_Y_INIT.
Wasora should be able to automatically detect which variables in phase-space are differential and
which are purely algebraic. However, the `DIFFERENTIAL` keyword may be used to explicitly define them.
See the (SUNDIALS documentation)[https://computation.llnl.gov/casc/sundials/documentation/ida_guide.pdf] for further information.

###  LOAD_PLUGIN

> Load a wasora plug-in from a dynamic shared object.

~~~wasora
LOAD_PLUGIN { <file_path> | <plugin_name> }
~~~


A wasora plugin in the form of a dynamic shared object (i.e. `.so`) can be loaded
either with the `LOAD_PLUGIN` keyword or from the command line with the `-p` option.
Either a file path or a plugin name can be given. The following locations are tried:

 * the current directory `./`
 * the parent directory `../`
 * the user’s `LD_LIBRARY_PATH`
 * the cache file `/etc/ld.so.cache`
 * the directories `/lib`, `/usr/lib`, `/usr/local/lib`

If a wasora plugin was compiled and installed following the `make install` procedure,
the plugin should be loaded by just passing the name to `LOAD_PLUGIN`.

###  LOAD_ROUTINE

> Load one or more routines from a dynamic shared object.

~~~wasora
LOAD_ROUTINE <file_path> <routine_1> [ <routine_2> ... <routine_n> ]
~~~



###  M4

> Call the `m4` macro processor with definitions from wasora variables or expressions.

~~~wasora
M4 { INPUT_FILE <file_id> | FILE_PATH <file_path> } { OUTPUT_FILE <file_id> | OUTPUT_FILE_PATH <file_path> } [ EXPAND <name> ] ... } [ MACRO <name> [ <format> ] <definition> ] ... }
~~~



###  MATRIX

> Define a matrix.

~~~wasora
MATRIX <name> ROWS <expr> COLS <expr> [ DATA num_expr_1 num_expr_2 ... num_expr_n ]
~~~



###  MINIMIZE

> Find the combination of arguments that give a (relative) minimum of a function, i.e. run an optimization problem.

~~~wasora
MINIMIZE <function> <function> [ METHOD { conjugate_fr | conjugate_pr | vector_bfgs2 | vector_bfgs | steepest_descent | nmsimplex2 | nmsimplex | nmsimplex2rand } [ GRADIENT <expr_1> <expr_2> ... <expr_n> ] [ GUESS <expr_1> <expr_2> ... <expr_n> ] [ MIN <expr_1> <expr_2> ... <expr_n> ] [ MAX <expr_1> <expr_2> ... <expr_n> ] [ STEP <expr_1> <expr_2> ... <expr_n> ] [ VERBOSE ] [ NORERUN ] [ MAX_ITER <expr> ] [ TOL <expr> ] [ GRADTOL <expr> ]
~~~



###  PARAMETRIC

> Systematically sweep a zone of the parameter space, i.e. perform a parametric run.

~~~wasora
PARAMETRIC <var_1> [ ... <var_n> ] [ TYPE { linear | logarithmic | random | gaussianrandom | sobol | niederreiter | halton | reversehalton } ] [ MIN <num_expr_1> ... <num_expr_n> ] [ MAX <num_expr_1> ... <num_expr_n> ] [ STEP <num_expr_1> ... <num_expr_n> ] [ NSTEPS <num_expr_1> ... <num_expr_n> ] [ OUTER_STEPS <num_expr> ] [ MAX_DAUGHTERS <num_expr> ] [ OFFSET <num_expr> ] [ ADIABATIC ]
~~~



###  PHASE_SPACE

> Define which variables, vectors and/or matrices belong to the phase space of the DAE system to be solved.

~~~wasora
PHASE_SPACE { <vars> | <vectors> | <matrices> }
~~~



###  PRINT

> Print plain-text and/or formatted data to the standard output or into an output file.

~~~wasora
PRINT [ FILE <file_id> | FILE_PATH <file_path> ] [ NONEWLINE ] [ SEP <string> ] [ NOSEP ] [ HEADER ] [ SKIP_STEP <expr> ] [ SKIP_STATIC_STEP <expr> ] [ SKIP_TIME <expr> ] [ SKIP_HEADER_STEP <expr> ] [ <object_1> <object_2> ... <object_n> ] [ TEXT <string_1> ... TEXT <string_n> ]
~~~


Each argument `object` that is not a keyword is expected to be part of the output, can be either a matrix, a vector, an scalar algebraic expression.
If the given object cannot be solved into a valid matrix, vector or expression, it is treated as a string literal if `IMPLICIT` is `ALLOWED`, otherwise a parser error is raised.
To explicitly interpret an object as a literal string even if it resolves to a valid numerical expression, it should be prefixed with the `TEXT` keyword.
Hashes `#` appearing literal in text strings have to be quoted to prevent the parser to treat them as comments within the wasora input file and thus ignoring the rest of the line.
Whenever an argument starts with a porcentage sign `%`, it is treated as a C `printf`-compatible format definition and all the objects that follow it are printed using the given format until a new format definition is found.
The objects are treated as double-precision floating point numbers, so only floating point formats should be given. The default format is `"%g"`.
Matrices, vectors, scalar expressions, format modifiers and string literals can be given in any desired order, and are processed from left to right.
Vectors are printed element-by-element in a single row. See `PRINT_VECTOR` to print vectors column-wise.
Matrices are printed element-by-element in a single line using row-major ordering if mixed with other objects but in the natural row and column fashion if it is the only given object.
If the `FILE` keyword is not provided, default is to write to stdout.
If the `NONEWLINE` keyword is not provided, default is to write a newline `\n` character after all the objects are processed.
The `SEP` keywords expects a string used to separate printed objects, the default is a tab 'DEFAULT_PRINT_SEPARATOR' character.
Use the `NOSEP` keyword to define an empty string as object separator.
If the `HEADER` keyword is given, a single line containing the literal text
given for each object is printed at the very first time the `PRINT` instruction is
processed, starting with a hash `#` character.
If the `SKIP_STEP` (`SKIP_STATIC_STEP`)keyword is given, the instruction is processed only every
the number of transient (static) steps that results in evaluating the expression,
which may not be constant. By default the `PRINT` instruction is evaluated every
step. The `SKIP_HEADER_STEP` keyword works similarly for the optional `HEADER` but
by default it is only printed once. The `SKIP_TIME` keyword use time advancements
to choose how to skip printing and may be useful for non-constant time-step problems.

###  PRINT_FUNCTION

> Print one or more functions as a table of values of dependent and independent variables.

~~~wasora
PRINT_FUNCTION <function_1> [ { function_2 | expr_1 } ... { function_n | expr_n-1 } ] [ FILE <file_id> | FILE_PATH <file_path> ] [ HEADER ] [ MIN <expr_1> <expr_2> ... <expr_m> ] [ MAX <expr_1> <expr_2> ... <expr_m> ] [ STEP <expr_1> <expr_2> ... <expr_m> ] [ NSTEPs <expr_1> <expr_2> ... <expr_m> ] [ FORMAT <print_format> ] [ PHYSICAL_ENTITY <name> ]
~~~



###  PRINT_VECTOR

> Print the elements of one or more vectors.

~~~wasora
PRINT_VECTOR [ FILE <file_id> ] FILE_PATH <file_path> ] [ { VERTICAL | HORIZONTAL } ] [ ELEMS_PER_LINE <expr> ] [ FORMAT <print_format> ] <vector_1> [ vector_2 ... vector_n ]
~~~



###  READ

> Read data (variables, vectors o matrices) from files or shared-memory segments.

~~~wasora
[ READ | WRITE ] [ SHM <name> ] [ { ASCII_FILE_PATH | BINARY_FILE_PATH } <file_path> ] [ { ASCII_FILE | BINARY_FILE } <identifier> ] [ IGNORE_NULL ] [ object_1 object_2 ... object_n ]
~~~



###  SEMAPHORE

> Perform either a wait or a post operation on a named shared semaphore.

~~~wasora
[ SEMAPHORE | SEM ] <name> { WAIT | POST }
~~~



###  SHELL

> Execute a shell command.

~~~wasora
SHELL <print_format> [ expr_1 expr_2 ... expr_n ]
~~~



###  SOLVE

> Solve a non-linear system of\ $n$ equations with\ $n$ unknowns.

~~~wasora
SOLVE <n> UNKNOWNS <var_1> <var_2> ... <var_n> RESIDUALS <expr_1> <expr_2> ... <expr_n> ] GUESS <expr_1> <expr_2> ... <expr_n> ] [ METHOD { dnewton | hybrid | hybrids | broyden } ] [ EPSABS <expr> ] [ EPSREL <expr> ] [ MAX_ITER <expr> ] [ VERBOSE ]
~~~



###  TIME_PATH

> Force transient problems to pass through specific instants of time.

~~~wasora
TIME_PATH <expr_1> [ <expr_2> [ ... <expr_n> ] ]
~~~


The time step `dt` will be reduced whenever the distance between
the current time `t` and the next expression in the list is greater
than `dt` so as to force `t` to coincide with the expressions given.
The list of expresssions should evaluate to a sorted list of values.

###  VAR

> Define one or more scalar variables.

~~~wasora
VAR <name_1> [ <name_2> ] ... [ <name_n> ]
~~~



###  VECTOR

> Define a vector.

~~~wasora
VECTOR <name> SIZE <expr> [ DATA <expr_1> <expr_2> ... <expr_n> | FUNCTION_DATA <function> ]
~~~



###  VECTOR_SORT

> Sort the elements of a vector using a specific numerical order,
potentially making the same rearrangement of another vector.

~~~wasora
VECTOR_SORT <vector> [ ASCENDING_ORDER | DESCENDING_ORDER ] [ <vector> ]
~~~



###  WRITE

> Write data (variables, vectors o matrices) to files or shared-memory segments.
See the `READ` keyword for usage details.







## Fino variables

###  delta_sigma_max

> 

The uncertainty of the maximum Von Mises stress\ $\sigma$ of the elastic problem.
Not to be confused with the maximum uncertainty of the Von Mises stress.

###  displ_max

> 

The module of the maximum displacement of the elastic problem.

###  displ_max_x

> 

The\ $x$ coordinate of the maximum displacement of the elastic problem.

###  displ_max_y

> 

The\ $y$ coordinate of the maximum displacement of the elastic problem.

###  displ_max_z

> 

The\ $z$ coordinate of the maximum displacement of the elastic problem.

###  fino_abstol

> 

Absolute tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http:
Default `1e-50`.

###  fino_divtol

> 

Divergence tolerance,
as passed to PETSc’s
[`KSPSetTolerances`](http:
Default `1e+4`.

###  fino_gamg_threshold

> 

Relative threshold to use for dropping edges in aggregation graph for the
[Geometric Algebraic Multigrid Preconditioner](http:
as passed to PETSc’s
[`PCGAMGSetThreshold`](http:
A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
Default `0.01`.

###  fino_iterations

> 

This variable contains the actual number of iterations used
by the solver. It is set after `FINO_STEP`.

###  fino_max_iterations

> 

Number of maximum iterations before diverging,
as passed to PETSc’s
[`KSPSetTolerances`](http:
Default `10000`.

###  fino_penalty_weight

> 

The weight $w$ used when setting multi-freedom boundary conditions.
Higher values mean better precision in the constrain but distort
the matrix condition number.
Default is `1e8`.

###  fino_reltol

> 

Relative tolerance of the linear solver,
as passed to PETSc’s
[`KSPSetTolerances`](http:
Default `1e-6`.

###  fino_residual_norm

> 

This variable contains the residual obtained
by the solver. It is set after `FINO_STEP`.

###  lambda

> 

Requested eigenvalue. It is equal to 1.0 until
`FINO_STEP` is executed.

###  memory

> 

Maximum resident set size (global memory used), in bytes.

###  memory_available

> 

Total available memory, in bytes.

###  memory_petsc

> 

Maximum resident set size (memory used by PETSc), in bytes.

###  nodes_rough

> 

The number of nodes of the mesh in `ROUGH` mode.

###  petsc_flops

> 

Number of floating point operations performed by PETSc/SLEPc.

###  sigma_max

> 

The maximum von Mises stress\ $\sigma$ of the elastic problem.

###  sigma_max_x

> 

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.

###  sigma_max_y

> 

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.

###  sigma_max_z

> 

The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.

###  strain_energy

> 

The strain energy stored in the solid, computed as
$1/2 \cdot \vec{u}^T K \vec{u}$
where $\vec{u}$ is the displacements vector and $K$ is the stiffness matrix.

###  time_cpu_build

> 

CPU time insumed to build the problem matrices, in seconds.

###  time_cpu_solve

> 

CPU time insumed to solve the problem, in seconds.

###  time_cpu_stress

> 

CPU time insumed to compute the stresses from the displacements, in seconds.

###  time_petsc_build

> 

CPU time insumed by PETSc to build the problem matrices, in seconds.

###  time_petsc_solve

> 

CPU time insumed by PETSc to solve the eigen-problem, in seconds.

###  time_petsc_stress

> 

CPU time insumed by PETSc to compute the stresses, in seconds.

###  time_wall_build

> 

Wall time insumed to build the problem matrices, in seconds.

###  time_wall_solve

> 

Wall time insumed to solve the problem, in seconds.

###  time_wall_stress

> 

Wall time insumed to compute the stresses, in seconds.

###  time_wall_total

> 

Wall time insumed to initialize, build and solve, in seconds.
CPU time insumed to initialize, build and solve, in seconds.
CPU time insumed by PETSc to initialize, build and solve, in seconds.

###  T_max

> 

The maximum temperature\ $T_\r{max}$ of the thermal problem.

###  T_min

> 

The minimum temperature\ $T_\r{min}$ of the thermal problem.

###  u_at_displ_max

> 

The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.

###  u_at_sigma_max

> 

The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.

###  v_at_displ_max

> 

The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.

###  v_at_sigma_max

> 

The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.

###  w_at_displ_max

> 

The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.

###  w_at_sigma_max

> 

The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.




