/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's parsing routines
 *
 *  Copyright (C) 2015--2020 Seamplex
 *
 *  This file is part of Fino <https://www.seamplex.com/fino>.
 *
 *  fino is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  fino is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with wasora.  If not, see <http://www.gnu.org/licenses/>.
 *------------------- ------------  ----    --------  --     -       -         -
 */
#define _GNU_SOURCE 
#include <stdio.h>
#ifndef _FINO_H
#include "fino.h"
#endif

int plugin_parse_line(char *line) {
  
  char *token;

  if ((token = wasora_get_next_token(line)) != NULL) {

// ---------------------------------------------------------------------
///kw+FINO_PROBLEM+usage FINO_PROBLEM
    if (strcasecmp(token, "FINO_PROBLEM") == 0) {
///kw+FINO_PROBLEM+desc Sets the problem type that Fino has to solve.      
      double xi;
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+FINO_PROBLEM+usage [
        
///kw+FINO_PROBLEM+usage mechanical
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+detail  * `mechanical` (or `elastic` or `break`, default) solves the mechanical elastic problem (default).        
        if (strcasecmp(token, "mechanical") == 0 || strcasecmp(token, "elastic") == 0 || strcasecmp(token, "break") == 0) {
          fino.problem_family = problem_family_mechanical;
///kw+FINO_PROBLEM+usage thermal
          
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+detail  * `thermal` (or `heat` or `bake`) solves the heat conduction problem.
        } else if (strcasecmp(token, "thermal") == 0 || strcasecmp(token, "heat") == 0 || strcasecmp(token, "bake") == 0) {
          fino.problem_family = problem_family_thermal;

///kw+FINO_PROBLEM+usage modal
///kw+FINO_PROBLEM+usage ]@
///kw+FINO_PROBLEM+detail  * `modal` (or `shake`) computes the natural frequencies and oscillation modes.        
        } else if (strcasecmp(token, "modal") == 0 || strcasecmp(token, "shake") == 0) {
#ifndef HAVE_SLEPC
          wasora_push_error_message("modal problems need a Fino binary linked against SLEPc.");
          return WASORA_PARSER_ERROR;
#endif
          fino.problem_family = problem_family_modal;
          if (fino.nev == 0) {
            fino.nev = DEFAULT_NMODES;
          }
          

///kw+FINO_PROBLEM+usage [
///kw+FINO_PROBLEM+usage AXISYMMETRIC
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+detail @
///kw+FINO_PROBLEM+detail If the `AXISYMMETRIC` keyword is given, the mesh is expected to be two-dimensional in the $x$-$y$ plane
///kw+FINO_PROBLEM+detail and the problem is assumed to be axi-symmetric around the axis given by `SYMMETRY_AXIS` (default is $y$). 
        } else if (strcasecmp(token, "AXISYMMETRIC") == 0) {
          fino.problem_kind = problem_kind_axisymmetric;
          if (fino.symmetry_axis == symmetry_axis_none) {
            fino.symmetry_axis = symmetry_axis_y;
          }

///kw+FINO_PROBLEM+usage PLANE_STRESS
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+detail If the problem type is mechanical and the mesh is two-dimensional on the $x$-$y$ plane and no
///kw+FINO_PROBLEM+detail axisymmetry is given, either `PLANE_STRESS` and `PLAIN_STRAIN` can be provided (default is plane stress). 
        } else if (strcasecmp(token, "PLANE_STRESS") == 0) {
          fino.problem_kind = problem_kind_plane_stress;

///kw+FINO_PROBLEM+usage PLANE_STRAIN
///kw+FINO_PROBLEM+usage ]
        } else if (strcasecmp(token, "PLANE_STRAIN") == 0) {
          fino.problem_kind = problem_kind_plane_strain;
          
///kw+FINO_PROBLEM+usage [ SYMMETRY_AXIS { x | y } ]
        } else if (strcasecmp(token, "SYMMETRY_AXIS") == 0) {
          char *keywords[] = { "x", "y" };
          int values[] = {symmetry_axis_x, symmetry_axis_y, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.symmetry_axis));

///kw+FINO_PROBLEM+usage [ LINEAR
        } else if (strcasecmp(token, "LINEAR") == 0) {
          fino.math_type = math_type_linear;

///kw+FINO_PROBLEM+usage | NON_LINEAR ]@
        } else if (strcasecmp(token, "NON_LINEAR") == 0) {
          fino.math_type = math_type_nonlinear;

///kw+FINO_PROBLEM+usage [ QUASISTATIC
        } else if (strcasecmp(token, "QUASISTATIC") == 0) {
          fino.transient_type = transient_type_quasistatic;

///kw+FINO_PROBLEM+usage | TRANSIENT ]@
        } else if (strcasecmp(token, "TRANSIENT") == 0) {
          fino.transient_type = transient_type_transient;
          
///kw+FINO_PROBLEM+detail By default Fino tries to detect wheter the computation should be linear or non-linear.
///kw+FINO_PROBLEM+detail An explicit mode can be set with either `LINEAR` on `NON_LINEAR`.
          
///kw+FINO_PROBLEM+usage [ DIMENSIONS <expr> ]
///kw+FINO_PROBLEM+detail The number of spatial dimensions of the problem needs to be given either with the keyword `DIMENSIONS`
///kw+FINO_PROBLEM+detail or by defining a `MESH` (with an explicit `DIMENSIONS` keyword) before `FINO_PROBLEM`.
        } else if (strcasecmp(token, "DIMENSIONS") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          fino.dimensions = (int)(xi);
          if (fino.dimensions < 1 || fino.dimensions > 3)  {
            wasora_push_error_message("either one, two or three dimensions should be selected instead of '%d'", fino.dimensions);
            return WASORA_PARSER_ERROR;
          }
          
          
///kw+FINO_PROBLEM+usage [ MESH <identifier> ] @
///kw+FINO_PROBLEM+detail If there are more than one `MESH`es define, the one over which the problem is to be solved
///kw+FINO_PROBLEM+detail can be defined by giving the explicit mesh name with `MESH`. By default, the first mesh to be
///kw+FINO_PROBLEM+detail defined in the input file is the one over which the problem is solved.
        } else if (strcasecmp(token, "MESH") == 0) {
          char *mesh_name;
          
          wasora_call(wasora_parser_string(&mesh_name));
          if ((fino.mesh = wasora_get_mesh_ptr(mesh_name)) == NULL) {
            wasora_push_error_message("unknown mesh '%s'", mesh_name);
            free(mesh_name);
            return WASORA_PARSER_ERROR;
          }
          free(mesh_name);
        
///kw+FINO_PROBLEM+usage [ N_MODES <expr> ] @
///kw+FINO_PROBLEM+detail The number of modes to be computed in the modal problem. The default is DEFAULT_NMODES.
        } else if (strcasecmp(token, "N_MODES") == 0 || strcasecmp(token, "N_EIGEN") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          fino.nev = (int)(xi);
          if (fino.nev < 1)  {
            wasora_push_error_message("a positive number of modes should be given instead of '%d'", fino.nev);
            return WASORA_PARSER_ERROR;
          }
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      } 

      // si no nos dieron explicitamente la malla, ponemos la principal
      if (fino.mesh == NULL && (fino.mesh = wasora_mesh.main_mesh) == NULL) {
        wasora_push_error_message("unknown mesh (no MESH keyword before FINO_PROBLEM)", token);
        return WASORA_PARSER_ERROR;
      }

      // por si ya nos dieron mesh, usamos las dimensiones de la malla si no nos las dieron aca
      if (fino.dimensions == 0 && fino.mesh->spatial_dimensions != 0) {
        fino.dimensions = fino.mesh->spatial_dimensions;
      }
      // al reves, si ya nos la dieron, se las damos a la malla
      if (fino.mesh != NULL && fino.mesh->spatial_dimensions == 0 && fino.dimensions != 0) {
        fino.mesh->spatial_dimensions = fino.dimensions;
      }
      
      // si no pusieron tipo, somos mecanicos
      if (fino.problem_family == problem_family_undefined) {
        fino.problem_family = problem_family_mechanical;
      }
      
      // si no pusieron familia, somos full 3d
      if (fino.problem_kind == problem_kind_undefined) {
        fino.problem_kind = problem_kind_full3d;
      }
      
      if (fino.problem_family == problem_family_mechanical ||
          fino.problem_family == problem_family_modal) {
        
        if (fino.problem_kind == problem_kind_full3d) {
          fino.dimensions = 3;
          fino.degrees = 3;
          
        } else if (fino.problem_kind == problem_kind_axisymmetric ||
                   fino.problem_kind == problem_kind_plane_stress ||
                   fino.problem_kind == problem_kind_plane_strain) {
          
          fino.dimensions = 2;
          fino.degrees = 2;
        }
        
        if (fino.problem_family == problem_family_modal) {
          fino.math_type = math_type_eigen;
        }
        
        fino.unknown_name = calloc(fino.degrees, sizeof(char *));
        fino.unknown_name[0] = strdup("u");
        fino.unknown_name[1] = strdup("v");
        if (fino.degrees == 3) {
          fino.unknown_name[2] = strdup("w");
        }
        
      } else if (fino.problem_family == problem_family_thermal) {
        
        fino.degrees = 1;
        fino.unknown_name = calloc(fino.degrees, sizeof(char *));
        fino.unknown_name[0] = strdup("T");
        
      }

      wasora_call(fino_define_functions());
      
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_SOLVER+usage FINO_SOLVER
///kw+FINO_SOLVER+desc Sets options related to the solver and the computation of gradients.
///kw+FINO_SOLVER+detail 
    } else if (strcasecmp(token, "FINO_SOLVER") == 0) {

      while ((token = wasora_get_next_token(NULL)) != NULL) {

        if (strcasecmp(token, "ROUTINE") == 0) {
          if ((token = wasora_get_next_token(NULL)) == NULL) {
            wasora_push_error_message("expected solver name");
            return WASORA_PARSER_ERROR;
          }

          if ((fino.user_provided_linearsolver = wasora_get_loadable_routine(token)) == NULL) {
            wasora_push_error_message("unknown routine '%s'", token);
            return WASORA_PARSER_ERROR;
          }

///kw+FINO_SOLVER+usage [ PROGRESS ]@
///kw+FINO_SOLVER+detail If the keyword `PROGRESS` is given, three ASCII lines will show in the terminal the
///kw+FINO_SOLVER+detail progress of the ensamble of the stiffness matrix (or matrices), the solution of the system of equations
///kw+FINO_SOLVER+detail and the computation of gradients (stresses).
        } else if (strcasecmp(token, "PROGRESS") == 0 || strcasecmp(token, "PROGRESS_ASCII") == 0) {
          fino.progress_ascii = PETSC_TRUE;

///kw+FINO_SOLVER+usage [ PC_TYPE { gamg | mumps | lu | hypre | sor | bjacobi | cholesky | ... } ]@
///kw+FINO_SOLVER+detail The preconditioner, linear and non-linear solver might be any of those available in PETSc:
///kw+FINO_SOLVER+detail @          
///kw+FINO_SOLVER+detail  * List of `PC_TYPE`s <http:/\/www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>.
        } else if (strcasecmp(token, "PC_TYPE") == 0) {
          wasora_call(wasora_parser_string((char **)&fino.pc_type));

///kw+FINO_SOLVER+usage [ KSP_TYPE { gmres | mumps | bcgs | bicg | richardson | chebyshev | ... } ]@
///kw+FINO_SOLVER+detail  * List of `KSP_TYPE`s <http:/\/www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>.
        } else if (strcasecmp(token, "KSP_TYPE") == 0) {
          wasora_call(wasora_parser_string((char **)&fino.ksp_type));
          
///kw+FINO_SOLVER+usage [ TS_TYPE { bdf | arkimex | rosw | glle | beuler | ... } ]@
///kw+FINO_SOLVER+detail  * List of `TS_TYPE`s <http:/\/www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSType.html>.
        } else if (strcasecmp(token, "TS_TYPE") == 0) {
          wasora_call(wasora_parser_string((char **)&fino.ts_type));

///kw+FINO_SOLVER+usage [ SNES_TYPE { newtonls | newtontr | nrichardson | ngmres | qn | ngs | ... } ]@
///kw+FINO_SOLVER+detail  * List of `SNES_TYPE`s <http:/\/www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESType.html>.
///kw+FINO_SOLVER+detail @
        } else if (strcasecmp(token, "SNES_TYPE") == 0) {
          wasora_call(wasora_parser_string((char **)&fino.snes_type));

///kw+FINO_SOLVER+detail If either `PC_TYPE` or `KSP_TYPE` is set to `mumps` (and PETSc is compiled with MUMPS support) then this direct solver is used.
///kw+FINO_SOLVER+detail For the mechanical problem, the default is to use GAMG as the preconditioner and PETSc’s default solver (GMRES).
///kw+FINO_SOLVER+detail For the thermal problem, the default is to use the default PETSc settings.
///kw+FINO_SOLVER+detail For the modal problem, the default is to use the default SLEPc settings.
          
//kw+FINO_SOLVER+usage [ SET_NEAR_NULLSPACE { rigidbody | fino | none } ]@
/*          
        } else if (strcasecmp(token, "SET_NEAR_NULLSPACE") == 0 || strcasecmp(token, "SET_NEAR_NULL_SPACE") == 0) {
          token = wasora_get_next_token(NULL);
          if (token != NULL) {
            if (strcmp(token, "rigidbody") == 0) {
              fino.set_near_nullspace = set_near_nullspace_rigidbody;
            } else if (strcmp(token, "fino") == 0) {
              fino.set_near_nullspace = set_near_nullspace_fino;
            } else if (strcmp(token, "none") == 0) {
              fino.set_near_nullspace = set_near_nullspace_none;
            } else {
              wasora_push_error_message("unknown nullspace method '%s'", token);
              return WASORA_PARSER_ERROR;
            }
          }

//kw+FINO_SOLVER+usage [ DO_NOT_SET_BLOCK_SIZE 
        } else if (strcasecmp(token, "DO_NOT_SET_BLOCK_SIZE") == 0) {
          fino.do_not_set_block_size = 1;
          
//kw+FINO_SOLVER+usage | SET_BLOCK_SIZE ]@
        } else if (strcasecmp(token, "SET_BLOCK_SIZE") == 0) {
          fino.do_not_set_block_size = 0;
*/
          
///kw+FINO_SOLVER+usage [ GRADIENT {
        } else if (strcasecmp(token, "GRADIENT") == 0) {
///kw+FINO_SOLVER+detail The `GRADIENT` keyword controls how the derivatives (i.e. strains) at the first-order nodes
///kw+FINO_SOLVER+detail are to be computed out of the primary unknowns (i.e. displacements).
///kw+FINO_SOLVER+detail @
          char *keywords[] = {
///kw+FINO_SOLVER+usage   gauss |
                         "gauss",
///kw+FINO_SOLVER+detail  * `gauss` (default) computes the derivatives at the gauss points and the extrapolates the values to the nodes
///kw+FINO_SOLVER+usage   nodes |
                         "nodes",
///kw+FINO_SOLVER+detail  * `nodes` computes the derivatives direcetly at the nodes
///kw+FINO_SOLVER+usage   none
                         "none",
///kw+FINO_SOLVER+detail  * `none` does not compute any derivative at all
                         ""};
          int values[] = {gradient_gauss_extrapolated, gradient_at_nodes, gradient_none, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.gradient_evaluation));
///kw+FINO_SOLVER+usage } ]@
///kw+FINO_SOLVER+detail @
          
///kw+FINO_SOLVER+usage [ GRADIENT_HIGHER {
        } else if (strcasecmp(token, "GRADIENT_HIGHER") == 0) {
///kw+FINO_SOLVER+detail The way derivatives are computed at high-order nodes (i.e. those at the middle of edges or faces)
///kw+FINO_SOLVER+detail is controlled with `GRADIENT_HIGHER`:
///kw+FINO_SOLVER+detail @
          char *keywords[] = {
///kw+FINO_SOLVER+usage   average |
                         "average",
///kw+FINO_SOLVER+detail  * `average` (default) assigns the plain average of the first-order nodes that surrond each high-order node
///kw+FINO_SOLVER+usage   nodes
                         "nodes",
///kw+FINO_SOLVER+detail  * `none` computes the derivatives at the location of the high-order nodes
                         ""};
          int values[] = {gradient_average, gradient_actual, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.gradient_highorder_nodes));
///kw+FINO_SOLVER+usage } ]@
///kw+FINO_SOLVER+detail @
          
///kw+FINO_SOLVER+usage [ SMOOTH {
        } else if (strcasecmp(token, "SMOOTH") == 0) {
///kw+FINO_SOLVER+detail The keyword `SMOOTH` controls how the gradient-based functions (i.e. strains, stresses, etc) are
///kw+FINO_SOLVER+detail smoothed---or not---to obtain nodal values out of data which primarily comes from element-wise evaluations at the Gauss points.
///kw+FINO_SOLVER+detail @
          char *keywords[] = {
///kw+FINO_SOLVER+usage   always |
                         "always",
///kw+FINO_SOLVER+detail  * `always` (default) computes a single value for each node by averaging the contributions of individual elements.
///kw+FINO_SOLVER+usage   never |
                         "never",
///kw+FINO_SOLVER+detail  * `never` keeps the contribution of each individual element separate.
///kw+FINO_SOLVER+detail This option implies that the output mesh is different from the input mesh as each element now
///kw+FINO_SOLVER+detail has a “copy” of the original shared nodes.
///kw+FINO_SOLVER+usage   material
                         "material",
///kw+FINO_SOLVER+detail  * `material` averages element contribution only for those elements that belong to the same material (i.e. physical group).
///kw+FINO_SOLVER+detail As with `never`, a new output mesh is created where the nodes are duplicated even for those elements which belong to the same physical group.
                         ""};
          int values[] = {0, 1, 2, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.rough));
          if (fino.rough == 2) {
            fino.rough = 1;
            fino.roughish = 1;
          }
///kw+FINO_SOLVER+detail @
///kw+FINO_SOLVER+usage } ]@
          
///kw+FINO_SOLVER+usage [ ELEMENT_WEIGHT {
///kw+FINO_SOLVER+detail The way individual contributions of different elements to the same node are averaged is controlled by `ELEMENT_WEIGHT`:
///kw+FINO_SOLVER+detail @
        } else if (strcasecmp(token, "ELEMENT_WEIGHT") == 0) {
          char *keywords[] = {
///kw+FINO_SOLVER+usage   volume_times_quality |
                         "volume_times_quality",
///kw+FINO_SOLVER+detail  * `volume_times_quality` (default) weights each element by the product of its volume times its quality
///kw+FINO_SOLVER+usage   volume |
                         "volume",
///kw+FINO_SOLVER+detail  * `volume` weights each element by the its volume
///kw+FINO_SOLVER+usage   quality |
                         "quality",
///kw+FINO_SOLVER+detail  * `quality` weights each element by the its quality
///kw+FINO_SOLVER+usage   flat
                         "flat",
///kw+FINO_SOLVER+detail  * `flat` performs plain averages (i.e. the same weight for all elements)
                         ""};
          int values[] = {gradient_weight_volume_times_quality, gradient_weight_volume, gradient_weight_quality, gradient_weight_flat, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.gradient_element_weight));
///kw+FINO_SOLVER+usage } ]@
///kw+FINO_SOLVER+detail @
          
///kw+FINO_SOLVER+detail Hourglass control can be triggered by providing a positive coefficient in `HOURGLASS`.
///kw+MESH+usage [ HOURGLASS <num_expr> ] @
        } else if (strcasecmp(token, "HOURGLASS") == 0) {
          wasora_call(wasora_parser_expression_in_string(&fino.hourglass_epsilon));
          
//kw+FINO_SOLVER+usage [ GRADIENT_QUALITY_THRESHOLD <expr> ]@
//kw+FINO_SOLVER+detail If the `GRADIENT_QUALITY_THRESHOLD` 
/*          
        } else if (strcasecmp(token, "GRADIENT_QUALITY_THRESHOLD") == 0) {
          wasora_call(wasora_parser_expression_in_string(&fino.gradient_quality_threshold));
*/
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_STEP+usage FINO_STEP
///kw+FINO_STEP+desc Ask Fino to solve the problem and advance one step.
///kw+FINO_STEP+detail The location of the `FINO_STEP` keyword within the input file marks the logical location where      
///kw+FINO_STEP+detail the problem is solved and the result functions (displacements, temperatures, stresses, etc.) are available
///kw+FINO_STEP+detail for output or further computation.
    } else if (strcasecmp(token, "FINO_STEP") == 0) {

      fino_step_t *fino_step = calloc(1, sizeof(fino_step_t));
      instruction_t *instruction;
      
      if (fino.mesh == NULL && (fino.mesh = wasora_mesh.main_mesh) == NULL) {
        wasora_push_error_message("no mesh found! (FINO_STEP before MESH)");
        return WASORA_PARSER_ERROR;
      }

      
      // chequeo de dimensiones
      if (fino.dimensions == 0 && fino.mesh->spatial_dimensions == 0) {
        // defaulteamos a tres
        fino.dimensions = 3;
      }
      // si alguna es cero, la rellenamos con la otra
      if (fino.dimensions == 0) {
        fino.dimensions = fino.mesh->spatial_dimensions;
      } else if (fino.mesh->spatial_dimensions == 0) {
        fino.mesh->spatial_dimensions = fino.dimensions;
      }
      
      // si son diferentes nos quejamos
      if (fino.dimensions != fino.mesh->spatial_dimensions) {
        wasora_push_error_message("inconsistent dimensions (FINO_PROBLEM = %d, MESH = %d)", fino.dimensions, fino.mesh->spatial_dimensions);
        return WASORA_PARSER_ERROR;
      }
      
      // defaulteamos a mechanical
      if (fino.problem_family == problem_family_undefined) {
        fino.problem_family = problem_family_mechanical;
        if (fino.dimensions == 3) {
          fino.problem_kind = problem_kind_full3d;
          fino.degrees = 3;
        } else if (fino.dimensions == 2) {
          fino.problem_kind = problem_kind_plane_stress;
          fino.degrees = 2;
        } else {
          wasora_push_error_message("no explicit problem type given with a one-dimensional mesh. Give me PROBLEM_TYPE.");
          return WASORA_PARSER_ERROR;
        }  
        
        fino.unknown_name = calloc(fino.degrees, sizeof(char *));
        fino.unknown_name[0] = strdup("u");
        fino.unknown_name[1] = strdup("v");
        if (fino.degrees == 3) {
          fino.unknown_name[2] = strdup("w");
        }  
      }
      
      // si nadie dijo nada tenemos un solo grado de libertado
      if (fino.degrees == 0) {
        wasora_push_error_message("zero degrees of freedom");
        return WASORA_PARSER_ERROR;
      }
      
      
      if (fino.solution == NULL) {
        wasora_call(fino_define_functions());
      }

      // chequeo de grados de libertad (despues de tener las phi definidas)
      if (fino.mesh->degrees_of_freedom == 0) {
        fino.mesh->degrees_of_freedom = fino.degrees;
      } else if (fino.mesh->degrees_of_freedom != fino.degrees) {
        wasora_push_error_message("inconsistent degrees of freedom (FINO_PROBLEM = %d, MESH = %d)", fino.degrees, fino.mesh->degrees_of_freedom);
        return WASORA_PARSER_ERROR;
      }


/*      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
//kw+FINO_STEP+usage [ JUST_BUILD ]
        if (strcasecmp(token, "JUST_BUILD") == 0) {
          fino_step->do_not_solve = 1;
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }

      }
 */
      
      if (fino.rough) {
        fino.mesh_rough = calloc(1, sizeof(mesh_t));
        fino.mesh_rough->name = strdup("rough");
        HASH_ADD_KEYPTR(hh, wasora_mesh.meshes, fino.mesh_rough->name, strlen(fino.mesh_rough->name), fino.mesh_rough);
      }
      
      instruction = wasora_define_instruction(fino_instruction_step, fino_step);
      // esto no me gusta pero es para callar al valgrind
      instruction->argument_alloced = 1;
      
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_REACTION+usage FINO_REACTION
///kw+FINO_REACTION+desc Computes the reaction at the selected physical group.
///kw+FINO_REACTION+detail The result is stored in the variable or vector provided, depending on the number of degrees of freedoms of the problem. 
///kw+FINO_REACTION+detail If the object passed as `RESULT` does not exist, an appropriate object (scalar variable or vector) is created.
///kw+FINO_REACTION+detail For the elastic problem, the components of the total reaction force are stored in the result vector.
///kw+FINO_REACTION+detail For the thermal problem, the total power passing through the entity is computed as an scalar.
      
    } else if ((strcasecmp(token, "FINO_REACTION") == 0)) {

      char *name;
      fino_reaction_t *reaction;
      
      if (fino.problem_family != problem_family_mechanical) {
        wasora_push_error_message("FINO_REACTION makes sense only in elastic problems");
        return WASORA_PARSER_ERROR;
      }
      
      reaction = calloc(1, sizeof(fino_reaction_t));
      LL_APPEND(fino.reactions, reaction);

      while ((token = wasora_get_next_token(NULL)) != NULL) {
///kw+FINO_REACTION+usage PHYSICAL_GROUP <physical_group>
        
        if (strcasecmp(token, "PHYSICAL_GROUP") == 0 || strcasecmp(token, "PHYSICAL_ENTITY") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((reaction->physical_entity = wasora_get_physical_entity_ptr(name, fino.mesh)) == NULL) {
            reaction->physical_entity = wasora_define_physical_entity(name, fino.mesh, 1);
          }
          free(name);
///kw+FINO_REACTION+usage RESULT { <variable> | <vector> }
        } else if (strcasecmp(token, "RESULT") == 0) {
          wasora_call(wasora_parser_string(&name));
          
          reaction->scalar = wasora_get_variable_ptr(name);
          reaction->vector = wasora_get_vector_ptr(name);
          
          if (reaction->scalar == NULL && reaction->vector == NULL) {
            // nos dieron algo que no existe, hay que crearlo
            if (fino.degrees == 1) {
              if ((reaction->scalar = wasora_define_variable(name)) == NULL) {
                return WASORA_PARSER_ERROR;
              }
            } else {
              if ((reaction->vector = wasora_define_vector(name, fino.degrees, NULL, NULL)) == NULL) {
                return WASORA_PARSER_ERROR;
              }
            }
          } else if (reaction->scalar != NULL && fino.degrees != 1) {
            wasora_push_error_message("RESULT should pass a vector of size %d not a variable", fino.degrees);
            return WASORA_PARSER_ERROR;
          } else if (reaction->vector != NULL && fino.degrees == 1) {
            wasora_push_error_message("RESULT should pass a variable not a vector");
            return WASORA_PARSER_ERROR;
          } else if (reaction->vector != NULL && fino.degrees != reaction->vector->size) {
            wasora_push_error_message("RESULT should pass a vector of size %d not of size %d", fino.degrees, reaction->vector->size);
            return WASORA_PARSER_ERROR;
          }
          
          free(name);
          
        }
      }
      
      wasora_define_instruction(fino_instruction_reaction, reaction);
      return WASORA_PARSER_OK;
      
// ---------------------------------------------------------------------
///kw+FINO_LINEARIZE+usage FINO_LINEARIZE
///kw+FINO_LINEARIZE+desc Performs stress linearization according to ASME VII-Sec 5 over a
///kw+FINO_LINEARIZE+desc Stress Classification Line
///kw+FINO_LINEARIZE+detail The Stress Classification Line (SCL) may be given either as a one-dimensional physical group
///kw+FINO_LINEARIZE+detail in the mesh or as the (continuous) spatial coordinates of two end-points.
    } else if ((strcasecmp(token, "FINO_LINEARIZE") == 0)) {
      
      char *name;
      fino_linearize_t *linearize, *tmp;
      int n_linearizes;
      
      if (fino.problem_family != problem_family_mechanical) {
        wasora_push_error_message("FINO_LINEARIZE makes sense only in elastic problems");
        return WASORA_PARSER_ERROR;
      }
      
      linearize = calloc(1, sizeof(fino_linearize_t));
      LL_APPEND(fino.linearizes, linearize);

      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+FINO_LINEARIZE+usage {
///kw+FINO_LINEARIZE+usage PHYSICAL_GROUP <physical_group>
///kw+FINO_LINEARIZE+detail If the SCL is given as a `PHYSICAL_GROUP`, the entity should be one-dimensional (i.e a line)
///kw+FINO_LINEARIZE+detail independently of the dimension of the problem.
        
        if (strcasecmp(token, "PHYSICAL_ENTITY") == 0 || strcasecmp(token, "PHYSICAL_GROUP") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->physical_entity = wasora_get_physical_entity_ptr(name, fino.mesh)) == NULL) {
            linearize->physical_entity = wasora_define_physical_entity(name, fino.mesh, 1);
          }
          free(name);

///kw+FINO_LINEARIZE+usage |
///kw+FINO_LINEARIZE+usage START_POINT <x1> <y1> <z1>
///kw+FINO_LINEARIZE+detail If the SCL is given with `START_POINT` and `END_POINT`, the number of coordinates given should
///kw+FINO_LINEARIZE+detail match the problem dimension (i.e three coordinates for full\ 3D problems and two coordinates for
///kw+FINO_LINEARIZE+detail axisymmetric or plane problems).
///kw+FINO_LINEARIZE+detail Coordinates can be given algebraic expressions that will be evaluated at the time of the linearization.
        } else if (strcasecmp(token, "START_POINT") == 0) {
          if (fino.dimensions == 0) {
            wasora_push_error_message("need to know the problem dimension before LINEARIZE START_POINT");
            return WASORA_PARSER_ERROR;
          }
          wasora_call(wasora_parser_expression(&linearize->x1));
          if (fino.dimensions > 1) {
            wasora_call(wasora_parser_expression(&linearize->y1));
          }
          if (fino.dimensions > 2) {
            wasora_call(wasora_parser_expression(&linearize->z1));
          }
          
///kw+FINO_LINEARIZE+usage END_POINT <x2> <y2> <z2>
///kw+FINO_LINEARIZE+usage }@
        } else if (strcasecmp(token, "END_POINT") == 0) {
          if (fino.dimensions == 0) {
            wasora_push_error_message("need to know the problem dimension before LINEARIZE END_POINT");
            return WASORA_PARSER_ERROR;
          }
          wasora_call(wasora_parser_expression(&linearize->x2));
          if (fino.dimensions > 1) {
            wasora_call(wasora_parser_expression(&linearize->y2));
          }
          if (fino.dimensions > 2) {
            wasora_call(wasora_parser_expression(&linearize->z2));
          }
          
          
///kw+FINO_LINEARIZE+detail If either a `FILE` or a `FILE_PATH` is given, the total, membrane and membrane plus bending
///kw+FINO_LINEARIZE+detail stresses are written as a function of a scalar $t \in [0,1]$.
///kw+FINO_LINEARIZE+detail Moreover, the individual elements of the membrane and bending stress tensors are written
///kw+FINO_LINEARIZE+detail within comments (i.e. lines starting with the hash symbol `#`).
//TODO: decir como se plotea y como se hace un PDF

///kw+FINO_LINEARIZE+usage [ FILE <file_id> | 
        } else if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&linearize->file));
     
///kw+FINO_LINEARIZE+usage FILE_PATH <file_path> ]@
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
            wasora_call(wasora_parser_file_path(&linearize->file, "w"));

///kw+FINO_LINEARIZE+detail By default, the linearization uses the Von\ Mises criterion for the composition of stresses.
///kw+FINO_LINEARIZE+detail The definition of what _total stress_ means can be changed using the `TOTAL` keyword.
///kw+FINO_LINEARIZE+usage [ TOTAL {
        } else if (strcasecmp(token, "TOTAL") == 0) {
          char *keywords[] = {
///kw+FINO_LINEARIZE+usage   vonmises
                         "vonmises",
///kw+FINO_LINEARIZE+usage   tresca |
                         "tresca",
///kw+FINO_LINEARIZE+usage   tresca |
                         "principal1",
///kw+FINO_LINEARIZE+usage   principal1 |
                         "principal2",
///kw+FINO_LINEARIZE+usage   principal2 |
                         "principal3",
///kw+FINO_LINEARIZE+usage   principal3 }@
                         ""};
          int values[] = {linearize_vonmises, linearize_tresca, linearize_principal1, linearize_principal2, linearize_principal3, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&linearize->total));
                 
///kw+FINO_LINEARIZE+detail The membrane, bending and peak stress tensor elements are combined using the
///kw+FINO_LINEARIZE+detail Von\  Mises criterion and stored as variables.
///kw+FINO_LINEARIZE+detail If no name for any of the variables is given, they are stored in
///kw+FINO_LINEARIZE+detail `M_group`, `B_group` and `P_group` respectively if there is a physical group.
///kw+FINO_LINEARIZE+detail Otherwise `M_1`, `B_1` and `P_1` for the first instruction, `M_2`... etc.
            
///kw+FINO_LINEARIZE+usage [ M <variable> ]@
        } else if (strcasecmp(token, "M") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->M = wasora_get_or_define_variable_ptr(name)) == NULL) {
            return WASORA_PARSER_ERROR;
          }
          free(name);
          
///kw+FINO_LINEARIZE+usage [ MB <variable> ]@
        } else if (strcasecmp(token, "MB") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->MB = wasora_get_or_define_variable_ptr(name)) == NULL) {
            return WASORA_PARSER_ERROR;
          }
          free(name);

///kw+FINO_LINEARIZE+usage [ PEAK <variable> ]@
        } else if (strcasecmp(token, "PEAK") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->P = wasora_get_or_define_variable_ptr(name)) == NULL) {
            return WASORA_PARSER_ERROR;
          }
          free(name);
          
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      n_linearizes = 0;
      LL_FOREACH(fino.linearizes, tmp) {
        n_linearizes++;
      }
      
      if (linearize->M == NULL) {
        if (linearize->physical_entity != NULL) {
          if (asprintf(&name, "M_%s", linearize->physical_entity->name) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        } else {
          if (asprintf(&name, "M_%d", n_linearizes) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        }
        if ((linearize->M = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
        
      }
      if (linearize->MB == NULL) {
        if (linearize->physical_entity != NULL) {
          if (asprintf(&name, "MB_%s", linearize->physical_entity->name) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        } else {
          if (asprintf(&name, "MB_%d", n_linearizes) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        }
        if ((linearize->MB = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
      }
      if (linearize->P == NULL) {
        if (linearize->physical_entity != NULL) {
          if (asprintf(&name, "P_%s", linearize->physical_entity->name) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        } else {
          if (asprintf(&name, "P_%d", n_linearizes) == -1) {
            wasora_push_error_message("memory allocation error");
            return WASORA_PARSER_ERROR;
          }
        }

        if ((linearize->P = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
      }      
      
      if (linearize->physical_entity == NULL && linearize->x1.n_tokens == 0 && linearize->x2.n_tokens == 0) {
        wasora_push_error_message("need to know what the SCL is either as a PHYSICAL_GROUP or START_POINT and END_POINT");
        return WASORA_PARSER_ERROR;
      }
      
      wasora_define_instruction(fino_instruction_linearize, linearize);

      return WASORA_PARSER_OK;
      
// ---------------------------------------------------------------------
//kw+FINO_DEBUG+usage FINO_DEBUG
//kw+FINO_DEBUG+desc Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.
// temporarily disabled
/*      
    } else if ((strcasecmp(token, "FINO_DEBUG") == 0)) {
      
      fino_debug_t *debug;
      debug = calloc(1, sizeof(fino_debug_t));
      LL_APPEND(fino.debugs, debug);

      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
//kw+FINO_DEBUG+usage [ FILE <file_id> | 
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&debug->file));
          
//kw+FINO_DEBUG+usage FILE_PATH <file_path> ]@
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
            wasora_call(wasora_parser_file_path(&debug->file, "w"));
          
//kw+FINO_DEBUG+usage [ MATRICES_ASCII ]@
        } else if (strcasecmp(token, "MATRICES_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII;
//kw+FINO_DEBUG+usage [ MATRICES_ASCII_STRUCTURE ]@
        } else if (strcasecmp(token, "MATRICES_ASCII_STRUCTURE") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII_STRUCT;
//kw+FINO_DEBUG+usage [ MATRICES_PETSC_BINARY ]@
        } else if (strcasecmp(token, "MATRICES_PETSC_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_BINARY;
//kw+FINO_DEBUG+usage [ MATRICES_PETSC_COMPRESSED_BINARY ]@
        } else if (strcasecmp(token, "MATRICES_PETSC_COMPRESSED_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_COMPRESSED_BINARY;
//kw+FINO_DEBUG+usage [ MATRICES_PETSC_ASCII ]@
        } else if (strcasecmp(token, "MATRICES_PETSC_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_ASCII;
//kw+FINO_DEBUG+usage [ MATRICES_PETSC_OCTAVE ]@
        } else if (strcasecmp(token, "MATRICES_PETSC_OCTAVE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_OCTAVE;
//kw+FINO_DEBUG+usage [ MATRICES_PETSC_DENSE ]@
        } else if (strcasecmp(token, "MATRICES_PETSC_DENSE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_DENSE;
//kw+FINO_DEBUG+usage [ MATRICES_X ]@
        } else if (strcasecmp(token, "MATRICES_X") == 0) {
          debug->matrices |= DEBUG_MATRICES_X;
//kw+FINO_DEBUG+usage [ MATRICES_SNG ]@
        } else if (strcasecmp(token, "MATRICES_SNG") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG;
//kw+FINO_DEBUG+usage [ MATRICES_SNG_STRUCT ]@
        } else if (strcasecmp(token, "MATRICES_SNG_STRUCT") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG_STRUCT;

//kw+FINO_DEBUG+usage [ MATRICES_SIZE <expr> ]@
        } else if (strcasecmp(token, "MATRICES_SIZE") == 0 || strcasecmp(token, "MATRICES_X_SIZE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_size));
          
//kw+FINO_DEBUG+usage [ MATRICES_STRIDE <expr> ]@
        } else if (strcasecmp(token, "MATRICES_STRIDE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_stride));

          
//kw+FINO_DEBUG+usage [ INCLUDE_INPUT ]@
        } else if (strcasecmp(token, "INCLUDE_INPUT") == 0) {
          debug->include_input = 1;
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      // si pidieron DEBUG, le pedimos a petsc que loguee cosas
      PetscMemorySetGetMaximumUsage();
#if PETSC_VERSION_LT(3,7,0)
      PetscLogBegin();
#else
      PetscLogDefaultBegin();
#endif
      wasora_define_instruction(fino_instruction_debug, debug);

      return WASORA_PARSER_OK;
*/
    }
  }
  
  // si no entendimos la linea, probamos pasarsela al parser de mallas
  strcpy(line, wasora.line);
  return wasora_mesh_parse_line(line);    
      
  return WASORA_PARSER_UNHANDLED;
  
  
}



int fino_define_functions(void) {
  
  char *name = NULL;
  char *gradname = NULL;
  char *modename = NULL;
  int i, g, m;
  
  // las definimos solo si ya sabemos cuantas dimensiones tiene el problema
  if (fino.dimensions == 0 || fino.degrees == 0) {
    wasora_push_error_message("do not know how many dimensions the problem has, tell me with DIMENSIONS in either FINO_PROBLEM or MESH");
    return WASORA_PARSER_ERROR;
  }

  fino.solution = calloc(fino.degrees, sizeof(function_t *));
  fino.gradient = calloc(fino.degrees, sizeof(function_t *));
  fino.delta_gradient = calloc(fino.degrees, sizeof(function_t *));
  if (fino.nev > 0) {
    fino.mode = calloc(fino.degrees, sizeof(function_t *));
  }

  for (g = 0; g < fino.degrees; g++) {
    if (fino.unknown_name == NULL) {
      if (asprintf(&name, "phi%d", g+1) == -1) {
        wasora_push_error_message("memory allocation error");
        return WASORA_RUNTIME_ERROR;
      }
    } else {
      if (asprintf(&name, "%s", fino.unknown_name[g]) == -1) {
        wasora_push_error_message("memory allocation error");
        return WASORA_RUNTIME_ERROR;
      }
    }
    if ((fino.solution[g] = wasora_define_function(name, fino.dimensions)) == NULL) {
      return WASORA_PARSER_ERROR;
    }
    
    // esto lo ponemos aca por si alguien quiere hacer un PRINT o algo sin pasar por el STEP
    fino.solution[g]->mesh = (fino.rough==0)?fino.mesh:fino.mesh_rough;
    fino.solution[g]->var_argument = calloc(fino.dimensions, sizeof(var_t *));
    fino.solution[g]->var_argument_alloced = 1;
    fino.solution[g]->type = type_pointwise_mesh_node;

    if (fino.nev == 0) {
      // las derivadas de las soluciones con respecto al espacio solo si no es modal
      fino.gradient[g] = calloc(fino.dimensions, sizeof(function_t *));
      fino.delta_gradient[g] = calloc(fino.dimensions, sizeof(function_t *));
      
      for (m = 0; m < fino.dimensions; m++) {
        fino.solution[g]->var_argument[m] = wasora_mesh.vars.arr_x[m];
      
        if (asprintf(&gradname, "d%sd%s", name, wasora_mesh.vars.arr_x[m]->name) == -1) {
          wasora_push_error_message("cannot asprintf");
          return WASORA_RUNTIME_ERROR;
        }
        if ((fino.gradient[g][m] = wasora_define_function(gradname, fino.dimensions)) == NULL) {
          return WASORA_PARSER_ERROR;
        }
        free(gradname);
        gradname = NULL;
        
        fino.gradient[g][m]->mesh = fino.solution[g]->mesh;
        fino.gradient[g][m]->var_argument = fino.solution[g]->var_argument;
        fino.gradient[g][m]->type = type_pointwise_mesh_node;
        fino.gradient[g][m]->spatial_derivative_of = fino.solution[g];
        fino.gradient[g][m]->spatial_derivative_with_respect_to = m;
        
        // lo mismo para la incerteza
        if (asprintf(&gradname, "delta_d%sd%s", name, wasora_mesh.vars.arr_x[m]->name) == -1) {
          wasora_push_error_message("cannot asprintf");
          return WASORA_RUNTIME_ERROR;
        }
        if ((fino.delta_gradient[g][m] = wasora_define_function(gradname, fino.dimensions)) == NULL) {
          return WASORA_PARSER_ERROR;
        }
        free(gradname);
        gradname = NULL;
        
        fino.delta_gradient[g][m]->mesh = fino.solution[g]->mesh;
        fino.delta_gradient[g][m]->var_argument = fino.solution[g]->var_argument;
        fino.delta_gradient[g][m]->type = type_pointwise_mesh_node;
      }
      
    } else {  
      // en modal tenemos muchas soluciones
      fino.mode[g] = calloc(fino.nev, sizeof(function_t *));
      for (i = 0; i < fino.nev; i++) {
        if (asprintf(&modename, "%s%d", name, i+1) == -1) {
          wasora_push_error_message("cannot asprintf");
          return WASORA_RUNTIME_ERROR;
        }
        wasora_call(fino_define_result_function(modename, &fino.mode[g][i]));
        free(modename);
        modename = NULL;
        
        fino.mode[g][i]->mesh = fino.solution[g]->mesh;
        fino.mode[g][i]->var_argument = fino.solution[g]->var_argument;
        fino.mode[g][i]->type = fino.solution[g]->type;
      }
    }
    free(name);
    name = NULL;
  }

  if (fino.problem_family == problem_family_mechanical) {

    // TODO: describir las funciones para reference
    wasora_call(fino_define_result_function("sigmax", &fino.sigmax));
    wasora_call(fino_define_result_function("sigmay", &fino.sigmay));
    wasora_call(fino_define_result_function("tauxy", &fino.tauxy));

    if (fino.dimensions == 3) {
      wasora_call(fino_define_result_function("sigmaz", &fino.sigmaz));
      wasora_call(fino_define_result_function("tauyz", &fino.tauyz));
      wasora_call(fino_define_result_function("tauzx", &fino.tauzx));
    }
    
    wasora_call(fino_define_result_function("sigma1", &fino.sigma1));
    wasora_call(fino_define_result_function("sigma2", &fino.sigma2));
    wasora_call(fino_define_result_function("sigma3", &fino.sigma3));
    wasora_call(fino_define_result_function("sigma", &fino.sigma));
    wasora_call(fino_define_result_function("delta_sigma", &fino.delta_sigma));
    wasora_call(fino_define_result_function("tresca", &fino.tresca));
        
  }
    
  if (fino.nev > 0) {
///va+M_T+name M_T
///va+M_T+desc Total mass\ $m$ computed from the mass matrix\ $M$ as
///va+M_T+desc 
///va+M_T+desc \[ M_T = \frac{1}{n_\text{DOFs}} \cdot \vec{1}^T \cdot M \cdot \vec{1} \]
///va+M_T+desc 
///va+M_T+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node.
///va+M_T+desc Note that this is only approximately equal to the actual mass, i.e. the integral of the density $\rho(x,y,z)$ over the problem domain.
    fino.vars.M_T = wasora_define_variable("M_T");
    
///ve+f+name f
///ve+f+desc _Size:_ number of requested eigen-pairs.
///ve+f+desc _Elements:_ The frequency $f_i$ of the $i$-th mode, in cycles per unit of time.
    fino.vectors.f = wasora_define_vector("f", fino.nev, NULL, NULL);    

///ve+omega+name omega
///ve+omega+desc _Size:_ number of requested eigen-pairs.
///ve+omega+desc _Elements:_ The angular frequency $\omega_i$ of the $i$-th mode, in radians per unit of time.
    fino.vectors.omega = wasora_define_vector("omega", fino.nev, NULL, NULL);    

    
///ve+m+name m
///ve+m+desc _Size:_ number of requested eigen-pairs.
///ve+m+desc _Elements:_ The generalized modal mass $M_i$ of the $i$-th mode computed as
///ve+m+desc
///ve+m+desc \[ \text{m}_i = \frac{1}{n_\text{DOFs}} \vec{\phi}_i^T \cdot M \cdot \vec{\phi}_i \]
///va+m+desc 
///va+m+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node, $M$ is the mass matrix
///va+m+desc and $\vec{\phi}_i$ is the $i$-th eigenvector normalized such that the largest element is equal to one.

    fino.vectors.m = wasora_define_vector("m", fino.nev, NULL, NULL);    

///ve+L+name L
///ve+L+desc _Size:_ number of requested eigen-pairs.
///ve+L+desc _Elements:_ The excitation factor $L_i$ of the $i$-th mode computed as
///ve+L+desc
///ve+L+desc \[ L_i = \frac{1}{n_\text{DOFs}} \cdot \vec{\phi}_i^T \cdot M \cdot \vec{1} \]
///va+L+desc 
///va+L+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node, $M$ is the mass matrix
///va+L+desc and $\vec{\phi}_i$ is the $i$-th eigenvector normalized such that the largest element is equal to one.
    fino.vectors.L = wasora_define_vector("L", fino.nev, NULL, NULL);    

///ve+Gamma+name Gamma
///ve+Gamma+desc _Size:_ number of requested eigen-pairs.
///ve+Gamma+desc _Elements:_ The participation factor $\Gamma_i$ of the $i$-th mode computed as
///ve+Gamma+desc
///ve+Gamma+desc \[ \Gamma_i = \frac{ \vec{\phi}_i^T \cdot M \cdot \vec{1} }{ \vec{\phi}_i^T \cdot M \cdot \vec{\phi}} \]
    fino.vectors.Gamma = wasora_define_vector("Gamma", fino.nev, NULL, NULL);    
    
///ve+mu+name mu
///ve+mu+desc _Size:_ number of requested eigen-pairs.
///ve+mu+desc _Elements:_ The relatve effective modal mass $\mu_i$ of the $i$-th mode computed as
///ve+mu+desc
///ve+mu+desc \[ \mu_i = \frac{L_i^2}{M_t \cdot n_\text{DOFs} \cdot m_i} \]
///ve+mu+desc
///ve+mu+desc Note that $\sum_{i=1}^N m_i = 1$, where $N$ is total number of degrees of freedom ($n_\text{DOFs}$ times the number of nodes).
    fino.vectors.mu = wasora_define_vector("mu", fino.nev, NULL, NULL);    
    
///ve+Mu+name Mu
///ve+Mu+desc _Size:_ number of requested eigen-pairs.
///ve+Mu+desc _Elements:_ The accumulated relative effective modal mass $\Mu_i$ up to the $i$-th mode computed as
///ve+Mu+desc
///ve+Mu+desc \[ \Mu_i = \sum_{j=1}^i \mu_i \]
///ve+Mu+desc
///ve+Mu+desc Note that $\Mu_N = 1$, where $N$ is total number of degrees of freedom ($n_\text{DOFs}$ times the number of nodes).
    fino.vectors.Mu = wasora_define_vector("Mu", fino.nev, NULL, NULL);    
    
    fino.vectors.phi = malloc(fino.nev * sizeof(vector_t *));
    for (i = 0; i < fino.nev; i++) {
      if (asprintf(&modename, "phi%d", i+1) == -1) {
        wasora_push_error_message("memory allocation error");
        return WASORA_RUNTIME_ERROR;
      }
      
      if ((fino.vectors.phi[i] = wasora_define_vector(modename, 0, NULL, NULL)) == NULL) {
        wasora_push_error_message("cannot define vector %s", modename);
        return WASORA_RUNTIME_ERROR;
      }
      free(modename);
    }
    
  }
  
  // TODO: heat flux
  
  return WASORA_PARSER_OK;
}
