/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's parsing routines
 *
 *  Copyright (C) 2015--2018 jeremy theler
 *
 *  This file is part of fino.
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
#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "plugin_parse_line"
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
///kw+FINO_PROBLEM+usage BAKE
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc 
///kw+FINO_PROBLEM+desc  * `BAKE` (or `HEAT` or `THERMAL`) solves the heat conduction problem.
        if (strcasecmp(token, "BAKE") == 0 || strcasecmp(token, "HEAT") == 0 || strcasecmp(token, "THERMAL") == 0) {
          fino.problem_family = problem_family_bake;
          fino.problem_kind = problem_kind_full3d;
          fino.math_type = math_type_linear;
          fino.degrees = 1;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("T");

///kw+FINO_PROBLEM+usage SHAKE
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc  * `SHAKE` (or `MODAL`) computes the natural frequencies and modes.        
        } else if (strcasecmp(token, "SHAKE") == 0 || strcasecmp(token, "MODAL") == 0) {
#ifndef HAVE_SLEPC
          wasora_push_error_message("MODAL needs a fino binary linked agains SLEPc.");
          return WASORA_PARSER_ERROR;
#endif
          fino.problem_family = problem_family_shake;
          fino.problem_kind = problem_kind_full3d;
          fino.math_type = math_type_eigen;
          fino.dimensions = 3;
          fino.degrees = 3;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          fino.unknown_name[2] = strdup("w");
          
///kw+FINO_PROBLEM+usage BREAK
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc  * `BREAK` (or `ELASTIC` or `MECHANICAL`) solves the elastic problem.        
        } else if (strcasecmp(token, "BREAK") == 0 || strcasecmp(token, "ELASTIC") == 0 || strcasecmp(token, "MECHANICAL") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_full3d;
          fino.math_type = math_type_linear;
          fino.dimensions = 3;
          fino.degrees = 3;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          fino.unknown_name[2] = strdup("w");

///kw+FINO_PROBLEM+usage HEAT_AXISYMMETRIC
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc  * `HEAT_AXISYMMETRIC` solves the heat conduction problem in an axysimmetric way.          
        } else if (strcasecmp(token, "HEAT_AXISYMMETRIC") == 0) {
          fino.problem_family = problem_family_bake;
          fino.problem_kind = problem_kind_axisymmetric;
          if (fino.symmetry_axis == symmetry_axis_none) {
            fino.symmetry_axis = symmetry_axis_y;
          }
          fino.math_type = math_type_linear;
          fino.dimensions = 2;
          fino.degrees = 1;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("T");

///kw+FINO_PROBLEM+usage PLANE_STRESS
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc  * `PLANE_STRESS` solves the plane stress elastic problem. 
        } else if (strcasecmp(token, "PLANE_STRESS") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_plane_stress;
          if (fino.symmetry_axis == symmetry_axis_none) {
            fino.symmetry_axis = symmetry_axis_y;
          }
          fino.math_type = math_type_linear;
          fino.dimensions = 2;
          fino.degrees = 2;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");

///kw+FINO_PROBLEM+usage PLANE_STRAIN
///kw+FINO_PROBLEM+usage |
///kw+FINO_PROBLEM+desc  * `PLANE_STRAIN` solves the plane strain elastic problem.
        } else if (strcasecmp(token, "PLANE_STRAIN") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_plane_strain;
          if (fino.symmetry_axis == symmetry_axis_none) {
            fino.symmetry_axis = symmetry_axis_y;
          }
          fino.math_type = math_type_linear;
          fino.dimensions = 2;
          fino.degrees = 2;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          
///kw+FINO_PROBLEM+usage ELASTIC_AXISYMMETRIC
///kw+FINO_PROBLEM+usage ]
///kw+FINO_PROBLEM+desc  * `ELASTIC_AXISYMMETRIC` solves the elastic problem in an axysimmetric way.
///kw+FINO_PROBLEM+desc 
        } else if (strcasecmp(token, "ELASTIC_AXISYMMETRIC") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_axisymmetric;
          if (fino.symmetry_axis == symmetry_axis_none) {
            fino.symmetry_axis = symmetry_axis_y;
          }
          fino.math_type = math_type_linear;
          fino.dimensions = 2;
          fino.degrees = 2;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");

///kw+FINO_PROBLEM+usage [ DIMENSIONS <expr> ]
///kw+FINO_PROBLEM+desc The number of dimensions needs to be given with `DIMENSIONS` for the heat conduction problem.
        } else if (strcasecmp(token, "DIMENSIONS") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          fino.dimensions = (int)(xi);
          if (fino.dimensions < 1 || fino.dimensions > 3)  {
            wasora_push_error_message("either one, two or three dimensions should be selected instead of '%d'", fino.dimensions);
            return WASORA_PARSER_ERROR;
          }
          
///kw+FINO_PROBLEM+usage [ DEGREES <expr> ]
        } else if (strcasecmp(token, "DEGREES") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          fino.degrees = (int)(xi);
          if (fino.degrees < 1)  {
            wasora_push_error_message("a positive number of degrees should be given instead of '%d'", fino.degrees);
            return WASORA_PARSER_ERROR;
          }
          
///kw+FINO_PROBLEM+usage [ SYMMETRY_AXIS { x | y } ]
        } else if (strcasecmp(token, "SYMMETRY_AXIS") == 0) {
          char *keywords[] = { "x", "y" };
          int values[] = {symmetry_axis_x, symmetry_axis_y, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.symmetry_axis));
          
///kw+FINO_PROBLEM+usage [ MESH <identifier> ]
        } else if (strcasecmp(token, "MESH") == 0) {
          char *mesh_name;
          
          wasora_call(wasora_parser_string(&mesh_name));
          if ((fino.mesh = wasora_get_mesh_ptr(mesh_name)) == NULL) {
            wasora_push_error_message("unknown mesh '%s'", mesh_name);
            free(mesh_name);
            return WASORA_PARSER_ERROR;
          }
          free(mesh_name);
        
#ifdef HAVE_SLEPC
///kw+FINO_PROBLEM+usage [ N_EIGEN <expr> ]
        } else if (strcasecmp(token, "N_EIGEN") == 0) {
          wasora_call(wasora_parser_expression_in_string(&xi));
          fino.nev = (int)(xi);
          if (fino.nev < 1)  {
            wasora_push_error_message("a positive number of eigenvalues should be given instead of '%d'", fino.nev);
            return WASORA_PARSER_ERROR;
          }
#endif

///kw+FINO_PROBLEM+usage [ UNKNOWNS <name1> <name2> ... <name_degrees> ]
        } else if (strcasecmp(token, "SOLUTION") == 0 || strcasecmp(token, "SOLUTIONS") == 0 ||
                   strcasecmp(token, "SOLUTION_NAMES") == 0 ||
                   strcasecmp(token, "UNKNOWNS") == 0 || strcasecmp(token, "UNKNOWN") == 0) {

          int g;
          
          if (fino.degrees == 0) {
            wasora_push_error_message("UNKNOWNS before DEGREES");
            return WASORA_PARSER_ERROR;
          }
          
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          for (g = 0; g < fino.degrees; g++) {
            wasora_call(wasora_parser_string(&fino.unknown_name[g]));
          }
          
        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      // si no nos dieron explicitamente la malla, ponemos la principal
      if (fino.mesh == NULL && (fino.mesh = wasora_mesh.main_mesh) == NULL) {
        wasora_push_error_message("unknown mesh for FINO_PROBLEM (no MESH keyword)", token);
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

     
      wasora_call(fino_define_functions());
      
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_SOLVER+usage FINO_SOLVER
///kw+FINO_SOLVER+desc Sets options related to the eigen-solver.
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

///kw+FINO_SOLVER+usage [ KSP_TYPE { gmres | bcgs | bicg | richardson | chebyshev | ... } ]
///kw+FINO_SOLVER+detail List of `KSP_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html>
///kw+FINO_SOLVER+detail          
        } else if (strcasecmp(token, "KSP_TYPE") == 0) {
          wasora_call(wasora_parser_string(&fino.ksp_type));

///kw+FINO_SOLVER+usage [ PC_TYPE { lu | gamg | hypre | sor | bjacobi | cholesky | ... } ]
///kw+FINO_SOLVER+detail List of `PC_TYPE`s <http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html>
///kw+FINO_SOLVER+detail          
        } else if (strcasecmp(token, "PC_TYPE") == 0) {
          wasora_call(wasora_parser_string(&fino.pc_type));

///kw+FINO_SOLVER+usage [ SET_NEAR_NULLSPACE { rigidbody | fino | none } ]
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
          

///kw+FINO_SOLVER+usage [ DO_NOT_SET_BLOCK_SIZE 
        } else if (strcasecmp(token, "DO_NOT_SET_BLOCK_SIZE") == 0) {
          fino.do_not_set_block_size = 1;
          
///kw+FINO_SOLVER+usage | SET_BLOCK_SIZE ]
        } else if (strcasecmp(token, "SET_BLOCK_SIZE") == 0) {
          fino.do_not_set_block_size = 0;

///kw+FINO_SOLVER+usage [ GRADIENT_EVALUATION {
        } else if (strcasecmp(token, "GRADIENT_EVALUATION") == 0) {
          char *keywords[] = {
///kw+FINO_SOLVER+usage   mass_matrix_consistent |
                         "mass_matrix_consistent",
///kw+FINO_SOLVER+usage   mass_matrix_row_sum |
                         "mass_matrix_row_sum",
///kw+FINO_SOLVER+usage   mass_matrix_lobatto |
                         "mass_matrix_lobatto",
///kw+FINO_SOLVER+usage   mass_matrix_diagonal |
                         "mass_matrix_diagonal",
///kw+FINO_SOLVER+usage   node_average_all |
                         "node_average_all",
///kw+FINO_SOLVER+usage   node_average_corner |
                         "node_average_corner",
///kw+FINO_SOLVER+usage   gauss_average |
                         "gauss_average",
///kw+FINO_SOLVER+usage   none
                         "none",
                         ""};
          int values[] = {gradient_mass_matrix_consistent, gradient_mass_matrix_row_sum, gradient_mass_matrix_lobatto, gradient_mass_matrix_diagonal, gradient_node_average_all, gradient_node_average_corner, gradient_gauss_average, gradient_none, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&fino.gradient_evaluation));

///kw+FINO_SOLVER+usage } ]

///kw+FINO_SOLVER+usage [ GRADIENT_JACOBIAN_THRESHOLD <expr> ]
        } else if (strcasecmp(token, "GRADIENT_JACOBIAN_THRESHOLD") == 0) {
          wasora_call(wasora_parser_expression_in_string(&fino.gradient_jacobian_threshold));

///kw+FINO_SOLVER+usage [ PROGRESS_ASCII ]
        } else if (strcasecmp(token, "PROGRESS_ASCII") == 0) {
          fino.progress_ascii = 1;

          ///kw+FINO_SOLVER+usage [ PROGRESS_BUILD_SHM <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_PROGRESS_BUILD") == 0 || strcasecmp(token, "PROGRESS_BUILD_SHMEM") == 0) {
          wasora_call(wasora_parser_string(&fino.progress_build_shname));

///kw+FINO_SOLVER+usage [ PROGRESS_SOLVE_SHMEM <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_PROGRESS_SOLVE") == 0 || strcasecmp(token, "PROGRESS_SOLVE_SHMEM") == 0) {
          wasora_call(wasora_parser_string(&fino.progress_solve_shname));

///kw+FINO_SOLVER+usage [ MEMORY_USAGE_SHMEM <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_MEMORY") == 0 || strcasecmp(token, "MEMORY_USAGE_SHMEM") == 0) {
          wasora_call(wasora_parser_string(&fino.memory_shname));

        } else {
          wasora_push_error_message("undefined keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }

      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_STEP+usage FINO_STEP
///kw+FINO_STEP+desc Solves the linear eigenvalue problem.
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
        
      // defaulteamos a break
      if (fino.problem_family == problem_family_undefined) {
        fino.problem_family = problem_family_break;
        fino.problem_kind = problem_kind_full3d;
        fino.math_type = math_type_linear;
        fino.degrees = 3;
        fino.unknown_name = calloc(fino.degrees, sizeof(char *));
        fino.unknown_name[0] = strdup("u");
        fino.unknown_name[1] = strdup("v");
        fino.unknown_name[2] = strdup("w");
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
      
      // si nadie dijo nada tenemos un solo grado de libertado
      if (fino.degrees == 0) {
        fino.degrees = 1;
      }
      
      
      if (fino.solution == NULL) {
        wasora_call(fino_define_functions());
      }

      // chequeo de grados de libertad (despues de tener las phi definidas)
      // por default la cantidad de grupos es uno, asi que no hay manera de que sea cero
      if (fino.mesh->degrees_of_freedom == 0) {
        fino.mesh->degrees_of_freedom = fino.degrees;
      } else if (fino.mesh->degrees_of_freedom != fino.degrees) {
        wasora_push_error_message("inconsistent degrees of freedom (FINO_PROBLEM = %d, MESH = %d)", fino.degrees, fino.mesh->degrees_of_freedom);
        return WASORA_PARSER_ERROR;
      }

      
      while ((token = wasora_get_next_token(NULL)) != NULL) {
///kw+FINO_STEP+usage [ JUST_BUILD |
        if (strcasecmp(token, "JUST_BUILD") == 0) {
          fino_step->do_not_build = 0;
          fino_step->do_not_solve = 1;
///kw+FINO_STEP+usage JUST_SOLVE ]
        } else if (strcasecmp(token, "JUST_SOLVE") == 0) {
          fino_step->do_not_build = 1;
          fino_step->do_not_solve = 0;
          
///kw+FINO_STEP+usage [ DUMP_FILE_PATH <filepath> ]
/*          
        } else if (strcasecmp(token, "DUMP_FILE_PATH") == 0) {
          char *filepath;
          wasora_parser_string(&filepath);
          fino_step.dump = fopen(filepath, "w");
          free(filepath);
*/        
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }

      }
      
      
      instruction = wasora_define_instruction(fino_instruction_step, fino_step);
      // esto no me gusta pero es para callar al valgrind
      instruction->argument_alloced = 1;
      
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_LINEARIZE+usage FINO_LINEARIZE
///kw+FINO_LINEARIZE+desc Performs stress linearization according to ASME VII-Sec 5 over a
///kw+FINO_LINEARIZE+desc Stress Classification Line given either as a one-dimensional physical entity in the
///kw+FINO_LINEARIZE+desc mesh or as the (continuous) spatial coordinates of two end-points.
      
    } else if ((strcasecmp(token, "FINO_LINEARIZE") == 0)) {
      
      char *name;
      fino_linearize_t *linearize, *tmp;
      int n_linearizes;
      
      if (fino.problem_family != problem_family_break) {
        wasora_push_error_message("FINO_LINEARIZE makes sense only in elastic problems");
        return WASORA_PARSER_ERROR;
      }
      
      linearize = calloc(1, sizeof(fino_linearize_t));
      LL_APPEND(fino.linearizes, linearize);

      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+FINO_LINEARIZE+usage {
///kw+FINO_LINEARIZE+usage PHYSICAL_ENTITY <physical_entity_name>
///kw+FINO_LINEARIZE+desc If the SCL is given as a `PHYSICAL_ENTITY`, the entity should be one-dimensional (i.e a line)
///kw+FINO_LINEARIZE+desc independently of the dimension of the problem.
        
        if (strcasecmp(token, "PHYSICAL_ENTITY") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->physical_entity = wasora_get_physical_entity_ptr(name, fino.mesh)) == NULL) {
            linearize->physical_entity = wasora_define_physical_entity(name, fino.mesh, 1);
          }
          free(name);

///kw+FINO_LINEARIZE+usage |
///kw+FINO_LINEARIZE+usage START_POINT <x1> <y1> <z1>
///kw+FINO_LINEARIZE+desc If the SCL is given with `START_POINT` and `END_POINT`, the number of coordinates given should
///kw+FINO_LINEARIZE+desc match the problem dimension (i.e three coordinates for full\ 3D problems and two coordinates for
///kw+FINO_LINEARIZE+desc axisymmetric or plane problems).
///kw+FINO_LINEARIZE+desc Coordinates can be given algebraic expressions that will be evaluated at the time of the linearization.
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
///kw+FINO_LINEARIZE+usage }
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
          
          
///kw+FINO_LINEARIZE+desc If either a `FILE` or a `FILE_PATH` is given, the total, membrane and membrane plus bending
///kw+FINO_LINEARIZE+desc stresses are written as a function of a scalar $t \in [0,1]$.
///kw+FINO_LINEARIZE+desc Moreover, the individual elements of the membrane and bending stress tensors are written
///kw+FINO_LINEARIZE+desc within comments (i.e. lines starting with the hash symbol `#`).
//TODO: decir como se plotea y como se hace un PDF

///kw+FINO_LINEARIZE+usage [ FILE <file_id> | 
        } else if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&linearize->file));
     
///kw+FINO_LINEARIZE+usage FILE_PATH <file_path> ]
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
            wasora_call(wasora_parser_file_path(&linearize->file, "w"));

///kw+FINO_LINEARIZE+desc By default, the linearization uses the Von\ Mises criterion for the composition of stresses.
///kw+FINO_LINEARIZE+desc The definition of what _total stress_ means can be changed using the `TOTAL` keyword.
///kw+FINO_SOLVER+usage [ TOTAL {
        } else if (strcasecmp(token, "TOTAL") == 0) {
          char *keywords[] = {
///kw+FINO_SOLVER+usage   vonmises",
                         "vonmises",
///kw+FINO_SOLVER+usage   tresca |
                         "tresca",
///kw+FINO_SOLVER+usage   tresca |
                         "principal1",
///kw+FINO_SOLVER+usage   principal1 |
                         "principal2",
///kw+FINO_SOLVER+usage   principal2 |
                         "principal3",
///kw+FINO_SOLVER+usage   principal3 |
                         ""};
          int values[] = {linearize_vonmises, linearize_tresca, linearize_principal1, linearize_principal2, linearize_principal3, 0};
          wasora_call(wasora_parser_keywords_ints(keywords, values, (int *)&linearize->total));
                 
///kw+FINO_LINEARIZE+desc The membrane, bending and peak stress tensor elements are combined using the
///kw+FINO_LINEARIZE+desc Von\  Mises criterion and stored as variables.
///kw+FINO_LINEARIZE+desc If no name for any of the variables is given, they are stored in
///kw+FINO_LINEARIZE+desc `M_entity`, `B_entity` and `P_entity` respectively if there is a physical entity.
///kw+FINO_LINEARIZE+desc Otherwise `M_1`, `B_1` and `P_1` for the first instruction, `M_2`... etc.
            
///kw+FINO_LINEARIZE+usage [ M <variable> ]
        } else if (strcasecmp(token, "M") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->M = wasora_get_or_define_variable_ptr(name)) == NULL) {
            return WASORA_PARSER_ERROR;
          }
          free(name);
          
///kw+FINO_LINEARIZE+usage [ MB <variable> ]
        } else if (strcasecmp(token, "MB") == 0) {
          wasora_call(wasora_parser_string(&name));
          if ((linearize->MB = wasora_get_or_define_variable_ptr(name)) == NULL) {
            return WASORA_PARSER_ERROR;
          }
          free(name);

///kw+FINO_LINEARIZE+usage [ PEAK <variable> ]
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
          asprintf(&name, "M_%s", linearize->physical_entity->name);
        } else {
          asprintf(&name, "M_%d", n_linearizes);
        }
        if ((linearize->M = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
        
      }
      if (linearize->MB == NULL) {
        if (linearize->physical_entity != NULL) {
          asprintf(&name, "MB_%s", linearize->physical_entity->name);
        } else {
          asprintf(&name, "MB_%d", n_linearizes);
        }
        if ((linearize->MB = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
      }
      if (linearize->P == NULL) {
        if (linearize->physical_entity != NULL) {
          asprintf(&name, "P_%s", linearize->physical_entity->name);
        } else {
          asprintf(&name, "P_%d", n_linearizes);
        }

        if ((linearize->P = wasora_define_variable(name)) == NULL) {
          free(name);
          return WASORA_PARSER_ERROR;
        }
        free(name);
      }      
      
      if (linearize->physical_entity == NULL && linearize->x1.n_tokens == 0 && linearize->x2.n_tokens == 0) {
        wasora_push_error_message("need to know what the SCL is either as a PHYSICAL_ENTITY or START_POINT and END_POINT");
        return WASORA_PARSER_ERROR;
      }
      
      wasora_define_instruction(fino_instruction_linearize, linearize);

      return WASORA_PARSER_OK;
      
// ---------------------------------------------------------------------
///kw+FINO_DEBUG+usage FINO_DEBUG
///kw+FINO_DEBUG+desc Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.
    } else if ((strcasecmp(token, "FINO_DEBUG") == 0)) {
      
      fino_debug_t *debug;
      debug = calloc(1, sizeof(fino_debug_t));
      LL_APPEND(fino.debugs, debug);

      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+FINO_DEBUG+usage [ FILE <file_id> | 
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&debug->file));
          
///kw+FINO_DEBUG+usage FILE_PATH <file_path> ]
        } else if (strcasecmp(token, "FILE_PATH") == 0) {
            wasora_call(wasora_parser_file_path(&debug->file, "w"));
          
///kw+FINO_DEBUG+usage [ MATRICES_ASCII ]
        } else if (strcasecmp(token, "MATRICES_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII;
///kw+FINO_DEBUG+usage [ MATRICES_ASCII_STRUCTURE ]
        } else if (strcasecmp(token, "MATRICES_ASCII_STRUCTURE") == 0) {
          debug->matrices |= DEBUG_MATRICES_ASCII_STRUCT;
///kw+FINO_DEBUG+usage [ MATRICES_PETSC_BINARY ]
        } else if (strcasecmp(token, "MATRICES_PETSC_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_BINARY;
///kw+FINO_DEBUG+usage [ MATRICES_PETSC_COMPRESSED_BINARY ]
        } else if (strcasecmp(token, "MATRICES_PETSC_COMPRESSED_BINARY") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_COMPRESSED_BINARY;
///kw+FINO_DEBUG+usage [ MATRICES_PETSC_ASCII ]
        } else if (strcasecmp(token, "MATRICES_PETSC_ASCII") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_ASCII;
///kw+FINO_DEBUG+usage [ MATRICES_PETSC_OCTAVE ]
        } else if (strcasecmp(token, "MATRICES_PETSC_OCTAVE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_OCTAVE;
///kw+FINO_DEBUG+usage [ MATRICES_PETSC_DENSE ]
        } else if (strcasecmp(token, "MATRICES_PETSC_DENSE") == 0) {
          debug->matrices |= DEBUG_MATRICES_PETSC_DENSE;
///kw+FINO_DEBUG+usage [ MATRICES_X ]
        } else if (strcasecmp(token, "MATRICES_X") == 0) {
          debug->matrices |= DEBUG_MATRICES_X;
///kw+FINO_DEBUG+usage [ MATRICES_SNG ]
        } else if (strcasecmp(token, "MATRICES_SNG") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG;
///kw+FINO_DEBUG+usage [ MATRICES_SNG_STRUCT ]
        } else if (strcasecmp(token, "MATRICES_SNG_STRUCT") == 0) {
          debug->matrices |= DEBUG_MATRICES_SNG_STRUCT;

///kw+FINO_DEBUG+usage [ MATRICES_SIZE <expr> ]
        } else if (strcasecmp(token, "MATRICES_SIZE") == 0 || strcasecmp(token, "MATRICES_X_SIZE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_size));
          
///kw+FINO_DEBUG+usage [ MATRICES_STRIDE <expr> ]
        } else if (strcasecmp(token, "MATRICES_STRIDE") == 0) {
          wasora_call(wasora_parser_expression(&debug->matrices_stride));

          
///kw+FINO_DEBUG+usage [ INCLUDE_INPUT ]
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

    }
  }
  
  // si no entendimos la linea, probamos pasarsela al parser de mallas
  strcpy(line, wasora.line);
  return wasora_mesh_parse_line(line);    
      
  return WASORA_PARSER_UNHANDLED;
  
  
}



#undef  __FUNCT__
#define __FUNCT__ "fino_define_functions"
int fino_define_functions(void) {
  
  char *name;
  char *gradname;
  char *vibname;
  int i, g, d;
  
  // las definimos solo si ya sabemos cuantas dimensiones tiene el problema
  if (fino.dimensions == 0) {
    wasora_push_error_message("do not know how many dimensions the problem has, tell me with DIMENSIONS in either FINO_PROBLEM or MESH");
    return WASORA_PARSER_ERROR;
  } else if (fino.degrees == 0) {
    return WASORA_PARSER_ERROR;
  }

  fino.solution = calloc(fino.degrees, sizeof(function_t *));
  fino.gradient = calloc(fino.degrees, sizeof(function_t *));
  if (fino.nev > 1) {
    fino.vibration = calloc(fino.degrees, sizeof(function_t *));
  }

  for (g = 0; g < fino.degrees; g++) {
    if (fino.unknown_name == NULL) {
      if (asprintf(&name, "phi%d", g+1) == -1) {
        wasora_push_error_message("cannot asprintf");
        return WASORA_RUNTIME_ERROR;
      }
    } else {
      if (asprintf(&name, "%s", fino.unknown_name[g]) == -1) {
        wasora_push_error_message("cannot asprintf");
        return WASORA_RUNTIME_ERROR;
      }
    }
    if ((fino.solution[g] = wasora_define_function(name, fino.dimensions)) == NULL) {
      return WASORA_PARSER_ERROR;
    }
    
    // esto lo ponemos aca por si alguien quiere hacer un PRINT o algo sin pasar por el STEP
    fino.solution[g]->mesh = fino.mesh;
    fino.solution[g]->var_argument = calloc(fino.dimensions, sizeof(var_t *));
    fino.solution[g]->var_argument_alloced = 1;
    fino.solution[g]->type = type_pointwise_mesh_node;

    // las derivadas de las soluciones con respecto al espacio
    fino.gradient[g] = calloc(fino.dimensions, sizeof(function_t *));
    
    for (d = 0; d < fino.dimensions; d++) {
      fino.solution[g]->var_argument[d] = wasora_mesh.vars.arr_x[d];
      
      if (asprintf(&gradname, "d%sd%s", name, wasora_mesh.vars.arr_x[d]->name) == -1) {
        wasora_push_error_message("cannot asprintf");
        return WASORA_RUNTIME_ERROR;
      }
      if ((fino.gradient[g][d] = wasora_define_function(gradname, fino.dimensions)) == NULL) {
        return WASORA_PARSER_ERROR;
      }
/*      
      fino.gradient[g][d]->mesh = fino.mesh;
      fino.gradient[g][d]->var_argument = fino.solution[g]->var_argument;
      fino.gradient[g][d]->type = type_pointwise_mesh_node;
 */
      free(gradname);
    }
    
    if (fino.nev > 1) {

      fino.vibration[g] = calloc(fino.nev, sizeof(function_t *));
      for (i = 0; i < fino.nev; i++) {
        if (asprintf(&vibname, "%s%d", name, i+1) == -1) {
          wasora_push_error_message("cannot asprintf");
          return WASORA_RUNTIME_ERROR;
        }
        wasora_call(fino_define_result_function(vibname, &fino.vibration[g][i]));
        free(vibname);
      }
    }
    
    free(name);
  }

  if (fino.problem_family == problem_family_break) {

    wasora_call(fino_define_result_function("sigmax", &fino.sigmax));
    wasora_call(fino_define_result_function("sigmay", &fino.sigmay));
    wasora_call(fino_define_result_function("sigmaz", &fino.sigmaz));
    wasora_call(fino_define_result_function("tauxy", &fino.tauxy));

    if (fino.dimensions == 3) {
      wasora_call(fino_define_result_function("tauyz", &fino.tauyz));
      wasora_call(fino_define_result_function("tauzx", &fino.tauzx));
    }
    
    wasora_call(fino_define_result_function("sigma1", &fino.sigma1));
    wasora_call(fino_define_result_function("sigma2", &fino.sigma2));
    wasora_call(fino_define_result_function("sigma3", &fino.sigma3));
    wasora_call(fino_define_result_function("sigma", &fino.sigma));
    wasora_call(fino_define_result_function("tresca", &fino.tresca));
        
  }
    
  if (fino.nev > 1) {
///va+mass+name mass
///va+mass+desc Total mass\ $m$ computed from the mass matrix\ $M$ as
///va+mass+desc 
///va+mass+desc \[ m = \frac{1}{n_\text{DOFs}} \cdot [1]^T \cdot M \cdot [1] \]
///va+mass+desc 
///va+mass+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node.
    fino.vars.mass = wasora_define_variable("mass");
    
///ve+f+name f
///ve+f+desc _Size:_ number of requested eigen-pairs.
///ve+f+desc _Elements:_ The frequency $f_i$ of the $i$-th mode, in cycles per unit of time.
    fino.vectors.f = wasora_define_vector("f", fino.nev, NULL, NULL);    

///ve+omega+name omega
///ve+omega+desc _Size:_ number of requested eigen-pairs.
///ve+omega+desc _Elements:_ The angular frequency $\omega_i$ of the $i$-th mode, in radians per unit of time.
    fino.vectors.omega = wasora_define_vector("omega", fino.nev, NULL, NULL);    

    
///ve+M+name M
///ve+M+desc _Size:_ number of requested eigen-pairs.
///ve+M+desc _Elements:_ The generalized modal mass $M_i$ of the $i$-th mode computed as
///ve+M+desc
///ve+M+desc \[ M_i = \frac{1}{n_\text{DOFs}} \cdot \vec{\phi}_i^T \cdot M \cdot \vec{\phi}_i \]
///va+M+desc 
///va+M+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node.

    fino.vectors.M = wasora_define_vector("M", fino.nev, NULL, NULL);    

///ve+L+name L
///ve+L+desc _Size:_ number of requested eigen-pairs.
///ve+L+desc _Elements:_ The excitation factor $L_i$ of the $i$-th mode computed as
///ve+L+desc
///ve+L+desc \[ L_i = \frac{1}{n_\text{DOFs}} \cdot \vec{\phi}_i^T \cdot M \cdot [\vec{1}] \]
///va+L+desc 
///va+K+desc where $n_\text{DOFs}$ is the number of degrees of freedoms per node.
    fino.vectors.L = wasora_define_vector("L", fino.nev, NULL, NULL);    

///ve+Gamma+name Gamma
///ve+Gamma+desc _Size:_ number of requested eigen-pairs.
///ve+Gamma+desc _Elements:_ The participation factor $\Gamma_i$ of the $i$-th mode computed as
///ve+Gamma+desc
///ve+Gamma+desc \[ \Gamma_i = \frac{L_i}{M_i} \]
    fino.vectors.Gamma = wasora_define_vector("Gamma", fino.nev, NULL, NULL);    
    
///ve+Me+name Me
///ve+Me+desc _Size:_ number of requested eigen-pairs.
///ve+Me+desc _Elements:_ The effective modal mass factor $M^e_i$ of the $i$-th mode computed as
///ve+Me+desc
///ve+Me+desc \[ M^e_i = \frac{L_i^2}{n_\text{DOFs} \cdot M_i} \]
///ve+Me+desc
///ve+Me+desc Note that $\sum_{i=1}^N M^e_i = m$, where $N$ is number of nodes.
    fino.vectors.Me = wasora_define_vector("Me", fino.nev, NULL, NULL);    
    
    fino.vectors.phi = malloc(fino.nev * sizeof(vector_t *));
    for (i = 0; i < fino.nev; i++) {
      if (asprintf(&vibname, "phi%d", i+1) == -1) {
        wasora_push_error_message("cannot asprintf");
        return WASORA_RUNTIME_ERROR;
      }
      
      if ((fino.vectors.phi[i] = wasora_define_vector(vibname, 0, NULL, NULL)) == NULL) {
        wasora_push_error_message("cannot define vector %s", vibname);
        return WASORA_RUNTIME_ERROR;
      }
      free(vibname);
      
    }
    
  }
  
  // TODO: heat flux
  
  return WASORA_PARSER_OK;
}
