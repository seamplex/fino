/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's parsing routines
 *
 *  Copyright (C) 2015--2017 jeremy theler
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
      
      double xi;
      
      while ((token = wasora_get_next_token(NULL)) != NULL) {

///kw+FINO_PROBLEM+usage [ BAKE |
        if (strcasecmp(token, "BAKE") == 0 || strcasecmp(token, "HEAT") == 0) {
          fino.problem_family = problem_family_bake;
          fino.math_type = math_linear;
          fino.degrees = 1;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("T");

///kw+FINO_PROBLEM+usage SHAKE |
        } else if (strcasecmp(token, "SHAKE") == 0) {
          fino.problem_family = problem_family_shake;
          fino.math_type = math_eigen;
          fino.dimensions = 3;
          fino.degrees = 3;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          fino.unknown_name[2] = strdup("w");
          
///kw+FINO_PROBLEM+usage BREAK |
        } else if (strcasecmp(token, "BREAK") == 0 || strcasecmp(token, "ELASTIC") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_full3d;
          fino.math_type = math_linear;
          fino.dimensions = 3;
          fino.degrees = 3;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          fino.unknown_name[2] = strdup("w");

///kw+FINO_PROBLEM+usage PLANE_STRESS |
        } else if (strcasecmp(token, "PLANE_STRESS") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_plane_stress;
          fino.math_type = math_linear;
          fino.dimensions = 2;
          fino.degrees = 2;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");

          ///kw+FINO_PROBLEM+usage PLANE_STRAIN ]
        } else if (strcasecmp(token, "PLANE_STRESS") == 0) {
          fino.problem_family = problem_family_break;
          fino.problem_kind = problem_kind_plane_strain;
          fino.dimensions = 2;
          fino.degrees = 2;
          fino.unknown_name = calloc(fino.degrees, sizeof(char *));
          fino.unknown_name[0] = strdup("u");
          fino.unknown_name[1] = strdup("v");
          
///kw+FINO_PROBLEM+usage [ DIMENSIONS <expr> ]
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
      
/*      
      switch (fino.problem_family) {
        case problem_shake:
          if (fino.dimensions != 0 && fino.dimensions != 3) {
            wasora_push_error_message("SHAKE works only for three-dimensional cases");
            return WASORA_PARSER_ERROR;
          } else if (fino.degrees != 0 && fino.degrees != 3) {
            wasora_push_error_message("SHAKE works only for three degrees of freedom");
            return WASORA_PARSER_ERROR;
          } else if (fino.math_type != math_eigen) {
            wasora_push_error_message("SHAKE works only for eigenvalue problem types");
            return WASORA_PARSER_ERROR;
          }
        break;
        case problem_break:
          if (fino.dimensions != 0 && fino.dimensions != 3) {
            wasora_push_error_message("BREAK works only for three-dimensional cases");
            return WASORA_PARSER_ERROR;
          } else if (fino.degrees != 0 && fino.degrees != 3) {
            wasora_push_error_message("BREAK works only for three degrees of freedom");
            return WASORA_PARSER_ERROR;
          } else if (fino.math_type != math_linear) {
            wasora_push_error_message("BREAK works only for linear problem types");
            return WASORA_PARSER_ERROR;
          }
        break;
        case problem_bake:
          if (fino.dimensions == 0) {
            wasora_push_error_message("BAKE needs a DIMENSION setting");
            return WASORA_PARSER_ERROR;
          } else if (fino.degrees != 0 && fino.degrees != 1) {
            wasora_push_error_message("BAKE works only for one degree of freedom");
            return WASORA_PARSER_ERROR;
          }
        break;
        default:
        break;
      }
*/

     
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
///kw+FINO_SOLVER+usage   mass_matrix_consistent
                         "mass_matrix_consistent",
///kw+FINO_SOLVER+usage   mass_matrix_row_sum
                         "mass_matrix_row_sum",
///kw+FINO_SOLVER+usage   mass_matrix_lobatto
                         "mass_matrix_lobatto",
///kw+FINO_SOLVER+usage   mass_matrix_diagonal
                         "mass_matrix_diagonal",
///kw+FINO_SOLVER+usage   node_average_all
                         "node_average_all",
///kw+FINO_SOLVER+usage   node_average_corner
                         "node_average_corner",
///kw+FINO_SOLVER+usage   gauss_average
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
          
///kw+FINO_SOLVER+usage [ SHMEM_PROGRESS_BUILD <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_PROGRESS_BUILD") == 0) {
          wasora_call(wasora_parser_string(&fino.shmem_progress_build_name));

///kw+FINO_SOLVER+usage [ SHMEM_PROGRESS_SOLVE <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_PROGRESS_SOLVE") == 0) {
          wasora_call(wasora_parser_string(&fino.shmem_progress_solve_name));

///kw+FINO_SOLVER+usage [ MEMORY_SHMEM <shmobject> ]
        } else if (strcasecmp(token, "SHMEM_MEMORY") == 0) {
          wasora_call(wasora_parser_string(&fino.shmem_memory_name));

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
        fino.math_type = math_linear;
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
///kw+FINO_STEP+usage [ DO_NOT_COMPUTE_GRADIENTS |
        } else if (strcasecmp(token, "DO_NOT_COMPUTE_GRADIENT") == 0 || strcasecmp(token, "DO_NOT_COMPUTE_GRADIENTS") == 0) {
          fino_step->do_not_compute_gradients = 1;
///kw+FINO_STEP+usage COMPUTE_GRADIENTS ]
        } else if (strcasecmp(token, "COMPUTE_GRADIENT") == 0 || strcasecmp(token, "COMPUTE_GRADIENTS") == 0) {
          fino_step->do_not_compute_gradients = 0;
          
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
      
      wasora_define_instruction(fino_instruction_step, fino_step);
      
      return WASORA_PARSER_OK;

// ---------------------------------------------------------------------
///kw+FINO_DEBUG+usage FINO_DEBUG
///kw+FINO_DEBUG+desc Generates debugging and benchmarking output and/or dumps the matrices into files or the screen.
    } else if ((strcasecmp(token, "FINO_DEBUG") == 0)) {
      
      debug_t *debug;
      debug = calloc(1, sizeof(debug_t));
      LL_APPEND(fino.debugs, debug);

      while ((token = wasora_get_next_token(NULL)) != NULL) {
        
///kw+FINO_DEBUG+usage [ FILE <file_id> | 
        if (strcasecmp(token, "FILE") == 0) {
          wasora_call(wasora_parser_file(&debug->file));
          
///kw+FINO_DEBUG+usage [ FILE_PATH <file_path> ]
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
      
// ---------------------------------------------------------------------
///kw+FINO_REACTION+usage FINO_REACTION
///kw+FINO_REACTION+desc Asks Fino to compute the reactions at physical entities with Dirichlet boundary conditions.
    } else if ((strcasecmp(token, "FINO_REACTION") == 0)) {

      char *buff;
      fino_reaction_t *reaction = calloc(1, sizeof(fino_reaction_t));      

      while ((token = wasora_get_next_token(NULL)) != NULL) {
///kw+FINO_REACTION+usage PHYSICAL_ENTITY <physical_entity>
        if (strcasecmp(token, "PHYSICAL_ENTITY") == 0) {
          char *name;
          wasora_call(wasora_parser_string(&name));
          if ((reaction->physical_entity = wasora_get_physical_entity_ptr(name)) == NULL) {
            wasora_push_error_message("unknown physical entity '%s'", name);
            free(name);
            return WASORA_PARSER_ERROR;
          }

///kw+FINO_REACTION+usage [ NAME_ROOT <name> ]
        } else if (strcasecmp(token, "NAME_ROOT") == 0) {
          wasora_call(wasora_parser_string(&reaction->name_root));
          
        } else {
          wasora_push_error_message("unknown keyword '%s'", token);
          return WASORA_PARSER_ERROR;
        }
      }
      
      if (reaction->physical_entity == NULL) {
        wasora_push_error_message("FINO_REACTION needs a PHYSICAL_ENTITY");
        return WASORA_PARSER_ERROR;
      }
      
      if (reaction->name_root == NULL) {
        reaction->name_root = malloc(strlen(reaction->physical_entity->name)+8);
        snprintf(reaction->name_root, strlen(reaction->physical_entity->name)+7, "R_%s", reaction->physical_entity->name);
      }
      
      if (strpbrk(reaction->name_root, factorseparators) != NULL) {
        wasora_push_error_message("FINO_REACTION cannot use '%s' as the name root", reaction->name_root);
        return WASORA_PARSER_ERROR;
      }
      
      buff = malloc(strlen(reaction->name_root)+8);

      snprintf(buff, strlen(reaction->name_root)+7, "%s_x", reaction->name_root);
      if ((reaction->R[0] = wasora_define_variable(buff)) == NULL) {
        wasora_push_error_message("cannot define variable '%s'", buff);
        return WASORA_PARSER_ERROR;
      }
      snprintf(buff, strlen(reaction->name_root)+7, "%s_y", reaction->name_root);
      if ((reaction->R[1] = wasora_define_variable(buff)) == NULL) {
        wasora_push_error_message("cannot define variable '%s'", buff);
        return WASORA_PARSER_ERROR;
      }
      snprintf(buff, strlen(reaction->name_root)+7, "%s_z", reaction->name_root);
      if ((reaction->R[2] = wasora_define_variable(buff)) == NULL) {
        wasora_push_error_message("cannot define variable '%s'", buff);
        return WASORA_PARSER_ERROR;
      }
      free(buff);
      
      LL_APPEND(fino.reactions, reaction);

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
  int g, d;
  
  // las definimos solo si ya sabemos cuantas dimensiones tiene el problema
   if (fino.dimensions == 0) {
    return WASORA_PARSER_ERROR;
  } else if (fino.degrees == 0) {
    return WASORA_PARSER_ERROR;
  }

  fino.solution = calloc(fino.degrees, sizeof(function_t *));
  fino.gradient = calloc(fino.degrees, sizeof(function_t *));

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
      fino.gradient[g][d]->mesh = fino.mesh;
      fino.gradient[g][d]->var_argument = fino.solution[g]->var_argument;
      fino.gradient[g][d]->type = type_pointwise_mesh_node;
      free(gradname);
      
    }
    
    free(name);
  }

  if (fino.problem_family == problem_family_break) {
    if ((fino.sigma = wasora_define_function("sigma", fino.dimensions)) == NULL) {
      wasora_push_error_message("sigma defined twice");
      return WASORA_RUNTIME_ERROR;
    }
    fino.sigma->mesh = fino.mesh;
    fino.sigma->var_argument = fino.solution[0]->var_argument;
    fino.sigma->type = type_pointwise_mesh_node;

    if ((fino.sigma1 = wasora_define_function("sigma1", fino.dimensions)) == NULL) {
      wasora_push_error_message("sigma1 defined twice");
      return WASORA_RUNTIME_ERROR;
    }
    fino.sigma1->mesh = fino.mesh;
    fino.sigma1->var_argument = fino.solution[0]->var_argument;
    fino.sigma1->type = type_pointwise_mesh_node;

    if ((fino.sigma2 = wasora_define_function("sigma2", fino.dimensions)) == NULL) {
      wasora_push_error_message("sigma2 defined twice");
      return WASORA_RUNTIME_ERROR;
    }
    fino.sigma2->mesh = fino.mesh;
    fino.sigma2->var_argument = fino.solution[0]->var_argument;
    fino.sigma2->type = type_pointwise_mesh_node;
    
    if ((fino.sigma3 = wasora_define_function("sigma3", fino.dimensions)) == NULL) {
      wasora_push_error_message("sigma3 defined twice");
      return WASORA_RUNTIME_ERROR;
    }
    fino.sigma3->mesh = fino.mesh;
    fino.sigma3->var_argument = fino.solution[0]->var_argument;
    fino.sigma3->type = type_pointwise_mesh_node;
    
    
  }
  // TODO: heat flux
  
  return WASORA_PARSER_OK;
}