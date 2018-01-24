/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's boundary conditions
 *
 *  Copyright (C) 2015--2019 jeremy theler
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
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "fino_read_bcs"
int fino_read_bcs(void) {

  physical_entity_t *physical_entity;
  bc_string_based_t *bc;
  
  // barremos los physical entities y mapeamos cadenas a valores enteros
  LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
    if (physical_entity->bc_strings != NULL) {
      
      char *equal_sign = NULL;
      char *name = NULL;
      char *expr = NULL;
      
      // TODO: ver https://scicomp.stackexchange.com/questions/3298/appropriate-space-for-weak-solutions-to-an-elliptical-pde-with-mixed-inhomogeneo/3300#3300
      LL_FOREACH(physical_entity->bc_strings, bc) {
        
        // si creamos una bc dummy entonces tiene string NULL
        // pero no importa porque ya la procesamos
        if (bc->string == NULL) {
          continue;
        }
        
        // si hay signo igual hay expresion, sino no
        name = bc->string;
        if ((equal_sign = strchr(bc->string, '=')) != NULL) {
          *equal_sign = '\0';
          expr = (equal_sign+1);
        } else {
          expr = NULL;
        }

        if (fino.problem_family == problem_family_break || fino.problem_family == problem_family_shake) {
          if (strcmp(name, "fixed") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_dirichlet;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_displacement_fixed;
            
          } else if (strncmp(name, "mimic(", 6) == 0) {
            char *closing_bracket;
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_dirichlet;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_displacement_mimic;
            
            if (name[6] == 'u') {
              bc->dof = 0;
            } else if (name[6] == 'v') {
              bc->dof = 1;
            } else if (name[6] == 'w') {
              bc->dof = 2;
            } else {
              wasora_push_error_message("expected either 'u_', 'v_' or 'w_' instead of '%s' in mimic()", name+5);
              return WASORA_PARSER_ERROR;
            }
            
            if (name[7] != '_') {
              wasora_push_error_message("expected underscore after '%c' instead of '%s' in mimic()", name[6], name+5);
              return WASORA_PARSER_ERROR;
            }
            
            if ((closing_bracket = strchr(name, ')')) == NULL) {
              wasora_push_error_message("cannot find closing bracket in '%s'", name);
              return WASORA_PARSER_ERROR;
            }
            *closing_bracket = '\0';
            
            if ((bc->mimic_to = wasora_get_physical_entity_ptr(name+8)) == NULL) {
              wasora_push_error_message("unknown phyisical entity '%s'", name+8);
              return WASORA_PARSER_ERROR;
            }
            
            
          } else if (strcmp(name, "u") == 0 ||
                     strcmp(name, "v") == 0 ||
                     strcmp(name, "w") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_dirichlet;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_displacement;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
            if (strcmp(name, "u") == 0) {
              bc->dof = 0;
            } else if (strcmp(name, "v") == 0) {
              bc->dof = 1;
            } else if (strcmp(name, "w") == 0) {
              bc->dof = 2;
            }

          } else if (strcmp(name, "0") == 0 || strcmp(name, "implicit") == 0) {
            char *dummy;
            
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_dirichlet;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_displacement_constrained;
            
            // el cuento es asi: aca quisieramos que el usuario escriba algo en funcion
            // de x,y,z pero tambien de u,v y w. Pero u,v,w ya son funciones, asi que no
            // se pueden usar como variables
            // mi solucion: definir variables U,V,W y reemplazar u,v,w por U,V,W en
            // esta expresion

            // TODO: ver que haya un separador antes y despues
            // TODO: derivadas            
            dummy = expr;
            while (*dummy != '\0') {
              if (*dummy == 'u') {
                *dummy = 'U';
              } else if (*dummy == 'v') {
                *dummy = 'V';
              } else if (*dummy == 'w') {
                *dummy = 'W';
              }
              dummy++;
            }
            
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
          } else if (strcasecmp(name, "tx") == 0 ||
                     strcasecmp(name, "ty") == 0 ||
                     strcasecmp(name, "tz") == 0) {

            if (strcmp(name+1, "x") == 0) {
              bc->dof = 0;
            } else if (strcmp(name+1, "y") == 0) {
              bc->dof = 1;
            } else if (strcmp(name+1, "z") == 0) {
              bc->dof = 2;
            }
            
            if (isupper(name[0])) {
              physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
              physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_force;
            } else {
              physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
              physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_stress;
            }

            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
            if (physical_entity->bc_args == NULL) {
              physical_entity->bc_args = calloc(3, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, &physical_entity->bc_args[bc->dof]));
            

            
          } else if (strcmp(name, "p") == 0 || strcmp(name, "P") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
            if (strcmp(name, "p") == 0) {
              physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_pressure_normal;
            } else if (strcmp(name, "P") == 0) {
              physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_pressure_real;
            }
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            if (physical_entity->bc_args == NULL) {
              physical_entity->bc_args = calloc(1, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, physical_entity->bc_args));

          } else if (strcasecmp(name, "Mx") == 0 ||
                     strcasecmp(name, "My") == 0 ||
                     strcasecmp(name, "Mz") == 0) {

            if (strcmp(name+1, "x") == 0) {
              bc->dof = bc_dof_moment_offset+0;
            } else if (strcmp(name+1, "y") == 0) {
              bc->dof = bc_dof_moment_offset+1;
            } else if (strcmp(name+1, "z") == 0) {
              bc->dof = bc_dof_moment_offset+2;
            }
            
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_moment;

            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
            if (physical_entity->bc_args == NULL) {
              // 3 para los momentos y 3 para el centro (opcional)
              physical_entity->bc_args = calloc(6, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, &physical_entity->bc_args[bc->dof - bc_dof_moment_offset]));

          } else if (strcasecmp(name, "x") == 0 ||
                     strcasecmp(name, "y") == 0 ||
                     strcasecmp(name, "z") == 0) {

            if (physical_entity->bc_type_phys != bc_phys_moment) {
              wasora_push_error_message("spatial data before moment in BC for '%s'", physical_entity->name);
              return WASORA_RUNTIME_ERROR;
            }
            
            if (strcmp(name, "x") == 0) {
              bc->dof = bc_dof_coordinates_offset+0;
            } else if (strcmp(name, "y") == 0) {
              bc->dof = bc_dof_coordinates_offset+1;
            } else if (strcmp(name, "z") == 0) {
              bc->dof = bc_dof_coordinates_offset+2;
            }
            
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_moment;

            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
            if (physical_entity->bc_args == NULL) {
              wasora_push_error_message("spatial data before moment in BC for '%s'", physical_entity->name);
              return WASORA_RUNTIME_ERROR;
            }
            
            wasora_call(wasora_parse_expression(expr, &physical_entity->bc_args[3 + bc->dof - bc_dof_coordinates_offset]));
            
          } else {
            wasora_push_error_message("unknown boundary condition type '%s'", name);
            PetscFunctionReturn(WASORA_PARSER_ERROR);
          }
         
        } else if (fino.problem_family == problem_family_bake) {
          if (strcmp(name, "T") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_dirichlet;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_temperature;
            // ojo! aca se la ponemos a la bc string, en el resto a la physical entity
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
          } else if (strcmp(name, "q") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_neumann;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_heat_flux;
            if (physical_entity->bc_args == NULL) {
              physical_entity->bc_args = calloc(1, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, physical_entity->bc_args));
            
          } else if (strcmp(name, "h") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_robin;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_convection;

            if (physical_entity->bc_args == NULL) {
              physical_entity->bc_args = calloc(2, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, &physical_entity->bc_args[0]));
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            

          } else if (strcmp(name, "Tref") == 0 || strcmp(name, "Tinf") == 0) {
            physical_entity->bc_type_math = bc->bc_type_math = bc_math_robin;
            physical_entity->bc_type_phys = bc->bc_type_phys = bc_phys_convection;

            if (physical_entity->bc_args == NULL) {
              physical_entity->bc_args = calloc(2, sizeof(expr_t));
            }
            wasora_call(wasora_parse_expression(expr, &physical_entity->bc_args[1]));
            wasora_call(wasora_parse_expression(expr, &bc->expr));

          } else {
            wasora_push_error_message("unknown boundary condition type '%s'", name);
            PetscFunctionReturn(WASORA_PARSER_ERROR);
          }
        }
        
        // restauramos el signo igual porque en parametrico en una epoca pasaba de nuevo por aca vamos a volver a pasar por aca
        // ahora ya no pero por si acaso
        if (equal_sign != NULL) {
          *equal_sign = '=';
        }

      }
    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_evaluate_bc_expressions"
int fino_evaluate_bc_expressions(physical_entity_t *physical_entity, node_t *node, int degrees, double multiplier, double *result) {

  int g;
  expr_t *bc_arg = physical_entity->bc_args;

  // las cc dependen de la posicion
  mesh_update_coord_vars(node->x);

  for (g = 0; g < degrees; g++) {
    if (bc_arg == NULL) {
      wasora_push_error_message("boundary condition '%s' needs %d expressions", physical_entity->name, degrees);
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }

    result[g] = multiplier * wasora_evaluate_expression(bc_arg);
    bc_arg = bc_arg->next;
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_set_essential_bc"
int fino_set_essential_bc(Mat A, Vec b) {

  // TODO: hacer esto mas limpio!!
  
  PetscScalar *local_b;
  PetscScalar *rhs;
  PetscInt *indexes;

  physical_entity_t *physical_entity;
  element_list_item_t *associated_element;
  
  const PetscInt *cols;
  const PetscScalar *vals;
  size_t n_bcs;
  size_t current_size;
  size_t current_threshold;
  PetscInt ncols;
  Vec vec_rhs;
  int k = 0;

  int i, j, d;
  
  bc_string_based_t *bc;
  

  double n[3];
  double y1, y2;
  double h = 1e-3;

/*  
  if (wasora_var(fino.vars.dirichlet_diagonal) == 0) {
    wasora_var(fino.vars.dirichlet_diagonal) = 1;
  }
 */

  n_bcs = (fino.problem_size>999)?ceil(BC_FACTOR*fino.problem_size):fino.problem_size;
  current_size = n_bcs;
  current_threshold = n_bcs - 2*fino.degrees;

  
  rhs = malloc(n_bcs * sizeof(PetscScalar));
  indexes = malloc(n_bcs * sizeof(PetscInt));
  fino.dirichlet_row = calloc(n_bcs, sizeof(dirichlet_row_t));
  
  for (j = 0; j < fino.mesh->n_nodes; j++) {
    LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element) {
      if (associated_element->element->type->dim != fino.dimensions &&
          (physical_entity = associated_element->element->physical_entity) != NULL &&
          physical_entity->bc_type_math == bc_math_dirichlet) {
          
        if (k >= current_threshold) {
          current_size += n_bcs;
          current_threshold = current_size - 2*fino.degrees;
          indexes = realloc(indexes, current_size * sizeof(PetscInt));
          rhs = realloc(rhs, current_size * sizeof(PetscScalar));
          fino.dirichlet_row = realloc(fino.dirichlet_row, current_size * sizeof(dirichlet_row_t));
        }

        if (associated_element->element->type->dim > 1) {
          wasora_call(mesh_compute_outward_normal(associated_element->element, n));
          wasora_var_value(fino.vars.nx) = n[0];
          wasora_var_value(fino.vars.ny) = n[1];
          wasora_var_value(fino.vars.nz) = n[2];
        }
              
        wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
        wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
        wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
        
        // empezamos a ver que nos dieron
        if (physical_entity->bc_type_phys == bc_phys_displacement_fixed) {
          for (d = 0; d < fino.degrees; d++) {
            fino.dirichlet_row[k].physical_entity = physical_entity;
            fino.dirichlet_row[k].dof = d;
            indexes[k] = fino.mesh->node[j].index[d];
            rhs[k] = 0;
            k++;
          }

        } else if (physical_entity->bc_type_phys == bc_phys_displacement_mimic) {
          
          gsl_matrix *c;
          gsl_matrix *K;
          int l[2];
          int dof;
          int i, target_index;
              
          c = gsl_matrix_calloc(1, 2);
          K = gsl_matrix_calloc(2, 2);

          gsl_matrix_set(c, 0, 0, +1);
          gsl_matrix_set(c, 0, 1, -1);
          
          dof = physical_entity->bc_strings->dof;

          // lo que hay que mimicar
          target_index = -1;
          for (i = 0; i < fino.mesh->n_elements; i++) {
            if (fino.mesh->element[i].physical_entity != NULL &&
                fino.mesh->element[i].physical_entity == physical_entity->bc_strings->mimic_to) {
              target_index = fino.mesh->element[i].node[0]->index[dof];
              break;
            }
          }
          
          if (target_index == -1) {
            wasora_push_error_message("cannot find who to mimic");
            return WASORA_RUNTIME_ERROR;
          }
            
          l[0] = fino.mesh->node[j].index[dof];
          l[1] = target_index;
          
          wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, K));
          // esto lo necesitamos porque en mimic ponemos cualquier otra estructura diferente a la que ya pusimos antes
          petsc_call(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
          MatSetValues(fino.K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(K, 0, 0), ADD_VALUES);
          petsc_call(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
              
          gsl_matrix_free(c);
          gsl_matrix_free(K);

        } else if (physical_entity->bc_strings != NULL) {

          LL_FOREACH(physical_entity->bc_strings, bc) {

            if (bc->bc_type_phys == bc_phys_displacement ||
                bc->bc_type_phys == bc_phys_temperature) {
              fino.dirichlet_row[k].physical_entity = physical_entity;
              fino.dirichlet_row[k].dof = bc->dof;

              indexes[k] = fino.mesh->node[j].index[bc->dof];

              if (fino.math_type == math_type_linear && (strcmp(bc->expr.string, "0") != 0)) {
                rhs[k] = wasora_evaluate_expression(&bc->expr);
              } else {
                rhs[k] = 0;
              }

              k++;
              
            } else if (bc->bc_type_phys == bc_phys_displacement_constrained) {

              gsl_matrix *c;
              gsl_matrix *K;
              int l[3];
              
              c = gsl_matrix_calloc(1, fino.degrees);
              K = gsl_matrix_calloc(fino.degrees, fino.degrees);
              
              wasora_var_value(fino.vars.U[0]) = 0;
              wasora_var_value(fino.vars.U[1]) = 0;
              wasora_var_value(fino.vars.U[2]) = 0;
              
              for (d = 0; d < fino.degrees; d++) {
                l[d] = fino.mesh->node[j].index[d];
                
                // TODO: evaluar las derivadas con GSL
                wasora_var_value(fino.vars.U[d]) = +h;
                y1 = wasora_evaluate_expression(&bc->expr);
                wasora_var_value(fino.vars.U[d]) = -h;
                y2 = wasora_evaluate_expression(&bc->expr);
                wasora_var_value(fino.vars.U[d]) = 0;
                  
                gsl_matrix_set(c, 0, d, -(y1-y2)/(2.0*h));
              }
              
              wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, K));
              MatSetValues(fino.K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(K, 0, 0), ADD_VALUES);
              
              // TODO: RHS
              
              gsl_matrix_free(c);
              gsl_matrix_free(K);
              
              

            }
          }
        }
      }
    }
  }

  fino.n_dirichlet_rows = k;

  // esto solo si hubo constrains
  wasora_call(fino_assembly());

    
  // antes de romper las filas de dirichlet, nos las acordamos para calcular las reacciones  
  // ojo! aca estamos contando varias veces el mismo nodo, porque un nodo pertenece a varios elementos
  // TODO: hacer lo que dijo barry
  if (fino.math_type != math_type_eigen) {
    petsc_call(VecGetArray(b, &local_b));
    for (i = 0; i < fino.n_dirichlet_rows; i++) {
      petsc_call(MatGetRow(A, indexes[i], &ncols, &cols, &vals));
      fino.dirichlet_row[i].ncols = ncols;
      if (ncols != 0) {
        fino.dirichlet_row[i].cols = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscInt *));
        fino.dirichlet_row[i].vals = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscScalar *));
        // por si acaso igualamos en lugar de hacer memcpy
        for (j = 0; j < ncols; j++) {
          fino.dirichlet_row[i].cols[j] = cols[j];
          fino.dirichlet_row[i].vals[j] = vals[j];
        }
        petsc_call(MatRestoreRow(A, indexes[i], &ncols, &cols, &vals));
      } else {
        wasora_push_error_message("topology error, please check the mesh connectivity in physical entity '%s' (do you have volumetric elements?)", fino.dirichlet_row->physical_entity->name);
        PetscFunctionReturn(WASORA_RUNTIME_ERROR);
      }
    
      // el cuento es asi: hay que poner un cero en b para romper las fuerzas volumetricas
      // despues con rhs le ponemos lo que va
      local_b[indexes[i]] = 0;
    }
    petsc_call(VecRestoreArray(b, &local_b));
  }
  
  
  petsc_call(MatCreateVecs(A, &vec_rhs, NULL));
  petsc_call(VecSetValues(vec_rhs, k, indexes, rhs, INSERT_VALUES));
  petsc_call(MatZeroRowsColumns(A, k, indexes, 1.0, vec_rhs, b));
  petsc_call(VecDestroy(&vec_rhs));
  
  if (fino.math_type == math_type_eigen) {
    petsc_call(MatZeroRowsColumns(fino.M, k, indexes, 0.0, PETSC_NULL, PETSC_NULL));
  }
    
  free(indexes);
  free(rhs);

  wasora_call(fino_assembly());

  PetscFunctionReturn(WASORA_RUNTIME_OK);

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_surface_objects"
// las condiciones de contorno son
// dphi/dn = a * phi + b
// para neumann, a = 0
int fino_build_surface_objects(element_t *element, expr_t *bc_a, expr_t *bc_b) {
  int v, g;
  double w_gauss;
  expr_t *current_bc_a;
  expr_t *current_bc_b;
  gsl_vector *Nb;
  gsl_matrix *Na;
  gsl_matrix *NaH;
  double *cc;

  if (element->type->dim == 0) {

    // si es un punto es un punto y chau
    cc = malloc(fino.degrees * sizeof(double));
    
    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[0], fino.degrees, 1, cc));
    VecSetValues(fino.b, fino.degrees, element->node[0]->index, cc, INSERT_VALUES);
    
    free(cc);
    
  } else if (element->type->dim == 1 && fino.dimensions == 3) {
    // las expresiones estan dadas en "lo que sea" por unidad de longitud, asi que le damos la mitad a cada nodo
    double halflength = 0.5 * gsl_hypot3(element->node[1]->x[0] - element->node[0]->x[0],
                                         element->node[1]->x[1] - element->node[0]->x[1],
                                         element->node[1]->x[2] - element->node[0]->x[2]);
    // si es una linea y el problema es 3 hacemos como si fuesen dos puntos
    cc = malloc(fino.degrees * sizeof(double));

    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[0], fino.degrees, halflength, cc));
    VecSetValues(fino.b, fino.degrees, element->node[0]->index, cc, INSERT_VALUES);

    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[1], fino.degrees, halflength, cc));
    VecSetValues(fino.b, fino.degrees, element->node[1]->index, cc, INSERT_VALUES);
    
    free(cc);
    
  } else {
    // sino hacemos la cuenta general
  
    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }
  
    Na = gsl_matrix_calloc(fino.degrees, fino.degrees);
    NaH = gsl_matrix_calloc(fino.degrees, fino.n_local_nodes);
    Nb = gsl_vector_calloc(fino.degrees);
  
    gsl_matrix_set_zero(fino.Ki);
    gsl_vector_set_zero(fino.bi);

    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);    
      current_bc_a = bc_a;
      current_bc_b = bc_b;
      for (g = 0; g < fino.degrees; g++) {
        if (current_bc_a != NULL) {
          gsl_matrix_set(Na, g, g, -wasora_evaluate_expression(current_bc_a));
        }
        gsl_vector_set(Nb, g, wasora_evaluate_expression(current_bc_b));

        if (current_bc_a != NULL) {
          current_bc_a = current_bc_a->next;
        }
        if (current_bc_b != NULL) {
          current_bc_b = current_bc_b->next;
        }
      }

      if (bc_a != NULL) {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Na, fino.mesh->fem.H, 0, NaH);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, fino.mesh->fem.H, NaH, 1, fino.Ki);
      }
    
      gsl_blas_dgemv(CblasTrans, -w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
    }
    
    MatSetValues(fino.K, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Ki, 0, 0), ADD_VALUES);
    VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    
    gsl_vector_free(Nb);
    gsl_matrix_free(Na);
    gsl_matrix_free(NaH);
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "fino_add_single_surface_term_to_rhs"
int fino_add_single_surface_term_to_rhs(element_t *element, bc_string_based_t *bc) {
  int v;
  double w_gauss;
  gsl_vector *Nb;
  double cc;
//  double n[3];
    
  if (element->type->dim == 0) {

    // si es un punto es un punto y chau
    // las cc dependen de la posicion
    wasora_var(wasora_mesh.vars.x) = element->node[0]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[0]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[0]->x[2];    
    cc = wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[0]->index[bc->dof], &cc, ADD_VALUES);
    
  } else if (element->type->dim == 1 && fino.dimensions == 3) {
    // las expresiones estan dadas en "lo que sea" por unidad de longitud, asi que le damos la mitad a cada nodo
    double halflength = 0.5 * gsl_hypot3(element->node[1]->x[0] - element->node[0]->x[0],
                                         element->node[1]->x[1] - element->node[0]->x[1],
                                         element->node[1]->x[2] - element->node[0]->x[2]);
    
    // si es una linea y el problema es 3 hacemos como si fuesen dos puntos
    wasora_var(wasora_mesh.vars.x) = element->node[0]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[0]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[0]->x[2];    
    cc = halflength * wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[0]->index[bc->dof], &cc, ADD_VALUES);

    wasora_var(wasora_mesh.vars.x) = element->node[1]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[1]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[1]->x[2];    
    cc = halflength * wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[1]->index[bc->dof], &cc, ADD_VALUES);
    
  } else {
    // sino hacemos la cuenta general
    
    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }
  
    Nb = gsl_vector_calloc(fino.degrees);
    gsl_vector_set_zero(fino.bi);

    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);    
      gsl_vector_set(Nb, bc->dof, wasora_evaluate_expression(&bc->expr));
      gsl_blas_dgemv(CblasTrans, w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
    }
    
    VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    
    gsl_vector_free(Nb);
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


