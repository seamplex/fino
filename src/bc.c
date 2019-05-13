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
#include <gsl/gsl_deriv.h>

#include "fino.h"

typedef struct {
  expr_t *expr;
  int dof;
} fino_gsl_function_of_uvw_params_t;


#undef  __FUNCT__
#define __FUNCT__ "fino_read_bcs"
int fino_bc_string2parsed(void) {

  physical_entity_t *physical_entity;
  bc_t *bc;
  bc_t *base_bc;
  bc_t *tmp;
  char *equal_sign;
  char *name;
  char *expr;
  int i;
  
  if (fino.mesh == NULL) {
    return WASORA_RUNTIME_OK;
  }
  
  // barremos los physical entities y mapeamos cadenas a valores enteros
  for (physical_entity = fino.mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
    // TODO: ver https://scicomp.stackexchange.com/questions/3298/appropriate-space-for-weak-solutions-to-an-elliptical-pde-with-mixed-inhomogeneo/3300#3300
    LL_FOREACH(physical_entity->bcs, bc) {
      
      fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);

      if (fino.problem_family == problem_family_break || fino.problem_family == problem_family_shake) {
        if (strcmp(name, "fixed") == 0) {
          bc->type_math = bc_math_dirichlet;
          bc->type_phys = bc_phys_displacement_fixed;

        } else if (strncmp(name, "mimic(", 6) == 0) {
          char *closing_bracket;
          bc->type_math = bc_math_dirichlet;
          bc->type_phys = bc_phys_displacement_mimic;

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

          if ((bc->slave = wasora_get_physical_entity_ptr(name+8, fino.mesh)) == NULL) {
            wasora_push_error_message("unknown phyisical entity '%s'", name+8);
            return WASORA_PARSER_ERROR;
          }


        } else if (strcmp(name, "u") == 0 ||
                   strcmp(name, "v") == 0 ||
                   strcmp(name, "w") == 0) {
          bc->type_math = bc_math_dirichlet;
          bc->type_phys = bc_phys_displacement;

          if (name[0] == 'u') bc->dof = 0;
          if (name[0] == 'v') bc->dof = 1;
          if (name[0] == 'w') bc->dof = 2;
          
          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if (strcmp(name, "0") == 0 || strcmp(name, "implicit") == 0) {
          char *dummy;

          bc->type_math = bc_math_dirichlet;
          bc->type_phys = bc_phys_displacement_constrained;

          // el cuento es asi: aca quisieramos que el usuario escriba algo en funcion
          // de x,y,z pero tambien de u,v y w. Pero u,v,w ya son funciones, asi que no
          // se pueden usar como variables
          // mi solucion: definir variables U,V,W y reemplazar u,v,w por U,V,W en
          // esta expresion

          // TODO: ver que haya un separador antes y despues
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

          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if (strcasecmp(name, "tx") == 0 ||
                   strcasecmp(name, "ty") == 0 ||
                   strcasecmp(name, "tz") == 0) {

          bc->type_math = bc_math_neumann;
          if (isupper(name[0])) {
            bc->type_phys = bc_phys_force;
          } else {
            bc->type_phys = bc_phys_stress;
          }

          if (name[1] == 'x') bc->dof = 0;
          if (name[1] == 'y') bc->dof = 1;
          if (name[1] == 'z') bc->dof = 2;

          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if (strcmp(name, "p") == 0 || strcmp(name, "P") == 0) {
          bc->type_math = bc_math_neumann;
          if (strcmp(name, "p") == 0) {
              bc->type_phys = bc_phys_pressure_normal;
          } else if (strcmp(name, "P") == 0) {
              bc->type_phys = bc_phys_pressure_real;
          }
          
          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if (strcasecmp(name, "Mx") == 0 ||
                   strcasecmp(name, "My") == 0 ||
                   strcasecmp(name, "Mz") == 0 ||
                   strcasecmp(name, "x0") == 0 ||
                   strcasecmp(name, "y0") == 0 ||
                   strcasecmp(name, "z0") == 0) {
          
          bc->type_math = bc_math_neumann;
          bc->type_phys = bc_phys_moment;
          
          // M necesita seis expresiones asi que las alocamos: Mx My Mz x0 y0 z0
          // las alocamos en la primera de las BCs
          base_bc = bc;
          base_bc->expr = calloc(6, sizeof(expr_t));

          // ahora el cuento es que barremos hasta el final de la linked list
          // y vamos parseando en expr, si alguna no aparece es cero
          do {
            // volvemos a poner el equal sign
            if (equal_sign != NULL) {
              *equal_sign = '=';
            }
            fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
            i = -1;
            if (name[1] == 'x') i = 0;
            if (name[1] == 'y') i = 1;
            if (name[1] == 'z') i = 2;
            if (name[0] == 'x') i = 3;
            if (name[0] == 'y') i = 4;
            if (name[0] == 'z') i = 5;
            if (i == -1) {
              wasora_push_error_message("expecting 'Mx', 'My', 'Mz', 'x0', 'y0' or 'z0' instead of '%s'", name);
              return WASORA_PARSER_ERROR;
            }
            wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
          } while ((bc = bc->next) != NULL);
          
        } else {
          wasora_push_error_message("unknown boundary condition type '%s'", name);
          PetscFunctionReturn(WASORA_PARSER_ERROR);
        }

      } else if (fino.problem_family == problem_family_bake) {
        if (strcmp(name, "T") == 0) {
          bc->type_math = bc_math_dirichlet;
          bc->type_phys = bc_phys_temperature;
          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if (strcmp(name, "q") == 0) {
          bc->type_math = bc_math_neumann;
          bc->type_phys = bc_phys_heat_flux;
          bc->expr = calloc(1, sizeof(expr_t));
          wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

        } else if ((strcmp(name, "h") == 0) ||
                   (strcmp(name, "Tref") == 0) ||
                   (strcmp(name, "Tinf") == 0)) {

          bc->type_math = bc_math_robin;
          bc->type_phys = bc_phys_convection;

          // conveccion necesita dos expresione
          // las alocamos en la primera de las BCs
          base_bc = bc;
          base_bc->expr = calloc(2, sizeof(expr_t));
          
          do {
            // volvemos a poner el equal sign, en la primera pasada
            // es para volver a parser, en las siguientes para no romper
            // la ultima se vuelve a arreglar fuera del loop grande
            if (equal_sign != NULL) {
              *equal_sign = '=';
            }
            fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
            i = -1;
            if (name[0] == 'h') i = 0;
            if (name[0] == 'T') i = 1;
            if (i == -1) {
              wasora_push_error_message("expecting 'h' or 'Tref' instead of '%s'", name);
              return WASORA_PARSER_ERROR;
            }             
            wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
            tmp = bc; // esto es para "volver para atras"
          } while ((bc = bc->next) != NULL);
          
          // ahora bc quedo apuntando a null, tenemos que volver para atras
          // sino el FOREACH de arriba palma
          bc = tmp;

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
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_bc_read_name_expr"
void fino_bc_read_name_expr(bc_t *bc, char **name, char **expr, char **equal_sign) {

  // si hay signo igual hay expresion, sino no
  *name = bc->string;
  if ((*equal_sign = strchr(bc->string, '=')) != NULL) {
    **equal_sign = '\0';
    *expr = (*equal_sign+1);
  } else {
    *expr = NULL;
  }
  return;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_set_essential_bc"
int fino_set_essential_bc(Mat A, Vec b) {

  PetscScalar *local_b;
  const PetscInt *cols;
  const PetscScalar *vals;
  
  physical_entity_t *physical_entity = NULL;
  physical_entity_t *physical_entity_last = NULL;
  element_list_item_t *associated_element = NULL;
  
  size_t n_bcs = 0;
  size_t current_size = 0;
  PetscInt ncols;
  Vec vec_rhs;
  int k = 0;

  int i, j, d;
  
  bc_t *bc;
  

  double n[3];

  if (fino.n_dirichlet_rows == 0) {
    n_bcs = (fino.problem_size>999)?ceil(BC_FACTOR*fino.problem_size):fino.problem_size;
    current_size = n_bcs;
    
    fino.dirichlet_indexes = calloc(n_bcs, sizeof(PetscInt));
    fino.dirichlet_rhs = calloc(n_bcs, sizeof(PetscScalar));
    fino.dirichlet_row = calloc(n_bcs, sizeof(dirichlet_row_t));
  }
  
  for (j = 0; j < fino.mesh->n_nodes; j++) {
    
//    physical_entity_last = NULL;
    
    LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element) {
      if (associated_element->element != NULL &&
          associated_element->element->type->dim < fino.dimensions &&
          associated_element->element->physical_entity != NULL &&
          associated_element->element->physical_entity != physical_entity_last) {
        LL_FOREACH(associated_element->element->physical_entity->bcs, bc) {
          if (bc->type_math == bc_math_dirichlet) {
//            physical_entity_last = physical_entity; // esto es para no pasar varias veces por lo mismo
            physical_entity = associated_element->element->physical_entity;

            if (fino.n_dirichlet_rows == 0 && k >= (current_size-16)) {
              current_size += n_bcs;
              fino.dirichlet_indexes = realloc(fino.dirichlet_indexes, current_size * sizeof(PetscInt));
              fino.dirichlet_rhs = realloc(fino.dirichlet_rhs, current_size * sizeof(PetscScalar));
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
            if (bc->type_phys == bc_phys_displacement_fixed) {
              for (d = 0; d < fino.degrees; d++) {
                fino.dirichlet_row[k].physical_entity = physical_entity;
                fino.dirichlet_row[k].dof = d;
                fino.dirichlet_indexes[k] = fino.mesh->node[j].index_dof[d];
                fino.dirichlet_rhs[k] = 0;
                k++;
              }

            } else if (bc->type_phys == bc_phys_displacement_mimic) {

              gsl_matrix *c;
              gsl_matrix *K;
              int l[2];
              int i, target_index;

              c = gsl_matrix_calloc(1, 2);
              K = gsl_matrix_calloc(2, 2);

              gsl_matrix_set(c, 0, 0, +1);
              gsl_matrix_set(c, 0, 1, -1);
          
              // lo que hay que mimicar
              target_index = -1;
              for (i = 0; i < fino.mesh->n_elements; i++) {
                if (fino.mesh->element[i].physical_entity != NULL &&
                fino.mesh->element[i].physical_entity == bc->slave) {
                  target_index = fino.mesh->element[i].node[0]->index_dof[bc->dof];
                  break;
                }
              }

              if (target_index == -1) {
                wasora_push_error_message("cannot find who to mimic");
                return WASORA_RUNTIME_ERROR;
              }

              l[0] = fino.mesh->node[j].index_dof[bc->dof];
              l[1] = target_index;

              if (l[0] != l[1]) {
                wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, K));
                // esto lo necesitamos porque en mimic ponemos cualquier otra estructura diferente a la que ya pusimos antes
                petsc_call(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
                MatSetValues(fino.K, 2, l, 2, l, gsl_matrix_ptr(K, 0, 0), ADD_VALUES);
                petsc_call(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
              }

              gsl_matrix_free(c);
              gsl_matrix_free(K);

            } else if (bc->type_phys == bc_phys_displacement ||
                       bc->type_phys == bc_phys_temperature) {
              fino.dirichlet_row[k].physical_entity = physical_entity;
              fino.dirichlet_row[k].dof = bc->dof;

              fino.dirichlet_indexes[k] = fino.mesh->node[j].index_dof[bc->dof];

              if (fino.math_type == math_type_linear && (strcmp(bc->expr[0].string, "0") != 0)) {
                fino.dirichlet_rhs[k] = wasora_evaluate_expression(&bc->expr[0]);
              } else {
                fino.dirichlet_rhs[k] = 0;
              }

              k++;

            } else if (bc->type_phys == bc_phys_displacement_constrained) {

              gsl_matrix *c;
              gsl_matrix *K;
              int l[3];
              
              gsl_function F;
              fino_gsl_function_of_uvw_params_t params;
              double result, abserr;
              double h = 1e-5;

              params.expr = bc->expr;
              
              F.function = fino_gsl_function_of_uvw;
              F.params = &params;

              c = gsl_matrix_calloc(1, fino.degrees);
              K = gsl_matrix_calloc(fino.degrees, fino.degrees);
              
              wasora_var_value(fino.vars.U[0]) = 0;
              wasora_var_value(fino.vars.U[1]) = 0;
              wasora_var_value(fino.vars.U[2]) = 0;
              

              for (d = 0; d < fino.degrees; d++) {
                l[d] = fino.mesh->node[j].index_dof[d];
                params.dof = d;
                gsl_deriv_central(&F, 0, h, &result, &abserr);
                gsl_matrix_set(c, 0, d, -result);
              }

              wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, K));
              MatSetValues(fino.K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(K, 0, 0), ADD_VALUES);

              // TODO: non-uniform RHS
              gsl_matrix_free(c);
              gsl_matrix_free(K);
            }
          }
        }
      }
    }
  }

  if (fino.n_dirichlet_rows == 0) {
    fino.n_dirichlet_rows = k;
    
    // si k == 0 esto es como hacer free
    fino.dirichlet_indexes = realloc(fino.dirichlet_indexes, fino.n_dirichlet_rows * sizeof(PetscInt));
    fino.dirichlet_rhs = realloc(fino.dirichlet_rhs, fino.n_dirichlet_rows * sizeof(PetscScalar));
    fino.dirichlet_row = realloc(fino.dirichlet_row, fino.n_dirichlet_rows * sizeof(dirichlet_row_t));
  }

  // esto se necesita solo si hubo constrains
  wasora_call(fino_assembly());
    
  // antes de romper las filas de dirichlet, nos las acordamos para calcular las reacciones  
  // ojo! aca estamos contando varias veces el mismo nodo, porque un nodo pertenece a varios elementos
  // TODO: hacer lo que dijo barry
  if (fino.problem_family == problem_family_break) {
    petsc_call(VecGetArray(b, &local_b));
    for (i = 0; i < fino.n_dirichlet_rows; i++) {
            
      petsc_call(MatGetRow(A, fino.dirichlet_indexes[i], &ncols, &cols, &vals));
//       printf("%d %d %d\n", i, ncols, fino.dirichlet_row[i].ncols);
      if (ncols != 0) {
        if (fino.dirichlet_row[i].ncols != ncols) {
          if (fino.dirichlet_row[i].ncols != 0) {
            free(fino.dirichlet_row[i].cols);
            free(fino.dirichlet_row[i].vals);
          }
          fino.dirichlet_row[i].ncols = ncols;
          fino.dirichlet_row[i].cols = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscInt *));
          fino.dirichlet_row[i].vals = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscScalar *));
        }
        
        // por si acaso igualamos en lugar de hacer memcpy
        for (j = 0; j < ncols; j++) {
          fino.dirichlet_row[i].cols[j] = cols[j];
          fino.dirichlet_row[i].vals[j] = vals[j];
        }
        
        petsc_call(MatRestoreRow(A, fino.dirichlet_indexes[i], &ncols, &cols, &vals));
        
      } else {
        wasora_push_error_message("topology error, please check the mesh connectivity in physical entity '%s' (do you have volumetric elements?)", fino.dirichlet_row->physical_entity->name);
        PetscFunctionReturn(WASORA_RUNTIME_ERROR);
      }
    
      // el cuento es asi: hay que poner un cero en b para romper las fuerzas volumetricas
      // despues con rhs le ponemos lo que va
//      printf("%d\n", fino.dirichlet_indexes[i]);
      local_b[fino.dirichlet_indexes[i]] = 0;
    }
    petsc_call(VecRestoreArray(b, &local_b));
  }
  
  
  petsc_call(MatCreateVecs(A, &vec_rhs, NULL));
  petsc_call(VecSetValues(vec_rhs, fino.n_dirichlet_rows, fino.dirichlet_indexes, fino.dirichlet_rhs, INSERT_VALUES));
  petsc_call(MatZeroRowsColumns(A, fino.n_dirichlet_rows, fino.dirichlet_indexes, 1.0, vec_rhs, b));
  petsc_call(VecDestroy(&vec_rhs));
  
  if (fino.math_type == math_type_eigen) {
    petsc_call(MatZeroRowsColumns(fino.M, fino.n_dirichlet_rows, fino.dirichlet_indexes, 0.0, PETSC_NULL, PETSC_NULL));
  }
    
  wasora_call(fino_assembly());

  PetscFunctionReturn(WASORA_RUNTIME_OK);

}



// esto sirve para calcular derivadas con GSL
double fino_gsl_function_of_uvw(double x, void *params) {

  double y;

  fino_gsl_function_of_uvw_params_t *dummy = (fino_gsl_function_of_uvw_params_t *)params;

  wasora_var_value(fino.vars.U[0]) = 0;
  wasora_var_value(fino.vars.U[1]) = 0;
  wasora_var_value(fino.vars.U[2]) = 0;
  wasora_var_value(fino.vars.U[dummy->dof]) = x;
  
  y = wasora_evaluate_expression(dummy->expr);

  if (gsl_isnan(y) || gsl_isinf(y)) {
    wasora_nan_error();
  }

  return y;

}
