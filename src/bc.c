/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's boundary conditions
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


int fino_bc_string2parsed(void) {

  physical_entity_t *physical_entity;
  bc_t *bc;
  char *equal_sign;
  char *name;
  char *expr;
  
  if (fino.mesh == NULL) {
    return WASORA_RUNTIME_OK;
  }
  
  // barremos los physical entities y mapeamos cadenas a valores enteros
  for (physical_entity = fino.mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
    // TODO: ver https://scicomp.stackexchange.com/questions/3298/appropriate-space-for-weak-solutions-to-an-elliptical-pde-with-mixed-inhomogeneo/3300#3300
    LL_FOREACH(physical_entity->bcs, bc) {

      // vemos si hay un "=>" que implica que la BC tiene una condicion
      if ((name = strstr(bc->string, "=>")) != NULL) {
        *name = '\0';
        name += 2;
        wasora_call(wasora_parse_expression(bc->string, &bc->condition));
      } else {
        name = bc->string;
      }
      
      // si hay signo igual hay expresion, sino no
      if ((equal_sign = strchr(name, '=')) != NULL) {
       *equal_sign = '\0';
       expr = equal_sign+1;
      } else {
        expr = NULL;
      } 

      if (fino.problem_family == problem_family_mechanical || fino.problem_family == problem_family_modal) {
        fino_bc_process_mechanical(bc, name, expr);
      } else if (fino.problem_family == problem_family_thermal) {
        fino_bc_process_thermal(bc, name, expr, equal_sign);
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

#undef  __func__
#define __func__ "fino_set_essential_bc"
// K - stiffness matrix: needs a one in the diagonal and the value in b and keep symmetry
// M - mass matrix: needs a zero in the diagonal and the same symmetry scheme that K
// J - Jacobian matrix: needs a one in the diagonal but does not need to keep the symmetry
// b - RHS: needs to be updated when modifying K
// x - solution: the BC values are set directly in order to be used as a initial condition or guess
// r - residual: the BC values are set to the difference between the value and the solution
int fino_set_essential_bcs(Mat K, Mat M, Mat J, Vec b, Vec phi, Vec r) {

  PetscScalar diag;
  
  physical_entity_t *physical_entity = NULL;
  physical_entity_t *physical_entity_last = NULL;
  element_list_item_t *associated_element = NULL;
  bc_t *bc = NULL;

  PetscScalar     *rhs;
  PetscScalar     *diff;
  PetscInt        *indexes;
  dirichlet_row_t *row;

  double n[3] = {0, 0, 0};

  size_t n_bcs = 0;
  size_t current_size = 0;
  
  int j, d;
  int k = 0;
  
  Vec vec_rhs;

  PetscFunctionBegin;  

  if (fino.n_dirichlet_rows != 0) {
    // if we are here then we know more or less the number of BCs we need
    n_bcs = fino.n_dirichlet_rows + 24;
  } else {
    n_bcs = (fino.global_size>999)?ceil(BC_FACTOR*fino.global_size):fino.global_size;
    if (n_bcs < 32) {
      n_bcs = 32;
    }
  }  
  current_size = n_bcs;

  indexes = calloc(n_bcs, sizeof(PetscInt));
  rhs = calloc(n_bcs, sizeof(PetscScalar));
  diff = calloc(n_bcs, sizeof(PetscScalar));
  row = calloc(n_bcs, sizeof(dirichlet_row_t));
  
  for (j = fino.first_node; j < fino.last_node; j++) {

// TODO: arreglar esta logica, los nodos de alto orden terminan con una constante de penalidad diferente que los de primer orden    
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

            if (k >= (current_size-16)) {            
              current_size += n_bcs;
              indexes = realloc(indexes, current_size * sizeof(PetscInt));
              rhs = realloc(rhs, current_size * sizeof(PetscScalar));
              row = realloc(row, current_size * sizeof(dirichlet_row_t));
            }

            // TODO: see if the outward normal is needed or not
            if ((fino.dimensions - associated_element->element->type->dim) == 1 && associated_element->element->type->dim > 0) {
              wasora_call(mesh_compute_normal(associated_element->element));
              n[0] = wasora_var_value(wasora_mesh.vars.nx);
              n[1] = wasora_var_value(wasora_mesh.vars.ny);
              n[2] = wasora_var_value(wasora_mesh.vars.nz);
            } else {
              n[0] = 0;
              n[1] = 0;
              n[2] = 0;
            }

            // TODO: see if the node coordinates are needed or not
            wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
            wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
            wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
            
            // if there is a condition we evaluated it now
            if (bc->condition.n_tokens == 0 || fabs(wasora_evaluate_expression(&bc->condition)) > 1e-3) {

              // let's see what the user asked for
              // TODO: call functions inside each physics implementation file
              if (bc->type_phys == bc_phys_displacement_fixed) {  

                for (d = 0; d < fino.degrees; d++) {
                  row[k].physical_entity = physical_entity;
                  row[k].dof = d;
                  indexes[k] = fino.mesh->node[j].index_dof[d];
                  rhs[k] = 0;
                  if (phi != NULL && r != NULL) {
                    // in this case the rhs is zero so the difference is just the solution
                    petsc_call(VecGetValues(phi, 1, &k, &diff[k]));
                  }
                  k++;
                }
                
              } else if (bc->type_phys == bc_phys_displacement ||
                         bc->type_phys == bc_phys_temperature) {

                row[k].physical_entity = physical_entity;
                row[k].dof = bc->dof;

                indexes[k] = fino.mesh->node[j].index_dof[bc->dof];

                if (fino.math_type != math_type_eigen && (strcmp(bc->expr[0].string, "0") != 0)) {
                  rhs[k] = wasora_evaluate_expression(&bc->expr[0]);
                } else {
                  rhs[k] = 0;
                }
                
                if (phi != NULL && r != NULL) {
                  petsc_call(VecGetValues(phi, 1, &k, &diff[k]));
                  diff[k] -= rhs[k];
                }

                k++;
                

              } else if (bc->type_phys == bc_phys_displacement_mimic) {

                if (K != NULL) {
                  // TODO: improve this, check that the master and slave have the same number of nodes
                  gsl_matrix *c;
                  gsl_matrix *P;
                  int l[2];
                  int i, target_index;

                  c = gsl_matrix_calloc(1, 2);
                  P = gsl_matrix_calloc(2, 2);

                  gsl_matrix_set(c, 0, 0, +1);
                  gsl_matrix_set(c, 0, 1, -1);

                  // what to mimic
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
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    // esto lo necesitamos porque en mimic ponemos cualquier otra estructura diferente a la que ya pusimos antes
                    petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
                    petsc_call(MatSetValues(K, 2, l, 2, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
                    petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
                  }

                  gsl_matrix_free(c);
                  gsl_matrix_free(P);
                  
                }  
                
              } else if (bc->type_phys == bc_phys_displacement_symmetry) {

                int coordinate_direction = -1;

                // check if the outward normal coincides with a coordinate axis
                for (d = 0; d < 3; d++) {
                  if (fabs(n[d]) > 1-1e-4) {
                    coordinate_direction = d;
                  }
                }

                if (coordinate_direction != -1) {
                  // traditional dirichlet on a single DOF
                  row[k].physical_entity = physical_entity;
                  row[k].dof = coordinate_direction;
                  indexes[k] = fino.mesh->node[j].index_dof[coordinate_direction];
                  rhs[k] = 0;
                  k++;

                } else {
                  
                  // generic multi-freedom
                  if (K != NULL) {
                    
                    int l[3];
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);

                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                      gsl_matrix_set(c, 0, d, n[d]);
                    }

                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));

                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
                  }

                }  

              } else if (bc->type_phys == bc_phys_displacement_radial) {

                if (K != NULL) {
                  int l[3];
                  double x[3];
                  double eps = 1e-2;

                  for (d = 0; d < 3; d++) {
                    if (bc->expr[d].n_tokens != 0) {
                      x[d] = fino.mesh->node[j].x[d]-wasora_evaluate_expression(&bc->expr[d]);
                    } else {
                      x[d] = fino.mesh->node[j].x[d]-physical_entity->cog[d];
                    }
                  }
  
                  // x-y
                  if (fabs(x[0]) > eps && fabs(x[1]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 0, +x[1]);
                    gsl_matrix_set(c, 0, 1, -x[0]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
  
                  }
  
                  // x-z
                  if (fabs(x[0]) > eps && fabs(x[2]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 0, +x[2]);
                    gsl_matrix_set(c, 0, 2, -x[0]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
  
                  }
  
                  // y-z
                  if (fabs(x[1]) > eps && fabs(x[2]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 1, +x[2]);
                    gsl_matrix_set(c, 0, 2, -x[1]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
  
                  }
  
                } else if (bc->type_phys == bc_phys_displacement_constrained) {
  
                  gsl_matrix *c;
                  gsl_matrix *P;
                  int l[3];
  
                  gsl_function F;
                  fino_gsl_function_of_uvw_params_t params;
                  double result, abserr;
                  double h = 1e-5;
  
                  params.expr = bc->expr;
  
                  F.function = fino_gsl_function_of_uvw;
                  F.params = &params;
  
                  c = gsl_matrix_calloc(1, fino.degrees);
                  P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                  wasora_var_value(fino.vars.U[0]) = 0;
                  wasora_var_value(fino.vars.U[1]) = 0;
                  wasora_var_value(fino.vars.U[2]) = 0;
  
                  for (d = 0; d < fino.degrees; d++) {
                    l[d] = fino.mesh->node[j].index_dof[d];
                    params.dof = d;
                    gsl_deriv_central(&F, 0, h, &result, &abserr);
                    gsl_matrix_set(c, 0, d, -result);
                  }
  
                  wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                  MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES);
  
                  // TODO: non-uniform RHS
                  gsl_matrix_free(c);
                  gsl_matrix_free(P);
                }
              }  
            }  
          }
        }
      }
    }
  }

  // now we know how many rows we need to change
  if (fino.n_dirichlet_rows != k) {
    fino.n_dirichlet_rows = k;
    
    // if k == 0 this like freeing
    indexes = realloc(indexes, fino.n_dirichlet_rows * sizeof(PetscInt));
    rhs = realloc(rhs, fino.n_dirichlet_rows * sizeof(PetscScalar));
    row = realloc(row, fino.n_dirichlet_rows * sizeof(dirichlet_row_t));
  }

  
  if (K != NULL) {
    // this is needed only if there were constrains
    // TODO: raise a flag before
    wasora_call(fino_assembly());
  
    // sometimes there are free nodes with no associated volumes
    // this can trigger zeros on the diagonal and MatZeroRowsColumns complains
    // TODO: change to MatGetDiagonal
    for (k = fino.first_row; k < fino.last_row; k++) {
      petsc_call(MatGetValues(K, 1, &k, 1, &k, &diag));
      if (diag == 0) {
        petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
        petsc_call(MatSetValue(K, k, k, 1.0, INSERT_VALUES));
        petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
        wasora_call(fino_assembly());
      }  
    }
    
    // this vector holds the dirichlet values
    petsc_call(MatCreateVecs(fino.K, &vec_rhs, NULL));
    petsc_call(VecSetValues(vec_rhs, fino.n_dirichlet_rows, indexes, rhs, INSERT_VALUES));
    // TODO: scale up the diagonal!
    petsc_call(MatZeroRowsColumns(K, fino.n_dirichlet_rows, indexes, 1.0, vec_rhs, b));
    petsc_call(VecDestroy(&vec_rhs));
  }  
  
  if (M != NULL) {
    // the mass matrix is like the stiffness one but with zero instead of one
    petsc_call(MatZeroRowsColumns(M, fino.n_dirichlet_rows, indexes, 0.0, NULL, NULL));
  }  

  if (J != NULL) {
    // the jacobian is exactly one for the dirichlet values and zero otherwise
    petsc_call(MatZeroRows(J, fino.n_dirichlet_rows, indexes, 1.0, NULL, NULL));
  }  
  
  if (phi != NULL) {
    
    PetscInt locked;
    
    petsc_call(VecLockGet(NULL, &locked));
    if (locked == 0) {
      // if the solution vector is not locked, put the desired condition in the solution
      petsc_call(VecSetValues(phi, fino.n_dirichlet_rows, indexes, rhs, INSERT_VALUES));
    }
   
    
    if (r != NULL) {
      // set the difference between the current solution and the dirichlet value as the residual
      petsc_call(VecSetValues(r, fino.n_dirichlet_rows, indexes, diff, INSERT_VALUES));
    }
  }  

  
  wasora_call(fino_assembly());
  
  
  free(indexes);
  free(rhs);
  free(row);


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
