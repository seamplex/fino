/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's computation of reactions
 *
 *  Copyright (C) 2016--2020 Seamplex
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
#include "fino.h"


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_reaction"
int fino_instruction_reaction(void *arg) {

  fino_reaction_t *reaction = (fino_reaction_t *)arg;
  element_t *element;
  double R[3] = {0, 0, 0};

  IS set_rows;
  IS set_cols;
  Mat K_row;
  Vec K_row_u;
  PetscScalar xi;
  
  PetscInt i, j, j_local, g, k;
  PetscInt add;
  PetscInt    **rows;

  if (reaction->vector->initialized == 0) {
    wasora_vector_init(reaction->vector);
  }
  
  rows = calloc(fino.degrees, sizeof(PetscInt *));
  for (g = 0; g < fino.degrees; g++) {
    rows[g] = calloc(fino.size_local, sizeof(PetscInt));
  }  
  
  k = 0;
  for (j = fino.first_node; j < fino.last_node; j++) {
//  for (j = 0; j < fino.mesh->n_nodes; j++) {      
    add = 0;
    for (i = 0; add == 0 && i < reaction->physical_entity->n_elements; i++) {
      element = &fino.mesh->element[reaction->physical_entity->element[i]];
      for (j_local = 0; add == 0 && j_local < element->type->nodes; j_local++) {
        if (element->node[j_local]->index_mesh == j) {
          add = 1;
          for (g = 0; g < fino.degrees; g++) {
            rows[g][k] = fino.mesh->node[j].index_dof[g];
          }
          k++;
        }
      }
    }
  }

  if (k > 0) {
    // el index set de las columnas es el mismo para todos los dofs
    printf("before %d %d %d %d\n", wasora.rank, fino.global_size, fino.first_row, fino.size_local);
    petsc_call(ISCreateStride(PETSC_COMM_WORLD, fino.size_local, fino.first_row, 1, &set_cols));
    printf("after %d %d %d %d\n", wasora.rank, fino.global_size, fino.first_row, fino.size_local);

    for (g = 0; g < fino.degrees; g++) {
      printf("%d %d\n", wasora.rank, g);
      ISCreateGeneral(PETSC_COMM_WORLD, k, rows[g], PETSC_USE_POINTER, &set_rows);
      MatCreateSubMatrix(fino.K_nobc, set_rows, set_cols, (g==0)?MAT_INITIAL_MATRIX:MAT_REUSE_MATRIX, &K_row);
//      MatSetOption(K_row, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      if (g == 0) {
        MatCreateVecs(K_row, PETSC_NULL, &K_row_u);
      }  
      MatMult(K_row, fino.phi, K_row_u);
      VecSum(K_row_u, &xi);
      R[g] += xi;
      ISDestroy(&set_rows);
    }  
    
    VecDestroy(&K_row_u);
    MatDestroy(&K_row);
    ISDestroy(&set_cols);
  }    
    
  if (fino.problem_kind == problem_kind_axisymmetric) {
    if (fino.symmetry_axis == symmetry_axis_y) {
      gsl_vector_set(reaction->vector->value, 0, 0);
      gsl_vector_set(reaction->vector->value, 1, 2*M_PI*R[1]);
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      gsl_vector_set(reaction->vector->value, 0, 2*M_PI*R[0]);
      gsl_vector_set(reaction->vector->value, 1, 0);
    }
  } else {  
    for (g = 0; g < fino.degrees; g++) {
      gsl_vector_set(reaction->vector->value, g, R[g]);
    }
  }
  
  // esto lo ponemos aca porque si k==0 hay leak
  for (g = 0; g < fino.degrees; g++) {
    free(rows[g]);
  }  
  free(rows);
  

  return WASORA_RUNTIME_OK;
}
