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
#ifndef _FINO_H
#include "fino.h"
#endif


int fino_instruction_reaction(void *arg) {

  fino_reaction_t *reaction = (fino_reaction_t *)arg;
  element_t *element;
  double R[3] = {0, 0, 0};

  IS set_cols;
  IS set_rows[3];
  Mat K_row[3];
  Vec K_row_u[3];
  
  PetscInt i, j, j_local, g, k;
  PetscInt add;
  PetscInt **rows;

  if (reaction->vector->initialized == 0) {
    wasora_vector_init(reaction->vector);
  }
  
  rows = calloc(fino.degrees, sizeof(PetscInt *));
  for (g = 0; g < fino.degrees; g++) {
    rows[g] = calloc(fino.size_local, sizeof(PetscInt));
  }  
  
  k = 0;
  for (j = fino.first_node; j < fino.last_node; j++) {
    add = 0;
    // esto se podria hacer barriendo los assigned_elements y viendo si pertenece a la entidad
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

  // el index set de las columnas es el mismo para todos los dofs
  petsc_call(ISCreateStride(PETSC_COMM_WORLD, fino.size_local, fino.first_row, 1, &set_cols));

  for (g = 0; g < fino.degrees; g++) {
    // el de las filas depende del dofs
    petsc_call(ISCreateGeneral(PETSC_COMM_WORLD, k, rows[g], PETSC_USE_POINTER, &set_rows[g]));
    // intente hacer un solo K_row y reusarlo para los G dofs pero daba mal en paralelo
#if PETSC_VERSION_GT(3,8,0)
    petsc_call(MatCreateSubMatrix(fino.K_nobc, set_rows[g], set_cols, MAT_INITIAL_MATRIX, &K_row[g]));
#else
    petsc_call(MatGetSubMatrix(fino.K_nobc, set_rows[g], set_cols, MAT_INITIAL_MATRIX, &K_row[g]));
#endif
    petsc_call(MatCreateVecs(K_row[g], PETSC_NULL, &K_row_u[g]));
    petsc_call(MatMult(K_row[g], fino.phi, K_row_u[g]));
    petsc_call(VecSum(K_row_u[g], &R[g]));
    gsl_vector_set(reaction->vector->value, g, R[g]);
    petsc_call(VecDestroy(&K_row_u[g]));
    petsc_call(MatDestroy(&K_row[g]));
    petsc_call(ISDestroy(&set_rows[g]));
    free(rows[g]);
  }  
  petsc_call(ISDestroy(&set_cols));
  free(rows);
    

  
  if (fino.problem_kind == problem_kind_axisymmetric) {
    if (fino.symmetry_axis == symmetry_axis_y) {
      gsl_vector_set(reaction->vector->value, 0, 0);
      gsl_vector_set(reaction->vector->value, 1, 2*M_PI*R[1]);
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      gsl_vector_set(reaction->vector->value, 0, 2*M_PI*R[0]);
      gsl_vector_set(reaction->vector->value, 1, 0);
    }
  }
  
  return WASORA_RUNTIME_OK;
}
