/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's computation of reactions
 *
 *  Copyright (C) 2016--2019 jeremy theler
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
//#include "wasora.h"
#include "fino.h"


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_reaction"
int fino_instruction_reaction(void *arg) {

  fino_reaction_t *reaction = (fino_reaction_t *)arg;
  
  node_t *node;
  element_t *element;
  
  PetscInt d, i, j, j_local, k;
  PetscInt add;
  PetscInt ncols;
  const PetscInt    *cols;
  const PetscScalar *vals;  
  PetscScalar    *u;
  double R[3] = {0, 0, 0};

  if (reaction->vector->initialized == 0) {
    wasora_vector_init(reaction->vector);
  }
  
  petsc_call(VecGetArray(fino.phi, &u));

  
  for (j = 0; j < fino.mesh->n_nodes; j++) {
    add = 0;
    for (i = 0; add == 0 && i < reaction->physical_entity->n_elements; i++) {
      element = &fino.mesh->element[reaction->physical_entity->element[i]];
      for (j_local = 0; add == 0 && j_local < element->type->nodes; j_local++) {
        if (element->node[j_local]->index_mesh == j) {
          add = 1;
        }
      }
    }
    
    if (add) {
      node = &fino.mesh->node[j];
      
      for (d = 0; d < fino.degrees; d++) {
        petsc_call(MatGetRow(fino.K_nobc, node->index_dof[d], &ncols, &cols, &vals));
        for (k = 0; k < ncols; k++) {
          R[d] += vals[k] * u[cols[k]];
        }
        petsc_call(MatRestoreRow(fino.K_nobc, node->index_dof[d], &ncols, &cols, &vals));
      }
    }
  }

  petsc_call(VecRestoreArray(fino.phi, &u));
  
  gsl_vector_set(reaction->vector->value, 0, R[0]);
  gsl_vector_set(reaction->vector->value, 1, R[1]);
  gsl_vector_set(reaction->vector->value, 2, R[2]);
  

  return WASORA_RUNTIME_OK;
}
