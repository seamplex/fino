/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's computation of reactions
 *
 *  Copyright (C) 2016 jeremy theler
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
#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_reactions"
int fino_break_compute_reactions(void) {

  // TODO: hacer lo que dijo barry de traer matgetsubmatrix

  int g, i, j, k;
  PetscScalar xi;
  physical_entity_t *physical_entity;
  
  LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
    if (physical_entity->bc_type_math == bc_math_dirichlet) {
      physical_entity->F[0] = 0;
      physical_entity->F[1] = 0;
      physical_entity->F[2] = 0;

      physical_entity->M[0] = 0;
      physical_entity->M[1] = 0;
      physical_entity->M[2] = 0;
      
      for (i = 0; i < fino.n_dirichlet_rows; i++) {
        if (fino.dirichlet_row[i].physical_entity == physical_entity) {
          g = fino.dirichlet_row[i].dof;
          for (k = 0; k < fino.dirichlet_row[i].ncols; k++) {
            petsc_call(VecGetValues(fino.phi, 1, &fino.dirichlet_row[i].cols[k], &xi));
            xi *= fino.dirichlet_row[i].vals[k];
            physical_entity->F[g] +=  xi;
            
            if (fino.mesh->ordering == ordering_node_based) {
              j = floor(fino.dirichlet_row[i].cols[k]/fino.degrees);
            } else {
              j = fino.dirichlet_row[i].cols[k] % fino.mesh->n_nodes;
            }

            switch (g) {
              case 0:
                physical_entity->M[1] += xi * (fino.mesh->node[j].x[2] - physical_entity->cog[2]);
                physical_entity->M[2] += xi * (fino.mesh->node[j].x[1] - physical_entity->cog[1]);            
                break;
              case 1:
                physical_entity->M[0] += xi * (fino.mesh->node[j].x[2] - physical_entity->cog[2]);
                physical_entity->M[2] += xi * (fino.mesh->node[j].x[0] - physical_entity->cog[0]);
                break;
              case 2:
                physical_entity->M[0] += xi * (fino.mesh->node[j].x[1] - physical_entity->cog[1]);
                physical_entity->M[1] += xi * (fino.mesh->node[j].x[0] - physical_entity->cog[0]);
                break;
            }
          }
        }
      }
      
      if (physical_entity->vector_R0 != NULL) {
        if (!physical_entity->vector_R0->initialized) {
          wasora_call(wasora_vector_init(physical_entity->vector_R0));
        }
        
        gsl_vector_set(physical_entity->vector_R0->value, 0, physical_entity->F[0]);
        gsl_vector_set(physical_entity->vector_R0->value, 1, physical_entity->F[1]);
        gsl_vector_set(physical_entity->vector_R0->value, 2, physical_entity->F[2]);
      }

      if (physical_entity->vector_R1 != NULL) {
        if (!physical_entity->vector_R1->initialized) {
          wasora_call(wasora_vector_init(physical_entity->vector_R1));
        }
        
        gsl_vector_set(physical_entity->vector_R1->value, 0, physical_entity->M[0]);
        gsl_vector_set(physical_entity->vector_R1->value, 1, physical_entity->M[1]);
        gsl_vector_set(physical_entity->vector_R1->value, 2, physical_entity->M[2]);
      }
      
    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}
