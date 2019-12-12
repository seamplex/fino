/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's computation of reactions
 *
 *  Copyright (C) 2016--2017 jeremy theler
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
  
  var_t result;
  mesh_integrate_t mesh_integrate;
  fino_reaction_t *reaction = (fino_reaction_t *)arg;

  result.initialized = 1;
  result.value = malloc(sizeof(double));
  
  mesh_integrate.mesh = fino.mesh;
  mesh_integrate.function = NULL;
  mesh_integrate.expr.n_tokens = 0;
  mesh_integrate.physical_entity = reaction->physical_entity;
  mesh_integrate.centering = centering_nodes;
  mesh_integrate.gauss_points = 2;
  mesh_integrate.result = &result;
  
  if (fino.problem_family == problem_family_break) {
    
    if (reaction->vector->initialized == 0) {
      wasora_vector_init(reaction->vector);
    }
    
    // Rx
    wasora_call(wasora_parse_expression("nx*sigmax(x,y,z)+ny*tauxy(x,y,z)+nz*tauzx(x,y,z)", &mesh_integrate.expr));
    wasora_call(wasora_instruction_mesh_integrate(&mesh_integrate));
    gsl_vector_set(reaction->vector->value, 0, *(result.value));
    
    // Ry
    wasora_call(wasora_parse_expression("nx*tauxy(x,y,z)+ny*sigmay(x,y,z)+nz*tauyz(x,y,z)", &mesh_integrate.expr));
    wasora_call(wasora_instruction_mesh_integrate(&mesh_integrate));
    gsl_vector_set(reaction->vector->value, 1, *(result.value));
    
    // Rz
    wasora_call(wasora_parse_expression("nx*tauzx(x,y,z)+ny*tauyz(x,y,z)+nz*sigmaz(x,y,z)", &mesh_integrate.expr));
    wasora_call(wasora_instruction_mesh_integrate(&mesh_integrate));
    gsl_vector_set(reaction->vector->value, 2, *(result.value));
    
  }
  
  free(result.value);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_reactions"
int fino_break_compute_reactions(void) {

  // TODO: hacer lo que dijo barry de traer matgetsubmatrix
/*
  int g, i, j, k;
  int dirichlet;
  PetscScalar xi;
  physical_entity_t *physical_entity;
  bc_t *bc;
  

  for (physical_entity = fino.mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
    dirichlet = 0;
    LL_FOREACH(physical_entity->bcs, bc) {
      if (bc->type_math == bc_math_dirichlet) {
        dirichlet = 1;
      }
    }
    if (dirichlet) {
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
                physical_entity->M[2] -= xi * (fino.mesh->node[j].x[1] - physical_entity->cog[1]);            
                break;
              case 1:
                physical_entity->M[0] -= xi * (fino.mesh->node[j].x[2] - physical_entity->cog[2]);
                physical_entity->M[2] += xi * (fino.mesh->node[j].x[0] - physical_entity->cog[0]);
                break;
              case 2:
                physical_entity->M[0] += xi * (fino.mesh->node[j].x[1] - physical_entity->cog[1]);
                physical_entity->M[1] -= xi * (fino.mesh->node[j].x[0] - physical_entity->cog[0]);
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
*/
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}
