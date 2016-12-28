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

  int i, k;
  PetscScalar xi;
  fino_reaction_t *reaction;
  
  LL_FOREACH(fino.reactions, reaction) {
    if (reaction->R[0] != NULL) {
      wasora_var_value(reaction->R[0]) = 0;
      wasora_var_value(reaction->R[1]) = 0;
      wasora_var_value(reaction->R[2]) = 0;

      for (i = 0; i < fino.n_dirichlet_rows; i++) {
        if (fino.dirichlet_row[i].physical_entity == reaction->physical_entity) {
          for (k = 0; k < fino.dirichlet_row[i].ncols; k++) {
            petsc_call(VecGetValues(fino.phi, 1, &fino.dirichlet_row[i].cols[k], &xi));
            wasora_var_value(reaction->R[fino.dirichlet_row[i].dof]) += fino.dirichlet_row[i].vals[k] * xi;
          }
        }
      }
    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}
