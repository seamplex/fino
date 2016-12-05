/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's handling of distributions
 *
 *  Copyright (C) 2015-2016 jeremy theler
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
#define __FUNCT__ "fino_distribution_init"
int fino_distribution_init(fino_distribution_t *distribution, const char *name) {

  PetscFunctionBegin;
  
  // primero intentamos con una variable
  if ((distribution->variable = wasora_get_variable_ptr(name)) != NULL) {
    distribution->defined = 1;
    PetscFunctionReturn(WASORA_RUNTIME_OK);
  }
  
  // si no camino, buscamos una funcion
  if ((distribution->function = wasora_get_function_ptr(name)) != NULL) {
    if (distribution->function->n_arguments != fino.dimensions) {
      wasora_push_error_message("function '%s' should have %d arguments instead of %d to be used as a distribution", distribution->function->name, fino.dimensions, distribution->function->n_arguments);
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    distribution->defined = 1;
    PetscFunctionReturn(WASORA_RUNTIME_OK);
  }

  // si no hay nada podria ser que sea opcional, dejamos que el caller maneje
  // el caso en el que distribution->defined == 0
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_distribution_evaluate"
double fino_distribution_evaluate(fino_distribution_t *distribution) {
  double x[3];
  
  PetscFunctionBegin;
  
  if (distribution->variable != NULL) {
    PetscFunctionReturn(wasora_var_value(distribution->variable));
  } else if (distribution->function != NULL) {
    // TODO: esto creo que ya lo hace alguien, y si no lo hace deberia hacerlo
    x[0] = wasora_var_value(wasora_mesh.vars.x);
    x[1] = wasora_var_value(wasora_mesh.vars.y);
    x[2] = wasora_var_value(wasora_mesh.vars.z);
    PetscFunctionReturn(wasora_evaluate_function(distribution->function, x));
  }
  
  PetscFunctionReturn(0);
}