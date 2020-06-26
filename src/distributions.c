/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's approach to handling distributions of properties
 *
 *  Copyright (C) 2015-2016 jeremy theler
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

int fino_distribution_init(fino_distribution_t *distribution, const char *name) {

  // primero intentamos con una variable
  if ((distribution->variable = wasora_get_variable_ptr(name)) != NULL) {
    distribution->defined = 1;
    return WASORA_RUNTIME_OK;
  }
  
  // despues con una propiedad (tiene que venir antes de chequear la funcion)
  HASH_FIND_STR(wasora_mesh.physical_properties, name, distribution->physical_property);
  if (distribution->physical_property != NULL) {
    distribution->defined = 1;
    return WASORA_RUNTIME_OK;
  }
  
  // si no camino, buscamos una funcion
  if ((distribution->function = wasora_get_function_ptr(name)) != NULL) {
    if (distribution->function->n_arguments != fino.dimensions) {
      wasora_push_error_message("function '%s' should have %d arguments instead of %d to be used as a distribution", distribution->function->name, fino.dimensions, distribution->function->n_arguments);
      return WASORA_RUNTIME_ERROR;
    }
    distribution->defined = 1;
    return WASORA_RUNTIME_OK;
  }

  // si no hay nada podria ser que sea opcional, dejamos que el caller maneje
  // el caso en el que distribution->defined == 0
  return WASORA_RUNTIME_OK;
}

double fino_distribution_evaluate(fino_distribution_t *distribution, material_t *material, double *x) {
  
  if (distribution->variable != NULL) {
    return wasora_var_value(distribution->variable);
    
  } else if (distribution->physical_property != NULL) {
    if (material != NULL) {
      property_data_t *property_data = NULL;
      HASH_FIND_STR(material->property_datums, distribution->physical_property->name, property_data);
      if (property_data != NULL) {
        // evaluamos la expresion del material, que es una expresion (no una funcion) de x,y,z
        wasora_var_value(wasora_mesh.vars.x) = x[0];
        if (fino.dimensions > 1) {
          wasora_var_value(wasora_mesh.vars.y) = x[1];
          if (fino.dimensions > 2) {
            wasora_var_value(wasora_mesh.vars.z) = x[2];
          }
        }
        return wasora_evaluate_expression(&property_data->expr);
      } else {
        wasora_push_error_message("cannot find property '%s' in material '%s'", distribution->physical_property->name, material->name);
        wasora_runtime_error();
      }
      
    } else {
      function_t *function;
      if ((function = wasora_get_function_ptr(distribution->physical_property->name)) == NULL) {
        wasora_push_error_message("cannot find neither property nor function '%s'", distribution->physical_property->name);
        wasora_runtime_error();
      }
    }
    
  } else if (distribution->function != NULL) {
    return wasora_evaluate_function(distribution->function, x);
    
  }
  
  return 0;
}