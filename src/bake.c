/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of heat conduction problem (bake)
 *
 *  Copyright (C) 2015 jeremy theler & ezequiel manavela chiapero
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

fino_distribution_t distribution_k;    // conductividad
fino_distribution_t distribution_Q;    // heat source


int fino_build_bake(element_t *element, int v) {
  // constantes para hacer mas rapida la milonga
//  static double k;
//  static double Q;
  double w_gauss;
  material_t *material;
  int j;

  if (distribution_k.defined == 0) {
    wasora_call(fino_distribution_init(&distribution_k, "k"));
    wasora_call(fino_distribution_init(&distribution_Q, "Q"));
  }
  if (distribution_k.defined == 0) {
    wasora_push_error_message("cannot find thermal conductivity 'k'");
    PetscFunctionReturn(WASORA_RUNTIME_ERROR);
  }
  
  if (element->physical_entity != NULL && element->physical_entity->material != NULL) {
    material =  element->physical_entity->material;
  } else {
    material = NULL;
  }
  
  
  w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v); 

  if (distribution_Q.defined != 0) {
    // el vector de fuente de calor volumetrica
    for (j = 0; j < element->type->nodes; j++) {
      gsl_vector_add_to_element(fino.bi, j, w_gauss * gsl_vector_get(fino.mesh->fem.h, j) * fino_distribution_evaluate(&distribution_Q, material, gsl_vector_ptr(fino.mesh->fem.x, 0)));
    }
  }

  // calculamos la matriz de stiffnes 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * fino_distribution_evaluate(&distribution_k, material, gsl_vector_ptr(fino.mesh->fem.x, 0)), fino.mesh->fem.B, fino.mesh->fem.B, 1.0, fino.Ai);
  
  return WASORA_RUNTIME_OK;
  
}