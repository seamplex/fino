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

int fino_build_bake(element_t *element, int v) {
  // constantes para hacer mas rapida la milonga
  static double k;
  static double Q;
  double w_gauss;
  int j;

  k = wasora_get_var_value("k");
  Q = wasora_get_var_value("Q");

  w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v); 

  // el vector de fuente de calor volumetrica
  for (j = 0; j < element->type->nodes; j++) {
    gsl_vector_add_to_element(fino.bi, j, w_gauss * gsl_vector_get(fino.mesh->fem.h, j) * Q);
  }  

  // calculamos la matriz de stiffnes 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * k, fino.mesh->fem.B, fino.mesh->fem.B, 1.0, fino.Ai);
  
  return WASORA_RUNTIME_OK;
  
}