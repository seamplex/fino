/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's routines for computing gradients of results
 *
 *  Copyright (C) 2015-2017 jeremy theler
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
#include <gsl/gsl_vector.h>
#include "fino.h"

//#define OLD

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_gradients"
int fino_compute_gradients(void) {
  
  double w_gauss, M_jj;
  double *den;
  double vol;
//  double k;
  int i, j, v;
  int m, g, l;
  int j_prime, l_prime;
  int local_j, local_jprime;
    
  element_list_item_t *associated_element;

  for (g = 0; g < fino.degrees; g++) {
    for (m = 0; m <fino.dimensions; m++) {
      free(fino.gradient[g][m]->data_value);
      fino.gradient[g][m]->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
      fino.gradient[g][m]->data_size = fino.mesh->n_nodes;
      fino.gradient[g][m]->data_argument = fino.solution[g]->data_argument;
    }
  }

  if (fino.mesh->order == 1) {

    // metodo lumped mass
    // calculamos las derivadas de todos los grados de libertad g con respecto a las dimensiones d
    // en los nodos con el metodo de matias rivero
    for (j = 0; j < fino.mesh->n_nodes; j++) {

      M_jj = 0;   // lumped mass
      LL_FOREACH (fino.mesh->node[j].associated_elements, associated_element) {
        if (associated_element->element->type->dim == fino.dimensions) {

          // buscamos a que indice local i corresponde el j global
          local_j = -1;
          for (local_jprime = 0; local_jprime < associated_element->element->type->nodes; local_jprime++) {
            if (associated_element->element->node[local_jprime]->id == j+1) {
              local_j = local_jprime;
              break;
            }
          }
          if (local_j == -1) {
            wasora_push_error_message("internal mismatch with local_j");
            return WASORA_RUNTIME_ERROR;
          }

          for (v = 0; v < associated_element->element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
            w_gauss = mesh_integration_weight(fino.mesh, associated_element->element, v);
            mesh_compute_x(associated_element->element, fino.mesh->fem.r, fino.mesh->fem.x);
            mesh_inverse(fino.mesh->bulk_dimensions, fino.mesh->fem.dxdr, fino.mesh->fem.drdx);
            mesh_compute_dhdx(associated_element->element, fino.mesh->fem.r, fino.mesh->fem.drdx, fino.mesh->fem.dhdx);

            for (local_jprime = 0; local_jprime < associated_element->element->type->nodes; local_jprime++) {
              M_jj += w_gauss * gsl_vector_get(fino.mesh->fem.h, local_jprime) * gsl_vector_get(fino.mesh->fem.h, local_j);

              for (g = 0; g < fino.degrees; g++) {
                for (m = 0; m <fino.dimensions; m++) {
                  fino.gradient[g][m]->data_value[j] += w_gauss * gsl_matrix_get(fino.mesh->fem.dhdx, local_jprime, m) * gsl_vector_get(fino.mesh->fem.h, local_j) * fino.solution[g]->data_value[associated_element->element->node[local_jprime]->id-1];
                }
              }
            }
          }
        }
      }

      for (g = 0; g < fino.degrees; g++) {
        for (m = 0; m <fino.dimensions; m++) {
          if (M_jj != 0) {
            fino.gradient[g][m]->data_value[j] /= M_jj;

            // si tenemos un gradiente base hay que sumarlo
            if (fino.base_gradient != NULL && fino.base_gradient[g] != NULL && fino.base_gradient[g][m]) {
              // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
              if (fino.base_gradient[g][m]->data_size == fino.spatial_unknowns) {
                fino.gradient[g][m]->data_value[j] += fino.base_gradient[g][m]->data_value[j];
              } else {
                fino.gradient[g][m]->data_value[j] += wasora_evaluate_function(fino.base_gradient[g][m], fino.mesh->node[j].x);
              }
            }
          }

        }
      }
    }    
    
  } else {
 
    // promediamos los valores nodales pesados con los volumenes de los elementos
    den = calloc(fino.mesh->n_nodes, sizeof(double));
//    k = 0.0*(5.0-sqrt(5))/20.0;

    for (i = 0; i < fino.mesh->n_elements; i++) {
      if (fino.mesh->element[i].type->dim == fino.dimensions) {
        vol = fino.mesh->element[i].type->element_volume(&fino.mesh->element[i]);
        for (j = 0; j < fino.mesh->element[i].type->nodes; j++) {

          gsl_vector_set(fino.mesh->fem.x, 0, fino.mesh->element[i].node[j]->x[0]);
          gsl_vector_set(fino.mesh->fem.x, 1, fino.mesh->element[i].node[j]->x[1]);
          gsl_vector_set(fino.mesh->fem.x, 2, fino.mesh->element[i].node[j]->x[2]);

          wasora_call(mesh_compute_r(&fino.mesh->element[i], fino.mesh->fem.x, fino.mesh->fem.r));

/*          
          // lo corremos un toque para no estar justo sobre el nodo
          for (m = 0; m < fino.dimensions; m++) {
            if (gsl_vector_get(fino.mesh->fem.r, m) > 0.99 || gsl_vector_get(fino.mesh->fem.r, m) < 0.01)  {
              gsl_vector_set(fino.mesh->fem.r, m, gsl_vector_get(fino.mesh->fem.r, m) + k*(0.5-gsl_vector_get(fino.mesh->fem.r, m)));
            }
          }
*/
/*          
          gsl_vector_set(fino.mesh->fem.r, 0, gsl_vector_get(fino.mesh->fem.r, 0) + k*(0.5-gsl_vector_get(fino.mesh->fem.r, 0)));
          gsl_vector_set(fino.mesh->fem.r, 1, gsl_vector_get(fino.mesh->fem.r, 1) + k*(0.5-gsl_vector_get(fino.mesh->fem.r, 1)));
          gsl_vector_set(fino.mesh->fem.r, 2, gsl_vector_get(fino.mesh->fem.r, 2) + k*(0.5-gsl_vector_get(fino.mesh->fem.r, 2)));
*/        
          
          mesh_compute_dxdr(&fino.mesh->element[i], fino.mesh->fem.r, fino.mesh->fem.dxdr);
          mesh_inverse(fino.mesh->spatial_dimensions, fino.mesh->fem.dxdr, fino.mesh->fem.drdx);
          mesh_compute_dhdx(&fino.mesh->element[i], fino.mesh->fem.r, fino.mesh->fem.drdx, fino.mesh->fem.dhdx);

          l = fino.mesh->element[i].node[j]->id - 1;
          den[l] += vol;
          for (j_prime = 0; j_prime < fino.mesh->element[i].type->nodes; j_prime++) {
            l_prime = fino.mesh->element[i].node[j_prime]->id - 1;
            for (g = 0; g < fino.degrees; g++) {
              for (m = 0; m < fino.dimensions; m++) {
                fino.gradient[g][m]->data_value[l] += vol * gsl_matrix_get(fino.mesh->fem.dhdx, j_prime, m) * fino.solution[g]->data_value[l_prime];
              }
            }
          }
        }
      }
    }

    for (j = 0; j < fino.mesh->n_nodes; j++) {
      if (den[j] != 0) {
        for (g = 0; g < fino.degrees; g++) {
          for (m = 0; m < fino.dimensions; m++) {
            fino.gradient[g][m]->data_value[j] /= den[j];
            
            // si tenemos un gradiente base hay que sumarlo
            if (fino.base_gradient != NULL && fino.base_gradient[g] != NULL && fino.base_gradient[g][m]) {
              // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
              if (fino.base_gradient[g][m]->data_size == fino.spatial_unknowns) {
                fino.gradient[g][m]->data_value[j] += fino.base_gradient[g][m]->data_value[j];
              } else {
                fino.gradient[g][m]->data_value[j] += wasora_evaluate_function(fino.base_gradient[g][m], fino.mesh->node[j].x);
             }
            } 
            
          }
        }
      }
    }
    free(den);
  }

  return WASORA_RUNTIME_OK;
}
