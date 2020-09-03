/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's computation of gradients routines
 *
 *  Copyright (C) 2020 Seamplex
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


int fino_compute_gradients_at_nodes(mesh_t *mesh, element_t *element) {
  
  int v, V;
  int j, J;
  int g, m;
  int j_global;
  int j_global_prime, j_local_prime;

  node_relative_t *parent;
  
  
  
  // TODO: choose full or actual
  V = element->type->gauss[mesh->integration].V;
  J = element->type->nodes;
      
  // if already alloced, no need to realloc
  if (element->dphidx_node == NULL) {
    element->dphidx_node = calloc(J, sizeof(gsl_matrix *));
    for (j = 0; j < J; j++) {
      element->dphidx_node[j] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
    }  
  }

  if (fino.rough == 0) {
    if (fino.gradient_element_weight == gradient_weight_volume) {
      // ok, this looks odd but still I’d rather use C than C++
      element->type->element_volume(element);
      element->weight = element->volume;
    } else if (fino.gradient_element_weight == gradient_weight_quality) {
      mesh_compute_quality(mesh, element);
      element->weight = element->quality;
    } else if (fino.gradient_element_weight == gradient_weight_volume_times_quality) {
      element->type->element_volume(element);
      mesh_compute_quality(mesh, element);
      element->weight = element->volume*GSL_MAX(element->quality, 1);;
    } else if (fino.gradient_element_weight == gradient_weight_flat) {
      element->weight = 1;
    }
  }
        
  // if we were asked to extrapolate from gauss, we compute all the nodal values
  // at once by left-multiplying the gauss values by the (possibly-rectangular) 
  // extrapolation matrix to get the nodal values
  if (fino.gradient_evaluation == gradient_gauss_extrapolated && element->type->gauss[mesh->integration].extrap != NULL) {
    
    gsl_vector *at_gauss = gsl_vector_alloc(V);
    gsl_vector *at_nodes = gsl_vector_alloc(J);
    if (element->dphidx_gauss == NULL) {
      element->dphidx_gauss = calloc(V, sizeof(gsl_matrix *));
    }  
    
    for (v = 0; v < V; v++) {
    
      if (element->dphidx_gauss[v] == NULL) {
        element->dphidx_gauss[v] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      } else {
        gsl_matrix_set_zero(element->dphidx_gauss[v]);
      }
      mesh_compute_dhdx_at_gauss(element, v, mesh->integration);

      // aca habria que hacer una matriz con los phi globales
      // (de j y g, que de paso no depende de v asi que se podria hacer afuera del for de v)
      // y ver como calcular la matriz dphidx como producto de dhdx y esta matriz
      for (g = 0; g < fino.degrees; g++) {
        for (m = 0; m < fino.dimensions; m++) {
          for (j = 0; j < element->type->nodes; j++) {
            j_global = element->node[j]->index_mesh;
            gsl_matrix_add_to_element(element->dphidx_gauss[v], g, m, gsl_matrix_get(element->dhdx[v], j, m) * mesh->node[j_global].phi[g]);
          }
        }
      }
    }
    
    // take the product of the extrapolation matrix times the values at the gauss points
    for (g = 0; g < fino.degrees; g++) {
      for (m = 0; m < fino.dimensions; m++) {
        for (v = 0; v < V; v++) {
          gsl_vector_set(at_gauss, v, gsl_matrix_get(element->dphidx_gauss[v], g, m));
        }  
        
        gsl_blas_dgemv(CblasNoTrans, 1.0, element->type->gauss[mesh->integration].extrap, at_gauss, 0, at_nodes);
        for (j = 0; j < J; j++) {
          gsl_matrix_set(element->dphidx_node[j], g, m, gsl_vector_get(at_nodes, j));
        }
      }
    }
    gsl_vector_free(at_gauss);
    gsl_vector_free(at_nodes);
    
  } else {

    for (j = 0; j < J; j++) {
      j_global = element->node[j]->index_mesh;
    
      if (element->type->node_parents != NULL && element->type->node_parents[j] != NULL && fino.gradient_highorder_nodes == gradient_average) {
        // average of parents
        double den = 0;
        LL_FOREACH(element->type->node_parents[j], parent) {
          den += 1.0;
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              gsl_matrix_add_to_element(element->dphidx_node[j], g, m, gsl_matrix_get(element->dphidx_node[parent->index], g, m));
            }  
          }
        }  
        gsl_matrix_scale(element->dphidx_node[j], 1.0/den);          
      
      } else {
        
        // direct evalution at the nodes
        gsl_matrix *dhdx = gsl_matrix_calloc(J, fino.dimensions);
        mesh_compute_dhdx(element, element->type->node_coords[j], NULL, dhdx);
      
        // las nueve derivadas (o menos)
        // TODO: como arriba, aunque hay que pelar ojo si hay menos DOFs
        for (g = 0; g < fino.degrees; g++) {
          for (m = 0; m < fino.dimensions; m++) {
            for (j_local_prime = 0; j_local_prime < J; j_local_prime++) {
              j_global_prime = element->node[j_local_prime]->index_mesh;
              gsl_matrix_add_to_element(element->dphidx_node[j], g, m, gsl_matrix_get(dhdx, j_local_prime, m) * mesh->node[j_global_prime].phi[g]);
            }
          }
        }
        gsl_matrix_free(dhdx);
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
  
}