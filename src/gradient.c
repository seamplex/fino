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
#include <petsc.h>
#include <petscsys.h>
#include <petscksp.h>

#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_gradients"
int fino_compute_gradients(void) {
  
  double w_gauss, w_lobatto;
  double *den;
  double total_mass;
  double det;
  double vol;
  double xi;

  Mat M;
  Vec b;
  Vec x;
  PetscScalar  *data;
  KSP ksp;
  
  gsl_vector *diag_M;
  double trace_M;
  element_t *element;
  element_list_item_t *item;

  int i, v;
  int m, g;
  int j_global, j_global_prime;
  int j_local,  j_local_prime;
  int j1, j2;
  

  for (g = 0; g < fino.degrees; g++) {
    for (m = 0; m <fino.dimensions; m++) {
      free(fino.gradient[g][m]->data_value);
      fino.gradient[g][m]->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
      fino.gradient[g][m]->data_size = fino.mesh->n_nodes;
      fino.gradient[g][m]->data_argument = fino.solution[g]->data_argument;
    }
  }
 
  // defaults
  if (fino.gradient_evaluation == gradient_undefined) {
    if (fino.mesh->order > 1) {
/*      
      if (wasora_mesh.materials != NULL) {
        fino.gradient_evaluation = gradient_node_average_corner; 
      } else {
        fino.gradient_evaluation = gradient_gauss_average; 
      }
 */
      fino.gradient_evaluation = gradient_gauss_average; 
    } else {
      fino.gradient_evaluation = gradient_mass_matrix_row_sum;
    }
  } else {
    if (wasora_mesh.materials != NULL &&
         (fino.gradient_evaluation == gradient_mass_matrix_diagonal ||
          fino.gradient_evaluation == gradient_mass_matrix_consistent ||
          fino.gradient_evaluation == gradient_mass_matrix_lobatto)) {
      wasora_push_error_message("neither the mass_matrix_diagonal nor mass_matrix_consistent nor the gradient_mass_matrix_lobatto methods for GRADIENT_EVALUATION does not work with multi-part geometries");
      return WASORA_RUNTIME_ERROR;
/*
    } else if (wasora_mesh.materials != NULL && fino.mesh->order > 1 &&
          fino.gradient_evaluation == gradient_gauss_average) {
      wasora_push_error_message("the gauss_average method for GRADIENT_EVALUATION does not work with multi-part geometries and high-order meshes");
      return WASORA_RUNTIME_ERROR;
*/    
    }
  }
  
  if (fino.gradient_jacobian_threshold == 0) {
    fino.gradient_jacobian_threshold = 1e-3;
  }

  // ahora viene la milonga  
  if (fino.gradient_evaluation == gradient_none) {
    return WASORA_RUNTIME_OK;
  } else if (fino.gradient_evaluation == gradient_mass_matrix_consistent ||
             fino.gradient_evaluation == gradient_mass_matrix_row_sum ||
             fino.gradient_evaluation == gradient_mass_matrix_lobatto ||
             fino.gradient_evaluation == gradient_mass_matrix_diagonal) {
    
    if (fino.gradient_evaluation == gradient_mass_matrix_consistent) {
      
      // la matriz
      petsc_call(MatCreate(PETSC_COMM_WORLD, &M));
      petsc_call(MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, fino.mesh->n_nodes, fino.mesh->n_nodes));
      petsc_call(MatSetFromOptions(M));
      petsc_call(MatMPIAIJSetPreallocation(M, fino.mesh->max_first_neighbor_nodes, PETSC_NULL, fino.mesh->max_first_neighbor_nodes, PETSC_NULL));
      petsc_call(MatSeqAIJSetPreallocation(M, fino.mesh->max_first_neighbor_nodes, PETSC_NULL));
  
      // los vectores
      petsc_call(MatCreateVecs(M, NULL, &x));
      petsc_call(MatCreateVecs(M, NULL, &b));
      
    } else {
      // en paralelo esto deberia ser petsc
      diag_M = gsl_vector_alloc(fino.mesh->n_nodes);
      trace_M = 0;
    }
    
    for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
      LL_FOREACH (fino.mesh->node[j_global].associated_elements, item) {
        element = item->element;
        if (element->type->dim == fino.dimensions &&
            (fino.mesh->node[j_global].master_material == NULL ||       // no hay interfaces 
             fino.mesh->node[j_global].materials_list == NULL ||        // no hay materiales diferentes
             fino.mesh->node[j_global].materials_list->next == NULL ||  // todos los elementos asociados tienen un solo material
             element->physical_entity->material == fino.mesh->node[j_global].master_material) // hay una interfaz pero el material es el master (TODO: falla con 3 materiales)
           ) {
          // calcular el j_local
          if ((j_local = mesh_compute_local_node_index(element, j_global)) == -1) {
            wasora_push_error_message("cannot find local index for global %d in element %d", j_global, element->id);
            return WASORA_RUNTIME_ERROR;
          }

          for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
            w_gauss = mesh_integration_weight(fino.mesh, element, v);
            mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
            mesh_inverse(fino.mesh->bulk_dimensions, fino.mesh->fem.dxdr, fino.mesh->fem.drdx);
            mesh_compute_dhdx(element, fino.mesh->fem.r, fino.mesh->fem.drdx, fino.mesh->fem.dhdx);

            for (j_local_prime = 0; j_local_prime < element->type->nodes; j_local_prime++) {

              j_global_prime = element->node[j_local_prime]->id - 1;
              
              for (g = 0; g < fino.degrees; g++) {
                for (m = 0; m <fino.dimensions; m++) {
                  fino.gradient[g][m]->data_value[j_global] += w_gauss * gsl_vector_get(fino.mesh->fem.h, j_local) * gsl_matrix_get(fino.mesh->fem.dhdx, j_local_prime, m) * fino.solution[g]->data_value[j_global_prime];
                }
              }

              xi = w_gauss * gsl_vector_get(fino.mesh->fem.h, j_local_prime) * gsl_vector_get(fino.mesh->fem.h, j_local);
              
              if (fino.gradient_evaluation == gradient_mass_matrix_consistent) {
                petsc_call(MatSetValue(M, j_global, j_global_prime, xi, ADD_VALUES));
              } else if (fino.gradient_evaluation == gradient_mass_matrix_row_sum) {
                gsl_vector_set(diag_M, j_global, gsl_vector_get(diag_M, j_global) + xi);
              } else if (fino.gradient_evaluation == gradient_mass_matrix_diagonal) {
                if (j_global_prime == j_global) {
                  gsl_vector_set(diag_M, j_global, gsl_vector_get(diag_M, j_global) + xi);
                  trace_M += xi;
                }
              }
            }
          }
            
          if (fino.gradient_evaluation == gradient_mass_matrix_lobatto) {
            wasora_call(mesh_compute_r_at_node(element, j_local_prime, fino.mesh->fem.r));
            mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
            mesh_compute_h(element, fino.mesh->fem.r, fino.mesh->fem.h);
            mesh_compute_dxdr(element, fino.mesh->fem.r, fino.mesh->fem.dxdr);
            det = mesh_determinant(element->type->dim, fino.mesh->fem.dxdr);
            w_lobatto = ((j_local < 4)?(-1.0/20.0):(+1.0/5.0));
//            w_lobatto = 1.0/10.0;
            gsl_vector_set(diag_M, j_global, gsl_vector_get(diag_M, j_global) + w_lobatto * 1.0/6.0 * det);
          }
        }
      }
    }
    
    if (fino.gradient_evaluation == gradient_mass_matrix_diagonal) {

      total_mass = 0;
      for (i = 0; i < fino.mesh->n_elements; i++) {
        // TODO: MAL!
        if (fino.mesh->element[i].type->dim == fino.dimensions) {
          total_mass += fino.mesh->element[i].type->element_volume(&fino.mesh->element[i]);
        }
      }
      gsl_vector_scale(diag_M, total_mass/trace_M);
      
    }

    if (fino.gradient_evaluation != gradient_mass_matrix_consistent) {
      
      total_mass = 0;
      for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
        total_mass += gsl_vector_get(diag_M, j_global);
      }
//      printf("%f\n", total_mass);
    //  gsl_vector_fprintf(stdout, diag_M, "%e");
      
      
      // para row_sum y lobatto no hace falta volver a loopear en j, se puede hacer en el anterior
      // pero aca tenemos lo del gradiente y eso
      for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
        for (g = 0; g < fino.degrees; g++) {
          for (m = 0; m <fino.dimensions; m++) {
            if (gsl_vector_get(diag_M, j_global) != 0) {
              fino.gradient[g][m]->data_value[j_global] /= gsl_vector_get(diag_M, j_global);

              // si tenemos un gradiente base hay que sumarlo
              if (fino.base_gradient != NULL && fino.base_gradient[g] != NULL && fino.base_gradient[g][m]) {
                // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
                if (fino.base_gradient[g][m]->data_size == fino.spatial_unknowns) {
                  fino.gradient[g][m]->data_value[j_global] += fino.base_gradient[g][m]->data_value[j_global];
                } else {
                  fino.gradient[g][m]->data_value[j_global] += wasora_evaluate_function(fino.base_gradient[g][m], fino.mesh->node[j_global].x);
                }
              }
            }
          }
        }
      }
      
      gsl_vector_free(diag_M);
      
    } else {

      // matriz consistente! hay que calcular...      
      petsc_call(MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY));
      petsc_call(MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY));
      
      petsc_call(KSPCreate(PETSC_COMM_WORLD, &ksp));
      petsc_call(KSPSetOperators(ksp, M, M));
      
      for (g = 0; g < fino.degrees; g++) {
        for (m = 0; m <fino.dimensions; m++) {
          VecGetArray(b, &data);
          for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
            data[j_global] = fino.gradient[g][m]->data_value[j_global];
          }
          VecRestoreArray(b, &data);
          
          petsc_call(KSPSolve(ksp, b, x));

          VecGetArray(x, &data);
          for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
            fino.gradient[g][m]->data_value[j_global] = data[j_global];
          }
          VecRestoreArray(x, &data);
          
        }
      }
      
      petsc_call(KSPDestroy(&ksp));
      petsc_call(VecDestroy(&x));
      petsc_call(VecDestroy(&b));
      petsc_call(MatDestroy(&M));
      
      
    }
        
  } else if (fino.gradient_evaluation == gradient_node_average_corner ||
             fino.gradient_evaluation == gradient_node_average_all) {
    
    // promediamos los valores nodales pesados con los volumenes de los elementos
    den = calloc(fino.mesh->n_nodes, sizeof(double));

    for (i = 0; i < fino.mesh->n_elements; i++) {
      element = &fino.mesh->element[i];
      if (element->type->dim == fino.dimensions) {
        vol = element->type->element_volume(element);
        for (j_local = 0; j_local < element->type->nodes; j_local++) {

          if (element->node[j_local]->master_material == NULL ||       // no hay interfaces 
               element->node[j_local]->materials_list == NULL ||        // no hay materiales diferentes
               element->node[j_local]->materials_list->next == NULL ||  // todos los elementos asociados tienen un solo material
               element->physical_entity->material == element->node[j_local]->master_material // hay una interfaz pero el material es el master
             ) {
            
            // esto da exactamente ceros o unos
            wasora_call(mesh_compute_r_at_node(element, j_local, fino.mesh->fem.r));
          
            // TODO: esto da lo mismo para todos los nodos en primer orden
            mesh_compute_dxdr(element, fino.mesh->fem.r, fino.mesh->fem.dxdr);
            det = mesh_determinant(element->type->dim, fino.mesh->fem.dxdr);
          
            if (det > fino.gradient_jacobian_threshold) {
//              printf("%d %g\n", element->id, det);
              mesh_inverse(fino.mesh->spatial_dimensions, fino.mesh->fem.dxdr, fino.mesh->fem.drdx);
              mesh_compute_dhdx(element, fino.mesh->fem.r, fino.mesh->fem.drdx, fino.mesh->fem.dhdx);

              j_global = element->node[j_local]->id - 1;
              den[j_global] += vol;
              for (j_local_prime = 0; j_local_prime < element->type->nodes; j_local_prime++) {
                j_global_prime = element->node[j_local_prime]->id - 1;
                for (g = 0; g < fino.degrees; g++) {
                  for (m = 0; m < fino.dimensions; m++) {
                    xi = gsl_matrix_get(fino.mesh->fem.dhdx, j_local_prime, m) * fino.solution[g]->data_value[j_global_prime];
                    fino.gradient[g][m]->data_value[j_global] += vol * xi;
                  }
                }
              }
            }
          }
        }
      }
    }
   
    for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
      if (den[j_global] != 0) {
        for (g = 0; g < fino.degrees; g++) {
          for (m = 0; m < fino.dimensions; m++) {
            fino.gradient[g][m]->data_value[j_global] /= den[j_global];
            
            // si tenemos un gradiente base hay que sumarlo
            if (fino.base_gradient != NULL && fino.base_gradient[g] != NULL && fino.base_gradient[g][m]) {
              // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
              if (fino.base_gradient[g][m]->data_size == fino.spatial_unknowns) {
                fino.gradient[g][m]->data_value[j_global] += fino.base_gradient[g][m]->data_value[j_global];
              } else {
                fino.gradient[g][m]->data_value[j_global] += wasora_evaluate_function(fino.base_gradient[g][m], fino.mesh->node[j_global].x);
              }
            } 
          }
        }
      }
    }

    free(den);    
    
  } else if (fino.gradient_evaluation == gradient_gauss_average) {


    // promediamos los valores de los puntos de gauss y se los damos al nodo que esta cerquita
    den = calloc(fino.mesh->n_nodes, sizeof(double));

    for (i = 0; i < fino.mesh->n_elements; i++) {
      element = &fino.mesh->element[i];
      if (element->type->dim == fino.dimensions) {
        vol = element->type->element_volume(element);
        for (j_local = 0; j_local < element->type->nodes; j_local++) {
          
          if (element->type->dim == fino.dimensions &&
              (element->node[j_local]->master_material == NULL ||       // no hay interfaces 
               element->node[j_local]->materials_list == NULL ||        // no hay materiales diferentes
               element->node[j_local]->materials_list->next == NULL ||  // todos los elementos asociados tienen un solo material
               element->physical_entity->material == element->node[j_local]->master_material) // hay una interfaz pero el material es el master (TODO: falla con 3 materiales)
             ) {
            
            if (j_local < element->type->gauss[GAUSS_POINTS_CANONICAL].V) {
              // para los nodos principales, hacemos gauss como siempre
              w_gauss = mesh_integration_weight(fino.mesh, element, j_local);
            } else {
              // para los nodos secundarios evaluamos en el nodo
              wasora_call(mesh_compute_r_at_node(element, j_local, fino.mesh->fem.r));
            }

            mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
            mesh_inverse(fino.mesh->bulk_dimensions, fino.mesh->fem.dxdr, fino.mesh->fem.drdx);
            mesh_compute_dhdx(element, fino.mesh->fem.r, fino.mesh->fem.drdx, fino.mesh->fem.dhdx);
            mesh_compute_dxdr(element, fino.mesh->fem.r, fino.mesh->fem.dxdr);
            det = mesh_determinant(element->type->dim, fino.mesh->fem.dxdr);

            if (det > fino.gradient_jacobian_threshold) {
              j_global = element->node[j_local]->id - 1;
              den[j_global] += vol;
              for (j_local_prime = 0; j_local_prime < element->type->nodes; j_local_prime++) {
                j_global_prime = element->node[j_local_prime]->id - 1;
                for (g = 0; g < fino.degrees; g++) {
                  for (m = 0; m < fino.dimensions; m++) {
                    xi = gsl_matrix_get(fino.mesh->fem.dhdx, j_local_prime, m) * fino.solution[g]->data_value[j_global_prime];
                    fino.gradient[g][m]->data_value[j_global] += vol * xi;
                  }
                }
              }
            }
          }
        }
      }
    }
   
    for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
      if (den[j_global] != 0) {
        for (g = 0; g < fino.degrees; g++) {
          for (m = 0; m < fino.dimensions; m++) {
            fino.gradient[g][m]->data_value[j_global] /= den[j_global];
            
            // si tenemos un gradiente base hay que sumarlo
            if (fino.base_gradient != NULL && fino.base_gradient[g] != NULL && fino.base_gradient[g][m]) {
              // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
              if (fino.base_gradient[g][m]->data_size == fino.spatial_unknowns) {
                fino.gradient[g][m]->data_value[j_global] += fino.base_gradient[g][m]->data_value[j_global];
              } else {
                fino.gradient[g][m]->data_value[j_global] += wasora_evaluate_function(fino.base_gradient[g][m], fino.mesh->node[j_global].x);
              }
            } 
          }
        }
      }
    }
    
    free(den);
    
  }
  
  if (fino.gradient_evaluation == gradient_node_average_corner ||
      fino.gradient_evaluation == gradient_gauss_average) {

    int interface_flag;
   
    for (i = 0; i < fino.mesh->n_elements; i++) {
      element = &fino.mesh->element[i];
      interface_flag = 0;
      
      for (j_local = 0; j_local < element->type->nodes; j_local++) {
        if (element->node[j_local]->master_material != NULL &&
            element->node[j_local]->materials_list != NULL &&
            element->node[j_local]->materials_list->next != NULL) {
          interface_flag = 1;
        }
      }

      if (interface_flag == 0) {
        // esto del promedio lo hacemos solo si no estamos en una interfaz, sino lo dejamos asi como estaba
      
        if (element->type->id == ELEMENT_TYPE_TETRAHEDRON10) {


          for (j_local = 4; j_local < 10; j_local++) {
            j_global = element->node[j_local]->id - 1;

            switch (j_local) {
              case 4:
                j1 = 0;
                j2 = 1;
              break;
              case 5:
                j1 = 1;
                j2 = 2;
              break;
              case 6:
                j1 = 0;
                j2 = 2;
              break;
              case 7:
                j1 = 0;
                j2 = 3;
              break;
              case 8:
                j1 = 2;
                j2 = 3;
              break;
              case 9:
                j1 = 1;
                j2 = 3;
              break;
            }

            for (g = 0; g < fino.degrees; g++) {
              for (m = 0; m < fino.dimensions; m++) {

                fino.gradient[g][m]->data_value[j_global] = 0.5*(fino.gradient[g][m]->data_value[element->node[j1]->id - 1] +
                                                                 fino.gradient[g][m]->data_value[element->node[j2]->id - 1]);

              }
            }
          }
        }
      }
    }
  }  

  return WASORA_RUNTIME_OK;
}
