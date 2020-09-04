/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of heat conduction problem (bake)
 *
 *  Copyright (C) 2015--2020 Seamplex, 2015 ezequiel manavela chiapero
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

#include "fino.h"

#define QX             0
#define QY             1
#define QZ             2
#define THERMAL_FUNCTIONS  3


fino_distribution_t distribution_k;     // conductivity 
fino_distribution_t distribution_Q;     // volumetric heat source
fino_distribution_t distribution_kappa; // thermal diffusivity
fino_distribution_t distribution_rho;   // density
fino_distribution_t distribution_cp;    // heat capacity

extern double hourglass_2d[];
extern double hourglass_3d[];


int fino_bc_process_thermal(bc_t **bc_pointer, char *name, char *expr, char *equal_sign) {

  int i;
  bc_t *bc = *(bc_pointer);
  bc_t *base_bc = NULL;
  bc_t *tmp = NULL;

  if (strcmp(name, "T") == 0) {
    // fixed temperature, dirichlet bondition
    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_temperature;
    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "q") == 0 || strcmp(name, "Q") == 0) {
    // fixed heat flux, neumann
    bc->type_math = bc_math_neumann;
    if (strcmp(name, "Q") == 0) {
      bc->type_phys = bc_phys_heat_total;
    } else {
      bc->type_phys = bc_phys_heat_flux;
    }
    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if ((strcmp(name, "h") == 0) ||
             (strcmp(name, "Tref") == 0) || (strcmp(name, "T_ref") == 0) ||
             (strcmp(name, "Tinf") == 0) || (strcmp(name, "T_inf") == 0)) {

    // convection, robin
    bc->type_math = bc_math_robin;
    bc->type_phys = bc_phys_convection;

    // convection needs two expressions
    base_bc = bc;
    base_bc->expr = calloc(2, sizeof(expr_t));

    do {
      // put back the equal sign, the first time is to parse again
      // the next one is not to break the string
      // the last time is fixed outside the large loop
      if (equal_sign != NULL) {
        *equal_sign = '=';
      }
      fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
      i = -1;
      if (name[0] == 'h') i = 0;
      if (name[0] == 'T') i = 1;
      if (i == -1) {
        wasora_push_error_message("expecting 'h' or 'Tref' instead of '%s'", name);
        return WASORA_PARSER_ERROR;
      }             
      wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
      tmp = bc; // esto es para "volver para atras"
    } while ((bc = bc->next) != NULL);

    // bc is now pointing to null, we need to put it back otherwise the foreach loop breaks
    *bc_pointer = tmp;

  } else {
    // TODO: radiation
    wasora_push_error_message("unknown boundary condition type '%s'", name);
    return WASORA_PARSER_ERROR;
  }
  
  
  return WASORA_RUNTIME_OK;
}




int fino_thermal_build_element(element_t *element, int v) {
  
  double k, rhocp;
  double r_for_axisymmetric;
  
  material_t *material = NULL;
  int j;

  // TODO: ver que se evaluen bien las distribuciones
  if (distribution_k.defined == 0) {
    wasora_call(fino_distribution_init(&distribution_k, "k"));
    wasora_call(fino_distribution_init(&distribution_Q, "Q"));
  }
  if (distribution_k.defined == 0) {
    wasora_push_error_message("cannot find thermal conductivity 'k'");
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.M != NULL) {
    if (distribution_kappa.defined == 0) {
      wasora_call(fino_distribution_init(&distribution_kappa, "kappa"));
      if (distribution_kappa.defined == 0) {
        wasora_call(fino_distribution_init(&distribution_rho, "rho"));
        wasora_call(fino_distribution_init(&distribution_cp, "cp"));
        if (distribution_cp.defined == 0) {
          wasora_push_error_message("cannot find neither thermal diffusivity 'kappa' nor heat capacity 'cp'");
          return WASORA_RUNTIME_ERROR;
        }
        if (distribution_rho.defined == 0) {
          wasora_push_error_message("cannot find neither thermal diffusivity 'kappa' nor density 'rho'");
          return WASORA_RUNTIME_ERROR;
        }
      }
    }
  }
  
  
  if (element->physical_entity != NULL && element->physical_entity->material != NULL) {
    material = element->physical_entity->material;
  }
  
  mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);
  mesh_compute_H_at_gauss(element, v, fino.degrees, fino.mesh->integration);
  mesh_compute_B_at_gauss(element, v, fino.degrees, fino.mesh->integration);
  mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
  r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);

  if (distribution_Q.defined != 0) {
    // the volumetric heat source term
    for (j = 0; j < element->type->nodes; j++) {
      gsl_vector_add_to_element(fino.bi, j,
        element->w[v] * r_for_axisymmetric * element->type->gauss[fino.mesh->integration].h[v][j] * fino_distribution_evaluate(&distribution_Q, material, element->x[v]));
    }
  }

  // thermal stiffness matrix
  k = fino_distribution_evaluate(&distribution_k, material, element->x[v]);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * k, element->B[v], element->B[v], 1.0, fino.Ki);

  // see if we need to do hourglass control
  if (fino.mesh->integration == integration_reduced && fino.hourglass_epsilon > 0) {  
    gsl_matrix *Gamma = NULL;  // matrix with the gamma vectors, one vector per row
    gsl_matrix_view H;         // hourglass vectors, one vector per row (not to confuse with shape functions)
    gsl_matrix *X = NULL;      // matrix with the coordinates of the nodes
    gsl_matrix *HX = NULL;     // product H*X
    gsl_matrix *BBt = NULL;    // product (B*B') (the elemental matrix is B'*B)
    gsl_matrix *Kstab;
    double *ptr_h = NULL;
    int n_h = 0;

    double trace = 0;
    double eps_tilde;
    int m;

    
    if (element->type->dim == 2 && element->type->nodes == 4) {
      // quad4
      n_h = 1;
      ptr_h = hourglass_2d;
    } else if (element->type->dim == 3 && element->type->nodes == 8 ) {
      n_h = 4;
      ptr_h = hourglass_3d;
    }
    
    if (n_h != 0) {
    
      Gamma = gsl_matrix_calloc(n_h, element->type->nodes);
      H = gsl_matrix_view_array(ptr_h, n_h, element->type->nodes);
      X = gsl_matrix_calloc(element->type->nodes, fino.dimensions);   // either this or the transpose
      HX = gsl_matrix_calloc(n_h, fino.dimensions);
      BBt = gsl_matrix_calloc(fino.dimensions, fino.dimensions);
      Kstab = gsl_matrix_alloc(element->type->nodes, element->type->nodes);

      // coordinates 
      for (j = 0; j < element->type->nodes; j++) {
        for (m = 0; m < fino.dimensions; m++) {
          gsl_matrix_set(X, j, m, element->node[j]->x[m]);
        }
      }
      printf("Ki %d =\n", element->tag);
      fino_print_gsl_matrix(fino.Ki, stdout);

      printf("H %d =\n", element->tag);
      fino_print_gsl_matrix(&H.matrix, stdout);

      printf("X %d =\n", element->tag);
      fino_print_gsl_matrix(X, stdout);

      printf("HX %d =\n", element->tag);
      fino_print_gsl_matrix(HX, stdout);
      
      gsl_matrix_memcpy(Gamma, &H.matrix);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &H.matrix, X, 0.0, HX);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, HX, element->B[v], 1.0, Gamma);


      printf("Gamma %d =\n", element->tag);
      fino_print_gsl_matrix(Gamma, stdout);
    
      // trace of B*B' for normalization
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, element->B[v], element->B[v], 1.0, BBt);
      for (m = 0; m < fino.dimensions; m++) {
        trace += gsl_matrix_get(BBt, m, m);
      }
    
      eps_tilde = 1.0/12.0 * fino.hourglass_epsilon * trace;
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, k * eps_tilde, Gamma, Gamma, 0.0, Kstab);
      
      printf("Ks %d =\n", element->tag);
      fino_print_gsl_matrix(Kstab, stdout);

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * k * eps_tilde, Gamma, Gamma, 1.0, fino.Ki);
  
      gsl_matrix_free(BBt);
      gsl_matrix_free(HX);
      gsl_matrix_free(X);

      gsl_matrix_free(Gamma);
//      gsl_matrix_free(H);
    }
  }
  

  if (fino.M != NULL) {
    // compute the mass matrix Ht*rho*cp*H
    if (distribution_kappa.defined)  {
      rhocp = k / fino_distribution_evaluate(&distribution_kappa, material, element->x[v]);
    } else {
      rhocp = fino_distribution_evaluate(&distribution_rho, material, element->x[v]);
    }
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * rhocp, element->H[v], element->H[v], 1.0, fino.Mi);
  } 
  
  return WASORA_RUNTIME_OK;
  
}


int fino_thermal_set_heat_flux(element_t *element, bc_t *bc) {
  double r_for_axisymmetric;
  double q;
  int v;

  if (bc->type_phys == bc_phys_heat_total && element->physical_entity->volume == 0) {
    wasora_push_error_message("physical group '%s' has zero volume", element->physical_entity->name);
    return WASORA_RUNTIME_ERROR;
  }
  
  for (v = 0; v < element->type->gauss[fino.mesh->integration].V; v++) {
    
    mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);
    mesh_compute_H_at_gauss(element, v, fino.degrees, fino.mesh->integration);
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    mesh_update_coord_vars(element->x[v]);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
    
    if (bc->type_phys == bc_phys_heat_total) {
      q = wasora_evaluate_expression(&bc->expr[0])/element->physical_entity->volume;
    } else {
      q = wasora_evaluate_expression(&bc->expr[0]);
    }
    
    gsl_vector_set(fino.Nb, 0, q);
    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*element->w[v], element->H[v], fino.Nb, 1.0, fino.bi); 
  }
  
  return WASORA_RUNTIME_OK;
}



int fino_thermal_set_convection(element_t *element, bc_t *bc) {
  double r_for_axisymmetric;
  double h = 0;
  double Tinf = 0;
  int v;
  gsl_matrix *Na;
  gsl_matrix *NaH;
  
  Na = gsl_matrix_calloc(fino.degrees, fino.degrees);
  NaH = gsl_matrix_calloc(fino.degrees, fino.n_local_nodes);
  
  gsl_matrix_set_zero(fino.Ki);
  gsl_vector_set_zero(fino.bi);

  for (v = 0; v < element->type->gauss[fino.mesh->integration].V; v++) {
    
    mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);
    mesh_compute_H_at_gauss(element, v, fino.degrees, fino.mesh->integration);
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    mesh_update_coord_vars(element->x[v]);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
    
    h = wasora_evaluate_expression(&bc->expr[0]);
    Tinf = wasora_evaluate_expression(&bc->expr[1]);
    
    gsl_matrix_set(Na, 0, 0, h);
    gsl_vector_set(fino.Nb, 0, h*Tinf);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Na, element->H[v], 0, NaH);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric, element->H[v], NaH, 1.0, fino.Ki);
    gsl_blas_dgemv(CblasTrans, element->w[v] * r_for_axisymmetric, element->H[v], fino.Nb, 1.0, fino.bi); 
  }

  mesh_compute_l(fino.mesh, element);
  MatSetValues(fino.K, fino.elemental_size, element->l, fino.elemental_size, element->l, gsl_matrix_ptr(fino.Ki, 0, 0), ADD_VALUES);
  // este lo hacemos afuera
//  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_matrix_free(Na);
  gsl_matrix_free(NaH);
  
  return WASORA_RUNTIME_OK;
}


int fino_thermal_compute_fluxes(void) {
  
//  double dTdx, dTdy, dTdz;
//  double delta_dTdx, delta_dTdy, delta_dTdz;
  
  double k;
  double T_max = -INFTY;
  double T_min = +INFTY;
  
  int i, j;
  int m;
  int j_global;
  int progress = 0;
  int step = 0; 
  int ascii_progress_chars = 0;

  element_t *element;  
  element_list_item_t *associated_element;
  node_t *node;
  
  
  // if k is a scalar this works for the whole routine
  if (distribution_k.variable != NULL) {
    k = fino_distribution_evaluate(&distribution_k, NULL, NULL);
    if (k < 0) {
      wasora_push_error_message("k is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
 
  // number of steps we need to make for progress bar
  if ((step = ceil((double)(fino.mesh->n_elements+fino.mesh->n_nodes)*wasora.nprocs/100.0)) < 1) {
    step = 1;
  }
  
  if (fino.gradient[0][0]->data_value == NULL) {
    // gradient vector
    for (m = 0; m < fino.dimensions; m++) {
      fino_fill_result_function(gradient[0][m]);
      fino_fill_result_function(delta_gradient[0][m]);
    }
  }
  
  
  // step 1. sweep elements and compute gradients at each node of each element
  for (i = 0; i < fino.mesh->n_elements; i++) {
    if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }

    element = &fino.mesh->element[i];
    if (element->type->dim == fino.dimensions) {
      wasora_call(fino_compute_gradients_at_nodes(fino.mesh, element));
    }
  }
  
  // step 2. sweep nodes of the output mesh
  for (j_global = 0; j_global < fino.mesh->n_nodes; j_global++) {
    if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }

    node = &fino.mesh->node[j_global];
    
    if (fino.rough == 0) {
      double *mean;
      double *current;
      double delta;
      double sum_weight = 0;
      double rel_weight = 0;
      gsl_matrix *m2 = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      int found = 0;
      
      // en iterative si no hacemos esto estamos leakando
      if (node->dphidx == NULL) {
        node->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      } else {
        gsl_matrix_set_zero(node->dphidx);
      }  
      if (node->delta_dphidx == NULL) {
        node->delta_dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      } else {
        gsl_matrix_set_zero(node->delta_dphidx);
      }
/*      
      if (node->f == NULL) {
        node->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));
      }
*/      
      
      int n = 0;
      LL_FOREACH(fino.mesh->node[j_global].associated_elements, associated_element) {
        element = associated_element->element;
        if (element->dphidx_node != NULL) {
          found = 0;
          for (j = 0; !found && j < element->type->nodes; j++) {
            if (element->node[j]->index_mesh == j_global) {

              n++;
              found = 1;
              sum_weight += element->weight;
              rel_weight = element->weight / sum_weight;
              
              // las derivadas con sus incertezas segun weldford
//              for (g = 0; g < fino.degrees; g++) {
                for (m = 0; m < fino.dimensions; m++) {
                  mean = gsl_matrix_ptr(node->dphidx, 0, m);
                  current = gsl_matrix_ptr(element->dphidx_node[j], 0, m);
                  delta = *current - *mean;
                  *mean += rel_weight * delta;
                  gsl_matrix_add_to_element(m2, 0, m, element->weight * delta * ((*current)-(*mean)));
                  gsl_matrix_set(node->delta_dphidx, 0, m, sqrt(gsl_matrix_get(m2, 0, m)/sum_weight));
                }
//              }

/*              
              if (element->property_node != NULL) {
                lambda = element->property_node[j][LAMBDA];
                mu = element->property_node[j][MU];
                alpha = element->property_node[j][ALPHA];
              }
              fino_break_compute_nodal_stresses(element, j, lambda, mu, alpha, &sigmax, &sigmay, &sigmaz, &tauxy, &tauyz, &tauzx);

              node->f[SIGMAX] += rel_weight*(sigmax - node->f[SIGMAX]);
              node->f[SIGMAY] += rel_weight*(sigmay - node->f[SIGMAY]);
              node->f[SIGMAZ] += rel_weight*(sigmaz - node->f[SIGMAZ]);
              node->f[TAUXY] += rel_weight*(tauxy - node->f[TAUXY]);
              node->f[TAUYZ] += rel_weight*(tauyz - node->f[TAUYZ]);
              node->f[TAUZX] += rel_weight*(tauzx - node->f[TAUZX]);
              
              // necesitamos el peor lambda y el peor mu para calcular el delta de von mises
              if (uniform_properties == 0) {
                if (lambda > lambda_max) {
                  lambda_max = lambda;
                }
                if (mu > mu_max) {
                  mu_max = mu;
                }
              }  
*/
              
            }
          }
        }
      }
      
      gsl_matrix_free(m2);
      
    }
    
    // ya tenemos los promedios ahora, rellenamos las funciones
/*    
    dTdx = gsl_matrix_get(node->dphidx, 0, 0);
    dTdy = (fino.dimensions > 1) ? gsl_matrix_get(node->dphidx, 0, 1) : 0;
    dTdz = (fino.dimensions > 2) ? gsl_matrix_get(node->dphidx, 0, 2) : 0;
 */
    
    fino.gradient[0][0]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 0);
    if (fino.dimensions > 1) {
      fino.gradient[0][1]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 1);
      if (fino.dimensions > 2) {
        fino.gradient[0][2]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 2);
      }
    }
    
    // incertezas
    // las incertezas
    if (fino.rough == 0) {
      fino.delta_gradient[0][0]->data_value[j_global] = gsl_matrix_get(node->delta_dphidx, 0, 0);
      if (fino.dimensions > 1) {
        fino.delta_gradient[0][1]->data_value[j_global] = gsl_matrix_get(node->delta_dphidx, 0, 1);
        if (fino.dimensions > 2) {
          fino.delta_gradient[0][2]->data_value[j_global] = gsl_matrix_get(node->delta_dphidx, 0, 2);
        }
      }  
    }
    // TODO: heat flux
    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice T_max
    if (fino.solution[0]->data_value[j_global] >= T_max) {
      T_max = fino.mesh->node[j_global].phi[0];
    }
    if (fino.solution[0]->data_value[j_global] <= T_min) {
      T_min = fino.mesh->node[j_global].phi[0];
    }
    
  }

  wasora_var(fino.vars.T_max) = T_max;
  wasora_var(fino.vars.T_min) = T_min;

  if (fino.progress_ascii == PETSC_TRUE) {
    if (wasora.nprocs == 1) {
      while (ascii_progress_chars++ < 100) {
        printf(CHAR_PROGRESS_GRADIENT);
      }
    }
    if (wasora.rank == 0) {
      printf("\n");  
      fflush(stdout);
    }  
  }
  
  
  return WASORA_RUNTIME_OK;
}

