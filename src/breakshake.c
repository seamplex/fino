/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of linear elastic problem (break) with optional vibration (shake)
 *  and evaluation of the stress tensor out of the gradients of the displacements
 *
 *  Copyright (C) 2015--2020 Seamplex
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#define LAMBDA             0
#define MU                 1
#define ALPHA              2
#define ELASTIC_PROPERTIES 3

#define SIGMAX             0
#define SIGMAY             1
#define SIGMAZ             2
#define TAUXY              3
#define TAUYZ              4
#define TAUZX              5
#define ELASTIC_FUNCTIONS  6


#define fino_fill_result_function(fun_nam) {\
        fino.fun_nam->mesh = fino.rough==0?fino.mesh:fino.mesh_rough; \
        fino.fun_nam->data_argument = fino.fun_nam->mesh->nodes_argument;   \
        fino.fun_nam->data_size = fino.fun_nam->mesh->n_nodes; \
        fino.fun_nam->data_value = calloc(fino.fun_nam->mesh->n_nodes, sizeof(double));}
/*
#define fino_fill_result_function(fun_nam) {\
        fino.fun_nam->mesh = fino.rough==0?fino.mesh:fino.mesh_rough; \
        fino.fun_nam->var_argument = fino.solution[0]->var_argument; \
        fino.fun_nam->type = type_pointwise_mesh_node; \
        fino.fun_nam->data_argument = fino.solution[0]->data_argument;   \
        fino.fun_nam->data_size = fino.fun_nam->mesh->n_nodes; \
        fino.fun_nam->data_value = calloc(fino.mesh->n_nodes, sizeof(double));}
*/
        
fino_distribution_t distribution_E;     // modulo de young
fino_distribution_t distribution_nu;    // coef de poisson
fino_distribution_t distribution_rho;   // densidad
fino_distribution_t distribution_fx;    // fuerza volumetrica en x
fino_distribution_t distribution_fy;    // fuerza volumetrica en y
fino_distribution_t distribution_fz;    // fuerza volumetrica en z
fino_distribution_t distribution_alpha; // coeficiente de expansion termica
fino_distribution_t distribution_T;     // temperatura

// este es un escalar pero lo ponemos como dist para ver si ya lo inicializamos
fino_distribution_t distribution_T0;    // temperatura de referencia (i.e. sin deformacion)
double T0;  // este es el escalar que usamos para evaluar


// estos son globales porque puede ser que sean definidos por variables
// (es decir uniformes) entonces los evaluamos una vez y ya
double nu;
double E;
double alpha;
double DT;
double lambda;
double mu;

#undef  __FUNCT__
#define __FUNCT__ "fino_break_build_element"
int fino_break_build_element(element_t *element, int v) {

  static size_t J;            // cantidad de nodos locales
  // TODO: hacer un campo descripcion en fino_distribution_t para documentar
 
  static size_t stress_strain_size = 0;
  // matrices de la formulacion del problema
  static gsl_matrix *C = NULL;
  static gsl_matrix *B = NULL;
  // vector para calcular las tensiones termicas
  static gsl_vector *et = NULL;

  // matriz intermedia
  static gsl_matrix *CB;
  // vector intermedio
  static gsl_vector *Cet;
  
  double *h;
  gsl_matrix *dhdx;
  
  material_t *material;

  double E, nu;
  double rho;
  double c;
  double alphaDT;
  
  double r_for_axisymmetric = 1.0;
  int j;

  PetscFunctionBegin;
  
  material = (element->physical_entity != NULL)?element->physical_entity->material:NULL;
  
  mesh_compute_integration_weight_at_gauss(element, v);
  mesh_compute_dhdx_at_gauss(element, v);

  dhdx = element->dhdx[v];
  h = element->type->gauss[GAUSS_POINTS_CANONICAL].h[v];
  
  // si la matriz C de la formulacion es null entonces allocamos y
  // buscamos las distribuciones espaciales de parametros
  if (C == NULL) {
    wasora_call(fino_distribution_init(&distribution_E, "E"));
    wasora_call(fino_distribution_init(&distribution_nu, "nu"));
    wasora_call(fino_distribution_init(&distribution_rho, "rho"));
    wasora_call(fino_distribution_init(&distribution_fx, "fx"));
    wasora_call(fino_distribution_init(&distribution_fy, "fy"));
    wasora_call(fino_distribution_init(&distribution_fz, "fz"));
    wasora_call(fino_distribution_init(&distribution_alpha, "alpha"));
    wasora_call(fino_distribution_init(&distribution_T, "T"));
    
    wasora_call(fino_distribution_init(&distribution_T0, "T0"));
    if (distribution_T0.defined) {
      T0 = fino_distribution_evaluate(&distribution_T0, NULL, NULL);
    } else {
      T0 = 0;
    }
    
    if (distribution_E.defined == 0) {
      wasora_push_error_message("cannot find Young modulus 'E'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    } else if (distribution_nu.defined == 0) {
      wasora_push_error_message("cannot find Poisson coefficient 'nu'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    if (fino.math_type == math_type_eigen && distribution_rho.defined == 0) {
      wasora_push_error_message("cannot find density 'rho'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    if (fino.problem_kind == problem_kind_full3d) {
      stress_strain_size = 6;
    } else if (fino.problem_kind == problem_kind_axisymmetric) {
      stress_strain_size = 4;
    } else {
      stress_strain_size = 3;
    }

    // matriz de stress-strain
    C = gsl_matrix_calloc(stress_strain_size, stress_strain_size);
    // expansion termica
    et = gsl_vector_calloc(stress_strain_size);
    
    // si E y nu son variables, calculamos C una sola vez y ya porque no dependen del espacio
    if (distribution_E.variable != NULL && distribution_nu.variable != NULL) {
      if ((E = fino_distribution_evaluate(&distribution_E, material, NULL)) <= 0) {
        wasora_push_error_message("E is not positive (%g)", E);
        return WASORA_RUNTIME_ERROR;
      }

      nu = fino_distribution_evaluate(&distribution_nu, material, NULL);
      if (nu > 0.499) {
        wasora_push_error_message("nu is greater than 1/2");
        return WASORA_RUNTIME_ERROR;
      } else if (nu < 0) {
        wasora_push_error_message("nu is negative");
        return WASORA_RUNTIME_ERROR;
      }
      wasora_call(fino_break_compute_C(C, E, nu));
    }
    
  }
  
  if (J != element->type->nodes) {
    J = element->type->nodes;
    gsl_matrix_free(B);
    B = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);
    
    gsl_matrix_free(CB);
    CB = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);
    
    // esto lo ponemos aca porque sino es mucho lio ponerlo en otro lado
    gsl_vector_free(Cet);
    Cet = gsl_vector_alloc(stress_strain_size);
  }  
  
  // la H es la del framework fem, pero la B no es la misma 
  // porque la formulacion es reducida, i.e hace un 6x6 (o 3x3) cuando deberia ser 9x9
  
  gsl_matrix_set_zero(B);

  for (j = 0; j < J; j++) {
    if (fino.problem_kind == problem_kind_full3d) {
      gsl_matrix_set(B, 0, 3*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 3*j+1, gsl_matrix_get(dhdx, j, 1));
      
      gsl_matrix_set(B, 2, 3*j+2, gsl_matrix_get(dhdx, j, 2));
    
      gsl_matrix_set(B, 3, 3*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 3, 3*j+1, gsl_matrix_get(dhdx, j, 0));

      gsl_matrix_set(B, 4, 3*j+1, gsl_matrix_get(dhdx, j, 2));
      gsl_matrix_set(B, 4, 3*j+2, gsl_matrix_get(dhdx, j, 1));

      gsl_matrix_set(B, 5, 3*j+0, gsl_matrix_get(dhdx, j, 2));
      gsl_matrix_set(B, 5, 3*j+2, gsl_matrix_get(dhdx, j, 0));
    
    } else if (fino.problem_kind == problem_kind_axisymmetric) {

      r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
      
      // ecuacion 3.5 AFEM CH.03 sec 3.3.2 pag 3.5
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(dhdx, j, 1));

      if (fino.symmetry_axis == symmetry_axis_y) {
        gsl_matrix_set(B, 2, 2*j+0, h[j]/r_for_axisymmetric);
      } else if (fino.symmetry_axis == symmetry_axis_x) {
        gsl_matrix_set(B, 2, 2*j+1, h[j]/r_for_axisymmetric);
      }
      
      gsl_matrix_set(B, 3, 2*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 3, 2*j+1, gsl_matrix_get(dhdx, j, 0));

    } else  {
      // plane stress y plane strain son iguales
      // ecuacion 14.18 IFEM CH.14 sec 14.4.1 pag 14-11
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(dhdx, j, 1));
    
      gsl_matrix_set(B, 2, 2*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 2, 2*j+1, gsl_matrix_get(dhdx, j, 0));
      
    }
    
    if ((fino.problem_family == problem_family_break) &&
        (distribution_fx.defined != 0 || distribution_fy.defined != 0 || distribution_fz.defined != 0)) {
      // el vector de fuerzas volumetricas
      c = r_for_axisymmetric * element->w[v] * h[j];
      mesh_compute_x_at_gauss(element, v);
      if (distribution_fx.defined) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+0, c * fino_distribution_evaluate(&distribution_fx, material, element->x[v]));
      }
      if (distribution_fy.defined) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+1, c * fino_distribution_evaluate(&distribution_fy, material, element->x[v]));
      }
      if (distribution_fz.defined && fino.degrees == 3) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+2, c * fino_distribution_evaluate(&distribution_fz, material, element->x[v]));
      }
    }
    
  }
  
  // si E y nu estan dadas por variables, C es constante y no la volvemos a evaluar
  // pero si alguna es una propiedad o una funcion, es otro cantar
  if (distribution_E.variable == NULL || distribution_nu.variable == NULL) {
    mesh_compute_x_at_gauss(element, v);
    wasora_call(fino_break_compute_C(C,
        fino_distribution_evaluate(&distribution_E,  material, element->x[v]),
        fino_distribution_evaluate(&distribution_nu, material, element->x[v])));
  }

  // calculamos Bt*C*B
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, r_for_axisymmetric, C, B, 0, CB);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v], B, CB, 1.0, fino.Ki);

  // expansion termica
  if (distribution_alpha.defined != 0) {
    // TODO: no volver a evaluar si no hace falta
    // este debe ser el medio!
    mesh_compute_x_at_gauss(element, v);
    alphaDT = fino_distribution_evaluate(&distribution_alpha, material, element->x[v]);
    if (alphaDT != 0) {
      alphaDT *= fino_distribution_evaluate(&distribution_T, material, element->x[v])-T0;
      gsl_vector_set(et, 0, alphaDT);
      gsl_vector_set(et, 1, alphaDT);
      gsl_vector_set(et, 2, alphaDT);
      gsl_blas_dgemv(CblasTrans, r_for_axisymmetric, C, et, 0, Cet);
      gsl_blas_dgemv(CblasTrans, element->w[v], B, Cet, 1.0, fino.bi);
    }
  }
  
  if (fino.has_mass) {
    // calculamos la matriz de masa Ht*rho*H
    mesh_compute_x_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    rho = fino_distribution_evaluate(&distribution_rho, material, element->x[v]);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * rho, element->H[v], element->H[v], 1.0, fino.Mi);
  } 
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_C"
int fino_break_compute_C(gsl_matrix *C, double E, double nu) {
  
  double lambda, mu, lambda2mu;
  
  PetscFunctionBegin;
  
  // esto es mas elegante pero la referencia es tabla 4.3 pag 194 Bathe
  lambda = E*nu/((1+nu)*(1-2*nu));
  mu = 0.5*E/(1+nu);
  lambda2mu = lambda + 2*mu;

  if (fino.problem_kind == problem_kind_full3d) {

    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    gsl_matrix_set(C, 0, 2, lambda);

    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);
    gsl_matrix_set(C, 1, 2, lambda);

    gsl_matrix_set(C, 2, 0, lambda);
    gsl_matrix_set(C, 2, 1, lambda);
    gsl_matrix_set(C, 2, 2, lambda2mu);
  
    gsl_matrix_set(C, 3, 3, mu);
    gsl_matrix_set(C, 4, 4, mu);
    gsl_matrix_set(C, 5, 5, mu);
    
  } else if (fino.problem_kind == problem_kind_plane_stress) {
    
    double c1, c2;
    
    c1 = E/(1-nu*nu);
    c2 = nu * c1;
    gsl_matrix_set(C, 0, 0, c1);
    gsl_matrix_set(C, 0, 1, c2);
    
    gsl_matrix_set(C, 1, 0, c2);
    gsl_matrix_set(C, 1, 1, c1);

    gsl_matrix_set(C, 2, 2, c1*0.5*(1-nu));
    
  } else if (fino.problem_kind == problem_kind_plane_strain) {
    
    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    
    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);

    gsl_matrix_set(C, 2, 2, mu);
    
  } else if (fino.problem_kind == problem_kind_axisymmetric) {
    
    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    gsl_matrix_set(C, 0, 2, lambda);
    
    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);
    gsl_matrix_set(C, 1, 2, lambda);

    gsl_matrix_set(C, 2, 0, lambda);
    gsl_matrix_set(C, 2, 1, lambda);
    gsl_matrix_set(C, 2, 2, lambda2mu);

    gsl_matrix_set(C, 3, 3, mu);
    
  }

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}    

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_nodal_stresses"
int fino_break_compute_nodal_stresses(element_t *element, int j, double *sigmax, double *sigmay, double *sigmaz, double *tauxy, double *tauyz, double *tauzx) {
  
  double dudx = 0;
  double dudy = 0;
  double dudz = 0;

  double dvdx = 0;
  double dvdy = 0;
  double dvdz = 0;
  
  double dwdx = 0;
  double dwdy = 0;
  double dwdz = 0;
  
  double ex = 0;
  double ey = 0;
  double ez = 0;
  double gammaxy = 0;
  double gammayz = 0;
  double gammazx = 0;
  
  double xi, DT;
  
  // nombres lindos
  dudx = gsl_matrix_get(element->dphidx_node[j], 0, 0);
  dudy = gsl_matrix_get(element->dphidx_node[j], 0, 1);

  dvdx = gsl_matrix_get(element->dphidx_node[j], 1, 0);
  dvdy = gsl_matrix_get(element->dphidx_node[j], 1, 1);

  if (fino.dimensions == 3) {
    dudz = gsl_matrix_get(element->dphidx_node[j], 0, 2);
    dvdz = gsl_matrix_get(element->dphidx_node[j], 1, 2);

    dwdx = gsl_matrix_get(element->dphidx_node[j], 2, 0);
    dwdy = gsl_matrix_get(element->dphidx_node[j], 2, 1);
    dwdz = gsl_matrix_get(element->dphidx_node[j], 2, 2);
  }

  // el tensor de deformaciones
  ex = dudx;
  ey = dvdy;

  if (fino.problem_kind == problem_kind_full3d) {
    ez = dwdz;
  } else if (fino.problem_kind == problem_kind_axisymmetric) {
    if (fino.symmetry_axis == symmetry_axis_y) {
      // etheta = u/r
      if (element->node[j]->x[0] > 1e-6) {
        ez = element->node[j]->phi[0]/element->node[j]->x[0];
      }
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      // etheta = v/r
      if (element->node[j]->x[1] > 1e-6) {
        ez = element->node[j]->phi[1]/element->node[j]->x[1];
      }
    }
  }

  gammaxy = dudy + dvdx;
  if (fino.problem_kind == problem_kind_full3d) {
    gammayz = dvdz + dwdy;
    gammazx = dwdx + dudz;
  }

  // ya tenemos derivadas y strains, ahora las tensiones
  // primero vemos si tenemos que recalcular E y/o nu
  if (element->property_node != NULL) {
    lambda = element->property_node[j][LAMBDA];
    mu = element->property_node[j][MU];
    alpha = element->property_node[j][ALPHA];
  }

  // tensiones normales
  xi = ex + ey + ez;
  *sigmax = lambda * xi + 2*mu * ex;
  *sigmay = lambda * xi + 2*mu * ey;
  *sigmaz = lambda * xi + 2*mu * ez;  // esta es sigmatheta en axi

  // restamos la contribucion termica porque nos interesan las tensiones mecanicas ver IFEM.Ch30
  if (alpha != 0) {
    DT = fino_distribution_evaluate(&distribution_T, element->physical_entity->material, fino.mesh->node[j].x) - T0;
    xi = E/(1-2*nu) * alpha * DT;

    *sigmax -= xi;
    *sigmay -= xi;
    *sigmaz -= xi;
  }

  // esfuerzos de corte
  *tauxy =  mu * gammaxy;
  if (fino.dimensions == 3) {
    *tauyz =  mu * gammayz;
    *tauzx =  mu * gammazx;
  } else {
    *tauyz = 0;
    *tauzx = 0;
  }
  
  
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_stresses"
int fino_break_compute_stresses(void) {
  
  double sigmax = 0;
  double sigmay = 0;
  double sigmaz = 0;
  double tauxy = 0;
  double tauyz = 0;
  double tauzx = 0;
  double sigma = 0;
  double sigma1 = 0;
  double sigma2 = 0;
  double sigma3 = 0;
  double tresca = 0;
  
  double displ2 = 0;
  double max_displ2 = 0;
  
  double weight_total, weight_normalized;
  
  int v, V;
  int i, j, g, m, n;
  int j_global, j_global_prime;
  int j_local_prime;
  int progress = 0;
  int step = 0; 
  int ascii_progress_chars = 0;

  mesh_t *mesh;
  node_relative_t *parent;
  element_t *element;  
  element_list_item_t *associated_element;
  node_t *node;
  
  PetscFunctionBegin;
  
  mesh = (fino.rough == 0) ? fino.mesh : fino.mesh_rough;
  
  
  // depende de si es gauss o node
  if (fino.gradient_evaluation == gradient_gauss_extrapolated) {
    step = ceil((double)(2*mesh->n_elements+mesh->n_nodes)/100.0);
  } else if (fino.gradient_evaluation == gradient_at_nodes) {
    step = ceil((double)(mesh->n_elements+mesh->n_nodes)/100.0);
  }  
  if (step < 1) {
    step = 1;
  }
  
  
  if (fino.sigma->data_value == NULL) {
    // gradiente
    for (g = 0; g < fino.degrees; g++) {
      for (m = 0; m < fino.dimensions; m++) {
        // derivada del grado de libertad g con respecto a la dimension m
        fino_fill_result_function(gradient[g][m]);
      }
    }
    
    // tensor de tensiones
    fino_fill_result_function(sigmax);
    fino_fill_result_function(sigmay);
    fino_fill_result_function(sigmaz);
    fino_fill_result_function(tauxy);
    
    if (fino.dimensions == 3) {
      fino_fill_result_function(tauyz);
      fino_fill_result_function(tauzx);
    }

    // tensiones principales
    fino_fill_result_function(sigma1);
    fino_fill_result_function(sigma2);
    fino_fill_result_function(sigma3);

    // von mises
    fino_fill_result_function(sigma);
    
    // tresca
    fino_fill_result_function(tresca);
  }
  
  // evaluamos nu, E y alpha, si son uniformes esto ya nos sirve para siempre
  if (distribution_nu.variable != NULL ) {
    nu = fino_distribution_evaluate(&distribution_nu, NULL, NULL);
    if (nu > 0.499) {
      wasora_push_error_message("nu is greater than 1/2");
      return WASORA_RUNTIME_ERROR;
    } else if (nu < 0) {
      wasora_push_error_message("nu is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
  if (distribution_E.variable != NULL) {
    if ((E = fino_distribution_evaluate(&distribution_E, NULL, NULL)) <= 0) {
      wasora_push_error_message("E is not positive (%g)", E);
      return WASORA_RUNTIME_ERROR;
    }
  }
  if (distribution_alpha.variable != NULL) {
    alpha = fino_distribution_evaluate(&distribution_alpha, NULL, NULL);
  }
  
  if (distribution_E.variable != NULL && distribution_nu.variable != NULL) {
    // esto ya sirve para toda la cosecha  
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = 0.5*E/(1+nu);
  }  


  // paso 0 (solo si es gauss extrapolate): calculamos las derivadas en los puntos de gauss
  if (fino.gradient_evaluation == gradient_gauss_extrapolated) {
    for (i = 0; i < mesh->n_elements; i++) {
      if (fino.progress_ascii && (progress++ % step) == 0) {
        printf(CHAR_PROGRESS_GRADIENT);  
        fflush(stdout);
        ascii_progress_chars++;
      }
      
      element = &mesh->element[i];
      if (element->type->dim == fino.dimensions){
        
        V = element->type->gauss[GAUSS_POINTS_CANONICAL].V;
        element->dphidx_gauss = calloc(V, sizeof(gsl_matrix *));
        
        for (v = 0; v < V; v++) {
        
          element->dphidx_gauss[v] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
          mesh_compute_dhdx_at_gauss(element, v);

          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              for (j = 0; j < element->type->nodes; j++) {
                j_global_prime = element->node[j]->index_mesh;
                gsl_matrix_add_to_element(element->dphidx_gauss[v], g, m, gsl_matrix_get(element->dhdx[v], j, m) * mesh->node[j_global_prime].phi[g]);
              }
            }
          }
        }
      }
    }  
  }  
  
  // paso 1. barremos elementos y calculamos los tensores en cada nodo de cada elemento
  for (i = 0; i < mesh->n_elements; i++) {
    if (fino.progress_ascii && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }
    
    element = &mesh->element[i];
    if (element->type->dim == fino.dimensions) {

      element->dphidx_node = calloc(element->type->nodes, sizeof(gsl_matrix *));
      V = element->type->gauss[GAUSS_POINTS_CANONICAL].V;

      if (fino.rough == 0) {
        if (fino.gradient_element_weight == gradient_weight_volume) {
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
      
      // si nu, E y/o alpha no son uniformes, los tenemos que evaluar en los nodos
      if (distribution_E.function != NULL || distribution_E.physical_property != NULL ||
          distribution_nu.function != NULL || distribution_nu.physical_property != NULL ||
          distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL) {
        element->property_node = calloc(element->type->nodes, sizeof(double *));
      }
      
      for (j = 0; j < element->type->nodes; j++) {
      
        element->dphidx_node[j] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
        j_global = element->node[j]->index_mesh;
        
        if (fino.gradient_evaluation == gradient_gauss_extrapolated && j < V && element->type->gauss[GAUSS_POINTS_CANONICAL].extrap != NULL) {
          gsl_vector *inner = gsl_vector_alloc(V);
          gsl_vector *outer = gsl_vector_alloc(V);
          
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              
              for (v = 0; v < V; v++) {
                gsl_vector_set(inner, v, gsl_matrix_get(element->dphidx_gauss[v], g, m));
              }  
                
              gsl_blas_dgemv(CblasNoTrans, 1.0, element->type->gauss[GAUSS_POINTS_CANONICAL].extrap, inner, 0, outer);
              gsl_matrix_set(element->dphidx_node[j], g, m, gsl_vector_get(outer, j));
              
            }
          }
          gsl_vector_free(inner);
          gsl_vector_free(outer);
          
        } else if (element->type->node_parents != NULL && element->type->node_parents[j] != NULL && fino.gradient_highorder_nodes == gradient_average) {
          // promedio de padres
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
          
          // evaluacion directa en los nodos
          gsl_matrix *dhdx = gsl_matrix_calloc(element->type->nodes, fino.dimensions); // esto esta al vesre
          mesh_compute_dhdx(element, element->type->node_coords[j], NULL, dhdx);
          
          // las nueve derivadas (o menos)
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              for (j_local_prime = 0; j_local_prime < element->type->nodes; j_local_prime++) {
                j_global_prime = element->node[j_local_prime]->index_mesh;
                gsl_matrix_add_to_element(element->dphidx_node[j], g, m, gsl_matrix_get(dhdx, j_local_prime, m) * mesh->node[j_global_prime].phi[g]);
              }
            }
          }
          
          gsl_matrix_free(dhdx);
          
        }
        

        
        // si nu, E y/o alpha no son uniformes, los tenemos que evaluar en los nodos
        if (distribution_E.function != NULL || distribution_E.physical_property != NULL ||
            distribution_nu.function != NULL || distribution_nu.physical_property != NULL ||
            distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL) {
          
          element->property_node[j] = calloc(ELASTIC_PROPERTIES, sizeof(double));
          wasora_var_value(wasora_mesh.vars.x) = mesh->node[j_global].x[0];
          wasora_var_value(wasora_mesh.vars.y) = mesh->node[j_global].x[1];
          wasora_var_value(wasora_mesh.vars.z) = mesh->node[j_global].x[2];
          
          if (distribution_E.variable == NULL) {
            if ((E = fino_distribution_evaluate(&distribution_E, element->physical_entity->material, element->node[j]->x)) <= 0) {
              wasora_push_error_message("E is not positive (%g)", E);
              return WASORA_RUNTIME_ERROR;
            }      
          }
          if (distribution_nu.variable == NULL ) {
            nu = fino_distribution_evaluate(&distribution_nu, element->physical_entity->material, element->node[j]->x);
            if (nu > 0.499) {
              wasora_push_error_message("nu is greater than 1/2");
              return WASORA_RUNTIME_ERROR;
            } else if (nu < 0) {
              wasora_push_error_message("nu is negative");
              return WASORA_RUNTIME_ERROR;
            }
          }
          if (distribution_alpha.variable == NULL) {
            alpha = fino_distribution_evaluate(&distribution_alpha, element->physical_entity->material, element->node[j]->x);
          }

          element->property_node[j][LAMBDA] = E*nu/((1+nu)*(1-2*nu));
          element->property_node[j][MU] = 0.5*E/(1+nu);
          element->property_node[j][ALPHA] = alpha;
          
        }
        
        if (fino.rough) {
          // si estamos en rough ya calculamos los valores nodales y ya
          element->node[j]->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
          element->node[j]->dphidx = element->dphidx_node[j];

          element->node[j]->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));
          fino_break_compute_nodal_stresses(element, j, &sigmax, &sigmay, &sigmaz, &tauxy, &tauyz, &tauzx);
          element->node[j]->f[SIGMAX] = sigmax;
          element->node[j]->f[SIGMAY] = sigmay;
          element->node[j]->f[SIGMAZ] = sigmaz;
          element->node[j]->f[TAUXY] = tauxy;
          element->node[j]->f[TAUYZ] = tauyz;
          element->node[j]->f[TAUZX] = tauzx;
        }
      }
    }
  }


  // paso 2. barremos nodos de la malla de salida (la misma en smooth, rough en rough)
  for (j_global = 0; j_global < mesh->n_nodes; j_global++) {
    if (fino.progress_ascii && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }

    node = &mesh->node[j_global];
    
    if (fino.rough == 0) {
      node->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      node->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));

      weight_total = 0;
      n = 0;
      LL_FOREACH(mesh->node[j_global].associated_elements, associated_element) {
        if (associated_element->element->dphidx_node != NULL) {
          weight_total += associated_element->element->weight;
          n++;
        }
      }

      LL_FOREACH(mesh->node[j_global].associated_elements, associated_element) {
        element = associated_element->element;
        if (element->dphidx_node != NULL) {
          if (weight_total != 0) {
            weight_normalized = element->weight / weight_total;
          } else {
            weight_normalized = 1.0/(double)n;
          }  
          for (j = 0; j < element->type->nodes; j++) {
            if (element->node[j]->index_mesh == j_global) {

              // las derivadas
              // TODO: blas level2
              for (g = 0; g < fino.degrees; g++) {
                for (m = 0; m < fino.dimensions; m++) {
                  gsl_matrix_add_to_element(node->dphidx, g, m, weight_normalized * gsl_matrix_get(element->dphidx_node[j], g, m));
                }
              }

              fino_break_compute_nodal_stresses(element, j, &sigmax, &sigmay, &sigmaz, &tauxy, &tauyz, &tauzx);

              node->f[SIGMAX] += weight_normalized * sigmax;
              node->f[SIGMAY] += weight_normalized * sigmay;
              node->f[SIGMAZ] += weight_normalized * sigmaz;
              node->f[TAUXY] += weight_normalized * tauxy;
              node->f[TAUYZ] += weight_normalized * tauyz;
              node->f[TAUZX] += weight_normalized * tauzx;

            }
          }
        }
      }
    }
    
    // ya tenemos los promedios ahora, rellenamos las funciones    
    fino.gradient[0][0]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 0);
    fino.gradient[0][1]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 1);
    
    fino.gradient[1][0]->data_value[j_global] = gsl_matrix_get(node->dphidx, 1, 0);
    fino.gradient[1][1]->data_value[j_global] = gsl_matrix_get(node->dphidx, 1, 1);

    if (fino.dimensions > 2) {
      fino.gradient[0][2]->data_value[j_global] = gsl_matrix_get(node->dphidx, 0, 2);
      fino.gradient[1][2]->data_value[j_global] = gsl_matrix_get(node->dphidx, 1, 2);

      fino.gradient[2][0]->data_value[j_global] = gsl_matrix_get(node->dphidx, 2, 0);
      fino.gradient[2][1]->data_value[j_global] = gsl_matrix_get(node->dphidx, 2, 1);
      fino.gradient[2][2]->data_value[j_global] = gsl_matrix_get(node->dphidx, 2, 2);
    }
    
    
    fino.sigmax->data_value[j_global] = node->f[SIGMAX];
    fino.sigmay->data_value[j_global] = node->f[SIGMAY];
    fino.sigmaz->data_value[j_global] = node->f[SIGMAZ];

    fino.tauxy->data_value[j_global] = node->f[TAUXY];

    if (fino.dimensions == 3) {
      fino.tauyz->data_value[j_global] = node->f[TAUYZ];
      fino.tauzx->data_value[j_global] = node->f[TAUZX];
    }

    wasora_call(fino_compute_principal_stress(node->f[SIGMAX], node->f[SIGMAY], node->f[SIGMAZ],
                                              node->f[TAUXY], node->f[TAUYZ], node->f[TAUZX],
                                              &sigma1, &sigma2, &sigma3));
    
    fino.sigma1->data_value[j_global] = sigma1;
    fino.sigma2->data_value[j_global] = sigma2;
    fino.sigma3->data_value[j_global] = sigma3;

    // tresca
    tresca = fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);
    fino.tresca->data_value[j_global] = tresca;

    // von mises
    sigma = fino_compute_vonmises_from_principal(sigma1, sigma2, sigma3);
    
    if ((fino.sigma->data_value[j_global] = sigma) > wasora_var(fino.vars.sigma_max)) {
      wasora_var(fino.vars.sigma_max) = fino.sigma->data_value[j_global];
      
      wasora_var(fino.vars.sigma_max_x) = mesh->node[j_global].x[0];
      wasora_var(fino.vars.sigma_max_y) = mesh->node[j_global].x[1];
      wasora_var(fino.vars.sigma_max_z) = mesh->node[j_global].x[2];
      
      wasora_var(fino.vars.u_at_sigma_max) = mesh->node[j_global].phi[0];
      wasora_var(fino.vars.v_at_sigma_max) = mesh->node[j_global].phi[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_sigma_max) = mesh->node[j_global].phi[2];
      }
    }
    
    displ2 = 0;
    for (g = 0; g < fino.degrees; g++) {
      displ2 += gsl_pow_2(mesh->node[j_global].phi[g]);
    }
    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice displ_max
    if (displ2 >= max_displ2) {
      max_displ2 = displ2;
      wasora_var(fino.vars.displ_max) = sqrt(displ2);
      wasora_var(fino.vars.displ_max_x) = mesh->node[j_global].x[0];
      wasora_var(fino.vars.displ_max_y) = mesh->node[j_global].x[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.displ_max_z) = mesh->node[j_global].x[2];
      }
      
      wasora_var(fino.vars.u_at_displ_max) = mesh->node[j_global].phi[0];
      wasora_var(fino.vars.v_at_displ_max) = mesh->node[j_global].phi[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_displ_max) = mesh->node[j_global].phi[2];
      }
    }
  }

  if (fino.progress_ascii) {
    while (ascii_progress_chars++ < 100) {
      printf(CHAR_PROGRESS_GRADIENT);
    }
    printf("\n");  
    fflush(stdout);
  }

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_neumann"
int fino_break_set_neumann(element_t *element, bc_t *bc) {
  
  double x0, y0, z0;
  double Ix, Iy, Iz;
  double r_for_axisymmetric;
  int i, j, v;
  

  if (bc->type_phys ==  bc_phys_force) {
    if (element->physical_entity->volume == 0) {
      wasora_push_error_message("physical entity '%s' has zero volume", element->physical_entity->name);
      return WASORA_RUNTIME_ERROR;
    }
  } else if (bc->type_phys == bc_phys_pressure_normal || bc->type_phys == bc_phys_pressure_real) {
    if ((fino.dimensions-element->type->dim != 1)) {
      wasora_push_error_message("pressure BCs can only be applied to surfaces");
      return WASORA_RUNTIME_ERROR;
    }
  } else if (bc->type_phys == bc_phys_moment)   {
    
    double xix, xiy, xiz;
    physical_entity_t *physical_entity = element->physical_entity;
    element_t *tmp_element;
  
    if ((fino.dimensions-element->type->dim != 1)) {
      wasora_push_error_message("moment BCs can only be applied to surfaces");
      return WASORA_RUNTIME_ERROR;
    }
    
    // el centro
    if (bc->expr[4].n_tokens != 0) {
      x0 = wasora_evaluate_expression(&bc->expr[3]);
    } else {
      x0 = physical_entity->cog[0];
    }
    if (bc->expr[5].n_tokens != 0) {
      y0 = wasora_evaluate_expression(&bc->expr[4]);
    } else {
      y0 = physical_entity->cog[1];
    }
    if (bc->expr[6].n_tokens != 0) {
      z0 = wasora_evaluate_expression(&bc->expr[5]);
    } else {
      z0 = physical_entity->cog[2];
    }

    // el momento polar de inercia
    Ix = 0;
    Iy = 0;
    Iz = 0;
    for (i = 0; i < physical_entity->n_elements; i++) {
      tmp_element = &fino.mesh->element[physical_entity->element[i]];
      for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
        mesh_compute_integration_weight_at_gauss(element, v);

        xix = 0;
        xiy = 0;
        xiz = 0;
        for (j = 0; j < element->type->nodes; j++) {
          xix += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[0]);
          xiy += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[1]);
          xiz += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[2]);
        }

        Ix += element->w[v] * xix;
        Iy += element->w[v] * xiy;
        Iz += element->w[v] * xiz;
      }
    }
  }  
  
  
  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    mesh_compute_integration_weight_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    mesh_compute_x_at_gauss(element, v);
    mesh_update_coord_vars(element->x[v]);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
    
    if (bc->type_phys == bc_phys_stress) {
      
      gsl_vector_set(fino.Nb, bc->dof, wasora_evaluate_expression(&bc->expr[0]));
      
    } else if (bc->type_phys ==  bc_phys_force) {
      
      gsl_vector_set(fino.Nb, bc->dof, wasora_evaluate_expression(&bc->expr[0])/element->physical_entity->volume);
      
    } else if (bc->type_phys == bc_phys_pressure_normal || bc->type_phys == bc_phys_pressure_real) {
      
      double p;
      
      // la p chica es la proyeccion del vector tension sobre la normal, lo que uno espera en matematica
      p = wasora_evaluate_expression(&bc->expr[0]);
      // la P grande es presion positiva cuando comprime, como lo que uno espera en ingenieria      
      if (bc->type_phys == bc_phys_pressure_real) {
        p = -p;
      }
      
      gsl_vector_set(fino.Nb, 0, wasora_var_value(wasora_mesh.vars.nx) * p);
      gsl_vector_set(fino.Nb, 1, wasora_var_value(wasora_mesh.vars.ny) * p);
      if (fino.dimensions == 3) {
        gsl_vector_set(fino.Nb, 2, wasora_var_value(wasora_mesh.vars.nz) * p);
      }
      
    } else if (bc->type_phys == bc_phys_moment)   {
      
      double dx, dy, dz;
      double d, theta, F;
      
      dx = element->x[v][0] - x0;
      dy = element->x[v][1] - y0;
      dz = element->x[v][2] - z0;
    
      // los tres primeros tienen las componentes Mx My y Mz
      if (bc->expr[0].n_tokens != 0) {
        d = gsl_hypot(dy, dz);
        theta = atan2(dy, dz);
        F = wasora_evaluate_expression(&bc->expr[0]) * 0.5 * d / (Iy+Iz);
        // ty = cos(theta)*dz
        gsl_vector_add_to_element(fino.Nb, 1, -F*cos(theta));
        // tz = -sin(theta)*dy
        gsl_vector_add_to_element(fino.Nb, 2, +F*sin(theta));
      }
    
      if (bc->expr[1].n_tokens != 0) {
        d = gsl_hypot(dx, dz);
        theta = atan2(dx, dz);
        F = wasora_evaluate_expression(&bc->expr[1]) * 0.5 * d / (Ix+Iz);
        gsl_vector_add_to_element(fino.Nb, 0, -F*cos(theta));
        gsl_vector_add_to_element(fino.Nb, 2, +F*sin(theta));
      }
    
      if (bc->expr[2].n_tokens != 0) {
        d = gsl_hypot(dx, dy);
        theta = atan2(dx, dy);
        F = wasora_evaluate_expression(&bc->expr[2]) * 0.5 * d / (Ix+Iy);
        gsl_vector_add_to_element(fino.Nb, 0, -F*cos(theta));
        gsl_vector_add_to_element(fino.Nb, 1, +F*sin(theta));
      }
    
    }  

    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*element->w[v], element->H[v], fino.Nb, 1.0, fino.bi); 
    
    
  }

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_principal_stress"
int fino_compute_principal_stress(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx, double *sigma1, double *sigma2, double *sigma3) {
  
  double I1, I2, I3;
  double c1, c2, c3;
  double phi;
  
  // stress invariants
  // https://en.wikiversity.org/wiki/Principal_stresses
  I1 = sigmax + sigmay + sigmaz;
  I2 = sigmax*sigmay + sigmay*sigmaz + sigmaz*sigmax - gsl_pow_2(tauxy) - gsl_pow_2(tauyz) - gsl_pow_2(tauzx);
  I3 = sigmax*sigmay*sigmaz - sigmax*gsl_pow_2(tauyz) - sigmay*gsl_pow_2(tauzx) - sigmaz*gsl_pow_2(tauxy) + 2*tauxy*tauyz*tauzx;

  // principal stresses
  c1 = sqrt(fabs(gsl_pow_2(I1) - 3*I2));
  phi = 1.0/3.0 * acos((2.0*gsl_pow_3(I1) - 9.0*I1*I2 + 27.0*I3)/(2.0*gsl_pow_3(c1)));
  if (isnan(phi)) {
    phi = 0;
  }

  c2 = I1/3.0;
  c3 = 2.0/3.0 * c1;
  if (sigma1 != NULL) {
    *sigma1 = c2 + c3 * cos(phi);
  }
  if (sigma2 != NULL) {
    *sigma2 = c2 + c3 * cos(phi - 2.0*M_PI/3.0);
  }
  if (sigma3 != NULL) {
    *sigma3 = c2 + c3 * cos(phi - 4.0*M_PI/3.0);
  }
  
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_vonmises_from_principal"
double fino_compute_vonmises_from_principal(double sigma1, double sigma2, double sigma3) {
  
  return sqrt(0.5*(gsl_pow_2(sigma1-sigma2) + gsl_pow_2(sigma2-sigma3) + gsl_pow_2(sigma3-sigma1)));
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_vonmises_from_tensor"
double fino_compute_vonmises_from_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {
  
  return sqrt(0.5*(gsl_pow_2(sigmax-sigmay) + gsl_pow_2(sigmay-sigmaz) + gsl_pow_2(sigmaz-sigmax) +
                       6.0 * (gsl_pow_2(tauxy) + gsl_pow_2(tauyz) + gsl_pow_2(tauzx))));
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_tresca_from_principal"
double fino_compute_tresca_from_principal(double sigma1, double sigma2, double sigma3) {

  double S12 = fabs(sigma1-sigma2);
  double S23 = fabs(sigma2-sigma3);
  double S31 = fabs(sigma3-sigma1);

  if (S31 >= S12 && S31 >= S23) {
    return S31;
  } else if (S12 >= S23 && S12 >= S31) {
    return S12;
  } else if (S23 >= S12 && S23 >= S31) {
    return S23;
  }

  return 0;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_tresca_from_tensor"
double fino_compute_tresca_from_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {

  double sigma1, sigma2, sigma3;

  wasora_call(fino_compute_principal_stress(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx, &sigma1, &sigma2, &sigma3));
  return fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);

}


#undef  __FUNCT__
#define __FUNCT__ "fino_compute_strain_energy"
int fino_compute_strain_energy(void) {

  PetscScalar e;
  Vec Kphi;
  
  petsc_call(VecDuplicate(fino.phi, &Kphi));
  petsc_call(MatMult(fino.K_nobc, fino.phi, Kphi));
  petsc_call(VecDot(fino.phi, Kphi, &e));
  wasora_var(fino.vars.strain_energy) = 0.5*e;
  
  if (fino.problem_kind == problem_kind_axisymmetric) {
    wasora_var(fino.vars.strain_energy) *= 2*M_PI;
  }
  petsc_call(VecDestroy(&Kphi));
  
  return WASORA_RUNTIME_OK;
}
