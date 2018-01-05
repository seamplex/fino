/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's method for stress linearization over SCLs according to ASME
 *
 *  Copyright (C) 2017--2018 jeremy theler
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

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>

#include "fino.h"

struct linearize_params_t {
  double x1, y1, z1;
  double x2, y2, z2;
  double t;
  
  function_t *function;
};

#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_linearize"
int fino_instruction_linearize(void *arg) {
  
  fino_linearize_t *linearize = (fino_linearize_t *)arg;
  
  element_t *element;
  
  double integrand_mx, integrand_my, integrand_mz, integrand_mxy, integrand_myz, integrand_mzx;
  double integrand_bx, integrand_by, integrand_bz, integrand_bxy, integrand_byz, integrand_bzx;
  double sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m;
  double sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b;
  double sigma1, sigma2, sigma3;
  double peak;
  double x,  y,  z;
  double x1, y1, z1;      // coordenadas del punto inicial de la SCL
  double w_gauss;
  double r_for_axisymmetric;
  double t, t_prime;
  double h, den;
  int i, j, k, v;
  
// http://www.eng-tips.com/faqs.cfm?fid=982
  // TODO: sacar un markdown con los detalles
  // TODO: opcion "ignore through-thickness bending stress"

  if (linearize->physical_entity != NULL) {
    if (linearize->physical_entity->dimension != 1) {
      wasora_push_error_message("physical entity for the SCL has to be one-dimensional, not %d-dimensional", linearize->physical_entity->dimension);
      return WASORA_PARSER_ERROR;
    }  

    if (linearize->physical_entity->n_elements == 0 || linearize->physical_entity->element == NULL) {
      wasora_push_error_message("physical entity '%s' has no elements to act as a SCL", linearize->physical_entity->name);
    }

    // el primer nodo del primer elemento de la SCL
    x1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[0];
    y1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[1];
    z1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[2];

    // la longitud de la SCL
    t = linearize->physical_entity->volume;

    sigmax_m = sigmay_m = sigmaz_m = tauxy_m = tauyz_m = tauzx_m = 0;
    sigmax_b = sigmay_b = sigmaz_b = tauxy_b = tauyz_b = tauzx_b = 0;
    peak = 0;
    for (i = 0; i < linearize->physical_entity->n_elements; i++) {
      element = &fino.mesh->element[linearize->physical_entity->element[i]];
      for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
        w_gauss = mesh_integration_weight(fino.mesh, element, v);
        r_for_axisymmetric = fino_compute_r_for_axisymmetric();
        mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);

        x = gsl_vector_get(fino.mesh->fem.x, 0);
        y = gsl_vector_get(fino.mesh->fem.x, 1);
        z = gsl_vector_get(fino.mesh->fem.x, 2);
        t_prime = sqrt(gsl_pow_2(x-x1) + gsl_pow_2(y-y1) + gsl_pow_2(z-z1));

        integrand_mx = integrand_my = integrand_mz = integrand_mxy = integrand_myz = integrand_mzx = 0;
        integrand_bx = integrand_by = integrand_bz = integrand_bxy = integrand_byz = integrand_bzx = 0;
        for (j = 0; j < element->type->nodes; j++) {
          h = gsl_vector_get(fino.mesh->fem.h, j);
          k = element->node[j]->id - 1;
          integrand_mx   += h * fino.sigmax->data_value[k];
          integrand_my   += h * fino.sigmay->data_value[k];
          integrand_mz   += h * fino.sigmaz->data_value[k];
          integrand_mxy  += h * fino.tauxy->data_value[k];
          integrand_myz  += h * fino.tauyz->data_value[k];
          integrand_mzx  += h * fino.tauzx->data_value[k];

          integrand_bx   += h * fino.sigmax->data_value[k] * (t/2 - t_prime);
          integrand_by   += h * fino.sigmay->data_value[k] * (t/2 - t_prime);
          integrand_bz   += h * fino.sigmaz->data_value[k] * (t/2 - t_prime);
          integrand_bxy  += h * fino.tauxy->data_value[k] * (t/2 - t_prime);
          integrand_byz  += h * fino.tauyz->data_value[k] * (t/2 - t_prime);
          integrand_bzx  += h * fino.tauzx->data_value[k] * (t/2 - t_prime);

          if (fino.sigma->data_value[k] > peak) {
            peak = fino.sigma->data_value[k];
          }
        }

        w_gauss *= r_for_axisymmetric;
        sigmax_m  += w_gauss * integrand_mx;
        sigmay_m  += w_gauss * integrand_my;
        sigmaz_m  += w_gauss * integrand_mz;
        tauxy_m   += w_gauss * integrand_mxy;
        tauyz_m   += w_gauss * integrand_myz;
        tauzx_m   += w_gauss * integrand_mzx;

        sigmax_b  += w_gauss * integrand_bx;
        sigmay_b  += w_gauss * integrand_by;
        sigmaz_b  += w_gauss * integrand_bz;
        tauxy_b   += w_gauss * integrand_bxy;
        tauyz_b   += w_gauss * integrand_byz;
        tauzx_b   += w_gauss * integrand_bzx;
      }
    }

  } else {
    
    // SCL definida por puntos
    double error;  
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DEFAULT_INTEGRATION_INTERVALS);
    gsl_function F;
    struct linearize_params_t params;
    
    params.x1 = wasora_evaluate_expression(&linearize->x1);
    params.y1 = wasora_evaluate_expression(&linearize->y1);
    params.z1 = wasora_evaluate_expression(&linearize->z1);

    params.x2 = wasora_evaluate_expression(&linearize->x2);
    params.y2 = wasora_evaluate_expression(&linearize->y2);
    params.z2 = wasora_evaluate_expression(&linearize->z2);
    
    params.t = gsl_hypot3(params.x2-params.x1, params.y2-params.y1, params.z2-params.z1);
    
    F.params = &params;
    
    // membrane
    F.function = &fino_linearization_integrand_membrane;
   
    params.function = fino.sigmax;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmax_m, &error);

    params.function = fino.sigmay;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmay_m, &error);

    params.function = fino.sigmaz;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmaz_m, &error);
    
    params.function = fino.tauxy;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauxy_m, &error);

    params.function = fino.tauyz;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauyz_m, &error);
    
    params.function = fino.tauzx;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauzx_m, &error);

    // membrane
    F.function = &fino_linearization_integrand_bending;
   
    params.function = fino.sigmax;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmax_b, &error);

    params.function = fino.sigmay;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmay_b, &error);

    params.function = fino.sigmaz;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmaz_b, &error);
    
    params.function = fino.tauxy;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauxy_b, &error);

    params.function = fino.tauyz;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauyz_b, &error);
    
    params.function = fino.tauzx;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauzx_b, &error);
        
  }
  

  sigmax_m /= t;
  sigmay_m /= t;
  sigmaz_m /= t;
  tauxy_m  /= t;
  tauyz_m  /= t;
  tauzx_m  /= t;
  
//  printf("membrane\n");
//  printf("%e\t%e\t%e\t%e\t%e\t%e\n", sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);

  den = gsl_pow_2(t)/6.0;
  sigmax_b /= den;
  sigmay_b /= den;
  sigmaz_b /= den;
  tauxy_b  /= den;
  tauyz_b  /= den;
  tauzx_b  /= den;

//  printf("bending\n");
//  printf("%e\t%e\t%e\t%e\t%e\t%e\n", sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b);

    
//  printf("LINEARIZATION\n");
  fino_compute_principal_stress(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_m, tauzx_m+tauzx_m, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m+tauyz_m, tauzx_m-tauzx_m, &sigma1, &sigma2, &sigma3);
//  printf("-------------\n\n");

  if (linearize->membrane != NULL) {
    wasora_var_value(linearize->membrane) = fino_compute_vonmises_from_tensor(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
  }
  if (linearize->bending != NULL) {
    wasora_var_value(linearize->bending) = fino_compute_vonmises_from_tensor(sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b);
  }
  if (linearize->peak != NULL) {
    wasora_var_value(linearize->peak) = peak;
  }

/*  
  printf("membrane %g\n", wasora_var_value(linearize->membrane));
  printf("bending  %g\n", wasora_var_value(linearize->bending));
  printf("m+b      %g\n", fino_compute_vonmises_from_tensor(integrand_mx+integrand_bx, integrand_my+integrand_by, integrand_mz+integrand_bz, integrand_mxy+integrand_bxy, integrand_myz+integrand_byz, integrand_mzx+integrand_bzx));
  printf("m-b      %g\n", fino_compute_vonmises_from_tensor(integrand_mx-integrand_bx, integrand_my-integrand_by, integrand_mz-integrand_bz, integrand_mxy-integrand_bxy, integrand_myz-integrand_byz, integrand_mzx-integrand_bzx));
*/
  return WASORA_RUNTIME_OK;
}


double fino_linearization_integrand_membrane(double t_prime, void *params) {
  double x[3];
  double f;
  struct linearize_params_t *p = (struct linearize_params_t *)params;
  
  x[0] = p->x1 + t_prime/p->t * (p->x2 - p->x1);
  x[1] = p->y1 + t_prime/p->t * (p->y2 - p->y1);
  x[2] = p->z1 + t_prime/p->t * (p->z2 - p->z1);

  f = wasora_evaluate_function(p->function, x);
  
  return f;
}

double fino_linearization_integrand_bending(double t_prime, void *params) {
  double x[3];
  double f;
  struct linearize_params_t *p = (struct linearize_params_t *)params;
  
  x[0] = p->x1 + t_prime/p->t * (p->x2 - p->x1);
  x[1] = p->y1 + t_prime/p->t * (p->y2 - p->y1);
  x[2] = p->z1 + t_prime/p->t * (p->z2 - p->z1);

  f = wasora_evaluate_function(p->function, x) * (0.5*p->t - t_prime);
  
  return f;
}
