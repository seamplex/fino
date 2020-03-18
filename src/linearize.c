/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's method for stress linearization over SCLs according to ASME
 *
 *  Copyright (C) 2017--2018 jeremy theler
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
  double M, B, MB, MBplus, MBminus, P, T1, T2, T;
  double x1[3], x2[3];
  double t_prime, t_over_two_minus_t_prime;
  double h, den;
  int i, j, k, v;

  struct linearize_params_t params;
  
  char *total_name = NULL;
  function_t *total_function = NULL;

  // http://www.eng-tips.com/faqs.cfm?fid=982

  if (linearize->physical_entity != NULL) {
    if (linearize->physical_entity->dimension != 1) {
      wasora_push_error_message("physical entity for the SCL has to be one-dimensional, not %d-dimensional", linearize->physical_entity->dimension);
      return WASORA_PARSER_ERROR;
    }  

    if (linearize->physical_entity->n_elements == 0 || linearize->physical_entity->element == NULL) {
      wasora_push_error_message("physical entity '%s' has no elements to act as a SCL", linearize->physical_entity->name);
      return WASORA_PARSER_ERROR;
    }

    // el primer nodo del primer elemento de la SCL
    params.x1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[0];
    params.y1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[1];
    params.z1 = fino.mesh->element[linearize->physical_entity->element[0]].node[0]->x[2];
    
    // el final es el doble con respecto al cog
    params.x2 = params.x1 + 2*(linearize->physical_entity->cog[0]-params.x1);
    params.y2 = params.y1 + 2*(linearize->physical_entity->cog[1]-params.y1);
    params.z2 = params.z1 + 2*(linearize->physical_entity->cog[2]-params.z1);

    // la longitud de la SCL
    params.t = linearize->physical_entity->volume;

    sigmax_m = sigmay_m = sigmaz_m = tauxy_m = tauyz_m = tauzx_m = 0;
    sigmax_b = sigmay_b = sigmaz_b = tauxy_b = tauyz_b = tauzx_b = 0;
    for (i = 0; i < linearize->physical_entity->n_elements; i++) {
      element = &fino.mesh->element[linearize->physical_entity->element[i]];
      for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
        mesh_compute_integration_weight_at_gauss(element, v);
        mesh_compute_x_at_gauss(element,v );

        t_prime = gsl_hypot3(element->x[v][0]-params.x1, element->x[v][1]-params.y1, element->x[v][2]-params.z1);
        t_over_two_minus_t_prime = params.t/2 - t_prime;

        integrand_mx = integrand_my = integrand_mz = integrand_mxy = integrand_myz = integrand_mzx = 0;
        integrand_bx = integrand_by = integrand_bz = integrand_bxy = integrand_byz = integrand_bzx = 0;
        for (j = 0; j < element->type->nodes; j++) {
          h = element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j];
          k = element->node[j]->index_mesh;
          integrand_mx   += h * fino.sigmax->data_value[k];
          integrand_my   += h * fino.sigmay->data_value[k];
          integrand_mz   += h * fino.sigmaz->data_value[k];
          integrand_mxy  += h * fino.tauxy->data_value[k];

          integrand_bx   += h * fino.sigmax->data_value[k] * t_over_two_minus_t_prime;
          integrand_by   += h * fino.sigmay->data_value[k] * t_over_two_minus_t_prime;
          integrand_bz   += h * fino.sigmaz->data_value[k] * t_over_two_minus_t_prime;
          integrand_bxy  += h * fino.tauxy->data_value[k] * t_over_two_minus_t_prime;

          if (fino.dimensions > 2) {
            integrand_myz  += h * fino.tauyz->data_value[k];
            integrand_mzx  += h * fino.tauzx->data_value[k];
            
            integrand_byz  += h * fino.tauyz->data_value[k] * t_over_two_minus_t_prime;
            integrand_bzx  += h * fino.tauzx->data_value[k] * t_over_two_minus_t_prime;
          }
        }

//        w_gauss *= r_for_axisymmetric;
        sigmax_m  += element->w[v] * integrand_mx;
        sigmay_m  += element->w[v] * integrand_my;
        sigmaz_m  += element->w[v] * integrand_mz;
        tauxy_m   += element->w[v] * integrand_mxy;

        sigmax_b  += element->w[v] * integrand_bx;
        sigmay_b  += element->w[v] * integrand_by;
        sigmaz_b  += element->w[v] * integrand_bz;
        tauxy_b   += element->w[v] * integrand_bxy;

        if (fino.dimensions > 2) {
          tauyz_m   += element->w[v] * integrand_myz;
          tauzx_m   += element->w[v] * integrand_mzx;
          tauyz_b   += element->w[v] * integrand_byz;
          tauzx_b   += element->w[v] * integrand_bzx;
        }
      }
    }

  } else {
    
    // SCL definida por puntos
    double error;  
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DEFAULT_INTEGRATION_INTERVALS);
    gsl_function F;
    
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

    if (fino.dimensions > 2) {
      params.function = fino.tauyz;
      gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauyz_m, &error);
    
      params.function = fino.tauzx;
      gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauzx_m, &error);
    }

    // bending
    F.function = &fino_linearization_integrand_bending;
   
    params.function = fino.sigmax;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmax_b, &error);

    params.function = fino.sigmay;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmay_b, &error);

    params.function = fino.sigmaz;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &sigmaz_b, &error);
    
    params.function = fino.tauxy;
    gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauxy_b, &error);

    if (fino.dimensions > 2) {
      params.function = fino.tauyz;
      gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauyz_b, &error);
    
      params.function = fino.tauzx;
      gsl_integration_qags (&F, 0, params.t, 0, DEFAULT_INTEGRATION_INTERVALS, DEFAULT_INTEGRATION_INTERVALS, w, &tauzx_b, &error);
    }
    
    gsl_integration_workspace_free(w);
  }
  

  sigmax_m /= params.t;
  sigmay_m /= params.t;
  sigmaz_m /= params.t;
  tauxy_m  /= params.t;
  
  den = gsl_pow_2(params.t)/6.0;
  sigmax_b /= den;
  sigmay_b /= den;
  sigmaz_b /= den;
  tauxy_b  /= den;

  if (fino.dimensions > 2) {
    tauyz_m  /= params.t;
    tauzx_m  /= params.t;
    tauyz_b  /= den;
    tauzx_b  /= den;
  }

  
  x1[0] = params.x1;
  x1[1] = params.y1;
  x1[2] = params.z1;
  x2[0] = params.x2;
  x2[1] = params.y2;
  x2[2] = params.z2;
 
  switch (linearize->total) {
    case linearize_vonmises:
      total_function = fino.sigma;
      total_name = strdup("Von Mises");
      M = fino_compute_vonmises_from_stress_tensor(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
      MBplus = fino_compute_vonmises_from_stress_tensor(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_b, tauzx_m+tauzx_b);
      MBminus = fino_compute_vonmises_from_stress_tensor(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m-tauyz_b, tauzx_m-tauzx_b);
      
      // OJO! esto no es asi, hay que hacer un tensor que sea la diferencia
      T1 = wasora_evaluate_function(fino.sigma, x1);
      T2 = wasora_evaluate_function(fino.sigma, x2);
    break;
      
    case linearize_tresca:
      total_function = fino.tresca;
      total_name = strdup("Tresca");
      M = fino_compute_tresca_from_stress_tensor(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
      MBplus = fino_compute_tresca_from_stress_tensor(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_b, tauzx_m+tauzx_b);
      MBminus = fino_compute_tresca_from_stress_tensor(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m-tauyz_b, tauzx_m-tauzx_b);
      
      T1 = wasora_evaluate_function(fino.tresca, x1);
      T2 = wasora_evaluate_function(fino.tresca, x2);
    break;

    case linearize_principal1:
      total_function = fino.sigma1;
      total_name = strdup("Principal 1");
      fino_compute_principal_stress(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m, &M, NULL, NULL);
      fino_compute_principal_stress(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_b, tauzx_m+tauzx_b, &MBplus, NULL, NULL);
      fino_compute_principal_stress(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m-tauyz_b, tauzx_m-tauzx_b, &MBminus, NULL, NULL);
      
      T1 = wasora_evaluate_function(fino.sigma1, x1);
      T2 = wasora_evaluate_function(fino.sigma1, x2);
    break;

    case linearize_principal2:
      total_function = fino.sigma2;
      total_name = strdup("Principal 2");
      fino_compute_principal_stress(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m, NULL, &M, NULL);
      fino_compute_principal_stress(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_b, tauzx_m+tauzx_b, NULL, &MBplus, NULL);
      fino_compute_principal_stress(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m-tauyz_b, tauzx_m-tauzx_b, NULL, &MBminus, NULL);
      
      T1 = wasora_evaluate_function(fino.sigma1, x1);
      T2 = wasora_evaluate_function(fino.sigma2, x2);
    break;
  
    case linearize_principal3:
      total_function = fino.sigma3;
      total_name = strdup("Principal 3");
      fino_compute_principal_stress(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m, NULL, NULL, &M);
      fino_compute_principal_stress(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_b, tauzx_m+tauzx_b, NULL, NULL, &MBplus);
      fino_compute_principal_stress(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m-tauyz_b, tauzx_m-tauzx_b, NULL, NULL, &MBminus);
      
      T1 = wasora_evaluate_function(fino.sigma3, x1);
      T2 = wasora_evaluate_function(fino.sigma3, x2);
    break;
  }

  if (linearize->total != linearize_principal3) {
    MB = (MBplus>MBminus)?MBplus:MBminus;
    T = (T1>T2)?T1:T2;
    B = MB - M;
  } else {
    MB = (MBplus<MBminus)?MBplus:MBminus;
    T = (T1<T2)?T1:T2;
    B = M - MB;
  }
  P = T - MB;
  
  wasora_var_value(linearize->M) = M;
  wasora_var_value(linearize->MB) = MB;
  wasora_var_value(linearize->P) = P;
  
  if (linearize->file != NULL) {
    if (linearize->file->pointer == NULL) {
      wasora_call(wasora_instruction_open_file(linearize->file));
    }

    fprintf(linearize->file->pointer, "# # Stress linearization\n");
    fprintf(linearize->file->pointer, "#\n");

    fprintf(linearize->file->pointer, "# Start point: (%g, %g, %g)\n", params.x1, params.y1, params.z1);
    fprintf(linearize->file->pointer, "#   End point: (%g, %g, %g)\n", params.x2, params.y2, params.z2);
    fprintf(linearize->file->pointer, "#       Total: %s\n", total_name);
    
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "# ## Membrane stress tensor\n");
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "#  $\\sigma_{x}$ = %g\n", sigmax_m);
    fprintf(linearize->file->pointer, "#  $\\sigma_{y}$ = %g\n", sigmay_m);
    fprintf(linearize->file->pointer, "#  $\\sigma_{z}$ = %g\n", sigmaz_m);
    fprintf(linearize->file->pointer, "#   $\\tau_{xy}$ = %g\n", tauxy_m);
    fprintf(linearize->file->pointer, "#   $\\tau_{yz}$ = %g\n", tauyz_m);
    fprintf(linearize->file->pointer, "#   $\\tau_{zx}$ = %g\n", tauzx_m);
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "# ## Bending stress tensor\n");
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "#  $\\sigma_{x}$ = %g\n", sigmax_b);
    fprintf(linearize->file->pointer, "#  $\\sigma_{y}$ = %g\n", sigmay_b);
    fprintf(linearize->file->pointer, "#  $\\sigma_{z}$ = %g\n", sigmaz_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{xy}$ = %g\n", tauxy_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{yz}$ = %g\n", tauyz_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{zx}$ = %g\n", tauzx_b);
    fprintf(linearize->file->pointer, "#\n");
/*
    fprintf(linearize->file->pointer, "# ## Membrane plus bending stress tensor\n");
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "#  $\\sigma_{x}$ = %g\n", sigmax_m+sigmax_b);
    fprintf(linearize->file->pointer, "#  $\\sigma_{y}$ = %g\n", sigmay_m+sigmay_b);
    fprintf(linearize->file->pointer, "#  $\\sigma_{z}$ = %g\n", sigmaz_m+sigmaz_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{xy}$ = %g\n", tauxy_m+tauxy_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{yz}$ = %g\n", tauyz_m+tauyz_b);
    fprintf(linearize->file->pointer, "#   $\\tau_{zx}$ = %g\n", tauzx_m+tauzx_b);
    fprintf(linearize->file->pointer, "#\n");
*/
    fprintf(linearize->file->pointer, "# ## Linearization results\n");
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "#                Membrane stress $M$ = %g\n", wasora_var_value(linearize->M));
    fprintf(linearize->file->pointer, "#  Membrane plus bending stress $MB$ = %g\n", wasora_var_value(linearize->MB));
    fprintf(linearize->file->pointer, "#                    Peak stress $P$ = %g\n", wasora_var_value(linearize->P));
    fprintf(linearize->file->pointer, "#\n");
    fprintf(linearize->file->pointer, "# t_prime\tM\tMB\tT\n");

    for (i = 0; i <= 50; i++) {
      t_prime = i/50.0;
      x1[0] = params.x1 + t_prime * (params.x2 - params.x1);
      x1[1] = params.y1 + t_prime * (params.y2 - params.y1);
      x1[2] = params.z1 + t_prime * (params.z2 - params.z1);
      
      if (t_prime < 0.5) {
        fprintf(linearize->file->pointer, "%.2f\t%e\t%e\t%e\n", t_prime, M, M+((T1>M)?(+1):(-1))*B*(1-2*t_prime), wasora_evaluate_function(total_function, x1));
      } else {
        fprintf(linearize->file->pointer, "%.2f\t%e\t%e\t%e\n", t_prime, M, M+((T2>M)?(-1):(+1))*B*(1-2*t_prime), wasora_evaluate_function(total_function, x1));
      }
    }
  }    
   
  free(total_name);
  
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
