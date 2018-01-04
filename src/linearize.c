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
#include "fino.h"

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
  double xi, yi, zi;      // coordenadas del punto inicial de la SCL
//  double x0, y0, z0;      // coordenadas del COG de la SCL
  double w_gauss;
  double r_for_axisymmetric;
  double t, t_prime;
  double h, den;
  int i, j, k, v;
  
// http://www.eng-tips.com/faqs.cfm?fid=982
  // TODO: sacar un markdown con los detalles

  if (linearize->scl->dimension != 1) {
    wasora_push_error_message("physical entity for the SCL has to be one-dimensional, not %d-dimensional", linearize->scl->dimension);
    return WASORA_PARSER_ERROR;
  }  
  
  if (linearize->scl->n_elements == 0 || linearize->scl->element == NULL) {
    wasora_push_error_message("physical entity '%s' has no elements to act as a SCL", linearize->scl->name);
  }
  
  // el primer nodo del primer elemento de la SCL
  xi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[0];
  yi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[1];
  zi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[2];
  
  // el centro de gravedad de la SCL
/*  
  x0 = linearize->scl->cog[0];
  y0 = linearize->scl->cog[1];
  z0 = linearize->scl->cog[2];
*/  
  // la longitud de la SCL
  t = linearize->scl->volume;
  
  sigmax_m = sigmay_m = sigmaz_m = tauxy_m = tauyz_m = tauzx_m = 0;
  sigmax_b = sigmay_b = sigmaz_b = tauxy_b = tauyz_b = tauzx_b = 0;
  peak = 0;
  for (i = 0; i < linearize->scl->n_elements; i++) {
    element = &fino.mesh->element[linearize->scl->element[i]];
    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_integration_weight(fino.mesh, element, v);
      r_for_axisymmetric = fino_compute_r_for_axisymmetric();
      mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);

      x = gsl_vector_get(fino.mesh->fem.x, 0);
      y = gsl_vector_get(fino.mesh->fem.x, 1);
      z = gsl_vector_get(fino.mesh->fem.x, 2);
      t_prime = sqrt(gsl_pow_2(x-xi) + gsl_pow_2(y-yi) + gsl_pow_2(z-zi));

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

  sigmax_m /= t;
  sigmay_m /= t;
  sigmaz_m /= t;
  tauxy_m  /= t;
  tauyz_m  /= t;
  tauzx_m  /= t;
  printf("membrane\n");
  printf("%e\t%e\t%e\t%e\t%e\t%e\n", sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
  
  den = gsl_pow_2(t)/6.0;
  sigmax_b /= den;
  sigmay_b /= den;
  sigmaz_b /= den;
  tauxy_b  /= den;
  tauyz_b  /= den;
  tauzx_b  /= den;

  printf("bending\n");
  printf("%e\t%e\t%e\t%e\t%e\t%e\n", sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b);
  
  printf("LINEARIZATION\n");
  fino_compute_principal_stress(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_m+sigmax_b, sigmay_m+sigmay_b, sigmaz_m+sigmaz_b, tauxy_m+tauxy_b, tauyz_m+tauyz_m, tauzx_m+tauzx_m, &sigma1, &sigma2, &sigma3);
  fino_compute_principal_stress(sigmax_m-sigmax_b, sigmay_m-sigmay_b, sigmaz_m-sigmaz_b, tauxy_m-tauxy_b, tauyz_m+tauyz_m, tauzx_m-tauzx_m, &sigma1, &sigma2, &sigma3);
  printf("-------------\n\n");

  if (linearize->membrane != NULL) {
    wasora_var_value(linearize->membrane) = fino_compute_vonmises_from_tensor(sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
  }
  if (linearize->bending != NULL) {
    wasora_var_value(linearize->bending) = fino_compute_vonmises_from_tensor(sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b);
  }
  if (linearize->peak != NULL) {
    wasora_var_value(linearize->peak) = peak;
  }

  printf("membrane %g\n", wasora_var_value(linearize->membrane));
  printf("bending  %g\n", wasora_var_value(linearize->bending));
  printf("m+b      %g\n", fino_compute_vonmises_from_tensor(integrand_mx+integrand_bx, integrand_my+integrand_by, integrand_mz+integrand_bz, integrand_mxy+integrand_bxy, integrand_myz+integrand_byz, integrand_mzx+integrand_bzx));
  printf("m-b      %g\n", fino_compute_vonmises_from_tensor(integrand_mx-integrand_bx, integrand_my-integrand_by, integrand_mz-integrand_bz, integrand_mxy-integrand_bxy, integrand_myz-integrand_byz, integrand_mzx-integrand_bzx));

  return WASORA_RUNTIME_OK;
}
