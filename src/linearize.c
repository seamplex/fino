/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's method for stress linearization over SCLs according to ASME
 *
 *  Copyright (C) 2017 jeremy theler
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
  
  double mx, my, mz, mxy, myz, mzx;
  double bx, by, bz, bxy, byz, bzx;
  double sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m;
  double sigmax_b, sigmay_b, sigmaz_b, tauxy_b, tauyz_b, tauzx_b;
  double m1, m2, m3;
  double b1, b2, b3;
  double x,  y,  z;
  double xi, yi, zi;      // coordenadas del punto inicial de la SCL
  double x0, y0, z0;      // coordenadas del COG de la SCL
  double w_gauss;
  double r_for_axisymmetric;
  double t, t_prime;
  double h, den;
  int i, j, k, v;
  
// http://www.eng-tips.com/faqs.cfm?fid=982

  if (linearize->scl->dimension != 1) {
    wasora_push_error_message("physical entity for the SCL has to be one-dimensional, not %d-dimensional", linearize->scl->dimension);
    return WASORA_PARSER_ERROR;
  }  
  
  if (linearize->scl->n_elements == 0 || linearize->scl->element == NULL) {
    wasora_push_error_message("physical entity '%s' as no elements to act as an SCL", linearize->scl->name);
  }
  
  // el primer nodo del primer elemento de la SCL
  xi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[0];
  yi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[1];
  zi = fino.mesh->element[linearize->scl->element[0]].node[0]->x[2];
  
  x0 = linearize->scl->cog[0];
  y0 = linearize->scl->cog[1];
  z0 = linearize->scl->cog[2];
  t = linearize->scl->volume;
  
  mx = my = mz = mxy = myz = mzx = 0;
  bx = by = bz = bxy = byz = bzx = 0;
  for (i = 0; i < linearize->scl->n_elements; i++) {
    element = &fino.mesh->element[linearize->scl->element[i]];
    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_integration_weight(fino.mesh, element, v);
      r_for_axisymmetric = fino_compute_r_for_axisymmetric();
      mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);

      x = gsl_vector_get(fino.mesh->fem.x, 0);
      y = gsl_vector_get(fino.mesh->fem.x, 1);
      z = gsl_vector_get(fino.mesh->fem.x, 2);
      t_prime = 2*sqrt((gsl_pow_2(x-xi) + gsl_pow_2(y-yi) + gsl_pow_2(z-zi))/(gsl_pow_2(x0-xi) + gsl_pow_2(y0-yi) + gsl_pow_2(z0-zi)));

      sigmax_m = sigmay_m = sigmaz_m = tauxy_m = tauyz_m = tauzx_m = 0;
      sigmax_b = sigmay_b = sigmaz_b = tauxy_b = tauyz_b = tauzx_b = 0;
      for (j = 0; j < element->type->nodes; j++) {
        h = gsl_vector_get(fino.mesh->fem.h, j);
        k = element->node[j]->id - 1;
//        printf("%d %d %d %d\n", i, v, j, k);
        sigmax_m += h * fino.sigmax->data_value[k];
        sigmay_m += h * fino.sigmay->data_value[k];
        sigmaz_m += h * fino.sigmaz->data_value[k];
        tauxy_m  += h * fino.tauxy->data_value[k];
        tauyz_m  += h * fino.tauyz->data_value[k];
        tauzx_m  += h * fino.tauzx->data_value[k];

        sigmax_b += h * fino.sigmax->data_value[k] * (t/2 - t_prime);
        sigmay_b += h * fino.sigmay->data_value[k] * (t/2 - t_prime);
        sigmaz_b += h * fino.sigmaz->data_value[k] * (t/2 - t_prime);
        tauxy_b  += h * fino.tauxy->data_value[k] * (t/2 - t_prime);
        tauyz_b  += h * fino.tauyz->data_value[k] * (t/2 - t_prime);
        tauzx_b  += h * fino.tauzx->data_value[k] * (t/2 - t_prime);
        
//        printf("%e\t%e\t%e\t%e\t%e\t%e\n", sigmax_m, sigmay_m, sigmaz_m, tauxy_m, tauyz_m, tauzx_m);
//        printf("%e\t%e\t%e\t%e\t%e\t%e\n", mx, my, mz, mxy, myz, mzx);
      }

      w_gauss *= r_for_axisymmetric;
      mx  += w_gauss * sigmax_m;
      my  += w_gauss * sigmay_m;
      mz  += w_gauss * sigmaz_m;
      mxy += w_gauss * tauxy_m;
      myz += w_gauss * tauyz_m;
      mzx += w_gauss * tauzx_m;

      bx  += w_gauss * sigmax_b;
      by  += w_gauss * sigmay_b;
      bz  += w_gauss * sigmaz_b;
      bxy += w_gauss * tauxy_b;
      byz += w_gauss * tauyz_b;
      bzx += w_gauss * tauzx_b;
    }
  }

//  printf("membrane\n");
//  printf("%e\t%e\t%e\t%e\t%e\t%e\n", mx, my, mz, mxy, myz, mzx);
  mx  /= t;
  my  /= t;
  mz  /= t;
  mxy /= t;
  myz /= t;
  mzx /= t;
  
  den = gsl_pow_2(t)/6.0;
  bx  /= den;
  by  /= den;
  bz  /= den;
  bxy /= den;
  byz /= den;
  bzx /= den;

//  printf("%e\t%e\t%e\t%e\t%e\t%e\n", mx, my, mz, mxy, myz, mzx);
  
  wasora_call(fino_compute_principal_stress(mx, my, mz, mxy, myz, mzx, &m1, &m2, &m3));
  wasora_call(fino_compute_principal_stress(bx, by, bz, bxy, byz, bzx, &b1, &b2, &b3));
  
  
  
  if (linearize->membrane != NULL) {
    wasora_var_value(linearize->membrane) = sqrt(0.5*(gsl_pow_2(m1-m2) + gsl_pow_2(m2-m3) + gsl_pow_2(m3-m1)));
  }
  if (linearize->bending != NULL) {
    wasora_var_value(linearize->bending) = sqrt(0.5*(gsl_pow_2(b1-b2) + gsl_pow_2(b2-b3) + gsl_pow_2(b3-b1)));
  }

  return WASORA_RUNTIME_OK;
}
