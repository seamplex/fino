/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of linear elastic problem (break) with optional vibration (shake)
 *
 *  Copyright (C) 2015--2016 jeremy theler
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

fino_distribution_t E;     // modulo de young
fino_distribution_t nu;    // coef de poisson
fino_distribution_t rho;   // densidad
fino_distribution_t fx;    // fuerza volumetrica en x
fino_distribution_t fy;    // fuerza volumetrica en y
fino_distribution_t fz;    // fuerza volumetrica en z
fino_distribution_t alpha; // coeficiente de expansion termica
fino_distribution_t T;     // temperatura

  
#undef  __FUNCT__
#define __FUNCT__ "fino_build_breakshake"
int fino_build_breakshake(element_t *element, int v) {

  static int J;            // cantidad de nodos locales
  // TODO: hacer un campo descripcion en finto_distribution_t para documentar
 
  // matrices de la formulacion del problema
  static gsl_matrix *C = NULL;
  static gsl_matrix *B = NULL;
  // vector para calcular las tensiones termicas
  static gsl_vector *et = NULL;

  // matriz intermedia
  static gsl_matrix *CB;
  // vector intermedio
  static gsl_vector *Cet;

  double c;
  double alphaT;
  
  double w_gauss;
  int j;

  PetscFunctionBegin;

  w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v); 
  
  // si la matriz C de la formulacion es null entonces allocamos y
  // buscamos las distribuciones espaciales de parametros
  if (C == NULL) {
    wasora_call(fino_distribution_init(&E, "E"));
    wasora_call(fino_distribution_init(&nu, "nu"));
    wasora_call(fino_distribution_init(&rho, "rho"));
    wasora_call(fino_distribution_init(&fx, "fx"));
    wasora_call(fino_distribution_init(&fy, "fy"));
    wasora_call(fino_distribution_init(&fz, "fz"));
    wasora_call(fino_distribution_init(&alpha, "alpha"));
    wasora_call(fino_distribution_init(&T, "T"));
    
    // TODO: allow lambda+mu
    if (E.defined == 0) {
      wasora_push_error_message("cannot find Young modulus 'E'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    } else if (nu.defined == 0) {
      wasora_push_error_message("cannot find Poisson coefficient 'nu'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    if (fino.math_type == math_eigen && rho.defined == 0) {
      wasora_push_error_message("cannot find density 'rho'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    C = gsl_matrix_calloc(6, 6);
    
    // si E y nu son variables, calculamos C una sola vez y ya
    if (E.variable != NULL && nu.variable != NULL) {
      wasora_call(fino_break_compute_C(C));
    }
    
    // expansion termica
    et = gsl_vector_calloc(6);
  }
  
  if (J != element->type->nodes) {
    J = element->type->nodes;
    gsl_matrix_free(B);
    B = gsl_matrix_alloc(6, 3*J);
    
    gsl_matrix_free(CB);
    CB = gsl_matrix_alloc(6, 3*J);
    
    gsl_vector_free(Cet);
    Cet = gsl_vector_alloc(6);
  }  
  
  // la H es la del framework fem, pero la B no es la misma 
  // porque la formulacion es reducida, i.e hace un 6x6 cuando deberia ser 9x9
  gsl_matrix_set_zero(B);

  for (j = 0; j < J; j++) {
    gsl_matrix_set(B, 0, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
    gsl_matrix_set(B, 1, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
    gsl_matrix_set(B, 2, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
    
    gsl_matrix_set(B, 3, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
    gsl_matrix_set(B, 3, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));

    gsl_matrix_set(B, 4, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
    gsl_matrix_set(B, 4, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));

    gsl_matrix_set(B, 5, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
    gsl_matrix_set(B, 5, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
    
    // TODO: ver si hace falta chequear que este definido
    if (fino.problem == problem_break) {
      // el vector de fuerzas volumetricas
      c = w_gauss * gsl_vector_get(fino.mesh->fem.h, j);
      gsl_vector_add_to_element(fino.bi, 3*j+0, c * fino_distribution_evaluate(&fx));
      gsl_vector_add_to_element(fino.bi, 3*j+1, c * fino_distribution_evaluate(&fy));
      gsl_vector_add_to_element(fino.bi, 3*j+2, c * fino_distribution_evaluate(&fz));
    }
    
  }
  
  // si E y nu estan dadas por variables, C es constante y no la volvemos a evaluar
  // pero si alguna es una funcion, es otro cantar
  if (E.function != NULL || nu.function != NULL) {
    wasora_call(fino_break_compute_C(C));
  }

  // calculamos Bt*C*B
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, B, 0, CB);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, B, CB, 1.0, fino.Ai);

  // expansion termica
  alphaT = fino_distribution_evaluate(&alpha)*fino_distribution_evaluate(&T);
  gsl_vector_set(et, 0, alphaT);
  gsl_vector_set(et, 1, alphaT);
  gsl_vector_set(et, 2, alphaT);
  gsl_blas_dgemv(CblasTrans, 1.0, C, et, 0, Cet);
  gsl_blas_dgemv(CblasTrans, w_gauss, B, Cet, 1.0, fino.bi);
  
  if (fino.problem == problem_shake) {
    // calculamos Ht*rho*H
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * fino_distribution_evaluate(&rho), fino.mesh->fem.H, fino.mesh->fem.H, 1.0, fino.Bi);
  } 
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_C"
int fino_break_compute_C(gsl_matrix *C) {
  
  double evaluated_nu;
  double c1, c2, c3, c4;
  
  PetscFunctionBegin;

  evaluated_nu = fino_distribution_evaluate(&nu);
  c1 = fino_distribution_evaluate(&E)/((1+evaluated_nu)*(1-2*evaluated_nu));
  c2 = c1 * evaluated_nu;
  c3 = c1 * (1-evaluated_nu);
  c4 = c1 * (1-2*evaluated_nu)/2;
    
  gsl_matrix_set(C, 0, 0, c3);
  gsl_matrix_set(C, 0, 1, c2);
  gsl_matrix_set(C, 0, 2, c2);

  gsl_matrix_set(C, 1, 0, c2);
  gsl_matrix_set(C, 1, 1, c3);
  gsl_matrix_set(C, 1, 2, c2);

  gsl_matrix_set(C, 2, 0, c2);
  gsl_matrix_set(C, 2, 1, c2);
  gsl_matrix_set(C, 2, 2, c3);
  
  gsl_matrix_set(C, 3, 3, c4);
  gsl_matrix_set(C, 4, 4, c4);
  gsl_matrix_set(C, 5, 5, c4);

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}    



#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_stresses"
int fino_break_compute_stresses(void) {
  
  double evaluatednu, evaluatedE;
  double ex, ey, ez;
  double gammaxy, gammayz, gammazx;
  double sigmax, sigmay, sigmaz;
  double tauxy, tauyz, tauzx;
  double I1, I2, I3, phi;
  double c1, c2, c3, c4, c5;
  double sigma, sigma1, sigma2, sigma3;
  double displ2;
  double max_displ2 = 0;
  int j;
  
  PetscFunctionBegin;

  // von misses  
  fino.sigma->data_argument = fino.gradient[0][0]->data_argument;
  free(fino.sigma->data_value);
  fino.sigma->data_size = fino.mesh->n_nodes;
  fino.sigma->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
  
  // principal
  fino.sigma1->data_argument = fino.gradient[0][0]->data_argument;
  free(fino.sigma1->data_value);
  fino.sigma1->data_size = fino.mesh->n_nodes;
  fino.sigma1->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
  
  fino.sigma2->data_argument = fino.gradient[0][0]->data_argument;
  free(fino.sigma2->data_value);
  fino.sigma2->data_size = fino.mesh->n_nodes;
  fino.sigma2->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

  fino.sigma3->data_argument = fino.gradient[0][0]->data_argument;
  free(fino.sigma3->data_value);
  fino.sigma3->data_size = fino.mesh->n_nodes;
  fino.sigma3->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

  
  // evaluamos nu y E, si son uniformes esto ya nos sirve para siempre
  evaluatednu = fino_distribution_evaluate(&nu);
  if (evaluatednu > 0.5) {
    wasora_push_error_message("nu is greater than 1/2");
    return WASORA_RUNTIME_ERROR;
  } else if (evaluatednu < 0) {
    wasora_push_error_message("nu is negative");
    return WASORA_RUNTIME_ERROR;
  }
  
  evaluatedE = fino_distribution_evaluate(&E);
  
  c1 = evaluatedE/((1+evaluatednu)*(1-2*evaluatednu));
  c2 = (1-2*evaluatednu)/2;
  wasora_var(fino.vars.sigma_max) = 0;
  
  for (j = 0; j < fino.mesh->n_nodes; j++) {


    if (nu.function != NULL) {
      evaluatednu = fino_distribution_evaluate(&nu);
      if (evaluatednu > 0.5) {
        wasora_push_error_message("nu is greater than 1/2");
        return WASORA_RUNTIME_ERROR;
      } else if (evaluatednu < 0) {
        wasora_push_error_message("nu is negative");
        return WASORA_RUNTIME_ERROR;
      }      
    }
    
    if (E.function != NULL) {
      evaluatedE = fino_distribution_evaluate(&nu);
    }
    
    if (nu.function != NULL || E.function != NULL) {
      c1 = fino_distribution_evaluate(&E)/((1+evaluatednu)*(1-2*evaluatednu));
      c2 = (1-2*evaluatednu)/2;
    }

    // deformaciones
/*  
e_x(x,y,z) := du_dx(x,y,z)
e_y(x,y,z) := dv_dy(x,y,z)
e_z(x,y,z) := dw_dz(x,y,z)
gamma_xy(x,y,z) := du_dy(x,y,z) + dv_dx(x,y,z)
gamma_yz(x,y,z) := dv_dz(x,y,z) + dw_dy(x,y,z)
gamma_zx(x,y,z) := dw_dx(x,y,z) + du_dz(x,y,z)
*/    
    ex = fino.gradient[0][0]->data_value[j];
    ey = fino.gradient[1][1]->data_value[j];
    ez = fino.gradient[2][2]->data_value[j];
    gammaxy = fino.gradient[0][1]->data_value[j] + fino.gradient[1][0]->data_value[j];
    gammayz = fino.gradient[1][2]->data_value[j] + fino.gradient[2][1]->data_value[j];
    gammazx = fino.gradient[2][0]->data_value[j] + fino.gradient[0][2]->data_value[j];

    // tensiones
/*
sigma_x(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_x(x,y,z) + nu*(e_y(x,y,z)+e_z(x,y,z)))
sigma_y(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_y(x,y,z) + nu*(e_x(x,y,z)+e_z(x,y,z)))
sigma_z(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_z(x,y,z) + nu*(e_x(x,y,z)+e_y(x,y,z)))
tau_xy(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_xy(x,y,z)
tau_yz(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_yz(x,y,z)
tau_zx(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_zx(x,y,z)
VM_stress(x,y,z) := sqrt(1/2*((sigma_x(x,y,z)-sigma_y(x,y,z))^2 + (sigma_y(x,y,z)-sigma_z(x,y,z))^2 + (sigma_z(x,y,z)-sigma_x(x,y,z))^2 + 6*(tau_xy(x,y,z)^2+tau_yz(x,y,z)^2+tau_zx(x,y,z)^2)))
*/
    
    sigmax = c1 * ((1-evaluatednu)*ex + evaluatednu*(ey+ez));
    sigmay = c1 * ((1-evaluatednu)*ey + evaluatednu*(ex+ez));
    sigmaz = c1 * ((1-evaluatednu)*ez + evaluatednu*(ex+ey));
    tauxy =  c1 * c2 * gammaxy;
    tauyz =  c1 * c2 * gammayz;
    tauzx =  c1 * c2 * gammazx;
    
    // stress invariants
    I1 = sigmax + sigmay + sigmaz;
    I2 = sigmax*sigmay + sigmay*sigmaz + sigmaz*sigmax - gsl_pow_2(tauxy) - gsl_pow_2(tauyz) - gsl_pow_2(tauzx);
    I3 = sigmax*sigmay*sigmaz - sigmax*gsl_pow_2(tauyz) - sigmay*gsl_pow_2(tauzx) - sigmaz*gsl_pow_2(tauxy) + 2*tauxy*tauyz*tauzx;

   // principal stresses
    c5 = sqrt(fabs(gsl_pow_2(I1) - 3*I2));
    phi = 1.0/3.0 * acos((2.0*gsl_pow_3(I1) - 9.0*I1*I2 + 27.0*I3)/(2.0*gsl_pow_3(c5)));
    if (isnan(phi)) {
      phi = 0;
    }
    
    c3 = I1/3.0;
    c4 = 2.0/3.0 * c5;
    sigma1 = c3 + c4 * cos(phi);
    sigma2 = c3 + c4 * cos(phi - 2.0*M_PI/3.0);
    sigma3 = c3 + c4 * cos(phi - 4.0*M_PI/3.0);

    fino.sigma1->data_value[j] = sigma1;
    fino.sigma2->data_value[j] = sigma2;
    fino.sigma3->data_value[j] = sigma3;
    

    
    // von misses
//    sigma = sqrt(0.5*(gsl_pow_2(sigmax-sigmay) + gsl_pow_2(sigmay-sigmaz) + gsl_pow_2(sigmaz-sigmax) +
//                                                    6.0 * (gsl_pow_2(tauxy) + gsl_pow_2(tauyz) + gsl_pow_2(tauzx))));
    sigma = sqrt(0.5*(gsl_pow_2(sigma1-sigma2) + gsl_pow_2(sigma2-sigma3) + gsl_pow_2(sigma3-sigma1)));
    
    if ((fino.sigma->data_value[j] = sigma) > wasora_var(fino.vars.sigma_max)) {
      wasora_var(fino.vars.sigma_max) = fino.sigma->data_value[j];
      
      wasora_var(fino.vars.sigma_max_x) = fino.mesh->node[j].x[0];
      wasora_var(fino.vars.sigma_max_y) = fino.mesh->node[j].x[1];
      wasora_var(fino.vars.sigma_max_z) = fino.mesh->node[j].x[2];
      
      wasora_var(fino.vars.u_at_sigma_max) = fino.solution[0]->data_value[j];
      wasora_var(fino.vars.v_at_sigma_max) = fino.solution[1]->data_value[j];
      wasora_var(fino.vars.w_at_sigma_max) = fino.solution[2]->data_value[j];
    }
    
    displ2 = gsl_pow_2(fino.solution[0]->data_value[j]) +
             gsl_pow_2(fino.solution[1]->data_value[j]) +        
             gsl_pow_2(fino.solution[2]->data_value[j]);
    
    if (displ2 > max_displ2) {
      max_displ2 = displ2;
      wasora_var(fino.vars.displ_max) = sqrt(displ2);
      wasora_var(fino.vars.displ_max_x) = fino.mesh->node[j].x[0];
      wasora_var(fino.vars.displ_max_y) = fino.mesh->node[j].x[1];
      wasora_var(fino.vars.displ_max_z) = fino.mesh->node[j].x[2];
      
      wasora_var(fino.vars.u_at_displ_max) = fino.solution[0]->data_value[j];
      wasora_var(fino.vars.v_at_displ_max) = fino.solution[1]->data_value[j];
      wasora_var(fino.vars.w_at_displ_max) = fino.solution[2]->data_value[j];
    }
    
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_reactions"
int fino_break_compute_reactions(void) {

  // TODO: hacer lo que dijo barry de traer matgetsubmatrix

  int i, k;
  PetscScalar xi;
  physical_entity_t *physical_entity;

  LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
    if (physical_entity->R[0] != NULL) {
      wasora_var_value(physical_entity->R[0]) = 0;
      wasora_var_value(physical_entity->R[1]) = 0;
      wasora_var_value(physical_entity->R[2]) = 0;
    }
  }

  for (i = 0; i < fino.n_dirichlet_rows; i++) {
    if (fino.dirichlet_row[i].physical_entity->R[0] != NULL) {
      for (k = 0; k < fino.dirichlet_row[i].ncols; k++) {
        petsc_call(VecGetValues(fino.phi, 1, &fino.dirichlet_row[i].cols[k], &xi));
        wasora_var_value(fino.dirichlet_row[i].physical_entity->R[fino.dirichlet_row[i].dof]) += fino.dirichlet_row[i].vals[k] * xi;
      }
    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}
