/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of linear elastic problem (break) with optional vibration (shake)
 *  and evaluation of the stress tensor out of the gradients of the displacements
 *
 *  Copyright (C) 2015--2017 jeremy theler
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

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

double T0;  // este es el escalar

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
  
  material_t *material;

  double c;
  double alphaDT;
  
  double w_gauss;
  double r_for_axisymmetric = 1.0;
  int j;

  PetscFunctionBegin;
  
  if (element->physical_entity != NULL && element->physical_entity->material != NULL) {
    material =  element->physical_entity->material;
  } else {
    material = NULL;
  }
  
  w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v); 
  
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
    
    // TODO: allow lambda+mu
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
      wasora_call(fino_break_compute_C(C, fino_distribution_evaluate(&distribution_E, material,NULL), fino_distribution_evaluate(&distribution_nu, material,NULL)));
    }
    
  }
  
  if (J != element->type->nodes) {
    J = element->type->nodes;
    gsl_matrix_free(B);
    B = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);
    
    gsl_matrix_free(CB);
    CB = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);
    
    gsl_vector_free(Cet);
    Cet = gsl_vector_alloc(stress_strain_size);
  }  
  
  // la H es la del framework fem, pero la B no es la misma 
  // porque la formulacion es reducida, i.e hace un 6x6 (o 3x3) cuando deberia ser 9x9
  
  gsl_matrix_set_zero(B);

  for (j = 0; j < J; j++) {
    if (fino.problem_kind == problem_kind_full3d) {
      gsl_matrix_set(B, 0, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
      
      gsl_matrix_set(B, 2, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
    
      gsl_matrix_set(B, 3, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
      gsl_matrix_set(B, 3, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));

      gsl_matrix_set(B, 4, 3*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
      gsl_matrix_set(B, 4, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));

      gsl_matrix_set(B, 5, 3*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 2));
      gsl_matrix_set(B, 5, 3*j+2, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
    
    } else if (fino.problem_kind == problem_kind_axisymmetric) {

      r_for_axisymmetric = fino_compute_r_for_axisymmetric();
      
      // ecuacion 3.5 AFEM CH.03 sec 3.3.2 pag 3.5
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));

      if (fino.symmetry_axis == symmetry_axis_y) {
        gsl_matrix_set(B, 2, 2*j+0, gsl_vector_get(fino.mesh->fem.h, j)/r_for_axisymmetric);
      } else if (fino.symmetry_axis == symmetry_axis_x) {
        gsl_matrix_set(B, 2, 2*j+1, gsl_vector_get(fino.mesh->fem.h, j)/r_for_axisymmetric);
      }
      
      gsl_matrix_set(B, 3, 2*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
      gsl_matrix_set(B, 3, 2*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));

    } else  {
      // plane stress y plane strain son iguales
      // ecuacion 14.18 IFEM CH.14 sec 14.4.1 pag 14-11
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
    
      gsl_matrix_set(B, 2, 2*j+0, gsl_matrix_get(fino.mesh->fem.dhdx, j, 1));
      gsl_matrix_set(B, 2, 2*j+1, gsl_matrix_get(fino.mesh->fem.dhdx, j, 0));
      
    }
    
    if ((fino.problem_family == problem_family_break) &&
        (distribution_fx.defined != 0 || distribution_fy.defined != 0 || distribution_fz.defined != 0)) {
      // el vector de fuerzas volumetricas
      c = r_for_axisymmetric * w_gauss * gsl_vector_get(fino.mesh->fem.h, j);
      gsl_vector_add_to_element(fino.bi, fino.degrees*j+0, c * fino_distribution_evaluate(&distribution_fx, material, gsl_vector_ptr(fino.mesh->fem.x, 0)));
      gsl_vector_add_to_element(fino.bi, fino.degrees*j+1, c * fino_distribution_evaluate(&distribution_fy, material, gsl_vector_ptr(fino.mesh->fem.x, 0)));
      if (fino.degrees == 3) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+2, c * fino_distribution_evaluate(&distribution_fz, material, gsl_vector_ptr(fino.mesh->fem.x, 0)));
      }
    }
    
  }
  
  // si E y nu estan dadas por variables, C es constante y no la volvemos a evaluar
  // pero si alguna es una propiedad o una funcion, es otro cantar
  if (distribution_E.variable == NULL || distribution_nu.variable == NULL) {
     if (material == NULL) {
       wasora_push_error_message("element %d does not have an associated material", element->id);
       PetscFunctionReturn(WASORA_RUNTIME_ERROR);
     }

    wasora_call(fino_break_compute_C(C, fino_distribution_evaluate(&distribution_E, material, gsl_vector_ptr(fino.mesh->fem.x, 0)), fino_distribution_evaluate(&distribution_nu, material, gsl_vector_ptr(fino.mesh->fem.x, 0))));
  }

  // calculamos Bt*C*B
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, r_for_axisymmetric, C, B, 0, CB);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, B, CB, 1.0, fino.Ki);

  // expansion termica
  if (distribution_alpha.defined != 0) {
    // este debe ser el medio!
    alphaDT = fino_distribution_evaluate(&distribution_alpha, material, gsl_vector_ptr(fino.mesh->fem.x, 0));
    if (alphaDT != 0) {
      alphaDT *= fino_distribution_evaluate(&distribution_T, material, gsl_vector_ptr(fino.mesh->fem.x, 0))-T0;
      gsl_vector_set(et, 0, alphaDT);
      gsl_vector_set(et, 1, alphaDT);
      gsl_vector_set(et, 2, alphaDT);
      gsl_blas_dgemv(CblasTrans, r_for_axisymmetric, C, et, 0, Cet);
      gsl_blas_dgemv(CblasTrans, w_gauss, B, Cet, 1.0, fino.bi);
    }
  }
  
  if (fino.has_mass) {
    // calculamos la matriz de masa Ht*rho*H
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * fino_distribution_evaluate(&distribution_rho, material, gsl_vector_ptr(fino.mesh->fem.x, 0)), fino.mesh->fem.H, fino.mesh->fem.H, 1.0, fino.Mi);
  } 
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_C"
int fino_break_compute_C(gsl_matrix *C, double E, double nu) {
  
  double c1, c2, c3, c4;
  
  PetscFunctionBegin;
  
  // tabla 4.3 pag 194 Bathe

  if (fino.problem_kind == problem_kind_full3d) {
  
    c1 = E/((1+nu)*(1-2*nu));
    c2 = c1 * nu;
    c3 = c1 * (1-nu);
    c4 = c1 * (1-2*nu)/2;
    
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
    
  } else if (fino.problem_kind == problem_kind_plane_stress) {
    
    c1 = E/(1-nu*nu);
    c2 = nu * c1;
    gsl_matrix_set(C, 0, 0, c1);
    gsl_matrix_set(C, 0, 1, c2);
    
    gsl_matrix_set(C, 1, 0, c2);
    gsl_matrix_set(C, 1, 1, c1);

    gsl_matrix_set(C, 2, 2, c1*0.5*(1-nu));
    
  } else if (fino.problem_kind == problem_kind_plane_strain) {
    
    c1 = E*(1-nu)/((1+nu)*(1-2*nu));
    c2 = nu/(1-nu) * c1;
    gsl_matrix_set(C, 0, 0, c1);
    gsl_matrix_set(C, 0, 1, c2);
    
    gsl_matrix_set(C, 1, 0, c2);
    gsl_matrix_set(C, 1, 1, c1);

    gsl_matrix_set(C, 2, 2, c1*0.5*(1-2*nu)/(1-nu));
    
  } else if (fino.problem_kind == problem_kind_axisymmetric) {
    
    c1 = E*(1-nu)/((1+nu)*(1-2*nu));
    c2 = nu/(1-nu) * c1;
    gsl_matrix_set(C, 0, 0, c1);
    gsl_matrix_set(C, 0, 1, c2);
    gsl_matrix_set(C, 0, 2, c2);
    
    gsl_matrix_set(C, 1, 0, c2);
    gsl_matrix_set(C, 1, 1, c1);
    gsl_matrix_set(C, 1, 2, c2);

    gsl_matrix_set(C, 2, 0, c2);
    gsl_matrix_set(C, 2, 1, c2);
    gsl_matrix_set(C, 2, 2, c1);

    gsl_matrix_set(C, 3, 3, c1*0.5*(1-2*nu)/(1-nu));
    
  }

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}    



#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_stresses"
int fino_break_compute_stresses(void) {
  
  double ex = 0;
  double ey = 0;
  double ez = 0;
  double gammaxy = 0;
  double gammayz = 0;
  double gammazx = 0;
  
  double sigmax = 0;
  double sigmay = 0;
  double sigmaz = 0;
  double tauxy = 0;
  double tauyz = 0;
  double tauzx = 0;
  
  double c1 = 1.0;
  double c1c2 = 1.0;
  double sigma = 0;
  
  double sigma1 = 0;
  double sigma2 = 0;
  double sigma3 = 0;
  double tresca = 0;
  
  double displ2 = 0;

  double max_displ2 = 0;
  double nu = 0;
  double E = 0;
  
  int j, m;
  
  PetscFunctionBegin;

  // von mises  
  if (fino.sigma->data_value == NULL) {
    // tensor de tensiones
    fino.sigmax->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigmax->data_size = fino.mesh->n_nodes;
    fino.sigmax->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

    fino.sigmay->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigmay->data_size = fino.mesh->n_nodes;
    fino.sigmay->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

    fino.sigmaz->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigmaz->data_size = fino.mesh->n_nodes;
    fino.sigmaz->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
      
    fino.tauxy->data_argument = fino.gradient[0][0]->data_argument;
    fino.tauxy->data_size = fino.mesh->n_nodes;
    fino.tauxy->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
    
    if (fino.dimensions == 3) {
      fino.tauyz->data_argument = fino.gradient[0][0]->data_argument;
      fino.tauyz->data_size = fino.mesh->n_nodes;
      fino.tauyz->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
      
      fino.tauzx->data_argument = fino.gradient[0][0]->data_argument;
      fino.tauzx->data_size = fino.mesh->n_nodes;
      fino.tauzx->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
    }

    // tensiones principales
    fino.sigma1->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigma1->data_size = fino.mesh->n_nodes;
    fino.sigma1->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

    fino.sigma2->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigma2->data_size = fino.mesh->n_nodes;
    fino.sigma2->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

    fino.sigma3->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigma3->data_size = fino.mesh->n_nodes;
    fino.sigma3->data_value = calloc(fino.mesh->n_nodes, sizeof(double));

    // von mises
    fino.sigma->data_argument = fino.gradient[0][0]->data_argument;
    fino.sigma->data_size = fino.mesh->n_nodes;
    fino.sigma->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
    
    // tresca
    fino.tresca->data_argument = fino.gradient[0][0]->data_argument;
    fino.tresca->data_size = fino.mesh->n_nodes;
    fino.tresca->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
  }
  
  // evaluamos nu y E, si son uniformes esto ya nos sirve para siempre
  if (distribution_nu.variable != NULL) {
    nu = fino_distribution_evaluate(&distribution_nu, NULL, NULL);
    if (nu > 0.5) {
      wasora_push_error_message("nu is greater than 1/2");
      return WASORA_RUNTIME_ERROR;
    } else if (nu < 0) {
      wasora_push_error_message("nu is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
  if (distribution_E.variable != NULL) {
    E = fino_distribution_evaluate(&distribution_E, NULL, NULL);
    if (E < 0) {
      wasora_push_error_message("E is negative (%g)", E);
      return WASORA_RUNTIME_ERROR;
    }
  }
  
 
  wasora_var(fino.vars.sigma_max) = 0;
  
  for (j = 0; j < fino.mesh->n_nodes; j++) {

    wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
    wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
    wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
    
    if (distribution_nu.physical_property != NULL) {
      nu = fino_distribution_evaluate(&distribution_nu, fino.mesh->node[j].master_material, fino.mesh->node[j].x);
      
      if (nu > 0.5) {
        wasora_push_error_message("nu is greater than 1/2 at node %d", j+1);
        return WASORA_RUNTIME_ERROR;
      } else if (nu < 0) {
        wasora_push_error_message("nu is negative at node %d", j+1);
        return WASORA_RUNTIME_ERROR;
      }      
    }
    
    if (distribution_E.physical_property != NULL) {
      E = fino_distribution_evaluate(&distribution_E, fino.mesh->node[j].master_material, fino.mesh->node[j].x);
      
      if (E < 0) {
        wasora_push_error_message("E is negative at node %d", j+1);
        return WASORA_RUNTIME_ERROR;
      }      
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
    if (fino.problem_kind == problem_kind_full3d) {
      ez = fino.gradient[2][2]->data_value[j];
    } else if (fino.problem_kind == problem_kind_axisymmetric) {
      if (fino.symmetry_axis == symmetry_axis_y) {
        // etheta = u/r
        if (fino.solution[0]->data_argument[0][j] > 1e-3) {
          ez = fino.solution[0]->data_value[j]/fino.solution[0]->data_argument[0][j];
        }
      } else if (fino.symmetry_axis == symmetry_axis_x) {
        // etheta = v/r
        if (fino.solution[1]->data_argument[1][j] > 1e-3) {
          ez = fino.solution[1]->data_value[j]/fino.solution[1]->data_argument[1][j];
        }
      }
    }
    
    gammaxy = fino.gradient[0][1]->data_value[j] + fino.gradient[1][0]->data_value[j];
    if (fino.problem_kind == problem_kind_full3d) {
      gammayz = fino.gradient[1][2]->data_value[j] + fino.gradient[2][1]->data_value[j];
      gammazx = fino.gradient[2][0]->data_value[j] + fino.gradient[0][2]->data_value[j];
    }
    
    
    // tensiones
    if (fino.problem_kind == problem_kind_axisymmetric || fino.problem_kind == problem_kind_full3d) {
  /*
  sigma_x(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_x(x,y,z) + nu*(e_y(x,y,z)+e_z(x,y,z)))
  sigma_y(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_y(x,y,z) + nu*(e_x(x,y,z)+e_z(x,y,z)))
  sigma_z(x,y,z) := E/((1+nu)*(1-2*nu))*((1-nu)*e_z(x,y,z) + nu*(e_x(x,y,z)+e_y(x,y,z)))
  tau_xy(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_xy(x,y,z)
  tau_yz(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_yz(x,y,z)
  tau_zx(x,y,z) :=  E/((1+nu)*(1-2*nu))*(1-2*nu)/2*gamma_zx(x,y,z)
  VM_stress(x,y,z) := sqrt(1/2*((sigma_x(x,y,z)-sigma_y(x,y,z))^2 + (sigma_y(x,y,z)-sigma_z(x,y,z))^2 + (sigma_z(x,y,z)-sigma_x(x,y,z))^2 + 6*(tau_xy(x,y,z)^2+tau_yz(x,y,z)^2+tau_zx(x,y,z)^2)))
  */
      // constantes    
      c1 = E/((1+nu)*(1-2*nu));
      c1c2 = c1 * 0.5*(1-2*nu);
    }

    sigmax = c1 * ((1-nu)*ex + nu*(ey+ez));
    sigmay = c1 * ((1-nu)*ey + nu*(ex+ez));
    sigmaz = c1 * ((1-nu)*ez + nu*(ex+ey));  // esta es sigmatheta en axi
    tauxy =  c1c2 * gammaxy;
    if (fino.dimensions == 3) {
      tauyz =  c1c2 * gammayz;
      tauzx =  c1c2 * gammazx;
    }

    wasora_call(fino_compute_principal_stress(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx, &sigma1, &sigma2, &sigma3));

    fino.sigmax->data_value[j] = sigmax;
    fino.sigmay->data_value[j] = sigmay;
    fino.tauxy->data_value[j] = tauxy;
    fino.sigmaz->data_value[j] = sigmaz;
    
    if (fino.dimensions == 3) {
      fino.tauyz->data_value[j] = tauyz;
      fino.tauzx->data_value[j] = tauzx;
    }

    fino.sigma1->data_value[j] = sigma1;
    fino.sigma2->data_value[j] = sigma2;
    fino.sigma3->data_value[j] = sigma3;

    // tresca
    tresca = fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);
    fino.tresca->data_value[j] = tresca;

    // von mises
    sigma = fino_compute_vonmises_from_principal(sigma1, sigma2, sigma3);
    //sigma = fino_compute_vonmises_from_tensor(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx);
      
    
    if ((fino.sigma->data_value[j] = sigma) > wasora_var(fino.vars.sigma_max)) {
      wasora_var(fino.vars.sigma_max) = fino.sigma->data_value[j];
      
      wasora_var(fino.vars.sigma_max_x) = fino.mesh->node[j].x[0];
      wasora_var(fino.vars.sigma_max_y) = fino.mesh->node[j].x[1];
      wasora_var(fino.vars.sigma_max_z) = fino.mesh->node[j].x[2];
      
      wasora_var(fino.vars.u_at_sigma_max) = fino.solution[0]->data_value[j];
      wasora_var(fino.vars.v_at_sigma_max) = fino.solution[1]->data_value[j];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_sigma_max) = fino.solution[2]->data_value[j];
      }
    }
    
    displ2 = 0;
    for (m = 0; m < fino.dimensions; m++) {
      displ2 += gsl_pow_2(fino.solution[m]->data_value[j]);
    }
    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice displ_max
    if (displ2 >= max_displ2) {
      max_displ2 = displ2;
      wasora_var(fino.vars.displ_max) = sqrt(displ2);
      wasora_var(fino.vars.displ_max_x) = fino.mesh->node[j].x[0];
      wasora_var(fino.vars.displ_max_y) = fino.mesh->node[j].x[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.displ_max_z) = fino.mesh->node[j].x[2];
      }
      
      wasora_var(fino.vars.u_at_displ_max) = fino.solution[0]->data_value[j];
      wasora_var(fino.vars.v_at_displ_max) = fino.solution[1]->data_value[j];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_displ_max) = fino.solution[2]->data_value[j];
      }
    }
    
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}



#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_stress"
int fino_break_set_stress(element_t *element) {
  int v, g;
  double w_gauss;
  double r_for_axisymmetric;  
  gsl_vector *Nb;
    
  if ((fino.dimensions == 3 && element->type->dim != 2) ||
      (fino.dimensions == 2 && element->type->dim != 1)) {
    wasora_push_error_message("stress BCs can only be applied to surfaces");
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }

  Nb = gsl_vector_calloc(fino.degrees);
  gsl_vector_set_zero(fino.bi);

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric();
    mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
    mesh_update_coord_vars(gsl_vector_ptr(fino.mesh->fem.x, 0));

    for (g = 0; g < fino.degrees; g++) {
      gsl_vector_set(Nb, g, wasora_evaluate_expression(&element->physical_entity->bc_args[g]));
    }
    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_moment"
int fino_break_set_moment(element_t *element) {
  int v;
  double w_gauss;
  double theta, dx, dy, dz, M;
  double r_for_axisymmetric;  
  gsl_vector *Nb;
    
  if ((fino.dimensions == 3 && element->type->dim != 2) ||
      (fino.dimensions == 2 && element->type->dim != 1)) {
    wasora_push_error_message("moment BCs can only be applied to surfaces");
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }

  Nb = gsl_vector_calloc(fino.degrees);
  gsl_vector_set_zero(fino.bi);

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric();
    mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
    mesh_update_coord_vars(gsl_vector_ptr(fino.mesh->fem.x, 0));

    dx = gsl_vector_get(fino.mesh->fem.x, 0) - element->physical_entity->cog[0];
    dy = gsl_vector_get(fino.mesh->fem.x, 1) - element->physical_entity->cog[1];
    dz = gsl_vector_get(fino.mesh->fem.x, 2) - element->physical_entity->cog[2];
    
    gsl_vector_set_zero(Nb);

    // los tres primeros tienen las componentes Mx My y Mz
    if (element->physical_entity->bc_args[0].n_tokens != 0) {
      M = wasora_evaluate_expression(&element->physical_entity->bc_args[0]);
      theta = atan2(dy, dz);
      // ty = cos(theta)*dz
      gsl_vector_add_to_element(Nb, 1, -M*cos(theta));
      // tz = -sin(theta)*dy
      gsl_vector_add_to_element(Nb, 2, +M*sin(theta));
    }
    
    if (element->physical_entity->bc_args[1].n_tokens != 0) {
      M = wasora_evaluate_expression(&element->physical_entity->bc_args[1]);
      theta = atan2(dx, dz);
      gsl_vector_add_to_element(Nb, 0, -M*cos(theta));
      gsl_vector_add_to_element(Nb, 2, +M*sin(theta));
    }
    
    if (element->physical_entity->bc_args[2].n_tokens != 0) {
      M = wasora_evaluate_expression(&element->physical_entity->bc_args[2]);
      theta = atan2(dx, dy);
      gsl_vector_add_to_element(Nb, 0, -M*cos(theta));
      gsl_vector_add_to_element(Nb, 1, +M*sin(theta));
    }
    
    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_force"
int fino_break_set_force(element_t *element) {
  int v, g;
  double w_gauss;
  double r_for_axisymmetric;
  gsl_vector *Nb;
    
  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }
  
  Nb = gsl_vector_calloc(fino.degrees);
  gsl_vector_set_zero(fino.bi);
  

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric();

    mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
    mesh_update_coord_vars(gsl_vector_ptr(fino.mesh->fem.x, 0));

    for (g = 0; g < fino.degrees; g++) {
      gsl_vector_set(Nb, g, wasora_evaluate_expression(&element->physical_entity->bc_args[g])/element->physical_entity->volume);
    }
    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_pressure"
int fino_break_set_pressure(element_t *element) {
  double w_gauss;
  double p;
  double r_for_axisymmetric;
  int v;
  gsl_vector *Nb;

  if ((fino.dimensions == 3 && element->type->dim != 2) ||
      (fino.dimensions == 2 && element->type->dim != 1)) {
    wasora_push_error_message("pressure BCs can only be applied to surfaces");
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }

  Nb = gsl_vector_calloc(fino.degrees);
  gsl_vector_set_zero(fino.bi);
  
  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric();
    mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
    mesh_update_coord_vars(gsl_vector_ptr(fino.mesh->fem.x, 0));
    
    // la p chica es la proyeccion del vector tension sobre la normal, lo que uno espera en matematica
    p = wasora_evaluate_expression(&  element->physical_entity->bc_args[0]);
    // la P grande es presion positiva cuando comprime, como lo que uno espera en ingenieria
    if (element->physical_entity->bc_type_phys == bc_phys_pressure_real) {
      p = -p;
    }
    gsl_vector_set(Nb, 0, wasora_var_value(fino.vars.nx) * p);
    gsl_vector_set(Nb, 1, wasora_var_value(fino.vars.ny) * p);
    if (fino.dimensions == 3) {
      gsl_vector_set(Nb, 2, wasora_var_value(fino.vars.nz) * p);
    }
    
    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  
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
  *sigma1 = c2 + c3 * cos(phi);
  *sigma2 = c2 + c3 * cos(phi - 2.0*M_PI/3.0);
  *sigma3 = c2 + c3 * cos(phi - 4.0*M_PI/3.0);

//  printf("%g %g %g %g %g %g %g\n", I1, I2, I3, phi, *sigma1, *sigma2, *sigma3);
  
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


