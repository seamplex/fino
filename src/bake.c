/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of heat conduction problem (bake)
 *
 *  Copyright (C) 2015--2018 jeremy theler & ezequiel manavela chiapero
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

fino_distribution_t distribution_k;     // conductividad
fino_distribution_t distribution_Q;     // heat source
fino_distribution_t distribution_kappa; // thermal diffusivity
fino_distribution_t distribution_rho;   // density
fino_distribution_t distribution_cp;    // heat capacity

#undef  __FUNCT__
#define __FUNCT__ "fino_bake_step_initial"
int fino_bake_step_initial(void) {
  
  int j;
  double xi;
  function_t *ic;

  if ((ic = wasora_get_function_ptr("T_0")) != NULL) {

    if (ic->n_arguments != fino.dimensions) {
      wasora_push_error_message("initial condition function T_0 ought to have %d arguments instead of %d", fino.dimensions, ic->n_arguments);
      return WASORA_RUNTIME_ERROR;
    }

    for (j = 0; j < fino.mesh->n_nodes; j++) {
      xi = wasora_evaluate_function(ic, fino.mesh->node[j].x);
      VecSetValue(fino.phi, fino.mesh->node[j].index[0], xi, INSERT_VALUES);
    }

  } else {
    
    wasora_call(fino_build_bulk());           // ensamblamos objetos elementales
    wasora_call(fino_set_essential_bc(fino.K, fino.b));     // condiciones de contorno esenciales
    wasora_call(fino_solve_linear_petsc(fino.K, fino.b));
    
    // este ksp ya no sirve mas, porque despues usamos otras matrice y demas
    petsc_call(KSPDestroy(&fino.ksp));
    fino.ksp = NULL;

  }
    
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_bake_step_transient"
int fino_bake_step_transient(void) {

  // resolvemos 
  //   A*T(n+1) = b
  // con
  //   A = K*theta + C/dt 
  //   b = (-K*(1-theta)+C/dt)*T(n) + b

  // TODO: elegir esto desde input
  double theta = 0.5;
  int nonlinear = 1;

  // TODO: ver como hacer esto mas eficiente
  petsc_call(MatZeroEntries(fino.K));
  petsc_call(MatZeroEntries(fino.M));
  petsc_call(VecZeroEntries(fino.b));
  
  wasora_call(fino_build_bulk());           // ensamblamos objetos elementales

  if (fino.has_transient == 0) {
    fino.has_transient = 1;
    MatDuplicate(fino.K, MAT_COPY_VALUES, &fino.A);
    MatDuplicate(fino.K, MAT_COPY_VALUES, &fino.B);
    VecDuplicate(fino.b, &fino.c);
    if (nonlinear) {
      MatDuplicate(fino.M, MAT_COPY_VALUES, &fino.lastM);
      MatDuplicate(fino.M, MAT_COPY_VALUES, &fino.dotM);
      VecDuplicate(fino.b, &fino.m);
    }
  } else {
    MatCopy(fino.K, fino.A, SAME_NONZERO_PATTERN);
    MatCopy(fino.K, fino.B, SAME_NONZERO_PATTERN);
    if (nonlinear) {
      MatCopy(fino.M, fino.dotM, SAME_NONZERO_PATTERN);
    }
  }

  if (nonlinear) {
    // correccion de Mpunto
    MatAXPY(fino.dotM, -1, fino.lastM, SAME_NONZERO_PATTERN);
    MatScale(fino.dotM, 1/wasora_var_value(wasora_special_var(dt)));
  
    // hacemos un K nuevo como K = K-Mdot
    //  MatAXPY(fino.A, -1, fino.dotM, SUBSET_NONZERO_PATTERN);
    //  MatAXPY(fino.B, -1, fino.dotM, SUBSET_NONZERO_PATTERN);
  }
 
  MatScale(fino.A, theta);
  MatAXPY(fino.A, 1/wasora_var_value(wasora_special_var(dt)), fino.M, SUBSET_NONZERO_PATTERN);

  MatScale(fino.B, -(1-theta));
  MatAXPY(fino.B, 1/wasora_var_value(wasora_special_var(dt)), fino.M, SUBSET_NONZERO_PATTERN);

  fino_assembly();
  
  MatMult(fino.B, fino.phi, fino.c);
  VecAXPY(fino.c, 1, fino.b);
/*  
  printf("phi\n");
  fino_print_petsc_vector(fino.phi, PETSC_VIEWER_STDOUT_SELF);
  printf("b\n");
  fino_print_petsc_vector(fino.b, PETSC_VIEWER_STDOUT_SELF);
  printf("B\n");
  fino_print_petsc_matrix(fino.B, PETSC_VIEWER_STDOUT_SELF);
  printf("c\n");
  fino_print_petsc_vector(fino.c, PETSC_VIEWER_STDOUT_SELF);
*/  
  if (nonlinear) {
    // o hacemos un c nuevo como c = c + Mdot*T
    MatMult(fino.dotM, fino.phi, fino.m);
    VecAXPY(fino.c, 1, fino.m);
  }

  fino_assembly();

  // hay que volver a poner esta
  wasora_call(fino_set_essential_bc(fino.A, fino.c));

/*  
  printf("A\n");
  fino_print_petsc_matrix(fino.A, PETSC_VIEWER_STDOUT_SELF);
  printf("c\n");
  fino_print_petsc_vector(fino.c, PETSC_VIEWER_STDOUT_SELF);
 */
  
  // y resolver
  wasora_call(fino_solve_linear_petsc(fino.A, fino.c));

/*  
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
  MatView(fino.A, PETSC_VIEWER_STDOUT_SELF);
  VecView(fino.c, PETSC_VIEWER_STDOUT_SELF);
  VecView(fino.phi, PETSC_VIEWER_STDOUT_SELF);
*/  
/*  
  printf("phi\n");
  fino_print_petsc_vector(fino.phi, PETSC_VIEWER_STDOUT_SELF);
*/
  
  if (nonlinear) {
    MatCopy(fino.M, fino.lastM, SAME_NONZERO_PATTERN);
  }

  // esto viene aca porque el ksp tiene que ser != null, capaz que se pueda meter en el solver  
  petsc_call(KSPSetInitialGuessNonzero(fino.ksp, PETSC_TRUE));
 
/*  
  MatDestroy(&fino.A);
  MatDestroy(&fino.B);
  VecDestroy(&fino.c);
  fino.has_transient = 0;
*/
        
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_bake"
int fino_build_bake(element_t *element, int v) {
  
  double w_gauss;
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
    PetscFunctionReturn(WASORA_RUNTIME_ERROR);
  }

  if (fino.has_mass) {
    if (distribution_kappa.defined == 0) {
      wasora_call(fino_distribution_init(&distribution_kappa, "kappa"));
      if (distribution_kappa.defined == 0) {
        wasora_call(fino_distribution_init(&distribution_rho, "rho"));
        wasora_call(fino_distribution_init(&distribution_cp, "cp"));
        if (distribution_cp.defined == 0) {
          wasora_push_error_message("cannot find neither thermal diffusivity 'kappa' nor heat capacity 'cp'");
          PetscFunctionReturn(WASORA_RUNTIME_ERROR);
        }
        if (distribution_rho.defined == 0) {
          wasora_push_error_message("cannot find neither thermal diffusivity 'kappa' nor density 'rho'");
          PetscFunctionReturn(WASORA_RUNTIME_ERROR);
        }
      }
    }
  }
  
  
  if (element->physical_entity != NULL && element->physical_entity->material != NULL) {
    material =  element->physical_entity->material;
  }
  
  w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v); 
  r_for_axisymmetric = fino_compute_r_for_axisymmetric();

  if (distribution_Q.defined != 0) {
    // el vector de fuente de calor volumetrica
    for (j = 0; j < element->type->nodes; j++) {
      gsl_vector_add_to_element(fino.bi, j, w_gauss * r_for_axisymmetric * gsl_vector_get(fino.mesh->fem.h, j) * fino_distribution_evaluate(&distribution_Q, material, gsl_vector_ptr(fino.mesh->fem.x, 0)));
    }
  }

  // calculamos la matriz de stiffness
  k = fino_distribution_evaluate(&distribution_k, material, gsl_vector_ptr(fino.mesh->fem.x, 0));
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * r_for_axisymmetric * k, fino.mesh->fem.B, fino.mesh->fem.B, 1.0, fino.Ki);

  if (fino.has_mass) {
    // calculamos la matriz de masa Ht*rho*cp*H
    if (distribution_kappa.defined)  {
      rhocp = k / fino_distribution_evaluate(&distribution_kappa, material, gsl_vector_ptr(fino.mesh->fem.x, 0));
    } else {
      rhocp = fino_distribution_evaluate(&distribution_rho, material, gsl_vector_ptr(fino.mesh->fem.x, 0)) * fino_distribution_evaluate(&distribution_cp, material, gsl_vector_ptr(fino.mesh->fem.x, 0));
    }
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * r_for_axisymmetric *rhocp, fino.mesh->fem.H, fino.mesh->fem.H, 1.0, fino.Mi);
  } 
  
  return WASORA_RUNTIME_OK;
  
}


#undef  __FUNCT__
#define __FUNCT__ "fino_bake_set_heat_flux"
int fino_bake_set_heat_flux(element_t *element) {
  double w_gauss;
  double r_for_axisymmetric;
  double q;
  int v;
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
    
    q = wasora_evaluate_expression(&element->physical_entity->bc_args[0]);
    gsl_vector_set(Nb, 0, -q);
    
    gsl_blas_dgemv(CblasTrans, w_gauss * r_for_axisymmetric, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  
  return WASORA_RUNTIME_OK;
}



#undef  __FUNCT__
#define __FUNCT__ "fino_bake_set_heat_flux"
int fino_bake_set_convection(element_t *element) {
  double w_gauss;
  double r_for_axisymmetric;
  double h = 0;
  double Tinf = 0;
  int v;
  gsl_matrix *Na;
  gsl_matrix *NaH;
  gsl_vector *Nb;
  
  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }
 
  Na = gsl_matrix_calloc(fino.degrees, fino.degrees);
  NaH = gsl_matrix_calloc(fino.degrees, fino.n_local_nodes);
  Nb = gsl_vector_calloc(fino.degrees);
  
  gsl_matrix_set_zero(fino.Ki);
  gsl_vector_set_zero(fino.bi);

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric();
    mesh_compute_x(element, fino.mesh->fem.r, fino.mesh->fem.x);
    mesh_update_coord_vars(gsl_vector_ptr(fino.mesh->fem.x, 0));
    
    if (element->physical_entity->bc_args != NULL) {
      h = wasora_evaluate_expression(&element->physical_entity->bc_args[0]);
      Tinf = wasora_evaluate_expression(&element->physical_entity->bc_args[1]);
    }
    
    gsl_matrix_set(Na, 0, 0, h);
    gsl_vector_set(Nb, 0, h*Tinf);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Na, fino.mesh->fem.H, 0, NaH);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss * r_for_axisymmetric, fino.mesh->fem.H, NaH, 1, fino.Ki);
    gsl_blas_dgemv(CblasTrans, w_gauss * r_for_axisymmetric, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
  }

  MatSetValues(fino.K, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Ki, 0, 0), ADD_VALUES);
  VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);

  gsl_vector_free(Nb);
  gsl_matrix_free(Na);
  gsl_matrix_free(NaH);
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_bake_compute_fluxes"
int fino_bake_compute_fluxes(void) {
  
  material_t *material = NULL;
  element_list_item_t *associated_element = NULL;
  
  double k;
  double T_max = -INFTY;
  double T_min = +INFTY;
  
  int j;
  
  PetscFunctionBegin;
  
  // evaluamos nu y E, si son uniformes esto ya nos sirve para siempre
  if (distribution_k.variable != NULL) {
    k = fino_distribution_evaluate(&distribution_k, NULL, NULL);
    if (k < 0) {
      wasora_push_error_message("k is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
 
  for (j = 0; j < fino.mesh->n_nodes; j++) {
    
    wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
    wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
    wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
    
    material = NULL;
    if (distribution_k.physical_property != NULL) {
      LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element)  {
        if (associated_element->element->type->dim == fino.dimensions &&
            associated_element->element->physical_entity != NULL) {
          material = associated_element->element->physical_entity->material;
        }
      }
      if (material == NULL) {
        wasora_push_error_message("cannot find a material for node %d", fino.mesh->node[j].id);
        return WASORA_RUNTIME_ERROR;
      }
      k = fino_distribution_evaluate(&distribution_k, material, fino.mesh->node[j].x);
      if (k < 0) {
        wasora_push_error_message("k is negative");
        return WASORA_RUNTIME_ERROR;
      }      
    }
    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice T_max
    if (fino.solution[0]->data_value[j] >= T_max) {
      T_max = fino.solution[0]->data_value[j];
    }
    if (fino.solution[0]->data_value[j] <= T_min) {
      T_min = fino.solution[0]->data_value[j];
    }
    
  }

  wasora_var(fino.vars.T_max) = T_max;
  wasora_var(fino.vars.T_min) = T_min;

  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

