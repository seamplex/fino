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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <petscts.h>

#include "fino.h"

fino_distribution_t distribution_k;     // conductivity 
fino_distribution_t distribution_Q;     // volumetric heat source
fino_distribution_t distribution_kappa; // thermal diffusivity
fino_distribution_t distribution_rho;   // density
fino_distribution_t distribution_cp;    // heat capacity

PetscErrorCode fino_thermal_IFunctionHeat(TS ts, PetscReal t, Vec T, Vec T_dot, Vec r, void *ctx);
PetscErrorCode fino_thermal_IJacobianHeat(TS ts, PetscReal t, Vec T, Vec T_dot, PetscReal s, Mat A, Mat B,void *ctx);

int fino_bc_process_thermal(bc_t *bc, char *name, char *expr, char *equal_sign) {

  int i;
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
             (strcmp(name, "Tref") == 0) ||
             (strcmp(name, "Tinf") == 0)) {

    // convection, robins
    
    bc->type_math = bc_math_robin;
    bc->type_phys = bc_phys_convection;

    // conveccion necesita dos expresione
    // las alocamos en la primera de las BCs
    base_bc = bc;
    base_bc->expr = calloc(2, sizeof(expr_t));

    do {
      // volvemos a poner el equal sign, en la primera pasada
      // es para volver a parser, en las siguientes para no romper
      // la ultima se vuelve a arreglar fuera del loop grande
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

    // ahora bc quedo apuntando a null, tenemos que volver para atras
    // sino el FOREACH de arriba palma
    bc = tmp;

  } else {
    // TODO: radiation
    wasora_push_error_message("unknown boundary condition type '%s'", name);
    return WASORA_PARSER_ERROR;
  }
  
  
  return WASORA_RUNTIME_OK;
}


int fino_thermal_step(void) {
  
  int j;
  double xi;
  function_t *ic = NULL;
  
  if (wasora_special_var(in_static)) {
    // check if we were given an initial temperature distribution
    if ((ic = wasora_get_function_ptr("T_0")) != NULL) {

      if (ic->n_arguments != fino.dimensions) {
        wasora_push_error_message("initial condition function T_0 ought to have %d arguments instead of %d", fino.dimensions, ic->n_arguments);
        return WASORA_RUNTIME_ERROR;
      }

      for (j = fino.first_node; j < fino.last_node; j++) {
        xi = wasora_evaluate_function(ic, fino.mesh->node[j].x);
        petsc_call(VecSetValue(fino.phi, fino.mesh->node[j].index_dof[0], xi, INSERT_VALUES));
      }

      petsc_call(VecAssemblyBegin(fino.phi));
      petsc_call(VecAssemblyEnd(fino.phi));
    
    }
  
    // now, if the problem is steady-state or we do not have an initial condition then we just solve it
    // we use T0 (if provided) plus Dirichlet BCs as an initial guess 
    if (wasora_special_var(end_time) == 0 || ic == NULL) {
    
      if (fino.math_type == math_type_linear) {
        wasora_call(fino_solve_petsc_linear());
      } else if (fino.math_type == math_type_nonlinear) {
        wasora_call(fino_solve_petsc_nonlinear());
      }
      
    }
  
    wasora_call(fino_phi_to_solution(fino.phi));
    
  } else {
    
    PetscInt ts_steps;

    if (fino.ksp != NULL) {
      petsc_call(KSPDestroy(&fino.ksp));
      fino.ksp = NULL;
    }
   
    if (fino.ts == NULL) {
      petsc_call(TSCreate(PETSC_COMM_WORLD, &fino.ts));
      petsc_call(TSSetProblemType(fino.ts, TS_NONLINEAR));
      
      // TODO: tener apuntadores a funcion en la estructura fino, en tiempo de parseo
      // hacerlos apuntar aca y meter esto en un file generico petsc_ts y no en thermal
      petsc_call(TSSetIFunction(fino.ts, NULL, fino_thermal_IFunctionHeat, NULL));
      
      // si nos dieron una condicion inicial entonces fino.K no existe
      // en paralelo esto falla porque fino.J tiene que estar ensamblada y toda la milonga
      if (wasora.nprocs > 1 && wasora_get_function_ptr("T_0") != NULL) {
        wasora_call(fino_build_bulk());
//        wasora_call(fino_set_essential_bcs());
      }  
      petsc_call(MatDuplicate(fino.K, MAT_DO_NOT_COPY_VALUES, &fino.J));
      petsc_call(TSSetIJacobian(fino.ts, fino.J, fino.J, fino_thermal_IJacobianHeat, NULL));
   
      petsc_call(TSSetTimeStep(fino.ts, wasora_var_value(wasora_special_var(dt))));
      if (fino.ts_type != NULL) {
        petsc_call(TSSetType(fino.ts, fino.ts_type));
      } else {
        petsc_call(TSSetType(fino.ts, TSBDF));
      }
   
      // TODO: choose
      petsc_call(TSSetMaxStepRejections(fino.ts, 10000));
      
      petsc_call(TSSetExactFinalTime(fino.ts, TS_EXACTFINALTIME_STEPOVER));
      petsc_call(TSSetFromOptions(fino.ts));    
    }
   
    petsc_call(TSGetStepNumber(fino.ts, &ts_steps));
    petsc_call(TSSetMaxSteps(fino.ts, ts_steps+1));
   
    petsc_call(TSSolve(fino.ts, fino.phi));
    petsc_call(fino_phi_to_solution(fino.phi));
   
    petsc_call(TSGetStepNumber(fino.ts, &ts_steps));
    petsc_call(TSGetTimeStep(fino.ts, wasora_value_ptr(wasora_special_var(dt))));
  }

    
  return WASORA_RUNTIME_OK;
  
}


PetscErrorCode fino_thermal_IFunctionHeat(TS ts, PetscReal t, Vec T, Vec T_dot, Vec r, void *ctx) {
  
  // resolvemos 
  //   K*T + M*T_dot - b = 0

  Vec r_tran;

  wasora_var_value(wasora_special_var(t)) = t;
  
  // TODO: ver como hacer esto mas eficiente
  // esto no camina porque el transient seguro esta en b
  // pero si b tiene conveccion ya hay que recalcular K
/*  
  if (fino.K_nobc == NULL) {
    petsc_call(MatZeroEntries(fino.K));
    petsc_call(MatZeroEntries(fino.M));
    petsc_call(VecZeroEntries(fino.b));
    wasora_call(fino_build_bulk());
  } else {
    wasora_call(MatCopy(fino.K_nobc, fino.K, SAME_NONZERO_PATTERN));
    wasora_call(VecCopy(fino.b_nobc, fino.b));
  }
*/
  
  wasora_call(fino_build_bulk());
//  wasora_call(fino_set_essential_bc());
    
  
/*  
  printf("t = %g\n", t);
  printf("T\n");
  fino_print_petsc_vector(T, PETSC_VIEWER_STDOUT_SELF);
  printf("T_dot\n");
  fino_print_petsc_vector(T_dot, PETSC_VIEWER_STDOUT_SELF);
  printf("K\n");
  fino_print_petsc_matrix(fino.K, PETSC_VIEWER_STDOUT_SELF);
  printf("M\n");
  fino_print_petsc_matrix(fino.M, PETSC_VIEWER_STDOUT_SELF);
  printf("b\n");
  fino_print_petsc_vector(fino.b, PETSC_VIEWER_STDOUT_SELF);
*/  
  // armamos el residuo
  petsc_call(MatMult(fino.K, T, r));
/*  
  printf("KT\n");
  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
*/
  petsc_call(VecDuplicate(r, &r_tran));
//  petsc_call(VecZeroEntries(r_tran));
  petsc_call(MatMult(fino.M, T_dot, r_tran));
/*  
  printf("MT_dot\n");
  fino_print_petsc_vector(r_tran, PETSC_VIEWER_STDOUT_SELF);
*/
  
  petsc_call(VecAXPY(r, +1.0, r_tran));
  petsc_call(VecAXPY(r, -1.0, fino.b));
/*  
  printf("r\n");
  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
*/  
  VecDestroy(&r_tran);
  
  return WASORA_RUNTIME_OK;
}

PetscErrorCode fino_thermal_IJacobianHeat(TS ts, PetscReal t, Vec T, Vec T_dot, PetscReal s, Mat A, Mat B,void *ctx) {

  petsc_call(MatCopy(fino.K, A, SUBSET_NONZERO_PATTERN));
  petsc_call(MatAXPY(A, s, fino.M, SAME_NONZERO_PATTERN));
  petsc_call(MatCopy(A, B, SAME_NONZERO_PATTERN));
  
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
  
  mesh_compute_integration_weight_at_gauss(element, v);
  mesh_compute_H_at_gauss(element, v, fino.degrees);
  mesh_compute_B_at_gauss(element, v, fino.degrees);
  mesh_compute_x_at_gauss(element, v);
  r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);

  if (distribution_Q.defined != 0) {
    // el vector de fuente de calor volumetrica
    for (j = 0; j < element->type->nodes; j++) {
      gsl_vector_add_to_element(fino.bi, j,
        element->w[v] * r_for_axisymmetric * element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * fino_distribution_evaluate(&distribution_Q, material, element->x[v]));
    }
  }

  // calculamos la matriz de stiffness
  k = fino_distribution_evaluate(&distribution_k, material, element->x[v]);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * k, element->B[v], element->B[v], 1.0, fino.Ki);

  if (fino.M != NULL) {
    // calculamos la matriz de masa Ht*rho*cp*H
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
  
  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    
    mesh_compute_integration_weight_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    mesh_compute_x_at_gauss(element, v);
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

  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    
    mesh_compute_integration_weight_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    mesh_compute_x_at_gauss(element, v);
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
  
//  material_t *material = NULL;
//  element_list_item_t *associated_element = NULL;
  
  double k;
  double T_max = -INFTY;
  double T_min = +INFTY;
  
  int j;
  
  // evaluamos k, si es uniformes esto ya nos sirve para siempre
  if (distribution_k.variable != NULL) {
    k = fino_distribution_evaluate(&distribution_k, NULL, NULL);
    if (k < 0) {
      wasora_push_error_message("k is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
 
  for (j = 0; j < fino.mesh->n_nodes; j++) {
/*    
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
      if (material != NULL) {
//        wasora_push_error_message("cannot find a material for node %d", fino.mesh->node[j].tag);
//        return WASORA_RUNTIME_ERROR;
      
        if (distribution_k.variable == NULL) {
          k = fino_distribution_evaluate(&distribution_k, material, fino.mesh->node[j].x);
        }  
        if (k < 0) {
          wasora_push_error_message("k is negative");
          return WASORA_RUNTIME_ERROR;
        }  
      }      
    }
*/    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice T_max
    if (fino.solution[0]->data_value[j] >= T_max) {
      T_max = fino.mesh->node[j].phi[0];
    }
    if (fino.solution[0]->data_value[j] <= T_min) {
      T_min = fino.mesh->node[j].phi[0];
    }
    
  }

  wasora_var(fino.vars.T_max) = T_max;
  wasora_var(fino.vars.T_min) = T_min;

  
  return WASORA_RUNTIME_OK;
}

