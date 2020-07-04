/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino is a free finite-element thermo-mechanical solver
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

#include <sys/time.h>

#ifndef _FINO_H
#include "fino.h"
#endif


int fino_instruction_step(void *arg) {
//  fino_step_t *fino_step = (fino_step_t *)arg;

  int j;
  double xi;
  function_t *ic = NULL;


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_step"
int fino_instruction_step(void *arg) {
  fino_step_t *fino_step = (fino_step_t *)arg;
  fino_times_t wall;
  fino_times_t cpu;
  fino_times_t petsc;
  
  PetscFunctionBegin;
  
  //---------------------------------
  // initialize only if we did not initialized before
  // TODO: how to handle changes in the mesh within steps?
  //---------------------------------
  if (fino.spatial_unknowns == 0) {
    wasora_call(fino_problem_init());
  }

  if (wasora_var_value(wasora_special_var(in_static)) || fino.transient_type == transient_type_quasistatic) {
    // check if we were given an initial temperature distribution
    // TODO: generalize to other problems
    if ((ic = wasora_get_function_ptr("T_0")) != NULL) {

      if (ic->n_arguments != fino.dimensions) {
        wasora_push_error_message("initial condition function T_0 ought to have %d arguments instead of %d", fino.dimensions, ic->n_arguments);
        return WASORA_RUNTIME_ERROR;
      }

      for (j = fino.first_node; j < fino.last_node; j++) {
        xi = wasora_evaluate_function(ic, fino.mesh->node[j].x);
        petsc_call(VecSetValue(fino.phi, fino.mesh->node[j].index_dof[0], xi, INSERT_VALUES));
      }

      // TODO: assembly routine
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
      } else if (fino.math_type == math_type_eigen) {
#ifdef HAVE_SLEPC
        wasora_call(fino_solve_eigen_slepc());
        wasora_call(fino_eigen_nev()); 
#else 
        wasora_push_error_message("Fino should be linked against SLEPc to be able to solve eigen-problems");
        return WASORA_RUNTIME_ERROR;
#endif      
      }
    }
  
    wasora_call(fino_phi_to_solution(fino.phi, 1));
    
  } else {
    
    PetscInt ts_steps;

    if (fino.ksp != NULL) {
      petsc_call(KSPDestroy(&fino.ksp));
      fino.ksp = NULL;
    } else if (fino.snes != NULL) {
      petsc_call(SNESDestroy(&fino.snes));
    }
   
     if (fino.ts == NULL) {
      petsc_call(TSCreate(PETSC_COMM_WORLD, &fino.ts));
      petsc_call(TSSetProblemType(fino.ts, TS_NONLINEAR));
      
      petsc_call(TSSetIFunction(fino.ts, NULL, fino_ts_residual, NULL));
      
      // si nos dieron una condicion inicial entonces fino.K no existe
      // en paralelo esto falla porque fino.J tiene que estar ensamblada y toda la milonga
      if (wasora_get_function_ptr("T_0") != NULL) {      
        wasora_call(fino_build_bulk());
      }  
      petsc_call(MatDuplicate(fino.K, MAT_DO_NOT_COPY_VALUES, &fino.J));
      petsc_call(TSSetIJacobian(fino.ts, fino.J, fino.J, fino_ts_jacobian, NULL));
   
      petsc_call(TSSetTimeStep(fino.ts, wasora_var_value(wasora_special_var(dt))));
      
      // if BCs depend on time we need DAEs
      petsc_call(TSSetEquationType(fino.ts, TS_EQ_IMPLICIT));      
//      petsc_call(TSSetEquationType(fino.ts, TS_EQ_DAE_IMPLICIT_INDEX1));
//      petsc_call(TSARKIMEXSetFullyImplicit(fino.ts, PETSC_TRUE)); 
      
      // TODO: the default depends on the problem type
      if (fino.ts_type != NULL) {
        petsc_call(TSSetType(fino.ts, fino.ts_type));
      } else {
        petsc_call(TSSetType(fino.ts, TSBDF));
//        petsc_call(TSSetType(fino.ts, TSARKIMEX));
      }
   
      // TODO: choose
      petsc_call(TSSetMaxStepRejections(fino.ts, 10000));
      petsc_call(TSSetMaxSNESFailures(fino.ts, 1000));
      petsc_call(TSSetExactFinalTime(fino.ts, TS_EXACTFINALTIME_STEPOVER));
      
      // options overwrite
      petsc_call(TSSetFromOptions(fino.ts));    
    }
   
    petsc_call(TSGetStepNumber(fino.ts, &ts_steps));
    petsc_call(TSSetMaxSteps(fino.ts, ts_steps+1));
   
    petsc_call(TSSolve(fino.ts, fino.phi));
    petsc_call(fino_phi_to_solution(fino.phi, 1));
   
    petsc_call(TSGetStepNumber(fino.ts, &ts_steps));
    petsc_call(TSGetTimeStep(fino.ts, wasora_value_ptr(wasora_special_var(dt))));
  }
  

  
/*  
  if (wasora_var_value(wasora_special_var(end_time)) == 0 || fino.problem_family != problem_family_thermal) {
      
      // if the problem is not transient heat we solve the (quasi) steady state here
      if (fino.math_type == math_type_linear) {
        
        wasora_call(fino_solve_petsc_linear());
        wasora_call(fino_phi_to_solution(fino.phi));
        
      } else if (fino.math_type == math_type_nonlinear) {
        
        wasora_call(fino_solve_petsc_nonlinear());
        
      } else if (fino.math_type == math_type_eigen) {
        
#ifdef HAVE_SLEPC
        wasora_call(fino_solve_eigen_slepc());
        wasora_call(fino_eigen_nev()); 
        wasora_call(fino_phi_to_solution(fino.phi));
#else 
        wasora_push_error_message("fino should be linked against SLEPc to be able to solve eigen-problems");
        return WASORA_RUNTIME_ERROR;
#endif      
      }
    } else {
      
      // the transient heat problem is different
      // TODO: homogenize the logic
      if (wasora_var_value(wasora_special_var(in_static))) {
        wasora_call(fino_thermal_step_initial());
      } else {
        wasora_call(fino_thermal_step_transient());
      }
    }
    time_checkpoint(solve_end);
    
  // ------------------------------------
  // stresses
  // ------------------------------------
    
    time_checkpoint(stress_begin);
    
    if (fino.problem_family == problem_family_mechanical) {
  
      wasora_call(fino_compute_strain_energy());
      if (fino.gradient_evaluation != gradient_none) {
        wasora_call(fino_break_compute_stresses());
      }  
      
    } else if (fino.problem_family == problem_family_thermal) {
        
      wasora_call(fino_thermal_compute_fluxes());
      
    }
    
    time_checkpoint(stress_end);
    
  }
*/  
  wasora_var(fino.vars.time_petsc_build) += fino.petsc.build_end - fino.petsc.build_begin;
  wasora_var(fino.vars.time_wall_build)  += fino.wall.build_end  - fino.wall.build_begin;
  wasora_var(fino.vars.time_cpu_build)   += fino.cpu.build_end   - fino.cpu.build_begin;
  
  wasora_var(fino.vars.time_petsc_solve) += fino.petsc.solve_end - fino.petsc.solve_begin;
  wasora_var(fino.vars.time_wall_solve)  += fino.wall.solve_end  - fino.wall.solve_begin;
  wasora_var(fino.vars.time_cpu_solve)   += fino.cpu.solve_end   - fino.cpu.solve_begin;

  wasora_var(fino.vars.time_petsc_stress) += fino.petsc.stress_end - fino.petsc.stress_begin;
  wasora_var(fino.vars.time_wall_stress)  += fino.wall.stress_end  - fino.wall.stress_begin;
  wasora_var(fino.vars.time_cpu_stress)   += fino.cpu.stress_end   - fino.cpu.stress_begin;
  
  wasora_var(fino.vars.time_petsc_total) += wasora_var(fino.vars.time_petsc_build) + wasora_var(fino.vars.time_petsc_solve) + wasora_var(fino.vars.time_petsc_stress);
  wasora_var(fino.vars.time_wall_total)  += wasora_var(fino.vars.time_wall_build)  + wasora_var(fino.vars.time_wall_solve)  + wasora_var(fino.vars.time_wall_stress);
  wasora_var(fino.vars.time_cpu_total)   += wasora_var(fino.vars.time_cpu_build)   + wasora_var(fino.vars.time_cpu_solve)   + wasora_var(fino.vars.time_cpu_stress);

  getrusage(RUSAGE_SELF, &fino.resource_usage);
  wasora_value(fino.vars.memory) = (double)(1024.0*fino.resource_usage.ru_maxrss);

  
  return WASORA_RUNTIME_OK;
}





int fino_assembly(void) {
  // which is better?
/*  
  petsc_call(MatAssemblyBegin(fino.K, MAT_FINAL_ASSEMBLY));
  petsc_call(MatAssemblyEnd(fino.K, MAT_FINAL_ASSEMBLY));

  petsc_call(VecAssemblyBegin(fino.phi));
  petsc_call(VecAssemblyEnd(fino.phi));
  
  if (fino.M != NULL) {
    petsc_call(MatAssemblyBegin(fino.M, MAT_FINAL_ASSEMBLY));
    petsc_call(MatAssemblyEnd(fino.M, MAT_FINAL_ASSEMBLY));
  }
  if (fino.J != NULL) {
    petsc_call(VecAssemblyBegin(fino.b));
    petsc_call(VecAssemblyEnd(fino.b));
  }
*/
  
  petsc_call(VecAssemblyBegin(fino.phi));
  if (fino.b != NULL) {
    petsc_call(VecAssemblyBegin(fino.b));
  }  
  petsc_call(MatAssemblyBegin(fino.K, MAT_FINAL_ASSEMBLY));
  if (fino.M != NULL) {
    petsc_call(MatAssemblyBegin(fino.M, MAT_FINAL_ASSEMBLY));
  }  


  petsc_call(VecAssemblyEnd(fino.phi));
  if (fino.b != NULL) {
    petsc_call(VecAssemblyEnd(fino.b));
  }  
  petsc_call(MatAssemblyEnd(fino.K, MAT_FINAL_ASSEMBLY));
  if (fino.M != NULL) {
    petsc_call(MatAssemblyEnd(fino.M, MAT_FINAL_ASSEMBLY));
  }
  
  return WASORA_RUNTIME_OK;
}


int fino_phi_to_solution(Vec phi, int compute_gradients) {

  double xi;
  int i, j, g;

  Vec                phi_full;
  VecScatter         vscat;

  petsc_call(VecScatterCreateToZero(phi, &vscat, &phi_full));
  petsc_call(VecScatterBegin(vscat, phi, phi_full, INSERT_VALUES,SCATTER_FORWARD););
  petsc_call(VecScatterEnd(vscat, phi, phi_full, INSERT_VALUES,SCATTER_FORWARD););
    
  // make up G functions with the solution
  if (wasora.rank == 0) {
    for (j = 0; j < fino.spatial_unknowns; j++) {
      for (g = 0; g < fino.degrees; g++) {
        petsc_call(VecGetValues(phi_full, 1, &fino.mesh->node[j].index_dof[g], &fino.mesh->node[j].phi[g]));

        // if there is a base solution
        if (fino.base_solution != NULL && fino.base_solution[g] != NULL) {
          if (fino.base_solution[g]->mesh == fino.mesh) {
            fino.mesh->node[j].phi[g] += fino.base_solution[g]->data_value[j];
          } else {
            fino.mesh->node[j].phi[g] += wasora_evaluate_function(fino.base_solution[g], fino.mesh->node[j].x);
          }
        }

        // if we are not in rough mode we fill the solution here
        // because it is esay, in rough mode we need to iterate over the elements
        if (fino.rough == 0) {
          fino.solution[g]->data_value[j] = fino.mesh->node[j].phi[g];
        }

        if (fino.nev > 1) {
          for (i = 0; i < fino.nev; i++) {
            // the values already have the excitation factor
            petsc_call(VecGetValues(fino.eigenvector[i], 1, &fino.mesh->node[j].index_dof[g], &xi));
            fino.mode[g][i]->data_value[j] = xi;
            wasora_vector_set(fino.vectors.phi[i], j, xi);
          }
        }
      }
    }
  }
    
  petsc_call(VecDestroy(&phi_full));
  petsc_call(VecScatterDestroy(&vscat));
    
    
  if (fino.rough) {
    node_t *node;  
    // in rough mode we need to iterate over the elements first and then over the nodes
    // TODO: see if the order of the loops is the optimal one
    for (g = 0; g < fino.degrees; g++) {
      for (i = 0; i < fino.mesh_rough->n_elements; i++) {
        for (j = 0; j < fino.mesh_rough->element[i].type->nodes; j++) {
          node = fino.mesh_rough->element[i].node[j];  
          fino.solution[g]->data_value[node->index_mesh] = node->phi[g];  
        }
      }
    }
  }

  
  if (compute_gradients) {
    // TODO: function pointers
    if (fino.problem_family == problem_family_mechanical) {
  
      time_checkpoint(stress_begin);
      wasora_call(fino_compute_strain_energy());
      if (fino.gradient_evaluation != gradient_none) {
        wasora_call(fino_break_compute_stresses());
      }  
      time_checkpoint(stress_end);
      
    } else if (fino.problem_family == problem_family_thermal) {
       
      wasora_call(fino_thermal_compute_fluxes());
    }  
  }
    
  
  
  return WASORA_RUNTIME_OK;
}
