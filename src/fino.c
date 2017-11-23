/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino plugin for wasora
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

#include <sys/time.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "petscsys.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include "wasora.h"
#include "fino.h"

#define time_checkpoint(which) \
  petsc_call(PetscTime(&wall.which)); \
  petsc_call(PetscGetCPUTime(&petsc.which)); \
  cpu.which = fino_get_cpu_time();


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_step"
int fino_instruction_step(void *arg) {
  fino_step_t *fino_step = (fino_step_t *)arg;
  fino_times_t wall;
  fino_times_t cpu;
  fino_times_t petsc;
  int k, g;

  PetscFunctionBegin;
  
  //---------------------------------
  // inicializamos si hace falta
  // TODO: ver si cambia la malla  
  //---------------------------------
  if (fino.spatial_unknowns == 0) {
    wasora_call(fino_problem_init());
  }
  

  // ------------------------------------
  // build
  // ------------------------------------
  if (fino_step->do_not_build == 0) {
    time_checkpoint(build_begin);
    wasora_call(fino_build_bulk());           // ensamblamos objetos elementales
    wasora_call(fino_set_essential_bc());     // condiciones de contorno esenciales
    time_checkpoint(build_end);
  }

  // ------------------------------------
  // solve
  // ------------------------------------
   if (fino_step->do_not_solve == 0) {
    time_checkpoint(solve_begin);
    
    if (wasora_var_value(wasora_special_var(in_static)) || fino.problem_family != problem_family_bake) {
      // resolvemos un steady state
      if (fino.math_type == math_linear) {
        wasora_call(fino_solve_linear_petsc(fino.K, fino.b));
      } else if (fino.math_type == math_eigen) {
#ifdef HAVE_SLEPC
        wasora_call(fino_solve_eigen_slepc());
        wasora_var(fino.vars.lambda) = fino.lambda;        // leemos el autovalor
#else
        wasora_push_error_message("fino should be linked against SLEPc to be able to solve eigen-problems");
        return WASORA_RUNTIME_ERROR;
#endif      
      }
    } else {
   
      // resolvemos 
      //   A*T(n+1) = b
      // con
      //   A = K*theta + C/dt 
      //   b = (-K*(1-theta)+C/dt)*T(n) + b
      double theta = 0.5;

      if (fino.has_transient == 0) {
        fino.has_transient = 1;
        MatDuplicate(fino.K, MAT_COPY_VALUES, &fino.A);
        MatDuplicate(fino.K, MAT_COPY_VALUES, &fino.B);
      } else {
        MatCopy(fino.K, fino.A, SAME_NONZERO_PATTERN);
        MatCopy(fino.K, fino.B, SAME_NONZERO_PATTERN);
      }
      
      MatScale(fino.A, theta);
      MatAXPY(fino.A, 1/wasora_var_value(wasora_special_var(dt)), fino.M, SUBSET_NONZERO_PATTERN);

      MatScale(fino.B, -(1-theta));
      MatAXPY(fino.B, 1/wasora_var_value(wasora_special_var(dt)), fino.M, DIFFERENT_NONZERO_PATTERN);

//      MatMult(fino.B, fino.phi, fino.b);
//      VecAXPY(problem.b, 1, problem.R);

//      wasora_call(fino_set_essential_bc());     // condiciones de contorno esenciales

      wasora_call(fino_solve_linear_petsc(fino.K, fino.b));
    }

    
    // fabricamos G funciones con la solucion
    for (k = 0; k < fino.spatial_unknowns; k++) {
      for (g = 0; g < fino.degrees; g++) {
        petsc_call(VecGetValues(fino.phi, 1, &fino.mesh->node[k].index[g], &fino.solution[g]->data_value[k]));
        // si tenemos una solucion base hay que sumarla
        if (fino.base_solution != NULL && fino.base_solution[g] != NULL) {
          // cuales son las chances de que estas sean iguales y no esten sobre la misma malla?
          if (fino.base_solution[g]->data_size == fino.spatial_unknowns) {
            fino.solution[g]->data_value[k] += fino.base_solution[g]->data_value[k];
          } else {
            fino.solution[g]->data_value[k] += wasora_evaluate_function(fino.base_solution[g], fino.mesh->node[k].x);
          }
        }
      }
    }
    
    if (fino_step->do_not_compute_gradients == 0) {
      wasora_call(fino_compute_gradients());
      if (fino.problem_family == problem_family_break) {
        wasora_call(fino_break_compute_reactions());
        wasora_call(fino_break_compute_stresses());
      } else if (fino.problem_family == problem_family_bake) {
        wasora_call(fino_bake_compute_fluxes());
      }
    }
    time_checkpoint(solve_end);
  }
  
  if (fino_step->do_not_build == 0) {
    wasora_var(fino.vars.time_petsc_build) = petsc.build_end - petsc.build_begin;
    wasora_var(fino.vars.time_wall_build)  = wall.build_end  - wall.build_begin;
    wasora_var(fino.vars.time_cpu_build)   = cpu.build_end   - cpu.build_begin;
  }
  
  if (fino_step->do_not_solve == 0) {
    wasora_var(fino.vars.time_petsc_solve) = petsc.solve_end - petsc.solve_begin;
    wasora_var(fino.vars.time_wall_solve)  = wall.solve_end  - wall.solve_begin;
    wasora_var(fino.vars.time_cpu_solve)   = cpu.solve_end   - cpu.solve_begin;
  }
  
  wasora_var(fino.vars.time_petsc_total) = wasora_var(fino.vars.time_petsc_build) + wasora_var(fino.vars.time_petsc_solve);
  wasora_var(fino.vars.time_wall_total)  = wasora_var(fino.vars.time_wall_build)  + wasora_var(fino.vars.time_wall_solve);
  wasora_var(fino.vars.time_cpu_total)   = wasora_var(fino.vars.time_cpu_build)   + wasora_var(fino.vars.time_cpu_solve);

  getrusage(RUSAGE_SELF, &fino.resource_usage);
  wasora_value(fino.vars.memory_usage_global) = (double)(1024.0*fino.resource_usage.ru_maxrss);
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}





int fino_assembly(void) {
  MatAssemblyBegin(fino.K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(fino.K, MAT_FINAL_ASSEMBLY);
  if (fino.has_mass) {
    MatAssemblyBegin(fino.M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fino.M, MAT_FINAL_ASSEMBLY);
  }
  if (fino.has_transient) {
    MatAssemblyBegin(fino.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fino.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fino.B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fino.B, MAT_FINAL_ASSEMBLY);
  }
  if (fino.has_rhs) {
    VecAssemblyBegin(fino.b);
    VecAssemblyEnd(fino.b);
  }
  
  return WASORA_RUNTIME_OK;
}
