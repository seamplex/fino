/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino plugin for wasora
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
  
  PetscFunctionBegin;
  
  //---------------------------------
  // initialize only if we did not initialized before
  // TODO: how to handle changes in the mesh within steps?
  //---------------------------------
  if (fino.spatial_unknowns == 0) {
    wasora_call(fino_problem_init());
  }
  

  // ------------------------------------
  // build
  // ------------------------------------
  // this logic ought to change to have a generalized way of building the stiffness matrix
  // independently of the problem type
  if (wasora_var_value(wasora_special_var(end_time)) == 0 || fino.problem_family != problem_family_thermal) {
    time_checkpoint(build_begin);
    wasora_call(fino_build_bulk());
    wasora_call(fino_set_essential_bc());
    time_checkpoint(build_end);
  }

  // ------------------------------------
  // solve
  // ------------------------------------
  // the do_not_solve flag can be used to debug stuff when the problem does not converge
  if (fino_step->do_not_solve == 0) {
    
    time_checkpoint(solve_begin);
    if (wasora_var_value(wasora_special_var(end_time)) == 0 || fino.problem_family != problem_family_thermal) {
      
      // if the problem is not transient heat we solve the (quasi) steady state here
      if (fino.math_type == math_type_linear) {
        
        wasora_call(fino_solve_petsc_linear());
        wasora_call(fino_phi_to_solution(fino.phi));
        
      } else if (fino.math_type == math_type_nonlinear) {
        
        wasora_call(fino_solve_nonlinear_petsc());
        
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
  
  wasora_var(fino.vars.time_petsc_build) = petsc.build_end - petsc.build_begin;
  wasora_var(fino.vars.time_wall_build)  = wall.build_end  - wall.build_begin;
  wasora_var(fino.vars.time_cpu_build)   = cpu.build_end   - cpu.build_begin;
  
  if (fino_step->do_not_solve == 0) {
    wasora_var(fino.vars.time_petsc_solve) = petsc.solve_end - petsc.solve_begin;
    wasora_var(fino.vars.time_wall_solve)  = wall.solve_end  - wall.solve_begin;
    wasora_var(fino.vars.time_cpu_solve)   = cpu.solve_end   - cpu.solve_begin;

    wasora_var(fino.vars.time_petsc_stress) = petsc.stress_end - petsc.stress_begin;
    wasora_var(fino.vars.time_wall_stress)  = wall.stress_end  - wall.stress_begin;
    wasora_var(fino.vars.time_cpu_stress)   = cpu.stress_end   - cpu.stress_begin;
  }
  
  wasora_var(fino.vars.time_petsc_total) = wasora_var(fino.vars.time_petsc_build) + wasora_var(fino.vars.time_petsc_solve) + wasora_var(fino.vars.time_petsc_stress);
  wasora_var(fino.vars.time_wall_total)  = wasora_var(fino.vars.time_wall_build)  + wasora_var(fino.vars.time_wall_solve)  + wasora_var(fino.vars.time_wall_stress);
  wasora_var(fino.vars.time_cpu_total)   = wasora_var(fino.vars.time_cpu_build)   + wasora_var(fino.vars.time_cpu_solve)   + wasora_var(fino.vars.time_cpu_stress);

  getrusage(RUSAGE_SELF, &fino.resource_usage);
  wasora_value(fino.vars.memory) = (double)(1024.0*fino.resource_usage.ru_maxrss);
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}





int fino_assembly(void) {
  petsc_call(MatAssemblyBegin(fino.K, MAT_FINAL_ASSEMBLY));
  petsc_call(MatAssemblyEnd(fino.K, MAT_FINAL_ASSEMBLY));
  if (fino.has_mass) {
    petsc_call(MatAssemblyBegin(fino.M, MAT_FINAL_ASSEMBLY));
    petsc_call(MatAssemblyEnd(fino.M, MAT_FINAL_ASSEMBLY));
  }
  if (fino.has_rhs) {
    petsc_call(VecAssemblyBegin(fino.b));
    petsc_call(VecAssemblyEnd(fino.b));
  }
  
  return WASORA_RUNTIME_OK;
}


int fino_phi_to_solution(Vec phi) {

  double xi;
  int i, j, g;

  Vec                phi0;
  VecScatter         vscat;

  petsc_call(VecScatterCreateToZero(phi, &vscat, &phi0));
  petsc_call(VecScatterBegin(vscat, phi, phi0, INSERT_VALUES,SCATTER_FORWARD););
  petsc_call(VecScatterEnd(vscat, phi, phi0, INSERT_VALUES,SCATTER_FORWARD););
    
  // fabricamos G funciones con la solucion
  if (wasora.rank == 0) {
    for (j = 0; j < fino.spatial_unknowns; j++) {
      for (g = 0; g < fino.degrees; g++) {
        petsc_call(VecGetValues(phi0, 1, &fino.mesh->node[j].index_dof[g], &fino.mesh->node[j].phi[g]));

        // si tenemos una solucion la sumamos 
        if (fino.base_solution != NULL && fino.base_solution[g] != NULL) {
          if (fino.base_solution[g]->mesh == fino.mesh) {
            fino.mesh->node[j].phi[g] += fino.base_solution[g]->data_value[j];
          } else {
            fino.mesh->node[j].phi[g] += wasora_evaluate_function(fino.base_solution[g], fino.mesh->node[j].x);
          }
        }

        // si no estamos en rough rellenamos la solucion de los desplazamietos
        // porque es facil, en rough hay que iterar sobre los elementos
        if (fino.rough == 0) {
          fino.solution[g]->data_value[j] = fino.mesh->node[j].phi[g];
        }

        if (fino.nev > 1) {
          for (i = 0; i < fino.nev; i++) {
            // las funciones ya vienen con el factor de excitacion
            petsc_call(VecGetValues(fino.eigenvector[i], 1, &fino.mesh->node[j].index_dof[g], &xi));
            fino.mode[g][i]->data_value[j] = xi;
            wasora_vector_set(fino.vectors.phi[i], j, xi);
          }
        }
      }
    }
  }
    
  petsc_call(VecDestroy(&phi0));
  petsc_call(VecScatterDestroy(&vscat));
    
    
  if (fino.rough) {
    node_t *node;  
    // si estamos en rough rellenamos los desplazamientos iterando sobre elements
    // y despues sobre cada nodo
    // TODO: ver orden de fors! capaz que convenga que el g este afuera  
    for (g = 0; g < fino.degrees; g++) {
      for (i = 0; i < fino.mesh_rough->n_elements; i++) {
        for (j = 0; j < fino.mesh_rough->element[i].type->nodes; j++) {
          node = fino.mesh_rough->element[i].node[j];  
          fino.solution[g]->data_value[node->index_mesh] = node->phi[g];  
        }
      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}
