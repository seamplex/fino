/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's non-linear solver using PETSc routines
 *
 *  Copyright (C) 2020 jeremy theler
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


#include <petsc.h>
#include <petscsys.h>
#include <petscsnes.h>

#include "fino.h"


PetscErrorCode fino_solve_residual(SNES snes, Vec phi, Vec r,void *ctx) {
  
  // pasamos phi a la solucion porque K puede depender de phi
  wasora_call(fino_assembly());
  wasora_call(fino_phi_to_solution(phi));
  if (fino.problem_family == problem_family_mechanical) {
    wasora_call(fino_break_compute_stresses());
  }  
  wasora_call(fino_build_bulk());
  wasora_call(fino_set_essential_bc(fino.K, fino.b));
  wasora_call(fino_assembly());
  
  MatMult(fino.K, phi, r);
  VecCopy(r, fino.r);
  printf("residual\n");
//  fino_print_petsc_matrix(fino.K, PETSC_VIEWER_STDOUT_SELF);
//  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
  

//  printf("phi\n");
//  fino_print_petsc_vector(phi, PETSC_VIEWER_STDOUT_SELF);
  
  return 0;
}

PetscErrorCode fino_solve_jacobian(SNES snes,Vec phi, Mat J, Mat P, void *ctx) {
  
  // pasamos phi a la solucion porque K puede depender de phi
  double xi;
  double K, dKdT;
  double ri, Tj;
  int i, j;
/*  
  wasora_call(fino_assembly());
  wasora_call(fino_phi_to_solution(phi));
  if (fino.problem_family == problem_family_mechanical) {
    wasora_call(fino_break_compute_stresses());
  }  
  wasora_call(fino_build_bulk());
  wasora_call(fino_set_essential_bc(fino.K, fino.b));
  wasora_call(fino_assembly());
*/  
  MatDuplicate(fino.K, MAT_COPY_VALUES, &J);
  MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


  for (i = 0; i < fino.mesh->n_nodes; i++) {
    for (j = 0; j < fino.mesh->n_nodes; j++) {
      VecGetValues(fino.r, 1, &i, &ri);
      Tj = fino.mesh->node[j].phi[0];
      K = 1+Tj;
      dKdT = 1;
      xi = 1/K * dKdT * ri;
      if (xi != 0) {
        MatSetValues(J, 1, &i, 1, &j, &xi, ADD_VALUES);
      }  
    }
  }

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
      
  printf("jacobiano\n");
  fino_print_petsc_matrix(J, PETSC_VIEWER_STDOUT_SELF);
  
  return 0;
   
}

#undef  __FUNCT__
#define __FUNCT__ "fino_solve_nonlinear_petsc"
int fino_solve_nonlinear_petsc(void) {

  SNESConvergedReason reason;
  Mat            J;    // jacobiano
  Vec            r;    // residuo
  PetscInt       its;    

  if (fino.snes == NULL) {
    petsc_call(SNESCreate(PETSC_COMM_WORLD, &fino.snes));
  }
  petsc_call(VecDuplicate(fino.phi, &fino.r));
  petsc_call(VecDuplicate(fino.phi, &r));
  petsc_call(MatDuplicate(fino.K, MAT_COPY_VALUES, &J));
  
  petsc_call(SNESSetFunction(fino.snes, r, fino_solve_residual, NULL));
  petsc_call(SNESSetJacobian(fino.snes, J, J, fino_solve_jacobian, NULL));
  petsc_call(SNESSetTolerances(fino.snes, wasora_var(fino.vars.abstol),
                                          wasora_var(fino.vars.reltol),
                                          PETSC_DEFAULT,
                                          (PetscInt)wasora_var(fino.vars.max_iterations),
                                          PETSC_DEFAULT));
  // el SNES 
  if (fino.snes_type != NULL) {
    // si nos dieron lo ponemos, sino dejamos el dafault
    petsc_call(SNESSetType(fino.snes, fino.snes_type));
  }
  
  // el KSP
  petsc_call(SNESGetKSP(fino.snes, &fino.ksp));
  if (fino.ksp_type != NULL) {
    // si nos dieron lo ponemos, sino dejamos el default
    petsc_call(KSPSetType(fino.ksp, fino.ksp_type));
  }
  
  wasora_call(fino_ksp_set_pc(fino.K));
  
  petsc_call(SNESSetFromOptions(fino.snes));

  // initial value
  petsc_call(VecSet(fino.phi, 0));

  // monitor
  petsc_call(SNESMonitorSet(fino.snes, fino_snes_monitor, NULL, 0));
  
  // solve
  petsc_call(SNESSolve(fino.snes, fino.b, fino.phi));
  
  // chequeamos que haya convergido
  petsc_call(SNESGetConvergedReason(fino.snes, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's non-linear solver did not converge with reason '%s' (%d)", SNESConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }
  
  petsc_call(SNESGetIterationNumber(fino.snes, &its));  
  wasora_value(fino.vars.iterations) = (double)its;
  
//  printf("residuo final\n");
//  MatMult(fino.K, fino.phi, r);
//  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
//  petsc_call(SNESGetResidualNorm(snes, wasora_value_ptr(fino.vars.residual_norm)));
  
  
  return WASORA_RUNTIME_OK;

}

PetscErrorCode fino_snes_monitor(SNES snes, PetscInt n, PetscReal rnorm, void *dummy) {
  int i;
  double current_progress;
  
  if (fino.progress_r0 == 0) {
    fino.progress_r0 = rnorm;
  }
  
  if (rnorm < 1e-20) {
    current_progress = 1;
  } else {
    current_progress = log((rnorm/fino.progress_r0))/log(wasora_var_value(fino.vars.reltol));
    if (current_progress > 1) {
      current_progress = 1;
    }
  } 

  printf("%d %e %.0f\n", n, rnorm/fino.progress_r0, 100*current_progress);
  
  if (fino.progress_ascii) {
    for (i = (int)(100*fino.progress_last); i < (int)(100*current_progress); i++) {
      printf(CHAR_PROGRESS_SOLVE);
      fflush(stdout);
    }
    fino.progress_last = current_progress;
  }

  return 0;
}
