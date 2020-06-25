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

  printf("phi\n");
  fino_print_petsc_vector(phi, PETSC_VIEWER_STDOUT_SELF);
  
  // pasamos phi a la solucion porque K puede depender de phi
  wasora_call(fino_phi_to_solution(phi));
  wasora_call(fino_build_bulk());
  
  petsc_call(MatMult(fino.K, phi, r));
  petsc_call(VecAYPX(r, -1.0, fino.b));
  wasora_call(fino_set_essential_bcs(NULL, NULL, NULL, NULL, fino.phi, NULL));
  
  printf("residual\n");
  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
  
  return 0;
}

PetscErrorCode fino_solve_jacobian(SNES snes,Vec phi, Mat J, Mat P, void *ctx) {
  
//  printf("jacobiano\n");
  petsc_call(MatCopy(fino.K, J, SAME_NONZERO_PATTERN));
  wasora_call(fino_set_essential_bcs(NULL, NULL, J, NULL, NULL, NULL));
  petsc_call(MatCopy(J, P, SAME_NONZERO_PATTERN));
  //MatDuplicate(fino.K, MAT_COPY_VALUES, &J);
  
  return 0;
   
}

#undef  __func__
#define __func__ "fino_solve_petsc_linear"
int fino_solve_petsc_nonlinear(void) {

  SNESConvergedReason reason;
  Mat            J;    // jacobiano
  Vec            r;    // residuo
  PetscInt       its;    

  PetscFunctionBegin;
      
  time_checkpoint(build_begin);
  wasora_call(fino_build_bulk());
  time_checkpoint(build_end);
  
  if (fino.snes == NULL) {
    petsc_call(SNESCreate(PETSC_COMM_WORLD, &fino.snes));
  }
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
  
  wasora_call(fino_set_ksp());
  wasora_call(fino_set_pc());
  
  petsc_call(SNESSetFromOptions(fino.snes));

  // monitor
  petsc_call(SNESMonitorSet(fino.snes, fino_snes_monitor, NULL, 0));
  
  // initial guess
  wasora_call(fino_set_essential_bcs(NULL, NULL, NULL, NULL, fino.phi, NULL));
  
  // solve
  petsc_call(SNESSolve(fino.snes, NULL, fino.phi));
  
  // chequeamos que haya convergido
  petsc_call(SNESGetConvergedReason(fino.snes, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's non-linear solver did not converge with reason '%s' (%d)", SNESConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }
  
  petsc_call(SNESGetIterationNumber(fino.snes, &its));
  wasora_value(fino.vars.iterations) = (double)its;
  
//  petsc_call(SNESGetResidualNorm(fino.snes, wasora_value_ptr(fino.vars.residual_norm)));
  
  
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

//  printf("%d %e %e\n", n, rnorm/fino.progress_r0, current_progress);
  
  if (fino.progress_ascii == PETSC_TRUE) {
    for (i = (int)(100*fino.progress_last); i < (int)(100*current_progress); i++) {
      printf(CHAR_PROGRESS_SOLVE);
      fflush(stdout);
    }
    fino.progress_last = current_progress;
  }

  return 0;
}
