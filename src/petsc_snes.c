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
#ifndef _FINO_H
#include "fino.h"
#endif


PetscErrorCode fino_snes_residual(SNES snes, Vec phi, Vec r,void *ctx) {


  // this check is only to avoid building the first time because we already
  // built K in order to create the SNES, but the rest of the step we need
  // to re-build always because we are sure the problem is non-linear
  if (fino.already_built == PETSC_FALSE) {
    // pass phi to the solution becuase K (and the BCs) might depend on phi
    wasora_call(fino_phi_to_solution(phi));
    wasora_call(fino_build_bulk());
  }  
  
  petsc_call(MatMult(fino.K, phi, r));
  petsc_call(VecAXPY(r, -1.0, fino.b));
  wasora_call(fino_dirichlet_set_r(r, phi));
  
  fino.already_built = PETSC_FALSE;
  
  return 0;
}

PetscErrorCode fino_snes_jacobian(SNES snes,Vec phi, Mat J, Mat P, void *ctx) {
  
  petsc_call(MatCopy(fino.K, J, SAME_NONZERO_PATTERN));
  wasora_call(fino_dirichlet_set_J(J));
  petsc_call(MatCopy(J, P, SAME_NONZERO_PATTERN));
  
  return 0;
}

int fino_solve_petsc_nonlinear(void) {

  KSP ksp;
  PC pc;
  Mat J;
  Vec r;
  SNESConvergedReason reason;
  PetscInt       its;
  
  time_checkpoint(build_begin);
  
  if (fino.snes == NULL) {
    petsc_call(SNESCreate(PETSC_COMM_WORLD, &fino.snes));

    if (fino.snes_type != NULL) {
      // if we have an explicit type, we set it
      petsc_call(SNESSetType(fino.snes, fino.snes_type));
    }
    
    // we build the matrices here and put a flag that it is already built
    // so we do not build it again in the first step of the SNES
    wasora_call(fino_build_bulk());
    wasora_call(fino_dirichlet_eval(fino.K, fino.b));
    fino.already_built = PETSC_TRUE;
    time_checkpoint(build_end);
    
    // monitor
    petsc_call(SNESMonitorSet(fino.snes, fino_snes_monitor, NULL, 0));

    // options
    petsc_call(SNESSetFromOptions(fino.snes));
  }  
    
  petsc_call(VecDuplicate(fino.phi, &r));
  petsc_call(MatDuplicate(fino.K, MAT_COPY_VALUES, &J));
  
  petsc_call(SNESSetFunction(fino.snes, r, fino_snes_residual, NULL));
  petsc_call(SNESSetJacobian(fino.snes, J, J, fino_snes_jacobian, NULL));
  
  petsc_call(SNESSetTolerances(fino.snes, wasora_var(fino.vars.abstol),
                                          wasora_var(fino.vars.reltol),
                                          PETSC_DEFAULT,
                                          (PetscInt)wasora_var(fino.vars.max_iterations),
                                          PETSC_DEFAULT));

  // customize ksp and pc (this needs to come after setting the jacobian)
  petsc_call(SNESGetKSP(fino.snes, &ksp));
  wasora_call(fino_set_ksp(ksp));

  petsc_call(KSPGetPC(ksp, &pc));
  wasora_call(fino_set_pc(pc));
  
  // initial guess
  wasora_call(fino_dirichlet_set_phi(fino.phi));
  
  // solve
  petsc_call(SNESSolve(fino.snes, NULL, fino.phi));
  
  // check convergence
  petsc_call(SNESGetConvergedReason(fino.snes, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's non-linear solver did not converge with reason '%s' (%d)", SNESConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }
  
  petsc_call(SNESGetIterationNumber(fino.snes, &its));
  wasora_value(fino.vars.iterations) = (double)its;
  
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
