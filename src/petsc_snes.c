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

extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);

PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *ctx)
{
  PetscErrorCode ierr;
  Vec            x;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&x);CHKERRQ(ierr);
  printf("phi\n");
  fino_print_petsc_vector(x, PETSC_VIEWER_STDOUT_SELF);
  printf("b\n");
  fino_print_petsc_vector(fino.b, PETSC_VIEWER_STDOUT_SELF);
  printf("K\n");
  fino_print_petsc_matrix(fino.K, PETSC_VIEWER_STDOUT_SELF);
  return 0;
}


PetscErrorCode fino_solve_residual(SNES snes, Vec phi, Vec r,void *ctx) {
  
  // pasamos phi a la solucion porque K puede depender de phi
  wasora_call(fino_phi_to_solution(phi));
  wasora_call(fino_build_bulk());
  wasora_call(fino_set_essential_bc(fino.K, fino.b));
  
  MatMult(fino.K, phi, r);
  printf("residual\n");
  fino_print_petsc_vector(r, PETSC_VIEWER_STDOUT_SELF);
  
  return 0;
}

PetscErrorCode fino_solve_jacobian(SNES snes,Vec phi, Mat J, Mat P, void *ctx) {
  
  printf("jacobiano\n");
  MatDuplicate(fino.K, MAT_COPY_VALUES, &J);
  
  return 0;
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_solve_nonlinear_petsc"
int fino_solve_nonlinear_petsc(void) {

  SNESConvergedReason reason;
  SNES           snes;
  Mat            J;            /* Jacobian matrix */
  Vec            r;
  PetscInt       its;    

  SNESCreate(PETSC_COMM_WORLD, &snes);
  VecDuplicate(fino.phi, &r);
  MatDuplicate(fino.K, MAT_COPY_VALUES, &J);
  
  SNESSetFunction(snes, r, fino_solve_residual, NULL);
  SNESSetJacobian(snes, J, J, fino_solve_jacobian, NULL);
  
  // extract ksp and set
  // TODO
  
  SNESSetFromOptions(snes);

  // initial value
  VecSet(fino.phi, 0);
      
  // monitor
  // TODO

//  SNESMonitorSet(snes, SNESMonitorResidual, NULL, 0);
//  PetscViewerDrawOpen(PETSC_COMM_WORLD,0,0,0,0,400,400,&monP.viewer);
  SNESMonitorSet(snes,Monitor,NULL,0);
  
  // solve
  SNESSolve(snes, fino.b, fino.phi);

  
  // chequeamos que haya convergido
  petsc_call(SNESGetConvergedReason(snes, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's non-linear solver did not converge with reason '%s' (%d)", SNESConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }
  
  SNESGetIterationNumber(snes, &its);  
  wasora_value(fino.vars.iterations) = (double)its;
  
//  petsc_call(SNESGetResidualNorm(snes, wasora_value_ptr(fino.vars.residual_norm)));
  
  
  return WASORA_RUNTIME_OK;

}

