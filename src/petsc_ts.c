/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's transient solver using PETSc routines
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


PetscErrorCode fino_ts_residual(TS ts, PetscReal t, Vec phi, Vec phi_dot, Vec r, void *ctx) {
  
  wasora_var_value(wasora_special_var(t)) = t;

//  if (fino.math_type == math_type_nonlinear) {
    // TODO: separate volumetric from surface elements
    // (in case only natural BCs change with time)
    wasora_call(fino_phi_to_solution(phi));
    wasora_call(fino_build_bulk());
    wasora_call(fino_dirichlet_eval(fino.K, fino.b));
//  }  
    
  // compute the residual R(t,phi,phi_dot) = K*phi + M*phi_dot - b
  petsc_call(MatMult(fino.K, phi, r));
  petsc_call(MatMultAdd(fino.M, phi_dot, r, r));
  petsc_call(VecAXPY(r, -1.0, fino.b));

  // set dirichlet bcs  
  wasora_call(fino_dirichlet_set_r(r, phi));
  
  return WASORA_RUNTIME_OK;
}

PetscErrorCode fino_ts_jacobian(TS ts, PetscReal t, Vec T, Vec T_dot, PetscReal s, Mat J, Mat P,void *ctx) {

  Mat M;
  
  petsc_call(MatCopy(fino.K, J, SUBSET_NONZERO_PATTERN));
//  printf("J = K\n"); fino_print_petsc_matrix(J, PETSC_VIEWER_STDOUT_SELF);
  wasora_call(fino_dirichlet_set_J(J));
//  printf("J con bc\n"); fino_print_petsc_matrix(J, PETSC_VIEWER_STDOUT_SELF);
  
  petsc_call(MatDuplicate(fino.M, MAT_COPY_VALUES, &M));
//  printf("M = M\n"); fino_print_petsc_matrix(M, PETSC_VIEWER_STDOUT_SELF);
  wasora_call(fino_dirichlet_set_dRdphi_dot(M));
//  printf("M con bc\n"); fino_print_petsc_matrix(M, PETSC_VIEWER_STDOUT_SELF);

  petsc_call(MatAXPY(J, s, M, SAME_NONZERO_PATTERN));
  petsc_call(MatCopy(J, P, SAME_NONZERO_PATTERN));
  petsc_call(MatDestroy(&M));
  
  return WASORA_RUNTIME_OK;
}
