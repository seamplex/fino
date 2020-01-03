/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's non-linear solver using PETSc routines
 *
 *  Copyright (C) 2015--2017 jeremy theler
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

/*  
PetscErrorCode MyComputeFunction(SNES snes,Vec x,Vec r,void *ctx) {
  Vec minusb;
  VecDuplicate(fino.b, &minusb);
  VecScale(minusb, -1);
  MatMultAdd(fino.K, x, minusb, r);
  MatMult(fino.K, x, r);
  
  return(0);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_solve_linear_petsc"
int fino_solve_linear_petsc(Mat A, Vec b) {

  SNES           snes;
//  Vec            x,r;
  PetscInt       its;    

  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetFunction(snes,NULL,MyComputeFunction,NULL);
//  SNESSetJacobian(snes,NULL,NULL,MyComputeJacobian,NULL);
  
  SNESSetFromOptions(snes);
  
  SNESSolve(snes, fino.b, fino.phi);
  SNESGetIterationNumber(snes,&its);  
  
  if (fino.shmem_progress_solve != NULL) {
    *fino.shmem_progress_solve = 1.0;
  }


  return WASORA_RUNTIME_OK;

}
*/
