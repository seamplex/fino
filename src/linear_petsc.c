/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's linear solver using PETSc routines
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

#include <petsc.h>
#include <petscsys.h>
#include <petscksp.h>

#include "fino.h"


#undef  __FUNCT__
#define __FUNCT__ "fino_solve_linear_petsc"
int fino_solve_linear_petsc(void) {

  KSPConvergedReason reason;
  PetscInt iterations;
  PetscInt i, j, d;
  PCType pc_type;
//  PetscViewerAndFormat *vf;  

  PetscInt nearnulldim = 0;
  MatNullSpace nullsp = NULL;
  PetscScalar  dots[5];
  Vec          *nullvec;  
  const Vec *vecs;
  Vec cero;
  PetscReal norm;
  PetscBool has_const;
  PetscInt  n;
  PetscReal *coords;
  Vec       vec_coords;

  // creamos un solver lineal
  if (fino.ksp == NULL) {
    petsc_call(KSPCreate(PETSC_COMM_WORLD, &fino.ksp));
    petsc_call(KSPSetOperators(fino.ksp, fino.A, fino.A));
    
  }

  petsc_call(KSPSetTolerances(fino.ksp, wasora_var(fino.vars.reltol),
                                        wasora_var(fino.vars.abstol),
                                        wasora_var(fino.vars.divtol),
                                        (PetscInt)wasora_var(fino.vars.max_iterations)));

  // el KSP
  if (fino.ksp_type != NULL) {
    petsc_call(KSPSetType(fino.ksp, fino.ksp_type));
    // el default es gmres, no esta mal
  }

  // precondicionador
  petsc_call(KSPGetPC(fino.ksp, &fino.pc));
  
  // el precondicionador
  if (fino.pc_type != NULL) {
    petsc_call(PCSetType(fino.pc, fino.pc_type));
  } else {
    if (fino.problem == problem_break) {
      petsc_call(PCSetType(fino.pc, "gamg"));
    } else {
      petsc_call(PCSetType(fino.pc, "lu"));
    }
  }
  
  petsc_call(PCGetType(fino.pc, &pc_type));
  if (strcmp(pc_type, PCGAMG) == 0) {
    PCGAMGSetThreshold(fino.pc, (PetscReal)wasora_var_value(fino.vars.gamg_threshold));
  }

  // las coordenadas (solo para break)
// http://computation.llnl.gov/casc/linear_solvers/pubs/Baker-2009-elasticity.pdf
  switch(fino.set_near_nullspace) {
    
    case set_near_nullspace_rigidbody:
      petsc_call(VecCreate(MPI_COMM_WORLD, &vec_coords));
      petsc_call(VecSetBlockSize(vec_coords, fino.degrees));
      petsc_call(VecSetSizes(vec_coords, PETSC_DECIDE, fino.dimensions * fino.mesh->n_nodes));
      petsc_call(VecSetUp(vec_coords));
      petsc_call(VecGetArray(vec_coords, &coords));
      
      for (j = 0; j < fino.mesh->n_nodes; j++) {
        for (d = 0; d < fino.dimensions; d++) {
          coords[fino.mesh->node[j].index[d]] = fino.mesh->node[j].x[d];
        }
      }
      
      petsc_call(VecRestoreArray(vec_coords, &coords));
      petsc_call(MatNullSpaceCreateRigidBody(vec_coords, &nullsp));
      petsc_call(MatSetNearNullSpace(fino.A, nullsp));
      petsc_call(MatNullSpaceDestroy(&nullsp));
      petsc_call(VecDestroy(&vec_coords));      
    break;
    
    case set_near_nullspace_fino:
      nearnulldim = 6; 
      petsc_call(PetscMalloc1(nearnulldim, &nullvec));
      for (i = 0; i < nearnulldim; i++) {
        petsc_call(MatCreateVecs(fino.A, &nullvec[i], NULL));
      }
      for (j = 0; j < fino.mesh->n_nodes; j++) {
        // traslaciones
        VecSetValue(nullvec[0], fino.mesh->node[j].index[0], 1.0, INSERT_VALUES);
        VecSetValue(nullvec[1], fino.mesh->node[j].index[1], 1.0, INSERT_VALUES);
        VecSetValue(nullvec[2], fino.mesh->node[j].index[2], 1.0, INSERT_VALUES);

        // rotaciones
        VecSetValue(nullvec[3], fino.mesh->node[j].index[0], +fino.mesh->node[j].x[1], INSERT_VALUES);
        VecSetValue(nullvec[3], fino.mesh->node[j].index[1], -fino.mesh->node[j].x[0], INSERT_VALUES);

        VecSetValue(nullvec[4], fino.mesh->node[j].index[1], -fino.mesh->node[j].x[2], INSERT_VALUES);
        VecSetValue(nullvec[4], fino.mesh->node[j].index[2], +fino.mesh->node[j].x[1], INSERT_VALUES);

        VecSetValue(nullvec[5], fino.mesh->node[j].index[0], +fino.mesh->node[j].x[2], INSERT_VALUES);
        VecSetValue(nullvec[5], fino.mesh->node[j].index[2], -fino.mesh->node[j].x[0], INSERT_VALUES);
      }

      for (i = 0; i < 3; i++) {
        VecNormalize(nullvec[i], PETSC_NULL);
      }

      // tomado de MatNullSpaceCreateRigidBody()
      for (i = 3; i < nearnulldim; i++) {
        // Orthonormalize vec[i] against vec[0:i-1]
        VecMDot(nullvec[i], i, nullvec, dots);
        for (j= 0; j < i; j++) {
          dots[j] *= -1.;
        }
        VecMAXPY(nullvec[i],i,dots,nullvec);
        VecNormalize(nullvec[i], PETSC_NULL);
      }

      petsc_call(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, nearnulldim, nullvec, &nullsp));
      petsc_call(MatSetNearNullSpace(fino.A, nullsp));
//      petsc_call(MatNullSpaceDestroy(&nullsp));
//      petsc_call(VecDestroy(&vec_coords));      
    break;  
   
/*    
    case set_near_nullspace_setcoordinates:
      petsc_call(PetscMalloc1(fino.dimensions * fino.mesh->n_nodes, &coords));
      for (j = 0; j < fino.mesh->n_nodes; j++) {
        for (d = 0; d < fino.dimensions; d++) {
          coords[fino.mesh->node[j].index[d]] = fino.mesh->node[j].x[d];
        }
      }
      
      // esto termina en PCSetCoordinates_AGG() en src/ksp/pc/impls/gamg/agg.c
      petsc_call(PCSetCoordinates(fino.pc, fino.dimensions, fino.mesh->n_nodes, coords));
      petsc_call(PetscFree(coords));
 */
    break;
    
    case set_near_nullspace_none:
      ;
    break;
  }

/*  
  MatGetNearNullSpace(fino.A, &nullsp);
  MatNullSpaceGetVecs(nullsp, &has_const, &n, &nullvec);
  petsc_call(MatCreateVecs(fino.A, NULL, &cero));
  
  petsc_call(MatMult(fino.A, nullvec[5], cero));
  fino_print_petsc_vector(cero, PETSC_VIEWER_STDOUT_WORLD);
  exit(0);
 */
//  petsc_call(VecNorm(cero, NORM_1, &norm));
//  printf("%g\n", norm);
  
  // el monitor
//  petsc_call(KSPMonitorSet(fino.ksp, fino_monitor_dots, NULL, 0));
  
//  petsc_call(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf));
//  petsc_call(KSPMonitorSet(fino.ksp, KSPMonitorDefault, vf, NULL));  
   
  
  
  // sobreescribimos con la linea de comandos
  petsc_call(KSPSetFromOptions(fino.ksp));
  
  // do the work!
  petsc_call(KSPSolve(fino.ksp, fino.b, fino.phi));

  
//  petsc_call(MatGetNearNullSpace(fino.A, &nullsp));
//  petsc_call(MatNullSpaceGetVecs(nullsp, &has_const, &n, &vecs));
  
  // chequeamos que haya convergido
  petsc_call(KSPGetConvergedReason(fino.ksp, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's linear solver did not converge with reason '%s' (%d)", KSPConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }

  petsc_call(KSPGetIterationNumber(fino.ksp, &iterations));
  wasora_value(fino.vars.iterations) = (double)iterations;


  return WASORA_RUNTIME_OK;

}

