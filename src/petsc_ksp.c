/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's linear solver using PETSc routines
 *
 *  Copyright (C) 2015--2020 jeremy theler
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

int fino_solve_petsc_linear(void) {

  KSPConvergedReason reason;
  PetscInt iterations;
  PC pc;
  
  time_checkpoint(build_begin);
  wasora_call(fino_build_bulk());
  wasora_call(fino_dirichlet_eval(fino.K, fino.b));
  wasora_call(fino_dirichlet_set_K(fino.K, fino.b));
  time_checkpoint(build_end);

  time_checkpoint(solve_begin);
  // create a KSP object if needed
  if (fino.ksp == NULL) {
    petsc_call(KSPCreate(PETSC_COMM_WORLD, &fino.ksp));
    
    // set the monitor for the ascii progress
    if (fino.progress_ascii == PETSC_TRUE) {  
      petsc_call(KSPMonitorSet(fino.ksp, fino_ksp_monitor, NULL, 0));
    }
  }
  
  petsc_call(KSPSetOperators(fino.ksp, fino.K, fino.K));
  petsc_call(KSPSetTolerances(fino.ksp, wasora_var(fino.vars.reltol),
                                        wasora_var(fino.vars.abstol),
                                        wasora_var(fino.vars.divtol),
                                        (PetscInt)wasora_var(fino.vars.max_iterations)));
  
  // TODO: this call sets up the nearnullspace, shouldn't that be in another place?
  petsc_call(KSPGetPC(fino.ksp, &pc));
  wasora_call(fino_set_pc(pc));
  wasora_call(fino_set_ksp(fino.ksp));

  // K is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner
  petsc_call(MatSetOption(fino.K, MAT_SYMMETRIC, PETSC_TRUE));  
  
  // try to use the solution as the initial guess (it already has Dirichlet BCs
  // but in quasi-static it has the previous solution which should be similar)
  if ((fino.ksp_type == NULL || strcasecmp(fino.ksp_type, "mumps") != 0) &&
      (fino.pc_type  == NULL || strcasecmp(fino.pc_type,  "mumps") != 0)) {
    // mumps cannot be used with a non-zero guess  
    petsc_call(KSPSetInitialGuessNonzero(fino.ksp, PETSC_TRUE));
  } 
  
  
  // do the work!
  petsc_call(KSPSolve(fino.ksp, fino.b, fino.phi));
  
  // check for convergence
  petsc_call(KSPGetConvergedReason(fino.ksp, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's linear solver did not converge with reason '%s' (%d)", KSPConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }

  // finish the progress line
  if (fino.progress_ascii == PETSC_TRUE) {
    int i;
    if (wasora.nprocs == 1) {
      for (i = (int)(100*fino.progress_last); i < 100; i++) {
        printf(CHAR_PROGRESS_SOLVE);
      }
    }
    if (wasora.rank == 0) {
      printf("\n");
      fflush(stdout);
    }  
  }

    
  petsc_call(KSPGetIterationNumber(fino.ksp, &iterations));
  wasora_value(fino.vars.iterations) = (double)iterations;
  
  petsc_call(KSPGetResidualNorm(fino.ksp, wasora_value_ptr(fino.vars.residual_norm)));
  time_checkpoint(solve_end);


  return WASORA_RUNTIME_OK;

}

PetscErrorCode fino_ksp_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy) {
//  wasora_value(fino.vars.iterations) = (double)n;
//  wasora_var_value(fino.vars.residual_norm) = rnorm;
  int i;
  double current_progress;
  
  if (wasora.rank == 0) {
  
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

//    printf("%d %e %.0f\n", n, rnorm/r0, 100*current_progress);
    
    if (fino.progress_ascii == PETSC_TRUE) {
      for (i = (int)(100*fino.progress_last); i < (int)(100*current_progress); i++) {
        printf(CHAR_PROGRESS_SOLVE);
        fflush(stdout);
      }
      fino.progress_last = current_progress;
    }
  }  

  return 0;
}



int fino_set_ksp(KSP ksp) {

  // the KSP type
  if ((fino.ksp_type != NULL && strcasecmp(fino.ksp_type, "mumps") == 0) ||
      (fino.pc_type  != NULL && strcasecmp(fino.pc_type,  "mumps") == 0)) {
    // mumps is a particular case, see fino_set_pc
    KSPSetType(ksp, KSPPREONLY);
  } else if (fino.ksp_type != NULL) {
    petsc_call(KSPSetType(ksp, fino.ksp_type));
  } else {
    // by default use whatever PETSc/SLEPc like
//    petsc_call(KSPSetType(ksp, KSPGMRES));
  }  

  // read command-line options
  petsc_call(KSPSetFromOptions(ksp));

  return WASORA_RUNTIME_OK;
}

int fino_set_pc(PC pc) {

  PetscInt i, j, d;
  PCType pc_type;

  PetscInt nearnulldim = 0;
  MatNullSpace nullsp = NULL;
  PetscScalar  dots[5];
  Vec          *nullvec;  
  PetscReal    *coords;
  Vec          vec_coords;

  // if we were asked for mumps, then either LU o cholesky needs to be used
  // and MatSolverType to mumps
  if ((fino.ksp_type != NULL && strcasecmp(fino.ksp_type, "mumps") == 0) ||
      (fino.pc_type != NULL && strcasecmp(fino.pc_type,  "mumps") == 0)) {
#if PETSC_VERSION_GT(3,9,0)
    petsc_call(MatSetOption(fino.K, MAT_SPD, PETSC_TRUE)); /* set MUMPS id%SYM=1 */
    petsc_call(PCSetType(pc, PCCHOLESKY));

    petsc_call(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    petsc_call(PCFactorSetUpMatSolverType(pc)); /* call MatGetFactor() to create F */
//  petsc_call(PCFactorGetMatrix(pc, &F));    
#else
    wasora_push_error_message("solver MUMPS needs at least PETSc 3.9.x");
    return WASORA_RUNTIME_ERROR;
#endif
  } else {

    if (fino.pc_type != NULL) {
      petsc_call(PCSetType(pc, fino.pc_type));
    } else {
      if (fino.problem_family == problem_family_mechanical) {
        petsc_call(PCSetType(pc, PCGAMG));
      } else {
        // leave default
          ;
//        petsc_call(PCSetType(pc, PCCHOLESKY));        
      }
    }
  }  
  
  
  petsc_call(PCGetType(pc, &pc_type));
  if (pc_type != NULL && strcmp(pc_type, PCGAMG) == 0) {
#if PETSC_VERSION_LT(3,8,0)
    PCGAMGSetThreshold(pc, (PetscReal)wasora_var_value(fino.vars.gamg_threshold));
#else
    PCGAMGSetThreshold(pc, (PetscReal *)(wasora_value_ptr(fino.vars.gamg_threshold)), 1);
#endif
    petsc_call(PCGAMGSetNSmooths(pc, 1));
  }

  if (fino.problem_family == problem_family_mechanical) {
    // las coordenadas (solo para break)
  // http://computation.llnl.gov/casc/linear_solvers/pubs/Baker-2009-elasticity.pdf
    // ojo que si estamos en node ordering no podemos usar set_near_nullspace_rigidbody
    if (fino.mesh->ordering == ordering_dof_major && fino.set_near_nullspace == set_near_nullspace_rigidbody) {
      fino.set_near_nullspace = set_near_nullspace_fino;
    }

    switch(fino.set_near_nullspace) {

      case set_near_nullspace_rigidbody:
        petsc_call(MatCreateVecs(fino.K, NULL, &vec_coords));
        petsc_call(VecSetBlockSize(vec_coords, fino.degrees));
        petsc_call(VecSetUp(vec_coords));
        petsc_call(VecGetArray(vec_coords, &coords));

        for (j = fino.first_node; j < fino.last_node; j++) {          
          for (d = 0; d < fino.dimensions; d++) {
            coords[fino.mesh->node[j].index_dof[d]-fino.first_row] = fino.mesh->node[j].x[d];
          }
        }

        petsc_call(VecRestoreArray(vec_coords, &coords));
        petsc_call(MatNullSpaceCreateRigidBody(vec_coords, &nullsp));
        petsc_call(MatSetNearNullSpace(fino.K, nullsp));
        petsc_call(MatNullSpaceDestroy(&nullsp));
        petsc_call(VecDestroy(&vec_coords));
      break;

      case set_near_nullspace_fino:
        nearnulldim = 6; 
        petsc_call(PetscMalloc1(nearnulldim, &nullvec));
        for (i = 0; i < nearnulldim; i++) {
          petsc_call(MatCreateVecs(fino.K, &nullvec[i], NULL));
        }
        for (j = 0; j < fino.mesh->n_nodes; j++) {
          // traslaciones
          VecSetValue(nullvec[0], fino.mesh->node[j].index_dof[0], 1.0, INSERT_VALUES);
          VecSetValue(nullvec[1], fino.mesh->node[j].index_dof[1], 1.0, INSERT_VALUES);
          VecSetValue(nullvec[2], fino.mesh->node[j].index_dof[2], 1.0, INSERT_VALUES);

          // rotaciones
          VecSetValue(nullvec[3], fino.mesh->node[j].index_dof[0], +fino.mesh->node[j].x[1], INSERT_VALUES);
          VecSetValue(nullvec[3], fino.mesh->node[j].index_dof[1], -fino.mesh->node[j].x[0], INSERT_VALUES);

          VecSetValue(nullvec[4], fino.mesh->node[j].index_dof[1], -fino.mesh->node[j].x[2], INSERT_VALUES);
          VecSetValue(nullvec[4], fino.mesh->node[j].index_dof[2], +fino.mesh->node[j].x[1], INSERT_VALUES);

          VecSetValue(nullvec[5], fino.mesh->node[j].index_dof[0], +fino.mesh->node[j].x[2], INSERT_VALUES);
          VecSetValue(nullvec[5], fino.mesh->node[j].index_dof[2], -fino.mesh->node[j].x[0], INSERT_VALUES);
        }

        for (i = 0; i < 3; i++) {
          VecNormalize(nullvec[i], PETSC_NULL);
        }

        // from MatNullSpaceCreateRigidBody()
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
        petsc_call(MatSetNearNullSpace(fino.K, nullsp));
      break;  

      case set_near_nullspace_none:
        ;
      break;
    }
  }
  
  return WASORA_RUNTIME_OK;
}
