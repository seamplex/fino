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

#include <petsc.h>
#include <petscsys.h>
#include <petscksp.h>

#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "fino_solve_petsc_linear"
int fino_solve_petsc_linear(void) {

  KSPConvergedReason reason;
  PetscInt iterations;
  int i;

  // creamos un solver lineal
  if (fino.ksp == NULL) {
    petsc_call(KSPCreate(PETSC_COMM_WORLD, &fino.ksp));
  }
  
  petsc_call(KSPSetOperators(fino.ksp, fino.K, fino.K));
  petsc_call(KSPSetTolerances(fino.ksp, wasora_var(fino.vars.reltol),
                                        wasora_var(fino.vars.abstol),
                                        wasora_var(fino.vars.divtol),
                                        (PetscInt)wasora_var(fino.vars.max_iterations)));
  wasora_call(fino_set_pc());
  wasora_call(fino_set_ksp());

  // el monitor
  if (fino.progress_ascii) {  
    petsc_call(KSPMonitorSet(fino.ksp, fino_ksp_monitor, NULL, 0));
  }
  
  // do the work!
  petsc_call(KSPSolve(fino.ksp, fino.b, fino.phi));
  
  // chequeamos que haya convergido
  petsc_call(KSPGetConvergedReason(fino.ksp, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's linear solver did not converge with reason '%s' (%d)", KSPConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.progress_ascii) {
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


  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_ksp_monitor"
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
    
    if (fino.progress_ascii) {
      for (i = (int)(100*fino.progress_last); i < (int)(100*current_progress); i++) {
        printf(CHAR_PROGRESS_SOLVE);
        fflush(stdout);
      }
      fino.progress_last = current_progress;
    }
  }  

  return 0;
}



#undef  __FUNCT__
#define __FUNCT__ "fino_ksp_set"
int fino_set_ksp(void) {

  // el KSP
  if ((fino.ksp_type != NULL && strcasecmp(fino.ksp_type, "mumps") == 0) ||
      (fino.pc_type != NULL && strcasecmp(fino.pc_type,  "mumps") == 0) ||
      fino.commandline_mumps == PETSC_TRUE) {
    // mumps es un caso particular, hay que poner pre-only aca y en el pc otras cosas
    KSPSetType(fino.ksp, KSPPREONLY);
  } else if (fino.ksp_type != NULL) {
    petsc_call(KSPSetType(fino.ksp, fino.ksp_type));
  }
  // si no nos dieron nada, dejamos el default gmres, que no esta mal

  // sobreescribimos con la linea de comandos
  petsc_call(KSPSetFromOptions(fino.ksp));

  // set up  
  petsc_call(KSPSetUp(fino.ksp));
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_set_pc"
int fino_set_pc(void) {

  PetscInt i, j, d;
  PCType pc_type;

  PetscInt nearnulldim = 0;
  MatNullSpace nullsp = NULL;
  PetscScalar  dots[5];
  Vec          *nullvec;  
  PetscReal    *coords;
  Vec          vec_coords;
//  PetscViewer viewer;
  
  petsc_call(KSPGetPC(fino.ksp, &fino.pc));
  
  // si nos pidieron mumps, hay que usar LU o cholesky y poner MatSolverType
  if ((fino.ksp_type != NULL && strcasecmp(fino.ksp_type, "mumps") == 0) ||
      (fino.pc_type != NULL && strcasecmp(fino.pc_type,  "mumps") == 0) ||
      fino.commandline_mumps == PETSC_TRUE) {
#if PETSC_VERSION_GT(3,9,0)
//  petsc_call(PCSetType(fino.pc, PCLU));
    petsc_call(MatSetOption(fino.K, MAT_SPD, PETSC_TRUE)); /* set MUMPS id%SYM=1 */
    petsc_call(PCSetType(fino.pc, PCCHOLESKY));

    petsc_call(PCFactorSetMatSolverType(fino.pc, MATSOLVERMUMPS));
    petsc_call(PCFactorSetUpMatSolverType(fino.pc)); /* call MatGetFactor() to create F */
//  petsc_call(PCFactorGetMatrix(pc, &F));    
#else
    wasora_push_error_message("solver MUMPS needs at least PETSc 3.9.x");
    return WASORA_RUNTIME_ERROR;
#endif
  } else {

    if (fino.pc_type != NULL) {
      // si nos dieron uno, lo ponemos
      petsc_call(PCSetType(fino.pc, fino.pc_type));
    } else {
      // si no nos dieron y estamos en mechanical, ponemos GAMG 
      if (fino.problem_family == problem_family_mechanical) {
        petsc_call(PCSetType(fino.pc, PCGAMG));
      }
      // sino que quede el default
    }
  }  
  
  
  petsc_call(PCGetType(fino.pc, &pc_type));
  if (pc_type != NULL && strcmp(pc_type, PCGAMG) == 0) {
#if PETSC_VERSION_LT(3,8,0)
    PCGAMGSetThreshold(fino.pc, (PetscReal)wasora_var_value(fino.vars.gamg_threshold));
#else
    PCGAMGSetThreshold(fino.pc, (PetscReal *)(wasora_value_ptr(fino.vars.gamg_threshold)), 1);
#endif
  }

  if (fino.problem_family == problem_family_mechanical) {
    // las coordenadas (solo para break)
  // http://computation.llnl.gov/casc/linear_solvers/pubs/Baker-2009-elasticity.pdf
    // ojo que si estamos en node ordering no podemos usar set_near_nullspace_rigidbody
    if (fino.mesh->ordering == ordering_unknown_based && fino.set_near_nullspace == set_near_nullspace_rigidbody) {
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
//        petsc_call(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "coords.bin", FILE_MODE_WRITE, &viewer));
//        petsc_call(VecView(vec_coords, viewer));
//        petsc_call(PetscViewerDestroy(&viewer));
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
        petsc_call(MatSetNearNullSpace(fino.K, nullsp));
      break;  

      case set_near_nullspace_none:
        ;
      break;
    }
  }
  
  // esta linea trae problemas con GAMG en versiones viejas de PETSc y la verdad no hace falta
#if PETSC_VERSION_GT(3,8,0)
    petsc_call(PCSetUp(fino.pc));
#endif
  
  return WASORA_RUNTIME_OK;
}
