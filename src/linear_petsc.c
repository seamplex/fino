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
int fino_solve_linear_petsc(Mat A, Vec b) {

  KSPConvergedReason reason;
  PetscInt iterations;
  PetscInt i, j, d;
  PCType pc_type;

  PetscInt nearnulldim = 0;
  MatNullSpace nullsp = NULL;
  PetscScalar  dots[5];
  Vec          *nullvec;  
  PetscReal *coords;
  Vec       vec_coords;

  // creamos un solver lineal
  if (fino.ksp == NULL) {
    petsc_call(KSPCreate(PETSC_COMM_WORLD, &fino.ksp));
    petsc_call(KSPSetOperators(fino.ksp, A, A));
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
    petsc_call(PCSetType(fino.pc, "gamg"));
  }
  
  petsc_call(PCGetType(fino.pc, &pc_type));
  if (strcmp(pc_type, PCGAMG) == 0) {
#if PETSC_VERSION_LT(3,8,0)
    PCGAMGSetThreshold(fino.pc, (PetscReal)wasora_var_value(fino.vars.gamg_threshold));
#else
    PCGAMGSetThreshold(fino.pc, (PetscReal *)(wasora_value_ptr(fino.vars.gamg_threshold)), 1);
#endif
  }

  if (fino.problem_family == problem_family_break) {
    // las coordenadas (solo para break)
  // http://computation.llnl.gov/casc/linear_solvers/pubs/Baker-2009-elasticity.pdf
    // ojo que si estamos en node ordering no podemos usar set_near_nullspace_rigidbody
    if (fino.mesh->ordering == ordering_unknown_based && fino.set_near_nullspace == set_near_nullspace_rigidbody) {
      fino.set_near_nullspace = set_near_nullspace_fino;
    }

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
        petsc_call(MatSetNearNullSpace(A, nullsp));
        petsc_call(MatNullSpaceDestroy(&nullsp));
        petsc_call(VecDestroy(&vec_coords));      
      break;

      case set_near_nullspace_fino:
        nearnulldim = 6; 
        petsc_call(PetscMalloc1(nearnulldim, &nullvec));
        for (i = 0; i < nearnulldim; i++) {
          petsc_call(MatCreateVecs(A, &nullvec[i], NULL));
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
        petsc_call(MatSetNearNullSpace(A, nullsp));
      break;  

      case set_near_nullspace_none:
        ;
      break;
    }
  }

  // el monitor
  if (fino.shmem_progress_build != NULL || fino.shmem_progress_solve != NULL || fino.shmem_memory != NULL) {
    petsc_call(KSPMonitorSet(fino.ksp, fino_ksp_monitor, NULL, 0));
  }
  
  // sobreescribimos con la linea de comandos
  petsc_call(KSPSetFromOptions(fino.ksp));
  
  // do the work!
  petsc_call(KSPSolve(fino.ksp, b, fino.phi));
  
  // chequeamos que haya convergido
  petsc_call(KSPGetConvergedReason(fino.ksp, &reason));
  if (reason < 0) {
    wasora_push_error_message("PETSc's linear solver did not converge with reason '%s' (%d)", KSPConvergedReasons[reason], reason);
    return WASORA_RUNTIME_ERROR;
  }

  if (fino.shmem_progress_solve != NULL) {
    *fino.shmem_progress_solve = 1.0;
  }

    
  petsc_call(KSPGetIterationNumber(fino.ksp, &iterations));
  wasora_value(fino.vars.iterations) = (double)iterations;
  
  petsc_call(KSPGetResidualNorm(fino.ksp, wasora_value_ptr(fino.vars.residual_norm)));


  return WASORA_RUNTIME_OK;

}

PetscErrorCode fino_ksp_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy) {
//  wasora_value(fino.vars.iterations) = (double)n;
//  wasora_var_value(fino.vars.residual_norm) = rnorm;

  if (fino.shmem_memory != NULL) {
    getrusage(RUSAGE_SELF, &fino.resource_usage);
    *fino.shmem_memory = (double)(1024.0*fino.resource_usage.ru_maxrss);
  }

  if (fino.shmem_progress_solve != NULL) {
    *fino.shmem_progress_solve = 1 - 0.1*(log10(rnorm/wasora_var_value(fino.vars.reltol)));
//    printf("%d %g\n", n, 1 - 0.1*(log10(rnorm/wasora_var_value(fino.vars.reltol))));
  }

  
  return 0;
}
