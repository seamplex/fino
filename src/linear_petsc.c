/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's linear solver using PETSc routines
 *
 *  Copyright (C) 2015 jeremy theler
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
//  PetscViewerAndFormat *vf;  

  PetscInt nearnulldim;
  MatNullSpace nullsp;
  PetscScalar  dots[5];
  Vec          *nullvecs;  
  PetscReal *coords;

  // creamos un solver lineal
  if (fino.ksp == NULL) {
    petsc_call(KSPCreate(PETSC_COMM_WORLD, &fino.ksp));
    petsc_call(KSPSetOperators(fino.ksp, fino.A, fino.A));
    
  }

  petsc_call(KSPSetTolerances(fino.ksp, wasora_var(fino.vars.rtol),
                                        wasora_var(fino.vars.atol),
                                        wasora_var(fino.vars.divtol),
                                        (PetscInt)wasora_var(fino.vars.max_iterations)));

  // el KSP
  if (fino.ksp_type != NULL) {
    petsc_call(KSPSetType(fino.ksp, fino.ksp_type));
//  } else {
//    petsc_call(KSPSetType(fino.ksp, "bcgs"));
  }

  // precondicionador
  petsc_call(KSPGetPC(fino.ksp, &fino.pc));
  
  // el precondicionador
  if (fino.pc_type != NULL) {
    petsc_call(PCSetType(fino.pc, fino.pc_type));
//  } else if (fino.problem == problem_break) {
//    petsc_call(PCSetType(fino.pc, "gamg"));
  } else {
    petsc_call(PCSetType(fino.pc, "lu"));
  }

  // las coordenadas (solo para break)
// http://computation.llnl.gov/casc/linear_solvers/pubs/Baker-2009-elasticity.pdf
  if (fino.use_pcsetcoordinates) {
    petsc_call(PetscMalloc1(fino.dimensions * fino.mesh->n_nodes, &coords));
    for (j = 0; j < fino.mesh->n_nodes; j++) {
      for (d = 0; d < fino.dimensions; d++) {
        coords[fino.mesh->node[j].index[d]] = fino.mesh->node[j].x[d];
      }
    }
    petsc_call(PCSetCoordinates(fino.pc, fino.dimensions, fino.mesh->n_nodes, coords));
    petsc_call(PetscFree(coords));
  }
  
  if (fino.set_near_null_space) {
    nearnulldim = 6; 
    petsc_call(PetscMalloc1(nearnulldim, &nullvecs));
    for (i = 0; i < nearnulldim; i++) {
      petsc_call(MatCreateVecs(fino.A, &nullvecs[i], NULL));
    }
    for (j = 0; j < fino.mesh->n_nodes; j++) {
      // traslaciones
      VecSetValue(nullvecs[0], fino.mesh->node[j].index[0], 1.0, INSERT_VALUES);
      VecSetValue(nullvecs[1], fino.mesh->node[j].index[1], 1.0, INSERT_VALUES);
      VecSetValue(nullvecs[2], fino.mesh->node[j].index[2], 1.0, INSERT_VALUES);

      // rotaciones
      VecSetValue(nullvecs[3], fino.mesh->node[j].index[0], +fino.mesh->node[j].x[1], INSERT_VALUES);
      VecSetValue(nullvecs[3], fino.mesh->node[j].index[1], -fino.mesh->node[j].x[0], INSERT_VALUES);

      VecSetValue(nullvecs[4], fino.mesh->node[j].index[0], -fino.mesh->node[j].x[2], INSERT_VALUES);
      VecSetValue(nullvecs[4], fino.mesh->node[j].index[2], +fino.mesh->node[j].x[1], INSERT_VALUES);

      VecSetValue(nullvecs[5], fino.mesh->node[j].index[1], +fino.mesh->node[j].x[2], INSERT_VALUES);
      VecSetValue(nullvecs[5], fino.mesh->node[j].index[2], -fino.mesh->node[j].x[0], INSERT_VALUES);
    }

    for (i = 0; i < 3; i++) {
      VecNormalize(nullvecs[i], PETSC_NULL);
    }

    for (i=3; i < nearnulldim; i++) {
      // Orthonormalize vec[i] against vec[0:i-1]
      VecMDot(nullvecs[i], i, nullvecs, dots);
      for (j= 0; j < i; j++) {
        dots[j] *= -1.;
      }
      VecMAXPY(nullvecs[i],i,dots,nullvecs);
      VecNormalize(nullvecs[i], PETSC_NULL);
    }
  
    petsc_call(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, nearnulldim, nullvecs, &nullsp));
    petsc_call(MatSetNearNullSpace(fino.A, nullsp));
  }
  
  // el monitor
//  petsc_call(KSPMonitorSet(fino.ksp, fino_monitor_dots, NULL, 0));
  
//  petsc_call(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf));
//  petsc_call(KSPMonitorSet(fino.ksp, KSPMonitorDefault, vf, NULL));  
   
  
  
  // sobreescribimos con la linea de comandos
  petsc_call(KSPSetFromOptions(fino.ksp));
  
  // do the work!
  petsc_call(KSPSolve(fino.ksp, fino.b, fino.phi));

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

