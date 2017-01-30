/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino eigenvalue problem solver using SLEPc routines
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
#ifdef HAVE_SLEPC
#include <slepceps.h>
#include <petscksp.h>

#include "fino.h"


#undef  __FUNCT__
#define __FUNCT__ "fino_solve_eigen_slepc"
int fino_solve_eigen_slepc(void) {

  PetscReal xi = 1.0;
  PetscInt nconv;

  // creamos el contexto del eigensolver
  if (fino.eps != NULL) {
    petsc_call(EPSDestroy(&fino.eps));
  }
  petsc_call(EPSCreate(PETSC_COMM_WORLD, &fino.eps));
  petsc_call(PetscObjectSetName((PetscObject)fino.eps, "eigenvalue-problem_solver"));
 
  // y obtenemos los contextos asociados
  petsc_call(EPSGetST(fino.eps, &fino.st));
  petsc_call(PetscObjectSetName((PetscObject)fino.st, "spectral_transformation"));
  petsc_call(STGetKSP(fino.st, &fino.ksp));
  petsc_call(PetscObjectSetName((PetscObject)fino.ksp, "linear_solver"));
  petsc_call(KSPGetPC(fino.ksp, &fino.pc));
  petsc_call(PetscObjectSetName((PetscObject)fino.pc, "preconditioner"));

  // TODO
  //petsc_call(EPSSetOperators(fino.eps, fino.A, fino.B));
  //petsc_call(EPSSetWhichEigenpairs(fino.eps, fino.eigen_spectrum));
  
//  petsc_call(EPSSetOperators(fino.eps, fino.B, fino.A));
//  petsc_call(EPSSetWhichEigenpairs(fino.eps, EPS_SMALLEST_MAGNITUDE));

  petsc_call(EPSSetOperators(fino.eps, fino.B, fino.A));
  petsc_call(EPSSetWhichEigenpairs(fino.eps, EPS_LARGEST_MAGNITUDE));
  
  // problema generalizado no hermitico (por las condiciones de contorno)
  // TODO: no romper simetria!
  petsc_call(EPSSetProblemType(fino.eps, EPS_GNHEP));  
  
  // TODO: ver bien esto del guess inicial
  //petsc_call(EPSSetInitialSpace(fino.eps, 1, &fino.guess));

  // elegimos el metodo de solucion del eps
  if (fino.eps_type != NULL) {
    petsc_call(EPSSetType(fino.eps, fino.eps_type));
  }
  
  // la transformada espectral
  if (fino.st_type != NULL) {
    petsc_call(STSetType(fino.st, fino.st_type));
  }
  // si no esta seteado el tipo se queja la SLEPc
  if (fino.st_shift.n_tokens != 0) {
    petsc_call(STSetShift(fino.st, wasora_evaluate_expression(&fino.st_shift)));
  }
  if (fino.st_anti_shift.n_tokens != 0) {
    petsc_call(STCayleySetAntishift(fino.st, wasora_evaluate_expression(&fino.st_anti_shift)));
  }

  // el KSP
  if (fino.ksp_type != NULL) {
    petsc_call(KSPSetType(fino.ksp, fino.ksp_type));
  }

  // el precondicionador
  if (fino.pc_type != NULL) {
    petsc_call(PCSetType(fino.pc, fino.pc_type));
  } else {
    // default gamg
  }  petsc_call(PCSetType(fino.pc, "lu"));

  // convergencia con respecto a la norma de las matrices
  petsc_call(EPSSetConvergenceTest(fino.eps, EPS_CONV_NORM));
  
  // tolerancia
  // TODO
//  if (wasora_var(fino.vars.reltol) != 0) {
//    petsc_call(EPSSetTolerances(fino.eps, wasora_var(fino.vars.reltol), PETSC_DECIDE));
//  }

  // si no nos pidieron que autovalor quieren, pedimos el primero
  if (fino.nev == 0) {
    fino.nev = 1;
  }
  
  // dimension del sub espacio
  if (fino.eps_ncv.n_tokens != 0) {
    petsc_call(EPSSetDimensions(fino.eps, fino.nev, (PetscInt)(wasora_evaluate_expression(&fino.eps_ncv)), PETSC_DEFAULT));
  } else {
    petsc_call(EPSSetDimensions(fino.eps, fino.nev, PETSC_DEFAULT, PETSC_DEFAULT));
  }
  
  // sobreescribimos con la linea de comandos
  petsc_call(EPSSetFromOptions(fino.eps));

  // do the work!
  petsc_call(EPSSolve(fino.eps));

  // chequeamos que haya convergido
  petsc_call(EPSGetConverged(fino.eps, &nconv));
  if (nconv < fino.nev) {
    wasora_push_error_message("eigen-solver obtained only %d converged eigen-pairs (%d requested)", nconv, fino.nev);
    return WASORA_RUNTIME_ERROR;
  }
  
  // leemos la solucion
  petsc_call(EPSGetEigenpair(fino.eps, fino.nev-1, &fino.lambda, &xi, fino.phi, PETSC_NULL));
  
  // chequeamos que el autovalor sea real
  if (xi != 0) {
    wasora_push_error_message("eigen-solver found a complex eigenvalue (%g + i %g)", fino.lambda, xi);
    return WASORA_RUNTIME_ERROR;
  }
    
  // lo normalizamos
  VecNormalize(fino.phi, PETSC_NULL);

  // obtenemos informacion auxiliar
/*  
  petsc_call(EPSGetErrorEstimate(fino.eps, 0, &xi));
  wasora_var(fino.vars.error_estimate) = (double)xi;

  petsc_call(EPSComputeResidualNorm(fino.eps, 0, &xi));
  wasora_var(fino.vars.residual_norm) = (double)xi;
  
  petsc_call(EPSComputeRelativeError(fino.eps, 0, &xi)); 
 wasora_var(fino.vars.rel_error) = (double)xi;
*/
  
  return WASORA_RUNTIME_OK;

}
#endif
