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
int fino_solve_eigen_slepc(Mat A, Mat B) {

  int i;
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

  petsc_call(EPSSetOperators(fino.eps, B, A));
  petsc_call(EPSSetWhichEigenpairs(fino.eps, EPS_LARGEST_MAGNITUDE));
  
  // problema generalizado no hermitico (por las condiciones de contorno)
  // TODO: no romper simetria!
  petsc_call(EPSSetProblemType(fino.eps, EPS_GNHEP));  
//  petsc_call(EPSSetProblemType(fino.eps, EPS_GHEP));    
  
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
  }
  
  // convergencia con respecto a la norma de las matrices
  petsc_call(EPSSetConvergenceTest(fino.eps, EPS_CONV_NORM));
  
  // tolerancia
  // TODO
  if (wasora_var(fino.vars.reltol) != 0) {
    petsc_call(EPSSetTolerances(fino.eps, wasora_var(fino.vars.reltol), PETSC_DECIDE));
  }

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
  
  // leemos la solucion posta
  petsc_call(EPSGetEigenpair(fino.eps, fino.nev-1, &fino.lambda, &xi, fino.phi, PETSC_NULL));
  
  // chequeamos que el autovalor sea real
  if (xi != 0) {
    wasora_push_error_message("eigen-solver found a complex eigenvalue (%g + i %g)", fino.lambda, xi);
    return WASORA_RUNTIME_ERROR;
  }
  
  // lo normalizamos
  // ya viene normalizado con alguna de las matrices
//  VecNormalize(fino.phi, PETSC_NULL);

  if (fino.nev != 1) {
    free(fino.eigenvalue);
    free(fino.eigenvector);
    fino.eigenvalue = calloc(fino.nev, sizeof(PetscScalar));
    fino.eigenvector = calloc(fino.nev, sizeof(Vec));
    
    for (i = 0; i < fino.nev; i++) {
      petsc_call(MatCreateVecs(fino.K, NULL, &fino.eigenvector[i]));
      petsc_call(EPSGetEigenpair(fino.eps, i, &fino.eigenvalue[i], &xi, fino.eigenvector[i], PETSC_NULL));
      if (xi != 0) {
        wasora_push_error_message("the eigenvalue %d is complex (%g + i %g)", i+1, fino.eigenvalue[i], xi);
        return WASORA_RUNTIME_ERROR;
      }
    }
  }
  
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

#undef  __FUNCT__
#define __FUNCT__ "fino_eigen_nev"
int fino_eigen_nev() {
  int i, j;
  double xi, fu;

  Vec tmp;
  Vec one;
  PetscScalar norm, M, L;

  VecDuplicate(fino.phi, &tmp);
  VecDuplicate(fino.phi, &one);
  VecSet(one, 1.0);

  // masa total
  MatMult(fino.M, one, tmp);
  VecDot(one, tmp, &xi);
  wasora_var_value(fino.vars.mass) = xi/(double)fino.degrees;


  // la fiesta de la ineficiencia
  for (i = 0; i < fino.nev; i++) {
    // autovalor convertido a frequencia
    fu = 4.0/1.0; // factor fumanchu
    wasora_vector_set(fino.vectors.f, i, sqrt(fu*2*M_PI/fino.eigenvalue[i]));
    wasora_vector_set(fino.vectors.omega, i, sqrt(fu/fino.eigenvalue[i]));

    // autovector i
    fino.vectors.phi[i]->size = fino.problem_size;

    // normalizado para que el maximo sea uno
    VecNorm(fino.eigenvector[i], NORM_INFINITY, &norm);
//            VecNorm(fino.eigenvector[i], NORM_1, &norm);
//            VecNorm(fino.eigenvector[i], NORM_2, &norm);            
    VecScale(fino.eigenvector[i], 1.0/norm);

    for (j = 0; j < fino.problem_size; j++) {
      petsc_call(VecGetValues(fino.eigenvector[i], 1, &j, &xi));
      wasora_vector_set(fino.vectors.phi[i], j, xi);
    }

    // masa modal
    MatMult(fino.Morig, fino.eigenvector[i], tmp);
    VecDot(fino.eigenvector[i], tmp, &M);
    wasora_vector_set(fino.vectors.M, i, M);

    // excitacion
    MatMult(fino.Morig, one, tmp);
    VecDot(fino.eigenvector[i], tmp, &L);
    wasora_vector_set(fino.vectors.L, i, L);

    wasora_vector_set(fino.vectors.Gamma, i, L/M);
    wasora_vector_set(fino.vectors.Me, i, gsl_pow_2(L)/(fino.degrees*M));
  }
  
  return WASORA_RUNTIME_OK;
}

#endif
