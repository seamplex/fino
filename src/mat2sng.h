#include <petscsys.h>

PetscErrorCode petsc_err;
#define petsc_call(s) {petsc_err = s; CHKERRQ(petsc_err);}

extern PetscErrorCode mat2sng(Mat A, PetscInt, PetscInt, PetscInt, PetscViewer);
