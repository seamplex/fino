#include <petscsys.h>
#include <petscmat.h>

#include "mat2sng.h"

#define NTONES 20
#define NTHRES 10

PetscErrorCode mat2sng(Mat A, PetscInt maxwidth, PetscInt stride, PetscInt structure, PetscViewer viewer) {
  
  MatType   t;
  PetscBool isseqaij;
  PetscInt n, m;
  PetscInt w, h;
  PetscScalar **xi;
  PetscScalar s, max;
  
  PetscInt           i,j,ncols;
  const PetscInt    *cols;
  const PetscScalar *vals;
  
  
  // esto es RFC 2045 pero no camina
//  char base64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";
  // esto es segun el manual (y el fuente ) de sng, pero falla porque ni + ni / son isaplha() ni isdigit())
    char base64[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";

  petsc_call(MatGetType(A,&t));
  petsc_call(PetscStrcmp(t, MATSEQAIJ, &isseqaij));
  if (!isseqaij) {
    SETERRQ(PetscObjectComm((PetscObject)A), PETSC_ERR_ORDER, "mat2sng only works with seqaij matrices");
  }
  
  petsc_call(MatGetSize(A, &m, &n));
  
  if (maxwidth != 0) {
    stride = 1+n/maxwidth;
  } else if (stride == 0) {
    stride = 1;
  }
  
  w = n/stride + 1;
  h = m/stride + 1;

  // TODO: object name
  petsc_call(PetscViewerASCIIPrintf(viewer, "#SNG: PETSc matrix\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "\n")); 
  petsc_call(PetscViewerASCIIPrintf(viewer, "IHDR: {\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "  width: %d;\n", w));
  petsc_call(PetscViewerASCIIPrintf(viewer, "  height: %d;\n", h));
  petsc_call(PetscViewerASCIIPrintf(viewer, "  bitdepth: 8;\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "  using color: palette;\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "}\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  petsc_call(PetscViewerASCIIPrintf(viewer, "PLTE: {\n"));
  
  // 0 es blanco
  petsc_call(PetscViewerASCIIPrintf(viewer, "  (255, 255, 255)\n"));
  
  // NTONES tonos de rojo
  for (i = 0; i < NTONES; i++) {
    s = 1 - (PetscScalar)(i+1)/(PetscScalar)(NTONES);
    petsc_call(PetscViewerASCIIPrintf(viewer, "  (%d, %d, %d)\n", 255, (int)(255 * s), (int)(255 * s)));
  }
  
  // NTONES tonos de azul
  for (i = 0; i < NTONES; i++) {
    s = 1 - (PetscScalar)(i+1)/(PetscScalar)(NTONES);
    petsc_call(PetscViewerASCIIPrintf(viewer, "  (%d, %d, %d)\n", (int)(255 * s), (int)(255 * s), 255));
  }
  petsc_call(PetscViewerASCIIPrintf(viewer, "}\n")); 
  petsc_call(PetscViewerASCIIPrintf(viewer, "\n")); 
  petsc_call(PetscViewerASCIIPrintf(viewer, "IMAGE: {\n")); 
  petsc_call(PetscViewerASCIIPrintf(viewer, "  pixels base64\n")); 

  // listo, basta de chacara
  xi = calloc(h, sizeof(PetscScalar *));
  for (i = 0; i < h; i++) {
    xi[i] = calloc(w, sizeof(PetscScalar));
  }
  max = 0;

  for (i = 0; i < n; i++) {
    petsc_call(MatGetRow(A,i,&ncols,&cols,&vals));
    for (j = 0; j < ncols; j++) {
      if ((xi[i/stride][cols[j]/stride] += vals[j]) > max) {
        max = xi[i/stride][cols[j]/stride];
      }
    }
  }
  
  // a comerla!
  if (structure) {
    for (i = 0; i < w; i++) {
      for (j = 0; j < h; j++) {
      
        if (xi[i][j] == 0) {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[0]));
        } else if (xi[i][j] > 0) {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[NTONES]));
        } else {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[2*NTONES]));
        }
      }
      petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "}\n")); 
  } else {
    for (i = 0; i < w; i++) {
      for (j = 0; j < h; j++) {
      
        s = fabs(xi[i][j]) / max;
        
        if (xi[i][j] == 0) {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[0]));
        } else if (xi[i][j] > 0) {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[1 + 0      + NTHRES + (int)((NTONES-NTHRES)*s)]));
        } else {
          petsc_call(PetscViewerASCIIPrintf(viewer, "%c", base64[1 + NTONES + NTHRES + (int)((NTONES-NTHRES)*s)]));
        }
      }
      petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "}\n")); 
  }    
    
  
  for (i = 0; i < h; i++) {
    free(xi[i]);
  }
  free(xi);
  
  return 0;
}
