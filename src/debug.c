/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino debugging and benchmarking routines
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>



#include <gsl/gsl_errno.h>

#include <petsc.h>
#include <slepcsys.h>
#include <petscdraw.h>

#include "fino.h"
#include "mat2sng.h"
#include "version.h"

#define fino_debug_insert_spaces(n)     for (j = 0; j < (n); j++) petsc_call(PetscViewerASCIIPrintf(debug->viewer, " "));
#define wasora_null_free(p) if (p!=NULL) {free(p); p=NULL;}

#undef  __FUNCT__
#define __FUNCT__ "fino_debug_open"
int fino_debug_open(debug_t *debug) {
  int i;
  time_t tm;
  struct utsname computer;

  char libversion[BUFFER_SIZE];

  uname(&computer);
  time(&tm);

  if (debug->file == NULL) {
    wasora_push_error_message("no FILE given to FINO_DEBUG");
    return WASORA_RUNTIME_ERROR;
  }

  // necesitamos el open file para evaluar a ver si es igual que antes pero no nos sirve el handler
  // porque se lo tenemos que pasar a la petsc en su formato
  debug->file->mode = strdup("w");
  wasora_call(wasora_instruction_open_file(debug->file));
  wasora_instruction_close_file(debug->file);
  debug->file_opened = 1;
  
  petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, debug->file->path, &debug->viewer));

  

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% fino debugging and benchmarking output\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% %s\n", getlogin()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% %s\n", ctime(&tm)));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "fino debugging and benchmarking output\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "=========================================\n\n"));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "code invocation\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------\n\n"));
   
// PetscGetUserName()
// PetscGetHostName()  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "    %s@%s$", getlogin(), computer.nodename));
  for (i = 0; i < wasora.argc; i++) {
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s", wasora.argv[i]));
  }
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\non %s\n\n", ctime(&tm)));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "code version\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------------\n\n"));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "[%s](https://bitbucket.org/gtheler/fino) %s  \n", plugin_name(), plugin_version()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n\n", plugin_description()));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s\n", plugin_longversion()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s\n", plugin_copyright()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n"));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "fino was compiled on %s  \n", COMPILATION_DATE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "by %s on %s (%s) using %s  \n", COMPILATION_USERNAME, COMPILATION_HOSTNAME, COMPILATION_ARCH, CCOMPILER_VERSION));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "and linked against  \n"));
  petsc_call(PetscGetVersion(libversion, BUFFER_SIZE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n", libversion));
  petsc_call(SlepcGetVersion(libversion, BUFFER_SIZE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n", libversion));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));


  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));

  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_debug"
int fino_instruction_debug(void *arg) {
  
  debug_t *debug = (debug_t *)arg;

  struct rusage resource_usage;

  PetscViewer viewer, viewer2;
  char *filename;

  petsc_call(PetscMemoryGetMaximumUsage(wasora_value_ptr(fino.vars.memory_usage_petsc)));
  petsc_call(PetscGetFlops(wasora_value_ptr(fino.vars.flops_petsc)));

  int size;
  if (debug->matrices_size.n_tokens != 0) {
    size = (int)(wasora_evaluate_expression(&debug->matrices_size));
  } else {
    size = DEFAULT_MATRICES_X_SIZE;
  }

    
  if (debug->file != NULL) {
    if (debug->file_opened == 0) {
      wasora_call(fino_debug_open(debug));
    }
    filename = malloc(strlen(debug->file->path)+32);

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "problem static step %d\n", (int)(wasora_var(wasora_special_var(step_static)))));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-----------------------\n", (int)(wasora_var(wasora_special_var(step_static)))));


    // eigencurrent solution
    if (fino.eps != NULL) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n### SLEPc's EPSView output\n\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
      petsc_call(EPSView(fino.eps, debug->viewer));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
        
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "### Eigenvalue problem result\n\n"));

      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------- ---------------                                              ------------------\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "requested eigenvalue   $\\lambda$                                                    %.10f\n", wasora_value(fino.vars.lambda)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "residual norm          $\\|R \\phi - \\lambda F \\phi\\|_2$                              %.4g\n", wasora_value(fino.vars.residual_norm)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "relative error         $\\|R \\phi - \\lambda F \\phi\\|_2 / \\| \\lambda F \\phi \\|_2$     %.4g\n", wasora_value(fino.vars.rel_error)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "error estimate         $\\|\\lambda - \\lambda_\\text{real}\\|$                          %.4g\n", wasora_value(fino.vars.error_estimate)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------- ---------------                                              ------------------\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
        
    } else if (fino.ksp != NULL) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n### PETSc's KSPView output\n\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
      petsc_call(KSPView(fino.ksp, debug->viewer));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
    }

    // system resources
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "### System resource usage\n\n"));
/* TODO: esto da errores en el valgrind   
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "user CPU time                          %.4g seconds\n", resource_usage.ru_utime.tv_sec + 1e-6*resource_usage.ru_utime.tv_usec));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "system CPU time                        %.4g seconds\n", resource_usage.ru_stime.tv_sec + 1e-6*resource_usage.ru_stime.tv_usec));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "total                                  %.4g seconds\n", resource_usage.ru_utime.tv_sec + resource_usage.ru_stime.tv_sec + 1e-6*(resource_usage.ru_utime.tv_usec + resource_usage.ru_stime.tv_usec)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "time insumed to                        cpu [secs]         wall [secs]\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "build the matrices                     %.4g               %.4g\n", wasora_value(fino.vars.time_cpu_build), wasora_value(fino.vars.time_wall_build)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "solve the problem                      %.4g               %.4g\n", wasora_value(fino.vars.time_cpu_solve), wasora_value(fino.vars.time_wall_solve)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "total                                  %.4g               %.4g\n", wasora_value(fino.vars.time_cpu_build)+wasora_value(fino.vars.time_cpu_solve), wasora_value(fino.vars.time_wall_build)+wasora_value(fino.vars.time_wall_solve)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "maximum memory resident set size       %ld kb\n", resource_usage.ru_maxrss));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "PETSc's maximum memory set size        %.0f kb\n", wasora_value(fino.vars.memory_usage_petsc)/1024.0));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of soft page faults             %ld\n", resource_usage.ru_minflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of hard page faults             %ld\n", resource_usage.ru_majflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of swaps                        %ld\n", resource_usage.ru_majflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of block input operations       %ld\n", resource_usage.ru_inblock));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of block output operations      %ld\n", resource_usage.ru_oublock));
 */
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
    

    if (debug->matrices & DEBUG_MATRICES_ASCII) {
      PetscViewer ascii_file;

      sprintf(filename, "%s-A.txt", debug->file->path);
      petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
      fino_print_petsc_matrix(fino.A, ascii_file);
      petsc_call(PetscViewerDestroy(&ascii_file));
      
      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.txt", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        fino_print_petsc_matrix(fino.B, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.txt", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        fino_print_petsc_vector(fino.b, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_ASCII_STRUCT) {
      PetscViewer ascii_file;

      sprintf(filename, "%s-A.str", debug->file->path);
      petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
      fino_print_petsc_matrix_struct(fino.A, ascii_file);
      petsc_call(PetscViewerDestroy(&ascii_file));

      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.str", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        fino_print_petsc_matrix_struct(fino.B, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.txt", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        fino_print_petsc_vector(fino.b, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_BINARY) {
      sprintf(filename, "%s-A.bin", debug->file->path);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
      MatView(fino.A, viewer);
      PetscViewerDestroy(&viewer);
      
      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.bin", debug->file->path);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
        MatView(fino.B, viewer);
        PetscViewerDestroy(&viewer);
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.bin", debug->file->path);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
        VecView(fino.b, viewer);
        PetscViewerDestroy(&viewer);
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_COMPRESSED_BINARY) {
      sprintf(filename, "%s-A.gz", debug->file->path);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
      MatView(fino.A, viewer);
      PetscViewerDestroy(&viewer);
      
      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.gz", debug->file->path);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
        MatView(fino.B, viewer);
        PetscViewerDestroy(&viewer);
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.gz", debug->file->path);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
        VecView(fino.b, viewer);
        PetscViewerDestroy(&viewer);
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_ASCII) {
      sprintf(filename, "%s-A.asc", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
      MatView(fino.A, viewer);
      PetscViewerDestroy(&viewer);

      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.asc", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
        MatView(fino.B, viewer);
        PetscViewerDestroy(&viewer);
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.asc", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
        VecView(fino.b, viewer);
        PetscViewerDestroy(&viewer);
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_OCTAVE) {
      sprintf(filename, "%s-A.m", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
      MatView(fino.A, viewer);
      PetscViewerDestroy(&viewer);

      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.m", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView(fino.B, viewer);
        PetscViewerDestroy(&viewer);
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.m", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        VecView(fino.b, viewer);
        PetscViewerDestroy(&viewer);
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_DENSE) {
      sprintf(filename, "%s-A.den", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
      MatView(fino.A, viewer);
      PetscViewerDestroy(&viewer);

      if (fino.math_type == math_eigen) {
        sprintf(filename, "%s-B.den", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
        MatView(fino.B, viewer);
        PetscViewerDestroy(&viewer);
      } else if (fino.math_type == math_linear) {
        sprintf(filename, "%s-b.den", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
        VecView(fino.b, viewer);
        PetscViewerDestroy(&viewer);
      }
    }

    if (debug->matrices & DEBUG_MATRICES_SNG) {
      sprintf(filename, "%s-A.sng", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
      mat2sng(fino.A, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 0, viewer);
      petsc_call(PetscViewerDestroy(&viewer));
    }    

    if (debug->matrices & DEBUG_MATRICES_SNG_STRUCT) {
      sprintf(filename, "%s-str-A.sng", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT);
      mat2sng(fino.A, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 1, viewer);
      petsc_call(PetscViewerDestroy(&viewer));
    }   
    
    free(filename);

  }

  if (debug->matrices & DEBUG_MATRICES_X) {

    PetscDraw draw;
    PetscDraw draw2;

    petsc_call(PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL, "A", 100, 100, size, size, &viewer));
    MatView(fino.A, viewer);
    PetscViewerDrawGetDraw(viewer, 0, &draw);
    PetscDrawSetPause(draw, -1);

    if (fino.math_type == math_eigen) {
      PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL, "B", size+100, size+100, size, size, &viewer2);
      MatView(fino.B, viewer2);
      PetscViewerDrawGetDraw(viewer2, 0, &draw2);
      PetscDrawSetPause(draw2, -1);
    } else {
      PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL, "b", size+100, size+100, size, size, &viewer2);
      VecView(fino.b, viewer2);
      PetscViewerDrawGetDraw(viewer2, 0, &draw2);
      PetscDrawSetPause(draw2, -1);
    }
    
    PetscDrawPause(draw);

    PetscViewerDestroy(&viewer);
    PetscViewerDestroy(&viewer2);  
    
  }
  
  if ((int)(wasora_var(wasora_special_var(static_steps))) == 1 || (int)(wasora_var(wasora_special_var(done))) == 1) {
    fino_debug_close(debug);
  }

  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_debug_close"
int fino_debug_close(debug_t *debug) {

  int c;
  FILE *finput;

  if (debug->file != NULL) {
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "PETSc's LogView output\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "----------------------\n"));
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
    petsc_call(PetscLogView(debug->viewer));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
  
    if (debug->include_input) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "transcription of input file\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------------\n\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
    finput = fopen(wasora.argv[1], "r");
      while ((c = fgetc(finput)) != EOF) {
        fputc(c, debug->file->pointer);
      }
      fclose(finput);
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
    }

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n*  *  *  *\n\n"));
  
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "[fino"));
#ifdef PLUGIN_VCS_BRANCH
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s.%d", PLUGIN_VCS_MAJOR, PLUGIN_VCS_REVNO));
#else
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s", PACKAGE_VERSION));
#endif
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "](https://bitbucket.org/gtheler/fino)\n\n"));

    PetscViewerDestroy(&debug->viewer);
  }

  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "fino_print_petsc_vector"
int fino_print_petsc_vector(Vec b, PetscViewer viewer) {

  double xi;
  int i;
  int m;

  VecGetSize(b, &m);

  for (i = 0; i < m; i++) {
    VecGetValues(b, 1, &i, &xi);
    if (xi != 0) {
      petsc_call(PetscViewerASCIIPrintf(viewer, "% .1e ", xi));
    } else {
      petsc_call(PetscViewerASCIIPrintf(viewer, "    0    "));
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_print_petsc_matrix"
int fino_print_petsc_matrix(Mat A, PetscViewer viewer) {

  double xi;
  int i, j;
  int m, n;

  MatGetSize(A, &m, &n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      MatGetValues(A, 1, &i, 1, &j, &xi);
      if (xi != 0) {
        petsc_call(PetscViewerASCIIPrintf(viewer, "% .1e ", xi));
      } else {
        petsc_call(PetscViewerASCIIPrintf(viewer, "    0    "));
      }
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_print_petsc_matrix_struct"
int fino_print_petsc_matrix_struct(Mat A, PetscViewer viewer) {

  double xi;
  int i, j;
  int m, n;

  MatGetSize(A, &m, &n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      MatGetValues(A, 1, &i, 1, &j, &xi);
      if (xi != 0) {
         petsc_call(PetscViewerASCIIPrintf(viewer, "#"));
      } else {
        petsc_call(PetscViewerASCIIPrintf(viewer, " "));
      }
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }

  return WASORA_RUNTIME_OK;
  
}
