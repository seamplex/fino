/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's initialization routines
 *
 *  Copyright (C) 2015-2018 jeremy theler
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
#include <unistd.h>
#include <signal.h>
#include "fino.h"

#define NAME_SIZE 32

#undef  __FUNCT__
#define __FUNCT__ "plugin_init_before_parser"
int plugin_init_before_parser(void) {

  char *dummy;
  int i;
  
  if (sizeof(PetscReal) != sizeof(double)) {
    wasora_push_error_message("PETSc should be compiled with double-precision real scalar types");
    return WASORA_PARSER_ERROR;
  }
  
  // amasamos la linea de comandos original (porque la que saca getopt puede tener un orden que no nos sirve)
  // el chiste es que por ejemplo "-log_summary" es atrapado por el getopt de wasora como "-l"
  // hay que re-escrbir eso como "--slepc_opt log_summary"
  // si alguna opcion tiene argumento hay que ponerlo como "--slepc_opt pc_type=sor"
  for (i = 0; i < wasora.argc_orig; i++) {
    if (strcmp(wasora.argv_orig[i], "--petsc") == 0) {
      if (i >= (wasora.argc_orig-1)) {
        wasora_push_error_message("commandline option --petsc needs an argument");
        return WASORA_PARSER_ERROR;
      } else if (wasora.argv_orig[i+1][0] == '-') {
        wasora_push_error_message("the argument of commandline option --petsc should not start with a dash (it is added automatically)");
        return WASORA_PARSER_ERROR;
      }
      
      if ((dummy = strchr(wasora.argv_orig[i+1], '=')) != NULL)  {
        char *tmp1, *tmp2;
        *dummy = '\0';
        tmp1 = strdup(wasora.argv_orig[i+1]);
        tmp2 = strdup(dummy+1);
        wasora.argv_orig[i]   = realloc(wasora.argv_orig[i],   strlen(wasora.argv_orig[i+1])+2);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(dummy)+1);
        sprintf(wasora.argv_orig[i],  "-%s", tmp1);
        sprintf(wasora.argv_orig[i+1], "%s", tmp2);
        free(tmp1);
        free(tmp2);
        
      } else {
        char *tmp1;
        tmp1 = strdup(wasora.argv_orig[i+1]);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(tmp1)+1);
        wasora.argv_orig[i][0] = '\0';
        sprintf(wasora.argv_orig[i+1],  "-%s", tmp1);
        free(tmp1);
      }
      i++;
    }
  }

#ifdef HAVE_SLEPC  
  // inicializamos la slepc (que a su vez inicializa la petsc)
  // le pasamos la linea de comandos que acabamos de amasar
  petsc_call(SlepcInitialize(&wasora.argc_orig, &wasora.argv_orig, (char*)0, PETSC_NULL));
#else
  // inicializamos la petsc
  // le pasamos la linea de comandos que acabamos de amasar
  petsc_call(PetscInitialize(&wasora.argc_orig, &wasora.argv_orig, (char*)0, PETSC_NULL));
#endif
  fino.petscinit_called = 1;
  
  // los segfaults son segfaults, no queremos que la petsc meta las narices
  signal(SIGSEGV, SIG_DFL);

  // esto lo vamos a usar despues cuando hagamos el chiste en paralelo
  MPI_Comm_rank(PETSC_COMM_WORLD, &fino.rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &fino.size);

  // instalamos nuestro error handler para errores la petsc 
  PetscPushErrorHandler(&fino_handler, NULL);

  // inicializamos mesh
  if (!wasora_mesh.initialized) {
    wasora_call(wasora_mesh_init_before_parser());
  }

  
  // variables especiales de fino
///va+fino_abstol+name fino_abstol
///va+fino_abstol+desc Absolute tolerance of the linear solver,
///va+fino_abstol+desc as passed to PETSc’s
///va+fino_abstol+desc [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html)
  fino.vars.abstol = wasora_define_variable("fino_abstol");
  // TODO: poner el default automaticamente
///va+fino_abstol+desc Default `1e-50`.
  wasora_var(fino.vars.abstol) = 1e-50;   // igual al de PETSc
 
///va+fino_reltol+name fino_reltol
///va+fino_reltol+desc Relative tolerance of the linear solver,
///va+fino_reltol+desc as passed to PETSc’s
///va+fino_reltol+desc [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
fino.vars.reltol = wasora_define_variable("fino_reltol");
///va+fino_reltol+desc Default `1e-6`.
  wasora_var(fino.vars.reltol) = 1e-6;    // el de PETSc es 1e-5
  
///va+fino_divtol+name fino_divtol
///va+fino_divtol+desc Divergence tolerance,
///va+fino_divtol+desc as passed to PETSc’s
///va+fino_divtol+desc [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
  fino.vars.divtol = wasora_define_variable("fino_divtol");
///va+fino_divtol+desc Default `1e+4`.  
  wasora_var(fino.vars.divtol) = 1e+4;  // igual al de PETSc
  
///va+fino_max_iterations+name fino_max_iterations
///va+fino_max_iterations+desc Number of maximum iterations before diverging,
///va+fino_max_iterations+desc as passed to PETSc’s
///va+fino_max_iterations+desc [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
  fino.vars.max_iterations = wasora_define_variable("fino_max_iterations");
///va+fino_max_iterations+desc Default `10000`.
  wasora_var(fino.vars.max_iterations) = 10000;   // igual al de PETSc

///va+fino_gamg_threshold+name fino_gamg_threshold
///va+fino_gamg_threshold+desc Relative threshold to use for dropping edges in aggregation graph for the
///va+fino_gamg_threshold+desc [Geometric Algebraic Multigrid Preconditioner](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html)
///va+fino_gamg_threshold+desc as passed to PETSc’s
///va+fino_gamg_threshold+desc [`PCGAMGSetThreshold`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMGSetThreshold.html).
///va+fino_gamg_threshold+desc A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
  fino.vars.gamg_threshold = wasora_define_variable("fino_gamg_threshold");
///va+fino_gamg_threshold+desc Default `0.01`.  
  wasora_var(fino.vars.gamg_threshold) = 0.01;
  
//va+fino_dirichlet_diagonal+name fino_dirichlet_diagonal
//va+fino_dirichlet_diagonal+desc Value that is inserted in the diagonal of the rows
//va+fino_dirichlet_diagonal+desc that correspond to Dirichlet boundary conditions.
//va+fino_dirichlet_diagonal+desc Default is one, but PETSc internally scales it up
//va+fino_dirichlet_diagonal+desc automatically to keep a good condition number.
//  fino.vars.dirichlet_diagonal = wasora_define_variable("fino_dirichlet_diagonal");
  
///va+fino_penalty_weight+name fino_penalty_weight
///va+fino_penalty_weight+desc The weight $w$ used when setting multi-freedom boundary conditions.
///va+fino_penalty_weight+desc Higher values mean better precision in the constrain but distort
///va+fino_penalty_weight+desc the matrix condition number. 
  fino.vars.penalty_weight = wasora_define_variable("fino_penalty_weight");
///va+fino_penalty_weight+desc Default is `1e7`.
  wasora_var(fino.vars.penalty_weight) = 1e7;  
  
///va+fino_iterations+name fino_iterations
///va+fino_iterations+desc This variable contains the actual number of iterations used
///va+fino_iterations+desc by the solver. It is set after `FINO_STEP`.
  fino.vars.iterations = wasora_define_variable("fino_iterations");
  
///va+fino_residual_norm+name fino_residual_norm
///va+fino_residual_norm+desc This variable contains the residual obtained
///va+fino_residual_norm+desc by the solver. It is set after `FINO_STEP`.
  fino.vars.residual_norm= wasora_define_variable("fino_residual_norm");

  // estas son para las expresiones algebraicas implicitamente
  // las definimos en mayusculas porque ya hay funciones que se llaman asi en minuscula
  // antes de parsear la expresion algebraica les cambiamos el case en bc.
  fino.vars.U[0]= wasora_define_variable("U");
  fino.vars.U[1]= wasora_define_variable("V");
  fino.vars.U[2]= wasora_define_variable("W");

  // estas se las dejamos para las condiciones de contorno
  fino.vars.nx= wasora_define_variable("nx");
  fino.vars.ny= wasora_define_variable("ny");
  fino.vars.nz= wasora_define_variable("nz");
  
///va+displ_max+name displ_max
///va+displ_max+desc The module of the maximum displacement of the elastic problem.
  fino.vars.displ_max = wasora_define_variable("displ_max");

///va+displ_max_x+name displ_max_x
///va+displ_max_x+desc The\ $x$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_x = wasora_define_variable("displ_max_x");
///va+displ_max_y+name displ_max_y
///va+displ_max_y+desc The\ $y$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_y = wasora_define_variable("displ_max_y");
///va+displ_max_z+name displ_max_z
///va+displ_max_z+desc The\ $z$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_z = wasora_define_variable("displ_max_z");

///va+u_at_displ_max+name u_at_displ_max
///va+u_at_displ_max+desc The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.
  fino.vars.u_at_displ_max = wasora_define_variable("u_at_displ_max");
///va+v_at_displ_max+name v_at_displ_max
///va+v_at_displ_max+desc The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.
  fino.vars.v_at_displ_max = wasora_define_variable("v_at_displ_max");
///va+w_at_displ_max+name w_at_displ_max
///va+w_at_displ_max+desc The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.
  fino.vars.w_at_displ_max = wasora_define_variable("w_at_displ_max");
  
///va+sigma_max+name sigma_max
///va+sigma_max+desc The maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max = wasora_define_variable("sigma_max");

///va+sigma_max_x+name sigma_max_x
///va+sigma_max_x+desc The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_x = wasora_define_variable("sigma_max_x");
///va+sigma_max_y+name sigma_max_y
///va+sigma_max_y+desc The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_y = wasora_define_variable("sigma_max_y");
///va+sigma_max_z+name sigma_max_z
///va+sigma_max_z+desc The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_z = wasora_define_variable("sigma_max_z");
  
///va+u_at_sigma_max+name u_at_sigma_max
///va+u_at_sigma_max+desc The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.u_at_sigma_max = wasora_define_variable("u_at_sigma_max");
///va+v_at_sigma_max+name v_at_sigma_max
///va+v_at_sigma_max+desc The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.v_at_sigma_max = wasora_define_variable("v_at_sigma_max");
///va+w_at_sigma_max+name w_at_sigma_max
///va+w_at_sigma_max+desc The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.w_at_sigma_max = wasora_define_variable("w_at_sigma_max");
  

///va+T_max+name T_max
///va+T_max+desc The maximum temperature\ $T_\text{max}$ of the thermal problem.
  fino.vars.T_max = wasora_define_variable("T_max");

///va+T_min+name T_min
///va+T_min+desc The minimum temperature\ $T_\text{min}$ of the thermal problem.
  fino.vars.T_min = wasora_define_variable("T_min");

  // variables internas
///va+lambda+name lambda
///va+lambda+desc 
///va+lambda+desc Requested eigenvalue. It is equal to 1.0 until
///va+lambda+desc `FINO_STEP` is executed.  
  fino.vars.lambda = wasora_define_variable("lambda");
  wasora_var(fino.vars.lambda) = 1.0;
  
///va+time_wall_build+name time_wall_build
///va+time_wall_build+desc Wall time insumed to build the problem matrices, in seconds.
  fino.vars.time_wall_build = wasora_define_variable("time_wall_build");

///va+time_wall_solve+name time_wall_solve
///va+time_wall_solve+desc Wall time insumed to solve the eigen-problem, in seconds.
  fino.vars.time_wall_solve = wasora_define_variable("time_wall_solve");

///va+time_wall_total+name time_wall_total
///va+time_wall_total+desc Wall time insumed to initialize, build and solve, in seconds.
  fino.vars.time_wall_total = wasora_define_variable("time_wall_total");
  
///va+time_cpu_build+name time_cpu_build
///va+time_cpu_build+desc CPU time insumed to build the problem matrices, in seconds.
  fino.vars.time_cpu_build = wasora_define_variable("time_cpu_build");

///va+time_cpu_solve+name time_cpu_solve
///va+time_cpu_solve+desc CPU time insumed to solve the eigen-problem, in seconds.
  fino.vars.time_cpu_solve = wasora_define_variable("time_cpu_solve");

///va+time_wall_total+name time_cpu_total
///va+time_wall_total+desc CPU time insumed to initialize, build and solve, in seconds.
  fino.vars.time_cpu_total = wasora_define_variable("time_cpu_total");
  
///va+time_petsc_build+name time_petsc_build
///va+time_petsc_build+desc CPU time insumed by PETSc to build the problem matrices, in seconds.
  fino.vars.time_petsc_build = wasora_define_variable("time_petsc_build");

///va+time_petsc_solve+name time_petsc_solve
///va+time_petsc_solve+desc CPU time insumed by PETSc to solve the eigen-problem, in seconds.
  fino.vars.time_petsc_solve = wasora_define_variable("time_petsc_solve");

///va+time_wall_total+name time_wall_total
///va+time_wall_total+desc CPU time insumed by PETSc to initialize, build and solve, in seconds.
  fino.vars.time_petsc_total = wasora_define_variable("time_petsc_total");

  ///va+petsc_flops+name petsc_flops
///va+petsc_flops+desc Number of floating point operations performed by PETSc/SLEPc.
  fino.vars.flops_petsc = wasora_define_variable("flops_petsc");
         
///va+memory_use+name available_memory
///va+memory_use+desc Total available memory, in bytes.
  fino.vars.available_memory = wasora_define_variable("available_memory");
  wasora_value(fino.vars.available_memory) = sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE);

///va+memory_usage_global+name global_memory_use
///va+memory_usage_global+desc Maximum resident set size (global memory used), in bytes.
  fino.vars.memory_usage_global = wasora_define_variable("memory_usage_global");
  
///va+memory_usage_petsc+name petsc_memory_use
///va+memory_usage_petsc+desc Maximum resident set size (memory used by PETSc), in bytes.
  fino.vars.memory_usage_petsc = wasora_define_variable("memory_usage_petsc");

/*
 * los sacamos por ahora habria que ver de definirlas solo si lo necesitamos  
  // vector con las funciones de forma
  fino.h = wasora_define_vector("h", 0, NULL, NULL);

  // jacobianos
  // eventualmente despues hacemos un free y lo definimos de nuevo
  fino.dhdx = wasora_define_matrix("dhdx", 0, NULL, 0, NULL, NULL);

  // nombres de objetos por default
  // TODO: ponerlo en otro lado poque si no se necesitan se usan nombres copados
  // que podrian estar disponibles para otra cosa (por ejemplo coef. de heat transfer)
  fino.shape_name = strdup("h");
  fino.lhs_matrix_name = strdup("A");
  fino.rhs_matrix_name = strdup("B");
  fino.rhs_vector_name = strdup("b");
*/
  return WASORA_PARSER_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "plugin_init_after_parser"
int plugin_init_after_parser(void) {

  int m, g;
  
  wasora_call(fino_read_bcs());  
  
  // desplazamientos (y derivadas) anteriores
  if (fino.problem_family == problem_family_break) {
    fino.base_solution = calloc(fino.degrees, sizeof(function_t *));
    fino.base_gradient = calloc(fino.degrees, sizeof(function_t **));
    
    if (fino.dimensions == 3) {
      
      fino.base_solution[0] = wasora_get_function_ptr("u0");
      fino.base_solution[1] = wasora_get_function_ptr("v0");
      fino.base_solution[2] = wasora_get_function_ptr("w0");

      fino.base_gradient[0] = calloc(fino.dimensions, sizeof(function_t *));
      fino.base_gradient[0][0] = wasora_get_function_ptr("du0dx");
      fino.base_gradient[0][1] = wasora_get_function_ptr("du0dy");
      fino.base_gradient[0][2] = wasora_get_function_ptr("du0dz");
    
      fino.base_gradient[1] = calloc(fino.dimensions, sizeof(function_t *));
      fino.base_gradient[1][0] = wasora_get_function_ptr("dv0dx");
      fino.base_gradient[1][1] = wasora_get_function_ptr("dv0dy");
      fino.base_gradient[1][2] = wasora_get_function_ptr("dv0dz");

      fino.base_gradient[2] = calloc(fino.dimensions, sizeof(function_t *));
      fino.base_gradient[2][0] = wasora_get_function_ptr("dw0dx");
      fino.base_gradient[2][1] = wasora_get_function_ptr("dw0dy");
      fino.base_gradient[2][2] = wasora_get_function_ptr("dw0dz");
      
    } else if (fino.dimensions == 2) {
      
      fino.base_solution[0] = wasora_get_function_ptr("u0");
      fino.base_solution[1] = wasora_get_function_ptr("v0");

      fino.base_gradient[0] = calloc(fino.dimensions, sizeof(function_t *));
      fino.base_gradient[0][0] = wasora_get_function_ptr("du0dx");
      fino.base_gradient[0][1] = wasora_get_function_ptr("du0dy");
    
      fino.base_gradient[1] = calloc(fino.dimensions, sizeof(function_t *));
      fino.base_gradient[1][0] = wasora_get_function_ptr("dv0dx");
      fino.base_gradient[1][1] = wasora_get_function_ptr("dv0dy");
      
    }
    

    for (g = 0; g < fino.degrees; g++) {
      if (fino.base_solution[g] != NULL && fino.base_solution[g]->n_arguments != fino.dimensions) {
        wasora_push_error_message("function '%s' should have %d arguments instead of %d", fino.base_solution[g]->name, fino.degrees, fino.base_solution[g]->n_arguments);
        return WASORA_PARSER_ERROR;
      }
      
      for (m = 0; m < fino.dimensions; m++) {
        if (fino.base_gradient[g][m] != NULL && fino.base_gradient[g][m]->n_arguments != fino.dimensions) {
          wasora_push_error_message("function '%s' should have %d arguments instead of %d", fino.base_gradient[g][m]->name, fino.dimensions, fino.base_gradient[g][m]->n_arguments);
          return WASORA_PARSER_ERROR;
        }
      }
      
    }
  }
  
  // los objetos para mostrar el progress
  if (fino.shmem_progress_build_name != NULL) {
    fino.shmem_progress_build = wasora_get_shared_pointer(fino.shmem_progress_build_name, sizeof(double));
  }
  if (fino.shmem_progress_solve_name != NULL) {
    fino.shmem_progress_solve = wasora_get_shared_pointer(fino.shmem_progress_solve_name, sizeof(double));
  }
  if (fino.shmem_memory_name != NULL) {
    fino.shmem_memory = wasora_get_shared_pointer(fino.shmem_memory_name, sizeof(double));
  }  

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "plugin_init_before_run"
int plugin_init_before_run(void) {

  fino.problem_size = 0;
  fino.spatial_unknowns = 0;

  wasora_call(fino_problem_free());
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "plugin_finalize"
int plugin_finalize(void) {

  wasora_call(fino_problem_free());

  if (fino.petscinit_called) {
#ifdef HAVE_SLEPC  
    petsc_call(SlepcFinalize());
#else
    petsc_call(PetscFinalize());
#endif
  }
  
  // los objetos para mostrar el progress
  if (fino.shmem_progress_build_name != NULL) {
    wasora_free_shared_pointer(fino.shmem_progress_build, fino.shmem_progress_build_name, sizeof(double));
  }
  if (fino.shmem_progress_solve_name != NULL) {
    wasora_free_shared_pointer(fino.shmem_progress_solve, fino.shmem_progress_solve_name, sizeof(double));
  }
  if (fino.shmem_memory_name != NULL) {
    wasora_free_shared_pointer(fino.shmem_memory, fino.shmem_memory_name, sizeof(double));
  }  
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_problem_init"
// esto viene despues de haber leido la malla
int fino_problem_init(void) {

  int i, g;
  int width;
  
  physical_entity_t *physical_entity;

//---------------------------------
// inicializamos parametros
//---------------------------------

  if ((fino.mesh = wasora_mesh.meshes) == NULL) {
    wasora_push_error_message("no mesh defined");
    return WASORA_RUNTIME_ERROR;
  }

  LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
    if (physical_entity->bc_type_math != bc_math_undefined && physical_entity->n_elements == 0) {
      wasora_push_error_message("physical entity '%s' (id %d) has a BC but no associated elements", physical_entity->name, physical_entity->id);
      return WASORA_RUNTIME_ERROR;
    }
    if (physical_entity->material != NULL && physical_entity->n_elements == 0) {
      wasora_push_error_message("physical entity '%s' (id %d) has a material but no associated elements", physical_entity->name, physical_entity->id);
      return WASORA_RUNTIME_ERROR;
    }
  }

  
  
  // ponemos esto para hacer explicito que somos FEM y no FVM
  fino.spatial_unknowns = fino.mesh->n_nodes;
  fino.mesh->data_type = data_type_node;
  fino.problem_size = fino.spatial_unknowns * fino.degrees;
  

//---------------------------------
// alocamos objetos globales
//---------------------------------

  width = fino.mesh->max_first_neighbor_nodes * fino.degrees;
  
  // la matriz de stiffnes global
  petsc_call(MatCreate(PETSC_COMM_WORLD, &fino.K));
  petsc_call(MatSetSizes(fino.K, PETSC_DECIDE, PETSC_DECIDE, fino.problem_size, fino.problem_size));
  petsc_call(MatSetFromOptions(fino.K));
  petsc_call(MatMPIAIJSetPreallocation(fino.K, width, PETSC_NULL, width, PETSC_NULL));
  petsc_call(MatSeqAIJSetPreallocation(fino.K, width, PETSC_NULL));
  petsc_call(MatSetOption(fino.K, MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE));
  if (fino.do_not_set_block_size == 0) {
    petsc_call(MatSetBlockSize(fino.K, fino.degrees));
  }
  petsc_call(PetscObjectSetName((PetscObject)fino.K, "K"));
  
  // el vector incognita
  petsc_call(MatCreateVecs(fino.K, NULL, &fino.phi));
  petsc_call(PetscObjectSetName((PetscObject)fino.phi, "phi"));

  if (fino.math_type == math_type_linear) {
    // el vector del miembro derecho
    fino.has_rhs = 1;
    petsc_call(MatCreateVecs(fino.K, NULL, &fino.b));
    petsc_call(PetscObjectSetName((PetscObject)fino.b, "b"));
  }
  
  if (fino.problem_family == problem_family_shake ||
      (fino.problem_family == problem_family_bake && wasora_var_value(wasora_special_var(end_time)) != 0)) {
    // la matriz de masa para autovalores del problema elastico o para 
    fino.has_mass = 1;
    petsc_call(MatCreate(PETSC_COMM_WORLD, &fino.M));
    petsc_call(MatSetSizes(fino.M, PETSC_DECIDE, PETSC_DECIDE, fino.problem_size, fino.problem_size));
    petsc_call(MatSetFromOptions(fino.M));
    petsc_call(MatMPIAIJSetPreallocation(fino.M, width, PETSC_NULL, width, PETSC_NULL));
    petsc_call(MatSeqAIJSetPreallocation(fino.M, width, PETSC_NULL));
    if (fino.do_not_set_block_size == 0) {
      petsc_call(MatSetBlockSize(fino.M, fino.degrees));
    }
    
    petsc_call(PetscObjectSetName((PetscObject)fino.M, "M"));
  }
  
  if (fino.mesh->structured) {
    wasora_mesh_struct_init_rectangular_for_nodes(fino.mesh);
  }

  // rellenamos holders las funciones continuas que van a tener la solucion
  for (g = 0; g < fino.degrees; g++) {
    fino.solution[g]->data_size = fino.spatial_unknowns;
    fino.solution[g]->data_argument = fino.mesh->nodes_argument;
    fino.solution[g]->data_value = calloc(fino.spatial_unknowns, sizeof(double));
    
    if (fino.nev > 1) {
      for (i = 0; i < fino.nev; i++) {
        fino.vibration[g][i]->data_argument = fino.gradient[0][0]->data_argument;
        fino.vibration[g][i]->data_size = fino.mesh->n_nodes;
        fino.vibration[g][i]->data_value = calloc(fino.mesh->n_nodes, sizeof(double));
      }
    }
  }

  wasora_call(mesh_node_indexes(fino.mesh, fino.degrees));
  
  return WASORA_PARSER_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_problem_free"
int fino_problem_free(void) {
  int g, d;

  if (fino.mesh != NULL && fino.mesh->n_elements != 0) {
    for (d = 0; d < fino.dimensions; d++) {
      if (fino.gradient[0][0]->data_argument != NULL) {
        free(fino.gradient[0][0]->data_argument[d]);
      }
    }
    free(fino.gradient[0][0]->data_argument);
    fino.gradient[0][0]->data_argument = NULL;

    
    for (g = 0; g < fino.degrees; g++) {
      for (d = 0; d < fino.dimensions; d++) {
        free(fino.gradient[g][d]->data_value);
        fino.gradient[g][d]->data_value = NULL;
      }
      
      free(fino.solution[g]->data_value);
      fino.solution[g]->data_value = NULL;
    }
    mesh_free(fino.mesh);
  }
  
  if (fino.problem_family == problem_family_break) {
    fino_function_clean_nodal_data(fino.sigma1);
    fino_function_clean_nodal_data(fino.sigma2);
    fino_function_clean_nodal_data(fino.sigma3);
    fino_function_clean_nodal_data(fino.sigma);
    fino_function_clean_nodal_data(fino.tresca);
  }
  
     
  if (fino.phi != PETSC_NULL) {
    petsc_call(VecDestroy(&fino.phi));
  }
  if (fino.K != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.K));
  }
  if (fino.M != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.M));
  }
  if (fino.b != PETSC_NULL) {
    petsc_call(VecDestroy(&fino.b));
  }
  if (fino.ksp != PETSC_NULL) {
    petsc_call(KSPDestroy(&fino.ksp));
  }
#ifdef HAVE_SLEPC  
  if (fino.eps != PETSC_NULL) {
//    petsc_call(EPSDestroy(&fino.eps));
  }
#endif
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_function_clean_nodal_data"
int fino_function_clean_nodal_data(function_t *function) {
 
  if (function->data_value != NULL) {  
    free(function->data_value);
    function->data_value = NULL;
  }
  
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_function_clean_nodal_data"
int fino_function_clean_nodal_arguments(function_t *function) {
 
  int d;

  if (function->data_argument != NULL) {
    for (d = 0; d < fino.dimensions; d++) {
      free(function->data_argument[d]);
    }
    free(function->data_argument);
  }
  
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_define_result_function"
int fino_define_result_function(char *name, function_t **function) {

  if ((*function = wasora_define_function(name, fino.dimensions)) == NULL) {
    wasora_push_error_message("result function '%s' defined twice", name);
    return WASORA_RUNTIME_ERROR;
  }
  (*function)->mesh = fino.mesh;
  fino_function_clean_nodal_arguments(*function);
  (*function)->var_argument = fino.solution[0]->var_argument;
  (*function)->type = type_pointwise_mesh_node;

  return 0;
}

