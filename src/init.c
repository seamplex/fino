/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's initialization routines
 *
 *  Copyright (C) 2015--2020 Seamplex
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
#include <unistd.h>
#include <signal.h>
#include "fino.h"

// this is a wrapper because PetscOptionsHasName change its arguments after 3.7.0
#if PETSC_VERSION_LT(3,7,0)
 #define PetscOptionsHasNameWrapper(a, b, c) PetscOptionsHasName(a, b, c)
#else
 #define PetscOptionsHasNameWrapper(a, b, c) PetscOptionsHasName(PETSC_NULL, a, b, c)
#endif


#define NAME_SIZE 32

int plugin_init_before_parser(void) {

  char *dummy;
  int i;
  
  if (sizeof(PetscReal) != sizeof(double)) {
    wasora_push_error_message("PETSc should be compiled with double-precision real scalar types");
    return WASORA_PARSER_ERROR;
  }
  
  // we process the original command line (beucase the one that remains after getopt might have a different order)
  // for instance, "-log_summary" is trapped by wasora's getoopt as "-l"
  // so that needs to be rewritten as "--petsc_opt log_summary"
  // if the option has an argument it has to be put as "--slepc_opt pc_type=sor"
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
        wasora.argv_orig[i]   = realloc(wasora.argv_orig[i],   strlen(wasora.argv_orig[i+1])+6);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(dummy+1)+6);
        snprintf(wasora.argv_orig[i], strlen(wasora.argv_orig[i+1])+4, "-%s", tmp1);
        snprintf(wasora.argv_orig[i+1], strlen(dummy+1)+4, "%s", tmp2);
        free(tmp1);
        free(tmp2);
        
      } else {
        char *tmp1;
        tmp1 = strdup(wasora.argv_orig[i+1]);
        wasora.argv_orig[i+1] = realloc(wasora.argv_orig[i+1], strlen(tmp1)+6);
        wasora.argv_orig[i][0] = '\0';
        snprintf(wasora.argv_orig[i+1], strlen(tmp1)+4, "-%s", tmp1);
        free(tmp1);
      }
      i++;
    }
  }
  
#ifdef HAVE_SLEPC  
  // initialize SLEPc (which in turn initalizes PETSc)
  // we pass the processed command line
  petsc_call(SlepcInitialize(&wasora.argc_orig, &wasora.argv_orig, (char*)0, PETSC_NULL));
#else
  // initialize PETSc
  // we pass the processed command line
  petsc_call(PetscInitialize(&wasora.argc_orig, &wasora.argv_orig, (char*)0, PETSC_NULL));
#endif
  fino.petscinit_called = 1;
  
  // segfaults are segfaults, try to leave PETSC out of them
  signal(SIGSEGV, SIG_DFL);

  // get the number of processes and the rank
  petsc_call(MPI_Comm_size(PETSC_COMM_WORLD, &wasora.nprocs));
  petsc_call(MPI_Comm_rank(MPI_COMM_WORLD, &wasora.rank));

  // install out error handler for PETSc
  petsc_call(PetscPushErrorHandler(&fino_handler, NULL));

  // register events
  petsc_call(PetscClassIdRegister("Fino", &fino.petsc_classid));

  petsc_call(PetscLogStageRegister("Assembly", &fino.petsc_stage_build));
  petsc_call(PetscLogStageRegister("Solution", &fino.petsc_stage_solve));
  petsc_call(PetscLogStageRegister("Stress", &fino.petsc_stage_solve));
  
  petsc_call(PetscLogEventRegister("fino_build", fino.petsc_classid, &fino.petsc_event_build));
  petsc_call(PetscLogEventRegister("fino_solve", fino.petsc_classid, &fino.petsc_event_solve));
  petsc_call(PetscLogEventRegister("fino_stress", fino.petsc_classid, &fino.petsc_event_solve));

  // initialize wasora's mesh framework
  if (!wasora_mesh.initialized) {
    wasora_call(wasora_mesh_init_before_parser());
  }

  
  // Fino's special variables
///va+fino_abstol+name fino_abstol
///va+fino_abstol+detail Absolute tolerance of the linear solver,
///va+fino_abstol+detail as passed to PETSc’s
///va+fino_abstol+detail [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html)
  fino.vars.abstol = wasora_define_variable("fino_abstol");
  // TODO: poner el default automaticamente
///va+fino_abstol+detail Default `1e-50`.
  wasora_var(fino.vars.abstol) = 1e-50;   // igual al de PETSc
 
///va+fino_reltol+name fino_reltol
///va+fino_reltol+detail Relative tolerance of the linear solver,
///va+fino_reltol+detail as passed to PETSc’s
///va+fino_reltol+detail [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
fino.vars.reltol = wasora_define_variable("fino_reltol");
///va+fino_reltol+detail Default `1e-6`.
  wasora_var(fino.vars.reltol) = 1e-6;    // el de PETSc es 1e-5
  
///va+fino_divtol+name fino_divtol
///va+fino_divtol+detail Divergence tolerance,
///va+fino_divtol+detail as passed to PETSc’s
///va+fino_divtol+detail [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
  fino.vars.divtol = wasora_define_variable("fino_divtol");
///va+fino_divtol+detail Default `1e+4`.  
  wasora_var(fino.vars.divtol) = 1e+4;  // igual al de PETSc
  
///va+fino_max_iterations+name fino_max_iterations
///va+fino_max_iterations+detail Number of maximum iterations before diverging,
///va+fino_max_iterations+detail as passed to PETSc’s
///va+fino_max_iterations+detail [`KSPSetTolerances`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html).
  fino.vars.max_iterations = wasora_define_variable("fino_max_iterations");
///va+fino_max_iterations+detail Default `10000`.
  wasora_var(fino.vars.max_iterations) = 10000;   // igual al de PETSc

///va+fino_gamg_threshold+name fino_gamg_threshold
///va+fino_gamg_threshold+detail Relative threshold to use for dropping edges in aggregation graph for the
///va+fino_gamg_threshold+detail [Geometric Algebraic Multigrid Preconditioner](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html)
///va+fino_gamg_threshold+detail as passed to PETSc’s
///va+fino_gamg_threshold+detail [`PCGAMGSetThreshold`](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMGSetThreshold.html).
///va+fino_gamg_threshold+detail A value of 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph.
  fino.vars.gamg_threshold = wasora_define_variable("fino_gamg_threshold");
///va+fino_gamg_threshold+detail Default `0.01`.  
  wasora_var(fino.vars.gamg_threshold) = 0.01;
  
///va+fino_penalty_weight+name fino_penalty_weight
///va+fino_penalty_weight+detail The weight $w$ used when setting multi-freedom boundary conditions.
///va+fino_penalty_weight+detail Higher values mean better precision in the constrain but distort
///va+fino_penalty_weight+detail the matrix condition number. 
  fino.vars.penalty_weight = wasora_define_variable("fino_penalty_weight");
///va+fino_penalty_weight+detail Default is `1e8`.
  wasora_var(fino.vars.penalty_weight) = 1e8;  
  
///va+fino_iterations+name fino_iterations
///va+fino_iterations+detail This variable contains the actual number of iterations used
///va+fino_iterations+detail by the solver. It is set after `FINO_STEP`.
  fino.vars.iterations = wasora_define_variable("fino_iterations");
  
///va+fino_residual_norm+name fino_residual_norm
///va+fino_residual_norm+detail This variable contains the residual obtained
///va+fino_residual_norm+detail by the solver. It is set after `FINO_STEP`.
  fino.vars.residual_norm= wasora_define_variable("fino_residual_norm");

///va+nodes_rough+name nodes_rough
///va+nodes_rough+detail The number of nodes of the mesh in `ROUGH` mode.
  fino.vars.nodes_rough = wasora_define_variable("nodes_rough");
  
  // these are for the algebraic expressions in the  which are implicitly-defined BCs
  // i.e. 0=u*nx+v*ny
  // here they are defined as uppercase because there already exist functions named u, v and w
  // but the parser changes their case when an implicit BC is read
  fino.vars.U[0]= wasora_define_variable("U");
  fino.vars.U[1]= wasora_define_variable("V");
  fino.vars.U[2]= wasora_define_variable("W");

///va+strain_energy+name strain_energy
///va+strain_energy+detail The strain energy stored in the solid, computed as
///va+strain_energy+detail $1/2 \cdot \vec{u}^T  K \vec{u}$
///va+strain_energy+detail where $\vec{u}$ is the displacements vector and $K$ is the stiffness matrix.
  fino.vars.strain_energy = wasora_define_variable("strain_energy");

  ///va+displ_max+name displ_max
///va+displ_max+detail The module of the maximum displacement of the elastic problem.
  fino.vars.displ_max = wasora_define_variable("displ_max");

///va+displ_max_x+name displ_max_x
///va+displ_max_x+detail The\ $x$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_x = wasora_define_variable("displ_max_x");
///va+displ_max_y+name displ_max_y
///va+displ_max_y+detail The\ $y$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_y = wasora_define_variable("displ_max_y");
///va+displ_max_z+name displ_max_z
///va+displ_max_z+detail The\ $z$ coordinate of the maximum displacement of the elastic problem.
  fino.vars.displ_max_z = wasora_define_variable("displ_max_z");

///va+u_at_displ_max+name u_at_displ_max
///va+u_at_displ_max+detail The\ $x$ component\ $u$ of the maximum displacement of the elastic problem.
  fino.vars.u_at_displ_max = wasora_define_variable("u_at_displ_max");
///va+v_at_displ_max+name v_at_displ_max
///va+v_at_displ_max+detail The\ $y$ component\ $v$ of the maximum displacement of the elastic problem.
  fino.vars.v_at_displ_max = wasora_define_variable("v_at_displ_max");
///va+w_at_displ_max+name w_at_displ_max
///va+w_at_displ_max+detail The\ $z$ component\ $w$ of the maximum displacement of the elastic problem.
  fino.vars.w_at_displ_max = wasora_define_variable("w_at_displ_max");
  
///va+sigma_max+name sigma_max
///va+sigma_max+detail The maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max = wasora_define_variable("sigma_max");

///va+delta_sigma_max+name delta_sigma_max
///va+delta_sigma_max+detail The uncertainty of the maximum Von Mises stress\ $\sigma$ of the elastic problem.
///va+delta_sigma_max+detail Not to be confused with the maximum uncertainty of the Von Mises stress.
  fino.vars.delta_sigma_max = wasora_define_variable("delta_sigma_max");
  
///va+sigma_max_x+name sigma_max_x
///va+sigma_max_x+detail The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_x = wasora_define_variable("sigma_max_x");
///va+sigma_max_y+name sigma_max_y
///va+sigma_max_y+detail The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_y = wasora_define_variable("sigma_max_y");
///va+sigma_max_z+name sigma_max_z
///va+sigma_max_z+detail The\ $x$ coordinate of the maximum von Mises stress\ $\sigma$ of the elastic problem.
  fino.vars.sigma_max_z = wasora_define_variable("sigma_max_z");
  
///va+u_at_sigma_max+name u_at_sigma_max
///va+u_at_sigma_max+detail The\ $x$ component\ $u$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.u_at_sigma_max = wasora_define_variable("u_at_sigma_max");
///va+v_at_sigma_max+name v_at_sigma_max
///va+v_at_sigma_max+detail The\ $y$ component\ $v$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.v_at_sigma_max = wasora_define_variable("v_at_sigma_max");
///va+w_at_sigma_max+name w_at_sigma_max
///va+w_at_sigma_max+detail The\ $z$ component\ $w$ of the displacement where the maximum von Mises stress\ $\sigma$ of the elastic problem is located.
  fino.vars.w_at_sigma_max = wasora_define_variable("w_at_sigma_max");
  

///va+T_max+name T_max
///va+T_max+detail The maximum temperature\ $T_\text{max}$ of the thermal problem.
  fino.vars.T_max = wasora_define_variable("T_max");

///va+T_min+name T_min
///va+T_min+detail The minimum temperature\ $T_\text{min}$ of the thermal problem.
  fino.vars.T_min = wasora_define_variable("T_min");

///va+lambda+name lambda
///va+lambda+detail 
///va+lambda+detail Requested eigenvalue. It is equal to 1.0 until
///va+lambda+detail `FINO_STEP` is executed.  
  fino.vars.lambda = wasora_define_variable("lambda");
  wasora_var(fino.vars.lambda) = 1.0;
  
///va+time_wall_build+name time_wall_build
///va+time_wall_build+detail Wall time insumed to build the problem matrices, in seconds.
  fino.vars.time_wall_build = wasora_define_variable("time_wall_build");

///va+time_wall_solve+name time_wall_solve
///va+time_wall_solve+detail Wall time insumed to solve the problem, in seconds.
  fino.vars.time_wall_solve = wasora_define_variable("time_wall_solve");

///va+time_wall_stress+name time_wall_stress
///va+time_wall_stress+detail Wall time insumed to compute the stresses, in seconds.
  fino.vars.time_wall_stress = wasora_define_variable("time_wall_stress");

///va+time_wall_total+name time_wall_total
///va+time_wall_total+detail Wall time insumed to initialize, build and solve, in seconds.
  fino.vars.time_wall_total = wasora_define_variable("time_wall_total");
  
///va+time_cpu_build+name time_cpu_build
///va+time_cpu_build+detail CPU time insumed to build the problem matrices, in seconds.
  fino.vars.time_cpu_build = wasora_define_variable("time_cpu_build");

///va+time_cpu_solve+name time_cpu_solve
///va+time_cpu_solve+detail CPU time insumed to solve the problem, in seconds.
  fino.vars.time_cpu_solve = wasora_define_variable("time_cpu_solve");

///va+time_cpu_stress+name time_cpu_stress
///va+time_cpu_stress+detail CPU time insumed to compute the stresses from the displacements, in seconds.
  fino.vars.time_cpu_stress = wasora_define_variable("time_cpu_stress");
  
///va+time_wall_total+name time_cpu_total
///va+time_wall_total+detail CPU time insumed to initialize, build and solve, in seconds.
  fino.vars.time_cpu_total = wasora_define_variable("time_cpu_total");
  
///va+time_petsc_build+name time_petsc_build
///va+time_petsc_build+detail CPU time insumed by PETSc to build the problem matrices, in seconds.
  fino.vars.time_petsc_build = wasora_define_variable("time_petsc_build");

///va+time_petsc_solve+name time_petsc_solve
///va+time_petsc_solve+detail CPU time insumed by PETSc to solve the eigen-problem, in seconds.
  fino.vars.time_petsc_solve = wasora_define_variable("time_petsc_solve");

///va+time_petsc_stress+name time_petsc_solve
///va+time_petsc_stress+detail CPU time insumed by PETSc to compute the stresses, in seconds.
  fino.vars.time_petsc_stress = wasora_define_variable("time_petsc_stress");
  
///va+time_wall_total+name time_wall_total
///va+time_wall_total+detail CPU time insumed by PETSc to initialize, build and solve, in seconds.
  fino.vars.time_petsc_total = wasora_define_variable("time_petsc_total");

  ///va+petsc_flops+name petsc_flops
///va+petsc_flops+detail Number of floating point operations performed by PETSc/SLEPc.
  fino.vars.flops_petsc = wasora_define_variable("flops_petsc");
         
///va+memory_available+name memory_available
///va+memory_available+detail Total available memory, in bytes.
  fino.vars.memory_available = wasora_define_variable("memory_available");
  wasora_value(fino.vars.memory_available) = sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE);

///va+memory+name memory
///va+memory+detail Maximum resident set size (global memory used), in bytes.
  fino.vars.memory = wasora_define_variable("memory");
  
///va+memory_petsc+name memory_petsc
///va+memory_petsc+detail Maximum resident set size (memory used by PETSc), in bytes.
  fino.vars.memory_petsc = wasora_define_variable("memory_petsc");
  
  
  // empezamos con un valor muy negativo, si nadie lo toca ni calculamos la calidad
//  fino.gradient_quality_threshold = DEFAULT_GRADIENT_JACOBIAN_THRESHOLD;

  return WASORA_PARSER_OK;
}

int plugin_init_after_parser(void) {

  int g;
  
  wasora_call(fino_bc_string2parsed());  
  
  // desplazamientos (y derivadas) anteriores
  if (fino.problem_family == problem_family_mechanical) {
    fino.base_solution = calloc(fino.degrees, sizeof(function_t *));

    fino.base_solution[0] = wasora_get_function_ptr("u0");
    fino.base_solution[1] = wasora_get_function_ptr("v0");
    if (fino.dimensions == 3) {
      fino.base_solution[2] = wasora_get_function_ptr("w0");
    }

    for (g = 0; g < fino.degrees; g++) {
      if (fino.base_solution[g] != NULL && fino.base_solution[g]->n_arguments != fino.dimensions) {
        wasora_push_error_message("function '%s' should have %d arguments instead of %d", fino.base_solution[g]->name, fino.degrees, fino.base_solution[g]->n_arguments);
        return WASORA_PARSER_ERROR;
      }
    }
  }
  
  // si no nos dijeron explicitamente si quieren lineal o no lienal, tratamos de adivinar
  if (fino.math_type == math_type_undefined) {
    fino.math_type = math_type_linear;
  }
  
  
  return WASORA_RUNTIME_OK;
}

int plugin_init_before_run(void) {

  fino.global_size = 0;
  fino.spatial_unknowns = 0;
  fino.progress_r0 = 0;
  fino.already_built = PETSC_FALSE;
  fino.first_build = PETSC_TRUE;

  wasora_call(fino_problem_free());
  
  return WASORA_RUNTIME_OK;
}


int plugin_finalize(void) {

  wasora_call(fino_problem_free());

  if (fino.petscinit_called) {
#ifdef HAVE_SLEPC  
    petsc_call(SlepcFinalize());
#else
    petsc_call(PetscFinalize());
#endif
  }
  
  return WASORA_RUNTIME_OK;
}


// esto viene despues de haber leido la malla
int fino_problem_init(void) {

  int i, g;
  int width;
  PetscBool flag;
  
  physical_entity_t *physical_entity;

//---------------------------------
// read command-line arguments that take precedence over the options in the input file
//---------------------------------
  
  // check for further commandline options
  // see if the user asked for mumps in the command line
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--mumps", &flag));
  if (flag == PETSC_TRUE) {
    fino.ksp_type = strdup("mumps");
    fino.pc_type = strdup("mumps");
  }

  // see if the user asked for progress in the command line
  if (fino.progress_ascii == PETSC_FALSE) {
    petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--progress", &fino.progress_ascii));
  }  

  // see if the user asked for a forced problem type
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--linear", &flag));
  if (flag == PETSC_TRUE) {
    fino.math_type = math_type_linear;
  }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--non-linear", &flag));
  if (flag == PETSC_TRUE) {
    fino.math_type = math_type_nonlinear;
  }
  petsc_call(PetscOptionsHasNameWrapper(PETSC_NULL, "--nonlinear", &flag));
  if (flag == PETSC_TRUE) {
    fino.math_type = math_type_nonlinear;
  }
  
  
//---------------------------------
// initialize parameters
//---------------------------------

  if ((fino.mesh = wasora_mesh.meshes) == NULL) {
    wasora_push_error_message("no mesh defined");
    return WASORA_RUNTIME_ERROR;
  }

  for (physical_entity = fino.mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
/*    
 * TODO: poner una variable para elegir esto
    if (physical_entity->bc_type_math != bc_math_undefined && physical_entity->n_elements == 0) {
      wasora_push_error_message("physical entity '%s' has a BC but no associated elements", physical_entity->name);
      return WASORA_RUNTIME_ERROR;
    }
*/
    if (physical_entity->material != NULL && physical_entity->n_elements == 0) {
      wasora_push_error_message("physical group '%s' has a material but no associated elements", physical_entity->name);
      return WASORA_RUNTIME_ERROR;
    }
  }


  // set this explicitly, we are FEM not FVM
  fino.spatial_unknowns = fino.mesh->n_nodes;
  fino.mesh->data_type = data_type_node;
  fino.global_size = fino.spatial_unknowns * fino.degrees;
  

//---------------------------------
// aloccate global objects
//---------------------------------

  width = GSL_MAX(fino.mesh->max_nodes_per_element, fino.mesh->max_first_neighbor_nodes) * fino.degrees;

  // ask how many local nodes we own
  fino.nodes_local = PETSC_DECIDE;
  petsc_call(PetscSplitOwnership(PETSC_COMM_WORLD, &fino.nodes_local, &fino.mesh->n_nodes));
  fino.size_local = fino.degrees * fino.nodes_local;
  
  // the global stiffnes matrix
  petsc_call(MatCreate(PETSC_COMM_WORLD, &fino.K));
  petsc_call(PetscObjectSetName((PetscObject)fino.K, "K"));
  petsc_call(MatSetSizes(fino.K, fino.size_local, fino.size_local, fino.global_size, fino.global_size));
  petsc_call(MatSetFromOptions(fino.K));
  petsc_call(MatMPIAIJSetPreallocation(fino.K, width, PETSC_NULL, width, PETSC_NULL));
  petsc_call(MatSeqAIJSetPreallocation(fino.K, width, PETSC_NULL));
  petsc_call(MatSetOption(fino.K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));

  // TODO: add an option
//  petsc_call(MatSetOption(fino.K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  
  if (fino.degrees > 1) {
    petsc_call(MatSetBlockSize(fino.K, fino.degrees));
  }

  // the solution (unknown) vector
  petsc_call(MatCreateVecs(fino.K, &fino.phi, NULL));
  petsc_call(PetscObjectSetName((PetscObject)fino.phi, "phi"));
  petsc_call(VecSetFromOptions(fino.phi));
  // explicit initial value
  petsc_call(VecSet(fino.phi, 0));

  
  if (fino.math_type != math_type_eigen) {
    // the right-hand-side vector
    petsc_call(MatCreateVecs(fino.K, NULL, &fino.b));
    petsc_call(PetscObjectSetName((PetscObject)fino.b, "b"));
    petsc_call(VecSetFromOptions(fino.b));
  }
  
  if (fino.problem_family == problem_family_modal ||
      (fino.problem_family == problem_family_thermal && wasora_var_value(wasora_special_var(end_time)) != 0)) {
    // the mass matrix for modal or heat transient
    petsc_call(MatCreate(PETSC_COMM_WORLD, &fino.M));
    petsc_call(PetscObjectSetName((PetscObject)fino.M, "M"));
    petsc_call(MatSetSizes(fino.M, fino.size_local, fino.size_local, fino.global_size, fino.global_size));
    petsc_call(MatSetFromOptions(fino.M));
    petsc_call(MatMPIAIJSetPreallocation(fino.M, width, PETSC_NULL, width, PETSC_NULL));
    petsc_call(MatSeqAIJSetPreallocation(fino.M, width, PETSC_NULL));
    petsc_call(MatSetOption(fino.M, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
    if (fino.degrees > 1) {
      petsc_call(MatSetBlockSize(fino.M, fino.degrees));
    }
  }
  
  
  // ask for the local ownership range
  petsc_call(MatGetOwnershipRange(fino.K, &fino.first_row, &fino.last_row));
  fino.first_node = fino.first_row / fino.degrees;
  fino.last_node = fino.last_row / fino.degrees;
  
  // TODO: partition mesh
  // https://lists.mcs.anl.gov/pipermail/petsc-users/2014-April/021433.html
  fino.first_element = (fino.mesh->n_elements / wasora.nprocs) * wasora.rank;
  if (fino.mesh->n_elements % wasora.nprocs > wasora.rank) {
    fino.first_element += wasora.rank;
    fino.last_element = fino.first_element + (fino.mesh->n_elements / wasora.nprocs) + 1;
  } else {  
    fino.first_element += fino.mesh->n_elements % wasora.nprocs;
    fino.last_element = fino.first_element + (fino.mesh->n_elements / wasora.nprocs);
  }  

  if (fino.mesh->structured) {
    wasora_mesh_struct_init_rectangular_for_nodes(fino.mesh);
  }

  // fill in the holders of the continuous functions that will hold the solution
  
  if (fino.rough == 0) {
    for (g = 0; g < fino.degrees; g++) {
      fino.solution[g]->mesh = fino.mesh;
      fino.solution[g]->data_size = fino.spatial_unknowns;
      fino.solution[g]->data_argument = fino.mesh->nodes_argument;
      fino.solution[g]->data_value = calloc(fino.spatial_unknowns, sizeof(double));
    
      if (fino.nev > 0) {
        for (i = 0; i < fino.nev; i++) {
          fino.mode[g][i]->mesh = fino.mesh;
          fino.mode[g][i]->data_argument = fino.solution[0]->data_argument;
          fino.mode[g][i]->data_size = fino.mesh->n_nodes;
          fino.mode[g][i]->data_value = calloc(fino.spatial_unknowns, sizeof(double));
        }
      }
    }
    
  } else {

    mesh_post_t *post; 
    
    fino_init_rough_mesh();
    // maybe the macros fill_* can be used
    for (g = 0; g < fino.degrees; g++) {
      fino.solution[g]->mesh = fino.mesh_rough;
      fino.solution[g]->data_size = fino.mesh_rough->n_nodes;
      fino.solution[g]->data_argument = fino.mesh_rough->nodes_argument;
      fino.solution[g]->data_value = calloc(fino.mesh_rough->n_nodes, sizeof(double));
    }
    
    // si estamos en rough tenemos que cambiar la malla de salida de los MESH_POSTs
    LL_FOREACH(wasora_mesh.posts, post) {
      if (post->mesh == fino.mesh) {
        post->mesh = fino.mesh_rough;
      }
    }
    
  }

  wasora_call(mesh_node_indexes(fino.mesh, fino.degrees));
  
  return WASORA_PARSER_OK;
}

  
int fino_init_rough_mesh(void) {
  
  int i, i_global;
  int j, j_global;
  int m;

  element_t *element;
  node_t *node;
  element_list_item_t *element_list;
  
  
  // primera pasada, copiamos lo que podemos, alocamos y contamos la cantidad de nodos
  // ya la alocamos en parser
  fino.mesh_rough->bulk_dimensions = fino.mesh->bulk_dimensions;
  fino.mesh_rough->spatial_dimensions = fino.mesh->spatial_dimensions;
  fino.mesh_rough->n_elements = fino.mesh->n_elements;
  fino.mesh_rough->element = calloc(fino.mesh_rough->n_elements, sizeof(element_t));
  fino.mesh_rough->n_nodes = 0;
  i_global = 0;
  for (i = 0; i < fino.mesh_rough->n_elements; i++) {
    if (fino.mesh->element[i].type->dim == fino.mesh_rough->bulk_dimensions) {
      element = &fino.mesh_rough->element[i_global];

      element->index = i_global;
      element->tag = i+1;
      element->type = fino.mesh->element[i].type;
      element->physical_entity = fino.mesh->element[i].physical_entity;

      fino.mesh_rough->n_nodes += element->type->nodes;
      i_global++;
    }  
  }
  
  fino.mesh_rough->n_elements = i_global;
  fino.mesh_rough->element = realloc(fino.mesh_rough->element, fino.mesh_rough->n_elements*sizeof(element_t));
  
  // segunda pasada, creamos los nodos
  j_global = 0;
  fino.mesh_rough->node = calloc(fino.mesh_rough->n_nodes, sizeof(node_t));
  for (i = 0; i < fino.mesh_rough->n_elements; i++) {
    element = &fino.mesh_rough->element[i];
    element->node = calloc(element->type->nodes, sizeof(node_t));
    
    for (j = 0; j < element->type->nodes; j++) {
      node = &fino.mesh_rough->node[j_global];
      node->tag = j_global+1;
      node->index_mesh = j_global;
      node->x[0] = fino.mesh->element[element->tag-1].node[j]->x[0];
      node->x[1] = fino.mesh->element[element->tag-1].node[j]->x[1];
      node->x[2] = fino.mesh->element[element->tag-1].node[j]->x[2];
      
      node->phi = fino.mesh->element[element->tag-1].node[j]->phi;
//      node->dphidx = fino.mesh->element[element->tag-1].dphidx_node[j];
      
      element_list = calloc(1, sizeof(element_list_item_t));
      element_list->element = element;
      LL_APPEND(node->associated_elements, element_list);
      
      element->node[j] = node;
      
      j_global++;
    }
  }

  // rellenamos un array de nodos que pueda ser usado como argumento de funciones
  fino.mesh_rough->nodes_argument = calloc(fino.mesh_rough->spatial_dimensions, sizeof(double *));
  for (m = 0; m < fino.mesh_rough->spatial_dimensions; m++) {
    fino.mesh_rough->nodes_argument[m] = calloc(fino.mesh_rough->n_nodes, sizeof(double));
    for (j = 0; j < fino.mesh_rough->n_nodes; j++) {
      fino.mesh_rough->nodes_argument[m][j] = fino.mesh_rough->node[j].x[m]; 
    }
  }
  
  // ponemos a disposicion la cantidad de nodos rough
  wasora_var_value(fino.vars.nodes_rough) = fino.mesh_rough->n_nodes;
  
  return WASORA_RUNTIME_OK;

}
  
  
int fino_problem_free(void) {
  int g, d;

  if (fino.mesh != NULL && fino.mesh->n_elements != 0) {
    for (g = 0; g < fino.degrees; g++) {
      for (d = 0; d < fino.dimensions; d++) {
        if (fino.gradient != NULL && fino.gradient[g] != NULL) {  
          free(fino.gradient[g][d]->data_value);
          fino.gradient[g][d]->data_value = NULL;
        }  
        if (fino.delta_gradient != NULL && fino.delta_gradient[g] != NULL) {  
          free(fino.delta_gradient[g][d]->data_value);
          fino.delta_gradient[g][d]->data_value = NULL;
        }  
      }
      
      free(fino.solution[g]->data_value);
      fino.solution[g]->data_value = NULL;      
    }
    
    if (fino.sigma != NULL) {
      free(fino.sigmax->data_value);
      if (fino.dimensions > 1) {
        free(fino.sigmay->data_value);
        free(fino.tauxy->data_value);
        if (fino.dimensions > 2) {
          free(fino.sigmaz->data_value);
          free(fino.tauyz->data_value);
          free(fino.tauzx->data_value);
        }  
      }
      free(fino.sigma1->data_value);
      free(fino.sigma2->data_value);
      free(fino.sigma3->data_value);
      free(fino.sigma->data_value);
      free(fino.tresca->data_value);
      
      fino.sigmax->data_value = NULL;
      if (fino.dimensions > 1) {
        fino.sigmay->data_value = NULL;
        fino.tauxy->data_value = NULL;
        if (fino.dimensions > 2) {
          fino.sigmaz->data_value = NULL;
          fino.tauyz->data_value = NULL;
          fino.tauzx->data_value = NULL;
        }
      }
      fino.sigma1->data_value = NULL;
      fino.sigma2->data_value = NULL;
      fino.sigma3->data_value = NULL;
      fino.sigma->data_value = NULL;
      fino.tresca->data_value = NULL;
    }
    
    mesh_free(fino.mesh);
    
  }

  if (fino.unknown_name != NULL) {
    if (fino.degrees != 0) {
      for (g = 0; g < fino.degrees; g++) {
        free(fino.unknown_name[g]);
        fino.unknown_name[g] = NULL;
      }
    }  
    free(fino.unknown_name);
    fino.unknown_name = NULL;
  }
  
  
  wasora_call(fino_free_elemental_objects());
  
  if (fino.problem_family == problem_family_mechanical) {
    
    fino_function_clean_nodal_data(fino.sigma1);
    fino_function_clean_nodal_data(fino.sigma2);
    fino_function_clean_nodal_data(fino.sigma3);
    fino_function_clean_nodal_data(fino.sigma);
    fino_function_clean_nodal_data(fino.tresca);
    
  }
  
  fino.n_dirichlet_rows = 0;
  free(fino.dirichlet_indexes);
  free(fino.dirichlet_values);
     
  if (fino.phi != PETSC_NULL) {
    petsc_call(VecDestroy(&fino.phi));
  }
  if (fino.K != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.K));
  }
  if (fino.K_nobc != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.K_nobc));
  }
  
  if (fino.J != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.J));
  }
  if (fino.M != PETSC_NULL) {
    petsc_call(MatDestroy(&fino.M));
  }
  if (fino.b != PETSC_NULL) {
    petsc_call(VecDestroy(&fino.b));
  }

  // mind the order!
  if (fino.ts != PETSC_NULL) {
    petsc_call(TSDestroy(&fino.ts));
  }
  if (fino.snes != PETSC_NULL) {
    petsc_call(SNESDestroy(&fino.snes));
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

int fino_function_clean_nodal_data(function_t *function) {
 
  if (function != NULL && function->data_value != NULL) {  
    free(function->data_value);
    function->data_value = NULL;
  }
  
  return 0;
}

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

int fino_define_result_function(char *name, function_t **function) {

  // aca la definimos para que este disponible en tiempo de parseo  
  if ((*function = wasora_define_function(name, fino.dimensions)) == NULL) {
    wasora_push_error_message("result function '%s' defined twice", name);
    return WASORA_RUNTIME_ERROR;
  }
  (*function)->mesh = fino.mesh; // esto puede cambiar a rough despues  
  fino_function_clean_nodal_arguments(*function);
  (*function)->var_argument = fino.solution[0]->var_argument;
  (*function)->type = type_pointwise_mesh_node;

  return 0;
}
