/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino main header
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

#ifndef _FINO_H_
#define _FINO_H_
#include <sys/time.h>
#include <sys/resource.h>

#include <wasora.h>

#include <petscksp.h>
#include <petscpc.h>
#include <petsctime.h>
#ifdef HAVE_SLEPC
#include <slepceps.h>
#endif


#if defined(PETSC_USE_COMPLEX)
#error "PETSc should be compiled with real scalars to be used with fino."
#endif

PetscErrorCode petsc_err;
#define petsc_call(s) {petsc_err = s; CHKERRQ(petsc_err);}

// maxima cantidad de dimensiones espaciales
#define M_MAX   3
// maxima cantidad de nodos (para ponerle un tamanio al vector)
#define J_MAX   8

#define DEBUG_MATRICES_ASCII                       1
#define DEBUG_MATRICES_ASCII_STRUCT                2
#define DEBUG_MATRICES_PETSC_BINARY                4
#define DEBUG_MATRICES_PETSC_COMPRESSED_BINARY     8
#define DEBUG_MATRICES_PETSC_ASCII                16
#define DEBUG_MATRICES_PETSC_OCTAVE               32
#define DEBUG_MATRICES_PETSC_DENSE                64
#define DEBUG_MATRICES_X                         128
#define DEBUG_MATRICES_SNG                       256
#define DEBUG_MATRICES_SNG_STRUCT                512

#define DEFAULT_MATRICES_X_SIZE     640

#define bc_math_undefined                              0
#define bc_math_dirichlet                              1
#define bc_math_neumann                                2
#define bc_math_robin                                  3

#define bc_phys_undefined                              5
#define bc_phys_displacement_fixed                     6
#define bc_phys_displacement                           7
#define bc_phys_displacement_constrained               8
#define bc_phys_displacement_mimic                     9
#define bc_phys_pressure                              16
#define bc_phys_stress                                17
#define bc_phys_force                                 18
#define bc_phys_moment                                19
#define bc_phys_temperature                           32
#define bc_phys_heat_flux                             33
#define bc_phys_convection                            34

#define bc_dof_moment_offset              64
#define bc_dof_coordinates_offset        128


#define BC_FACTOR 0.1

// forward definitions
typedef struct fino_distribution_t fino_distribution_t;
typedef struct fino_step_t fino_step_t;
typedef struct fino_times_t fino_times_t;
typedef struct fino_linearize_t fino_linearize_t;
typedef struct fino_debug_t fino_debug_t;



typedef struct {
  physical_entity_t *physical_entity;
  
  int n_cols;
  PetscScalar *alg_val;
  PetscInt *alg_col;
  
  int dof;
  PetscInt ncols;
  PetscInt *cols;
  PetscScalar *vals;
} dirichlet_row_t;


// estructura admnistrativa
struct {

  enum {
    math_type_linear,
    math_type_nonlinear,
    math_type_eigen,
  } math_type;
  
  enum {
    problem_family_undefined,
    problem_family_bake,
    problem_family_break,
    problem_family_shake,
  } problem_family;
  
  enum {
    problem_kind_undefined,
    problem_kind_full3d,
    problem_kind_plane_stress,
    problem_kind_plane_strain,
    problem_kind_axisymmetric    
  } problem_kind;
  
  enum {
    symmetry_axis_none,
    symmetry_axis_x,
    symmetry_axis_y
  } symmetry_axis;

  int spatial_unknowns;  // cant de incognitas espaciales (= celdas o nodos)
  int degrees;
  int dimensions;
  
  int problem_size;
  
  mesh_t *mesh;
  fino_linearize_t *linearizes;
  fino_debug_t *debugs;
  
  // variables internas
  struct {
    var_t *abstol;
    var_t *reltol;
    var_t *divtol;
    var_t *max_iterations;
    var_t *gamg_threshold;
    var_t *iterations;
    var_t *residual_norm;
    
    var_t *dirichlet_diagonal;
  
    var_t *unknowns;
    
//    var_t *error_estimate;
//    var_t *rel_error;
    
    var_t *nx;
    var_t *ny;
    var_t *nz;
    
    var_t *U[3];

    var_t *displ_max;
    var_t *displ_max_x;
    var_t *displ_max_y;
    var_t *displ_max_z;

    var_t *u_at_displ_max;
    var_t *v_at_displ_max;
    var_t *w_at_displ_max;

    
    var_t *sigma_max;
    var_t *sigma_max_x;
    var_t *sigma_max_y;
    var_t *sigma_max_z;

    var_t *u_at_sigma_max;
    var_t *v_at_sigma_max;
    var_t *w_at_sigma_max;

    var_t *T_max;
    var_t *T_min;
    
    var_t *lambda;

    var_t *time_wall_build;
    var_t *time_wall_solve;
    var_t *time_wall_total;

    var_t *time_cpu_build;
    var_t *time_cpu_solve;
    var_t *time_cpu_total;

    var_t *time_petsc_build;
    var_t *time_petsc_solve;
    var_t *time_petsc_total;

    var_t *flops_petsc;
    
    var_t *available_memory;
    var_t *memory_usage_global;
    var_t *memory_usage_petsc;
    
  } vars;

  // vectores
  struct {
    vector_t *omega;
  } vectors;

  // flag
  PetscInt petscinit_called;
  
  // cosas para paralelizacion
  PetscInt rank;
  PetscInt size;

  // objetos globales
  Vec phi;       // el vector incognita
  Vec b;         // el vector del miembro derecho para el steady-state
  Mat K;         // la matriz de rigidez (con E para elastico y k para calor)
  Mat M;         // la matriz de masa (con rho para elastico y rho*cp para calor)
  PetscScalar lambda; // el autovalor
  
  PetscScalar *eigenvalue;    // los autovalores
  Vec *eigenvector;           // los autovectores

  Mat A;         // las matrices para el transient de calor
  Mat B;
  Vec c;         // el vector del miembro derecho para el transient

  Mat lastM;
  Mat dotM;
  Vec m;
  
  // contexto del solver de krylov
  KSP ksp;
  PC pc;
  
  int has_mass;
  int has_rhs;
  int has_transient;

  loadable_routine_t *user_provided_linearsolver;
  
  // strings con tipos
  char *ksp_type;
  char *pc_type;

  char *eps_type;
  char *st_type;
  
  enum {
    set_near_nullspace_rigidbody,
    set_near_nullspace_fino,
    set_near_nullspace_none
  } set_near_nullspace;
  
  int do_not_set_block_size;
  
  char *shmem_progress_build_name;  
  char *shmem_progress_solve_name;
  char *shmem_memory_name;
  
  double *shmem_progress_build;
  double *shmem_progress_solve;
  double *shmem_memory;
  
  
  
  expr_t eps_ncv;
  expr_t st_shift;
  expr_t st_anti_shift;  

  // para la memoria  
  struct rusage resource_usage;
  
  // objetos intermedios que evalua fino y se lo deja disponible a wasora
  vector_t *h;       // funciones de forma

  // jacobiano de las funciones de forma con respecto a las coordenadas reales
  matrix_t *dhdx;

  // objectos locales
  int n_local_nodes;            // cantidad de nodos locales actual
  int elemental_size;           // tamanio actual del elemento
  gsl_matrix *Ki;               // la matriz de rigidez elemental
  gsl_matrix *Mi;               // la matriz de masa elemental
  gsl_vector *bi;               // el vector del miembro derecho elemental

  // holder para calcular las reacciones de vinculo de BCs dirichlet
  int n_dirichlet_rows;
  dirichlet_row_t *dirichlet_row;
  int n_algebraic_rows;
  dirichlet_row_t *algebraic_row;
  
  // user-provided functions para los objetos elementales, las linkeamos
  // a las que dio el usuario en el input en init
  function_t ***Ai_function;
  function_t ***Bi_function;
  function_t **bi_function;
  
  
  // nombres custom 
  char **unknown_name;          // uno para cada grado de libertad
  
  // las funciones con la solucion (una para cada grado de libertad)
  function_t **solution;
  // las derivadas de la solucion con respecto al espacio;
  function_t ***gradient;
  // los modos de vibracion
  function_t ***vibration;
  
  // soluciones anteriores (por ejemplos desplazamientos)
  function_t **base_solution;
  function_t ***base_gradient;
  
  enum {
    gradient_undefined,
    gradient_mass_matrix_consistent,
    gradient_mass_matrix_row_sum,
    gradient_mass_matrix_lobatto,
    gradient_mass_matrix_diagonal,
    gradient_node_average_corner,
    gradient_node_average_all,
    gradient_gauss_average,
    gradient_none,
  } gradient_evaluation;
  
  double gradient_jacobian_threshold;
  
  // tensor de tensiones
  function_t *sigmax;
  function_t *sigmay;
  function_t *sigmaz;
  function_t *tauxy;
  function_t *tauyz;
  function_t *tauzx;

  function_t *sigma1;    // principales
  function_t *sigma2;
  function_t *sigma3;
  function_t *sigma;     // von misses
  function_t *tresca;

  // cosas del eigensolver
  // las pongo al final por si acaso (mezcla de plugins compilados con difentes libs, no se)
  int nev;      // el numero del autovalor pedido
#ifdef HAVE_SLEPC
  EPS eps;      // contexto eigensolver (SLEPc)
  ST st;        // contexto de la transformacion espectral asociada
  EPSWhich eigen_spectrum;
#endif
  
} fino;

/*
struct fino_reaction_t {
  physical_entity_t *physical_entity;
  char *name_root;
  var_t *R[3];
  
  fino_reaction_t *next;
};
*/

// se rellena value a partir o bien de una variable o bien de una funcion
struct fino_distribution_t {
  int defined;
  physical_property_t *physical_property;
  var_t *variable;
  function_t *function;
};

struct fino_step_t {
  int do_not_build;
  int do_not_solve;
  int do_not_compute_gradients;
};

struct fino_linearize_t {
  physical_entity_t *scl;

  var_t *membrane;
  var_t *bending;
  var_t *peak;
//  var_t *membrane_plus_bending;  
  
  fino_linearize_t *next;
};

struct fino_debug_t {
  file_t *file;
  PetscViewer viewer;
  
  int matrices;   // bitmap con flags de que hay que exportar
  expr_t matrices_size;
  expr_t matrices_stride;

  int include_input;

  int file_opened;
  
  fino_debug_t *next;
};

// para medir tiempos (wall y cpu)
struct fino_times_t {
  PetscLogDouble init_begin;
  PetscLogDouble init_end;
  PetscLogDouble build_begin;
  PetscLogDouble build_end;
  PetscLogDouble solve_begin;
  PetscLogDouble solve_end;
};

// fino.c
extern int fino_instruction_step(void *);
extern int fino_assembly(void);
extern PetscErrorCode fino_ksp_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);

// bc.c
extern int fino_read_bcs(void);
extern int fino_count_bc_expressions(expr_t *bc_args);
extern int fino_evaluate_bc_expressions(physical_entity_t *, node_t *, int, double, double *);
extern int fino_set_essential_bc(Mat, Vec);
extern int fino_build_surface_objects(element_t *, expr_t *, expr_t *);
extern int fino_add_single_surface_term_to_rhs(element_t *, bc_string_based_t *);

// bulk.c
extern int fino_allocate_elemental_objects(element_t *);
extern int fino_build_bulk(void);
extern int fino_build_element_volumetric(element_t *);
extern int fino_build_element_bc(element_t *);

extern double fino_compute_r_for_axisymmetric(void);

extern int fino_print_gsl_vector(gsl_vector *, FILE *);
extern int fino_print_gsl_matrix(gsl_matrix *, FILE *);

// debug.c
extern int fino_debug_open(fino_debug_t *);
extern int fino_debug_initial(fino_debug_t *);
extern int fino_instruction_debug(void *);
extern int fino_debug_close(fino_debug_t *);
extern int fino_print_petsc_vector(Vec, PetscViewer);
extern int fino_print_petsc_matrix(Mat, PetscViewer);
extern int fino_print_petsc_matrix_struct(Mat, PetscViewer);

// init.c
extern int plugin_init_before_parser(void);
extern int plugin_init_after_parser(void);
extern int plugin_init_before_run(void);
extern int plugin_finalize(void);
extern int fino_problem_init(void);
extern int fino_problem_free(void);
extern int fino_function_clean_nodal_data(function_t *);
extern int fino_function_clean_nodal_arguments(function_t *);
extern int fino_define_result_function(char *, function_t **);

// linearize.c
extern int fino_instruction_linearize(void *);

// parser.c
extern int fino_parse_line(char *);
extern int fino_define_functions(void);

// version.c
extern void fino_usage(char *);
extern void fino_version(FILE *, int, int);
extern void fino_license(FILE *);

// linear_petsc.c
extern int fino_solve_linear_petsc(Mat, Vec);

// eigen_slepc.c
extern int fino_solve_eigen_slepc(Mat, Mat);

// bulk.c
extern int fino_build_bulk(void);
extern int fino_build_element_volumetric(element_t *);
extern int fino_build_element_bc(element_t *);

// gradient.c
extern int fino_compute_gradients(void);

// petschandler.c
PetscErrorCode fino_handler(MPI_Comm comm, int, const char *, const char *, PetscErrorCode, PetscErrorType, const char *, void *);

// times.c
extern double fino_get_cpu_time(void);

// breakshake.c
extern int fino_break_build_element(element_t *, int);
extern int fino_break_compute_C(gsl_matrix *, double, double);
extern int fino_break_compute_stresses(void);
extern int fino_break_compute_reactions(void);
extern int fino_break_set_stress(element_t *);
extern int fino_break_set_force(element_t *);
extern int fino_break_set_pressure(element_t *);
extern int fino_break_set_moment(element_t *);
extern int fino_compute_principal_stress(double, double, double, double, double, double, double *, double *, double *);
extern double fino_compute_vonmises_from_principal(double, double, double);
extern double fino_compute_vonmises_from_tensor(double, double, double, double, double, double);
extern double fino_compute_tresca_from_principal(double, double, double);
extern double fino_compute_tresca_from_tensor(double, double, double, double, double, double);


// bake.c
extern int fino_bake_step_initial();
extern int fino_bake_step_transient();
extern int fino_build_bake(element_t *, int);
extern int fino_bake_set_heat_flux(element_t *);
extern int fino_bake_set_convection(element_t *element);
extern int fino_bake_compute_fluxes(void);

extern const char *plugin_name(void);
extern const char *plugin_version(void);
extern const char *plugin_description(void);
extern const char *plugin_longversion(void);
extern const char *plugin_copyright(void);

// distributions.c
extern int fino_distribution_init(fino_distribution_t *, const char *);
extern double fino_distribution_evaluate(fino_distribution_t *, material_t *, double *);

#endif  /* _FINO_H_ */