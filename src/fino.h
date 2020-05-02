/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino main header
 *
 *  Copyright (C) 2015--2020 jeremy theler
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

#ifndef _FINO_H
#define _FINO_H
#include <sys/time.h>
#include <sys/resource.h>

#include <wasora.h>

#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>
#include <petsctime.h>
#ifdef HAVE_SLEPC
#include <slepceps.h>
#endif


#if defined(PETSC_USE_COMPLEX)
 #error "PETSc should be compiled with real scalars to be used with fino."
#endif

PetscErrorCode petsc_err;
#define petsc_call(s) {petsc_err = s; CHKERRQ(petsc_err);}

#define CHAR_PROGRESS_BUILD    "."
#define CHAR_PROGRESS_SOLVE    "-"
#define CHAR_PROGRESS_GRADIENT "="

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

#define DEFAULT_MATRICES_X_SIZE                  640
#define DEFAULT_GRADIENT_JACOBIAN_THRESHOLD       -1

#define DEFAULT_NMODES        10

#define bc_math_undefined                              0
#define bc_math_dirichlet                              1
#define bc_math_neumann                                2
#define bc_math_robin                                  3

#define bc_phys_undefined                              5

#define bc_phys_displacement_fixed                     6
#define bc_phys_displacement                           7
#define bc_phys_displacement_constrained               8
#define bc_phys_displacement_symmetry                  9
#define bc_phys_displacement_radial                   10
#define bc_phys_displacement_mimic                    11

#define bc_phys_pressure_normal                       16   // normal stress positive in the outward normal direction, i.e. traction
#define bc_phys_stress                                17
#define bc_phys_force                                 18
#define bc_phys_moment                                19
#define bc_phys_pressure_real                         20   // normal stress positive in the inward normal direction, i.e. compression

#define bc_phys_temperature                           32
#define bc_phys_heat_flux                             33
#define bc_phys_heat_total                            34
#define bc_phys_convection                            35

#define bc_dof_moment_offset              64
#define bc_dof_coordinates_offset        128


#define BC_FACTOR 1.00  // TODO: si esto es != 1 da cosas raras con los reallocs, pensar mejor


// forward definitions
typedef struct fino_distribution_t fino_distribution_t;
typedef struct fino_step_t fino_step_t;
typedef struct fino_times_t fino_times_t;
typedef struct fino_reaction_t fino_reaction_t;
typedef struct fino_linearize_t fino_linearize_t;
typedef struct fino_debug_t fino_debug_t;
typedef struct fino_roughish_avg_t fino_roughish_avg_t;


typedef struct {
  physical_entity_t *physical_entity;
  
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
    math_type_undefined,
    math_type_linear,
    math_type_nonlinear,
    math_type_eigen,
  } math_type;
  
  enum {
    problem_family_undefined,
    problem_family_mechanical,
    problem_family_thermal,
    problem_family_modal,
  } problem_family;
  
  enum {
    problem_kind_undefined,
    problem_kind_full3d,
    problem_kind_axisymmetric,
    problem_kind_plane_stress,
    problem_kind_plane_strain
  } problem_kind;
  
  enum {
    symmetry_axis_none,
    symmetry_axis_x,
    symmetry_axis_y
  } symmetry_axis;

  int spatial_unknowns;  // cant de incognitas espaciales (= celdas o nodos)
  int degrees;
  int dimensions;
  
  int global_size;
  int rough;             // con esto se mantienen las contribuciones de cada elemento a las derivadas
  int roughish;          // con esto se promedian por entidad fisica 
  
  mesh_t *mesh;
  mesh_t *mesh_rough;   // esta es una malla donde cada node pertenece solo a un elemento
  
  fino_reaction_t *reactions;
  fino_linearize_t *linearizes;
  
  fino_debug_t *debugs;

  // esto capaz que deberia ir en otro lado
  PetscClassId petsc_classid;

  PetscLogStage petsc_stage_build;
  PetscLogStage petsc_stage_solve;
  PetscLogStage petsc_stage_stress;
  
  PetscLogEvent petsc_event_build;
  PetscLogEvent petsc_event_solve;
  PetscLogEvent petsc_event_stress;
  
  PetscLogDouble petsc_flops_build;
  PetscLogDouble petsc_flops_solve;
  PetscLogDouble petsc_flops_stress;

  
  // variables internas
  struct {
    var_t *abstol;
    var_t *reltol;
    var_t *divtol;
    var_t *max_iterations;
    var_t *gamg_threshold;
    var_t *iterations;
    var_t *residual_norm;
    
    var_t *penalty_weight;
    var_t *nodes_rough;
    var_t *unknowns;
    
//    var_t *error_estimate;
//    var_t *rel_error;

    var_t *U[3];

    var_t *strain_energy;
    
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
    var_t *delta_sigma_max;

    var_t *u_at_sigma_max;
    var_t *v_at_sigma_max;
    var_t *w_at_sigma_max;

    var_t *T_max;
    var_t *T_min;
    
    var_t *lambda;

    var_t *time_wall_build;
    var_t *time_wall_solve;
    var_t *time_wall_stress;
    var_t *time_wall_total;

    var_t *time_cpu_build;
    var_t *time_cpu_solve;
    var_t *time_cpu_stress;
    var_t *time_cpu_total;

    var_t *time_petsc_build;
    var_t *time_petsc_solve;
    var_t *time_petsc_stress;
    var_t *time_petsc_total;

    var_t *flops_petsc;
    
    var_t *memory_available;
    var_t *memory;
    var_t *memory_petsc;
    
    var_t *M_T;
    
  } vars;

  // vectores
  struct {
    vector_t *f;
    vector_t *omega;
    vector_t *m;
    vector_t *L;
    vector_t *Gamma;
    vector_t *mu;
    vector_t *Mu;
    
    vector_t **phi;
  } vectors;

  // flag
  PetscInt petscinit_called;
  
  // cosas para paralelizacion
//  PetscInt rank;     // estos ahora estan en wasora
//  PetscInt nprocs;
  PetscInt nodes_local, size_local;
  PetscInt first_row, last_row;
  PetscInt first_node, last_node;
  PetscInt first_element, last_element;

  // objetos globales
  Vec phi;       // el vector incognita
  Vec b;         // el vector del miembro derecho para el steady-state
  Mat K;         // la matriz de rigidez (con E para elastico y k para calor)
  Mat K_nobc;    // la matriz de rigidez antes de poner las condiciones de dirichlet (para reacciones)
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
  
  // contexto de solvers de PETSc
  TS ts;
  SNES snes;
  KSP ksp;
  PC pc;
  
  int has_mass;
  int has_rhs;
  int has_transient;

  loadable_routine_t *user_provided_linearsolver;
  
  // strings con tipos
  PetscBool commandline_mumps; // esta sobre-escribe todo
  SNESType snes_type;
  KSPType ksp_type;
  PCType pc_type;

  enum {
    set_near_nullspace_rigidbody,
    set_near_nullspace_fino,
    set_near_nullspace_none
  } set_near_nullspace;
  
  int do_not_set_block_size;
  
  int progress_ascii;
  double progress_r0;
  double progress_last;

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
  gsl_vector *Nb;               // para las BCs de neumann

  // holder para poner las BCs dirichlet (y calcular las reacciones de vinculo)
  int n_dirichlet_rows;
  PetscScalar     *dirichlet_rhs;
  PetscInt        *dirichlet_indexes;
  dirichlet_row_t *dirichlet_row;
  
 
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
  // la incerteza (i.e la desviacion estandar de la contribucion de cada elemento)
  function_t ***delta_gradient;  
  // los modos de vibracion
  function_t ***mode;
  
  // soluciones anteriores (por ejemplos desplazamientos)
  function_t **base_solution;
  
  enum {
    gradient_gauss_extrapolated,
    gradient_at_nodes,
    gradient_none
  } gradient_evaluation;

  enum {
    gradient_weight_volume,
    gradient_weight_quality,
    gradient_weight_volume_times_quality,
    gradient_weight_flat,
  } gradient_element_weight;

  enum {
    gradient_average,
    gradient_actual
  } gradient_highorder_nodes;
  
//  double gradient_quality_threshold;
  
  // tensor de tensiones
  function_t *sigmax;
  function_t *sigmay;
  function_t *sigmaz;
  function_t *tauxy;
  function_t *tauyz;
  function_t *tauzx;

  function_t *sigma1;      // principales
  function_t *sigma2;
  function_t *sigma3;
  function_t *sigma;       // von misses
  function_t *delta_sigma; // incerteza
  function_t *tresca;

  // cosas del eigensolver
  // las pongo al final por si acaso (mezcla de plugins compilados con difentes libs, no se)
  int nev;      // el numero del autovalor pedido
#ifdef HAVE_SLEPC
  EPSType eps_type;
  STType st_type;
  
  EPS eps;      // contexto eigensolver (SLEPc)
  ST st;        // contexto de la transformacion espectral asociada
  EPSWhich eigen_spectrum;
#endif
  
} fino;


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
};

struct fino_reaction_t {
  physical_entity_t *physical_entity;
  var_t *scalar;
  vector_t *vector;
  
  fino_reaction_t *next;
};

struct fino_linearize_t {
  physical_entity_t *physical_entity;
  expr_t x1;
  expr_t y1;
  expr_t z1;
  expr_t x2;
  expr_t y2;
  expr_t z2;

  file_t *file;
  
  enum {
    linearize_vonmises,
    linearize_tresca,
    linearize_principal1,
    linearize_principal2,
    linearize_principal3
  } total;
  
  int ignore_through_thickness;
  
  
  var_t *M;
  var_t *MB;  
  var_t *P;
  
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
  PetscLogDouble stress_begin;
  PetscLogDouble stress_end;
  PetscLogDouble solve_begin;
  PetscLogDouble solve_end;
};

struct fino_roughish_avg_t {
  int smooth_element;
  int local_node;
  fino_roughish_avg_t *next;
};


// fino.c
extern int fino_instruction_step(void *);
extern int fino_assembly(void);
extern int fino_phi_to_solution(Vec phi);

// bc.c
extern int fino_bc_string2parsed(void);
void fino_bc_read_name_expr(bc_t *, char **, char **, char **);
extern int fino_bc_process_mechanical(bc_t *, char *, char *);
extern int fino_bc_process_thermal(bc_t *, char *, char *, char *);
extern int fino_set_essential_bc(Mat, Vec);
extern double fino_gsl_function_of_uvw(double, void *);

// bulk.c
extern int fino_allocate_elemental_objects(element_t *);
extern int fino_free_elemental_objects(void);
extern int fino_build_bulk(void);
extern int fino_build_element_volumetric(element_t *);
extern int fino_build_element_bc(element_t *, bc_t *);

extern double fino_compute_r_for_axisymmetric(element_t *, int);

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
extern int fino_init_rough_mesh(void);
extern int fino_problem_free(void);
extern int fino_function_clean_nodal_data(function_t *);
extern int fino_function_clean_nodal_arguments(function_t *);
extern int fino_define_result_function(char *, function_t **);

// linearize.c
extern int fino_instruction_linearize(void *);
extern double fino_linearization_integrand_membrane(double, void *);
extern double fino_linearization_integrand_bending(double, void *);
  
// parser.c
extern int fino_parse_line(char *);
extern int fino_define_functions(void);

// reaction.c
extern int fino_instruction_reaction(void *);


// version.c
extern void fino_usage(char *);
extern void fino_version(FILE *, int, int);
extern void fino_license(FILE *);

// petsc_ksp.c
extern int fino_solve_petsc_linear(Mat, Vec);
extern PetscErrorCode fino_ksp_monitor(KSP, PetscInt, PetscReal, void *);
extern int fino_set_ksp(void);
extern int fino_set_pc(Mat);

// petsc_snes.c
extern int fino_solve_nonlinear_petsc();
extern PetscErrorCode fino_snes_monitor(SNES, PetscInt, PetscReal, void *);
// slepc_eigen.c
extern int fino_solve_eigen_slepc();
extern int fino_eigen_nev(void);

// gradient.c
extern int fino_compute_gradients(void);

// petschandler.c
PetscErrorCode fino_handler(MPI_Comm comm, int, const char *, const char *, PetscErrorCode, PetscErrorType, const char *, void *);

// times.c
extern double fino_get_cpu_time(void);

// breakshake.c
extern int fino_break_build_element(element_t *, int);
extern int fino_break_set_neumann(element_t *, bc_t *);
extern int fino_break_compute_C(gsl_matrix *, double, double);
extern int fino_break_compute_nodal_stresses(element_t *, int, double, double, double, double *, double *, double *, double *, double *, double *);
extern int fino_break_compute_stresses(void);
extern int fino_break_compute_reactions(void);
extern int fino_break_set_moment(element_t *, bc_t *);
extern int fino_compute_principal_stress(double, double, double, double, double, double, double *, double *, double *);
extern double fino_compute_vonmises_from_principal(double, double, double);
extern double fino_compute_vonmises_from_stress_tensor(double, double, double, double, double, double);
extern double fino_compute_vonmises_from_strains(double, double,
                                          double, double, double,
                                          double, double, double,
                                          double, double, double);
extern double fino_compute_tresca_from_principal(double, double, double);
extern double fino_compute_tresca_from_stress_tensor(double, double, double, double, double, double);
extern int fino_compute_strain_energy(void);

// bake.c
extern int fino_bake_step_initial();
extern int fino_bake_step_transient();
extern int fino_bake_build_element(element_t *, int);
extern int fino_bake_set_heat_flux(element_t *, bc_t *);
extern int fino_bake_set_convection(element_t *, bc_t *);
extern int fino_bake_compute_fluxes(void);

extern const char *plugin_name(void);
extern const char *plugin_version(void);
extern const char *plugin_description(void);
extern const char *plugin_longversion(void);
extern const char *plugin_copyright(void);

// distributions.c
extern int fino_distribution_init(fino_distribution_t *, const char *);
extern double fino_distribution_evaluate(fino_distribution_t *, material_t *, double *);

// post.c
extern int fino_instruction_post(void *);

#endif  /* _FINO_H_ */
