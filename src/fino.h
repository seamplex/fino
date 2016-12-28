/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino main header
 *
 *  Copyright (C) 2015--2016 jeremy theler
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

#define BC_UNDEFINED       0
#define BC_DIRICHLET       1
#define BC_NEUMANN         2
#define BC_ROBIN           3
#define BC_POINT_FORCE     4
#define BC_DIRICHLET_NULL  5
#define BC_DIRICHLET_ALG   6

#define BC_FACTOR 1.0

// forward definitions
typedef struct fino_reaction_t fino_reaction_t;
typedef struct fino_distribution_t fino_distribution_t;
typedef struct fino_step_t fino_step_t;
typedef struct fino_times_t fino_times_t;
typedef struct debug_t debug_t;


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
    math_automatic,
    math_linear,
    math_nonlinear,
    math_eigen,
  } math_type;
  
  enum {
    problem_undefined,
    problem_generic,
    problem_bake,
    problem_break,
    problem_shake,
    problem_break_plane_stress,
    problem_shake_plane_stress,
    problem_break_plane_strain,
    problem_shake_plane_strain,
  } problem;

  int spatial_unknowns;  // cant de incognitas espaciales (= celdas o nodos)
  int degrees;
  int dimensions;
  
  int problem_size;
  
  mesh_t *mesh;
  debug_t *debugs;
    
  // variables internas
  struct {
    var_t *atol;
    var_t *rtol;
    var_t *divtol;
    var_t *max_iterations;
    var_t *iterations;
    var_t *residual;
    
    var_t *dirichlet_diagonal;
  
    var_t *nodes;
    var_t *elements;
    var_t *unknowns;
    
    var_t *residual_norm;
    var_t *error_estimate;
    var_t *rel_error;
    
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
  
  // flag
  PetscInt petscinit_called;
  
  // cosas para paralelizacion
  PetscInt rank;
  PetscInt size;

  // objetos globales
  Vec phi;       // el vector incognita
  Vec b;         // el vector del miembro derecho
  Mat A;         // la matriz
  Mat B;         // la matriz del miembro derecho (autovalores)
  PetscScalar lambda; // el autovalor

  // contexto del solver de krylov
  KSP ksp;
  PC pc;

  loadable_routine_t *user_provided_linearsolver;
  
  // strings con tipos
  char *ksp_type;
  char *pc_type;

  char *eps_type;
  char *st_type;
  
  int set_near_null_space;
  int use_pcsetcoordinates;
  int set_block_size;
  
  expr_t eps_ncv;
  expr_t st_shift;
  expr_t st_anti_shift;  
  
  // objetos intermedios que evalua fino y se lo deja disponible a wasora
  vector_t *h;       // funciones de forma

  // jacobiano de las funciones de forma con respecto a las coordenadas reales
  matrix_t *dhdx;

  // objectos locales
  int n_local_nodes;            // cantidad de nodos locales actual
  int elemental_size;           // tamanio actual del elemento
  gsl_matrix *Ai;               // la matriz elemental
  gsl_matrix *Bi;               // la matriz del miembro derecho elemental (para autovalores)
  gsl_vector *bi;               // el vector del miembro derecho elemental

  // holder para calcular las reacciones de vinculo de BCs dirichlet
  int n_dirichlet_rows;
  dirichlet_row_t *dirichlet_row;
  
  fino_reaction_t *reactions;

  
  // user-provided functions para los objetos elementales, las linkeamos
  // a las que dio el usuario en el input en init
  function_t ***Ai_function;
  function_t ***Bi_function;
  function_t **bi_function;
  
  
  // nombres custom 
  char **unknown_name;          // uno para cada grado de libertad
  char *shape_name;             // despues le agregamos el numero (i.e. h1)
  char *lhs_matrix_name;        // despues le agregamos los indices (i.e A1.1)
  char *rhs_matrix_name;        // despues le agregamos los indices (i.e B1.1)
  char *rhs_vector_name;        // despues le agregamos los indices (i.e b1)
  
  // las funciones con la solucion (una para cada grado de libertad)
  function_t **solution;
  // las derivadas de la solucion con respecto al espacio;
  function_t ***gradient;
  
  // soluciones anteriores (por ejemplos desplazamientos)
  function_t **base_solution;
  function_t ***base_gradient;
  
  function_t *sigma;     // von misses
  function_t *sigma1;    // principal
  function_t *sigma2;
  function_t *sigma3;

  // cosas del eigensolver
  // las pongo al final por si acaso (mezcla de plugins compilados con difentes libs, no se)
#ifdef HAVE_SLEPC
  int nev;      // el numero del autovalor pedido
  EPS eps;      // contexto eigensolver (SLEPc)
  ST st;        // contexto de la transformacion espectral asociada
  EPSWhich eigen_spectrum;
#endif
  
} fino;

struct fino_reaction_t {
  physical_entity_t *physical_entity;
  char *name_root;
  var_t *R[3];
  
  fino_reaction_t *next;
};

// se rellena value a partir o bien de una variable o bien de una funcion
struct fino_distribution_t {
  int defined;
  var_t *variable;
  function_t *function;
};

struct fino_step_t {
  int do_not_build;
  int do_not_solve;
  int do_not_compute_gradients;
};

struct debug_t {
  file_t *file;
  PetscViewer viewer;
  
  int matrices;   // bitmap con flags de que hay que exportar
  expr_t matrices_size;
  expr_t matrices_stride;

  int include_input;

  int file_opened;
  
  debug_t *next;
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
extern PetscErrorCode fino_monitor_dots(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);

// bc.c
extern int fino_read_bcs(void);
extern int fino_count_bc_expressions(expr_t *bc_args);
extern int fino_evaluate_bc_expressions(physical_entity_t *, node_t *, int, double, double *);
extern int fino_set_essential_bc(void);
extern int fino_build_surface_objects(element_t *, expr_t *, expr_t *);
extern int fino_add_single_surface_term_to_rhs(element_t *, bc_string_based_t *);

// bulk.c
extern int fino_allocate_elemental_objects(element_t *);

extern int fino_print_gsl_vector(gsl_vector *, FILE *);
extern int fino_print_gsl_matrix(gsl_matrix *, FILE *);

// debug.c
extern int fino_debug_open(debug_t *);
extern int fino_debug_initial(debug_t *);
extern int fino_instruction_debug(void *);
extern int fino_debug_close(debug_t *);
extern int fino_print_petsc_vector(Vec, PetscViewer);
extern int fino_print_petsc_matrix(Mat, PetscViewer);
extern int fino_print_petsc_matrix_struct(Mat, PetscViewer);

// init.c
extern int plugin_init_before_parser(void);
extern int plugin_init_after_parser(void);
extern int plugin_init_before_run(void);
extern int plugin_finalize(void);
extern int fino_problem_free(void);

// parser.c
extern int fino_parse_line(char *);
extern int fino_define_functions(void);

// version.c
extern void fino_usage(char *);
extern void fino_version(FILE *, int, int);
extern void fino_license(FILE *);

// eigen_slepc.c
extern int fino_solve_linear_petsc(void);

// linear_petsc.c
extern int fino_solve_eigen_slepc(void);

// bulk.c
extern int fino_build_bulk(void);
extern int fino_build_element(element_t *);

// bc.c
extern int fino_set_essential_bc(void);
extern int fino_problem_init(void);

// gradient.c
extern int fino_compute_gradients(void);

// petschandler.c
PetscErrorCode fino_handler(MPI_Comm comm, int, const char *, const char *, PetscErrorCode, PetscErrorType, const char *, void *);

// times.c
extern double fino_get_cpu_time(void);

// breakshake.c
extern int fino_build_breakshake(element_t *, int);
extern int fino_break_compute_C(gsl_matrix *);
extern int fino_break_compute_stresses(void);
extern int fino_break_compute_reactions(void);


// bake.c
extern int fino_build_bake(element_t *, int);

extern const char *plugin_name(void);
extern const char *plugin_version(void);
extern const char *plugin_description(void);
extern const char *plugin_longversion(void);
extern const char *plugin_copyright(void);

// distributions.c
extern int fino_distribution_init(fino_distribution_t *, const char *);
extern double fino_distribution_evaluate(fino_distribution_t *);

#endif  /* _FINO_H_ */
