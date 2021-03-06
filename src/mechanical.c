/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's construction of linear elastic problem (break) with optional vibration (shake)
 *  and evaluation of the stress tensor out of the gradients of the displacements
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
#include <ctype.h>
#include <math.h>

#include "fino.h"

#define LAMBDA                          0
#define MU                              1
#define ALPHA                           2
#define ISOTROPIC_ELASTIC_PROPERTIES    3

#define EX                              0
#define EY                              1
#define EZ                              2
#define NUXY                            3
#define NUYZ                            4
#define NUZX                            5
#define GXY                             6
#define GYZ                             7
#define GZX                             8
#define ALPHAX                          9
#define ALPHAY                         10
#define ALPHAZ                         11
#define ORTHOTROPIC_ELASTIC_PROPERTIES 12

#define SIGMAX             0
#define SIGMAY             1
#define SIGMAZ             2
#define TAUXY              3
#define TAUYZ              4
#define TAUZX              5
#define ELASTIC_FUNCTIONS  6


// isotropic properties
fino_distribution_t distribution_E;     // Young's modulus
fino_distribution_t distribution_nu;    // Poisson's ratio

// orthotropic
fino_distribution_t distribution_Ex;    // Young's modulus along axis x
fino_distribution_t distribution_Ey;
fino_distribution_t distribution_Ez;

fino_distribution_t distribution_nuxy;  // Poisson's ratio corresponding to a contraction in y when an extension in x is applied
fino_distribution_t distribution_nuyz;
fino_distribution_t distribution_nuzx;

fino_distribution_t distribution_Gxy;   // shear modulus in direction y on the plane whose normal is x
fino_distribution_t distribution_Gyz;
fino_distribution_t distribution_Gzx;


fino_distribution_t distribution_rho;   // density
fino_distribution_t distribution_fx;    // volumetric load in x
fino_distribution_t distribution_fy;    // volumetric load in y
fino_distribution_t distribution_fz;    // volumetric load in z
fino_distribution_t distribution_alpha; // thermal expansion coefficient
fino_distribution_t distribution_T;     // temperature

// this should be a scalar but we put it as a distribution to check if it is already initialized
fino_distribution_t distribution_T0;    // reference temperature (i.e. no deformations)
double T0;  // this is the scalar used to evaluate it

int uniform_properties;  // flag global para saber si hay que evaluar cada vez o no

// the 6x6 strain-stress matrix
gsl_matrix *C;


double hourglass_2d[] = {+1, -1, +1, -1};
double hourglass_3d[] = {+1, +1, -1, -1, -1, -1, +1, +1,
                         +1, -1, -1, +1, -1, +1, +1, -1,
                         +1, -1, +1, -1, +1, -1, +1, -1,
                         -1, +1, -1, +1, +1, -1, +1, -1};


int fino_bc_process_mechanical(bc_t **bc_pointer, char *name, char *expr, char *equal_sign) {

  int i;
  bc_t *bc = *(bc_pointer);
  bc_t *base_bc = NULL;
  bc_t *tmp = NULL;

  if (strcmp(name, "fixed") == 0) {
    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_fixed;

  } else if (strncmp(name, "mimic(", 6) == 0) {
    char *closing_bracket;
    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_mimic;

    if (name[6] == 'u') {
      bc->dof = 0;
    } else if (name[6] == 'v') {
      bc->dof = 1;
    } else if (name[6] == 'w') {
      if (fino.dimensions < 3) {
        wasora_push_error_message("cannot set 'w' in a two-dimensional problem");
        return WASORA_PARSER_ERROR;
      }
      bc->dof = 2;
    } else {
      wasora_push_error_message("expected either 'u_', 'v_' or 'w_' instead of '%s' in mimic()", name+5);
      return WASORA_PARSER_ERROR;
    }

    if (name[7] != '_') {
      wasora_push_error_message("expected underscore after '%c' instead of '%s' in mimic()", name[6], name+5);
      return WASORA_PARSER_ERROR;
    }

    if ((closing_bracket = strchr(name, ')')) == NULL) {
      wasora_push_error_message("cannot find closing bracket in '%s'", name);
      return WASORA_PARSER_ERROR;
    }
    *closing_bracket = '\0';

    if ((bc->slave = wasora_get_physical_entity_ptr(name+8, fino.mesh)) == NULL) {
      wasora_push_error_message("unknown phyisical entity '%s'", name+8);
      return WASORA_PARSER_ERROR;
    }


  } else if (strcmp(name, "u") == 0 ||
             strcmp(name, "v") == 0 ||
             strcmp(name, "w") == 0) {
    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement;

    bc->dof = name[0]-'u';
//    if (name[0] == 'u') bc->dof = 0;
//    if (name[0] == 'v') bc->dof = 1;
//    if (name[0] == 'w') bc->dof = 2;
    if (bc->dof == 2 && fino.dimensions < 3) {
      wasora_push_error_message("cannot set 'w' in a two-dimensional problem");
      return WASORA_PARSER_ERROR;
    }
    

    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "symmetry") == 0 || strcmp(name, "tangential") == 0) {

    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_symmetry;

  } else if (strcmp(name, "radial") == 0) {

    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_radial;

    // radial might have three extra expressions
    base_bc = bc;
    base_bc->expr = calloc(3, sizeof(expr_t));
    
    tmp = bc;
    while ((bc = bc->next) != NULL) {
      // put back the equal sign, the first time is to parse again
      // the next one is not to break the string
      // the last time is fixed outside the large loop
      if (equal_sign != NULL) {
        *equal_sign = '=';
      }
      fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
      i = -1;
      if (name[0] == 'x') i = 0;
      if (name[0] == 'y') i = 1;
      if (name[0] == 'z') i = 2;
      if (i == -1) {
        wasora_push_error_message("expecting 'x0', 'y0' or 'z0' instead of '%s'", name);
        return WASORA_PARSER_ERROR;
      }             
      wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
      tmp = bc; // esto es para "volver para atras"
    } 

    // bc is now pointing to null, we need to put it back otherwise the foreach loop breaks
    *bc_pointer = tmp;

  } else if (strcmp(name, "0") == 0 || strcmp(name, "implicit") == 0) {
    char *dummy;

    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_constrained;

    // el cuento es asi: aca quisieramos que el usuario escriba algo en funcion
    // de x,y,z pero tambien de u,v y w. Pero u,v,w ya son funciones, asi que no
    // se pueden usar como variables
    // mi solucion: definir variables U,V,W y reemplazar u,v,w por U,V,W en
    // esta expresion

    // TODO: ver que haya un separador (i.e. un operador) antes y despues
    dummy = expr;
    while (*dummy != '\0') {
      if (*dummy == 'u') {
        *dummy = 'U';
      } else if (*dummy == 'v') {
        *dummy = 'V';
      } else if (*dummy == 'w') {
        *dummy = 'W';
      }
      dummy++;
    }

    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcasecmp(name, "tx") == 0 || strcasecmp(name, "ty") == 0 || strcasecmp(name, "tz") == 0 ||
             strcasecmp(name, "fx") == 0 || strcasecmp(name, "fy") == 0 || strcasecmp(name, "fz") == 0) {

    bc->type_math = bc_math_neumann;

    if (isupper(name[0])) {
      // Tx/Fx means force
      bc->type_phys = bc_phys_force;
    } else {
      // tx/fx means stress
      bc->type_phys = bc_phys_stress;
    }

    if (name[1] == 'x') bc->dof = 0;
    if (name[1] == 'y') bc->dof = 1;
    if (name[1] == 'z') bc->dof = 2;

    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "traction") == 0 || strcmp(name, "p") == 0 ||
             strcmp(name, "compression") == 0 || strcmp(name, "P") == 0) {

    bc->type_math = bc_math_neumann;
    if (strcmp(name, "traction") == 0 || strcmp(name, "p") == 0) {
        bc->type_phys = bc_phys_pressure_normal;
    } else if (strcmp(name, "compression") == 0 || strcmp(name, "P") == 0) {
        bc->type_phys = bc_phys_pressure_real;
    }

    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "Mx") == 0 ||
             strcmp(name, "My") == 0 ||
             strcmp(name, "Mz") == 0) {

   
    bc->type_math = bc_math_neumann;
    bc->type_phys = bc_phys_moment;

    // M can take up to 6 expressions
    // asi que las alocamos: Mx My Mz x0 y0 z0 en la primera de las BCs
    // solo ponemos el tipo a la base, las otras no hay que procesarlas en bulk
    base_bc = bc;
    base_bc->expr = calloc(6, sizeof(expr_t));
    
    tmp = bc;
    do {
      // put back the equal sign, the first time is to parse again
      // the next one is not to break the string
      // the last time is fixed outside the large loop
      if (equal_sign != NULL) {
        *equal_sign = '=';
      }
      fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
      i = -1;
      if (name[1] == 'x') i = 0;
      if (name[1] == 'y') i = 1;
      if (name[1] == 'z') i = 2;
      if (name[0] == 'x') i = 3;
      if (name[0] == 'y') i = 4;
      if (name[0] == 'z') i = 5;
      if (i == -1) {
        wasora_push_error_message("expecting 'Mx', 'My', 'Mz', 'x0', 'y0' or 'z0' instead of '%s'", name);
        return WASORA_PARSER_ERROR;
      }             
      wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
      tmp = bc; // esto es para "volver para atras"
    } while ((bc = bc->next) != NULL);

    // bc is now pointing to null, we need to put it back otherwise the foreach loop breaks
    *bc_pointer = tmp;

  } else {
    wasora_push_error_message("unknown boundary condition type '%s'", name);
    return WASORA_PARSER_ERROR;
  }
  
  
  return WASORA_RUNTIME_OK;
}



int fino_break_build_element(element_t *element, int v) {

  static size_t J;            // number of local nodes
  // TODO: hacer un campo descripcion en fino_distribution_t para documentar
 
  static size_t stress_strain_size = 0;
  // matrices needed for the problem formulation
  static gsl_matrix *B = NULL;
  // this vector is needed to compute thermal expansion
  static gsl_vector *et = NULL;

  // intermediate objects
  static gsl_matrix *CB;
  static gsl_vector *Cet;
  
  double *h;
  gsl_matrix *dhdx;
  
  material_t *material;

  double rho;
  double c;
  double alphaDT;
  
  double r_for_axisymmetric = 1.0;
  int j;

  material = (element->physical_entity != NULL)?element->physical_entity->material:NULL;
  
  mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);
  mesh_compute_dhdx_at_gauss(element, v, fino.mesh->integration);

  // convenience variables
  dhdx = element->dhdx[v];
  h = element->type->gauss[fino.mesh->integration].h[v];
  
  // if the stress-strain matrix C is null then we have to allocate some stuff
  // and resolve the special distributions with material properties
  if (C == NULL) {
    wasora_call(fino_distribution_init(&distribution_E,     "E"));
    wasora_call(fino_distribution_init(&distribution_nu,    "nu"));
    
    wasora_call(fino_distribution_init(&distribution_Ex,    "Ex"));
    wasora_call(fino_distribution_init(&distribution_Ey,    "Ey"));
    wasora_call(fino_distribution_init(&distribution_Ez,    "Ez"));

    wasora_call(fino_distribution_init(&distribution_nuxy,  "nuxy"));
    wasora_call(fino_distribution_init(&distribution_nuyz,  "nuyz"));
    wasora_call(fino_distribution_init(&distribution_nuzx,  "nuzx"));

    wasora_call(fino_distribution_init(&distribution_Gxy,   "Gxy"));
    wasora_call(fino_distribution_init(&distribution_Gyz,   "Gyz"));
    wasora_call(fino_distribution_init(&distribution_Gzx,   "Gzx"));    
    
    wasora_call(fino_distribution_init(&distribution_rho,   "rho"));
    wasora_call(fino_distribution_init(&distribution_fx,    "fx"));
    wasora_call(fino_distribution_init(&distribution_fy,    "fy"));
    wasora_call(fino_distribution_init(&distribution_fz,    "fz"));
    wasora_call(fino_distribution_init(&distribution_alpha, "alpha"));
    wasora_call(fino_distribution_init(&distribution_T,     "T"));
    
    // T0 (the temperature at which no expansion occurs) should be a scalar constant
    wasora_call(fino_distribution_init(&distribution_T0,    "T0"));
    if (distribution_T0.defined) {
      if (distribution_T0.variable != NULL) {
        T0 = fino_distribution_evaluate(&distribution_T0, NULL, NULL);
      } else {
        wasora_push_error_message("reference temperature 'T0' has to be a scalar variable");
        return WASORA_RUNTIME_ERROR;
      }  
    } else {
      T0 = 0;
    }
    
    if (distribution_E.defined) {
      
      if (distribution_nu.defined == 0) {
        wasora_push_error_message("cannot find Poisson coefficient 'nu'");
        return WASORA_RUNTIME_ERROR;
      }
      
      fino.material_type = material_type_linear_isotropic;
      
      if (distribution_E.variable != NULL && distribution_nu.variable != NULL) {
        uniform_properties = 1;
      }
      
    } else if (distribution_Ex.defined || distribution_Ey.defined || distribution_Ez.defined) {
      if (distribution_Ex.defined == 0 || distribution_Ey.defined == 0 || distribution_Ez.defined == 0) {
        wasora_push_error_message("Three Young modulus 'Ex', 'Ey' and 'Ez' needed and (at least) one missing");
        return WASORA_RUNTIME_ERROR;
      }

      if (distribution_nuxy.defined == 0 || distribution_nuyz.defined == 0 || distribution_nuzx.defined == 0) {
        wasora_push_error_message("Three Poisson’s ratios 'nuxy', 'nuyz' and 'nuzx' needed and (at least) one missing");
        return WASORA_RUNTIME_ERROR;
      }

      fino.material_type = material_type_linear_orthotropic;

      if (distribution_Ex.variable != NULL && distribution_nuxy.variable != NULL &&
          distribution_Ey.variable != NULL && distribution_nuyz.variable != NULL &&
          distribution_Ez.variable != NULL && distribution_nuzx.variable != NULL) {
        uniform_properties = 1;
      } else {
        wasora_push_error_message("non-uniform orthotropic properties are not supported (yet)");
        return WASORA_RUNTIME_ERROR;
      }
      
    } else {
      wasora_push_error_message("cannot find Young modulus 'E'");
      return WASORA_RUNTIME_ERROR;
    }
    
    if (fino.problem_family == problem_family_modal && distribution_rho.defined == 0) {
      wasora_push_error_message("cannot find density 'rho'");
      return WASORA_RUNTIME_ERROR;
    }
    
    if (fino.problem_kind == problem_kind_full3d) {
      stress_strain_size = 6;
    } else if (fino.problem_kind == problem_kind_axisymmetric) {
      stress_strain_size = 4;
    } else {
      stress_strain_size = 3;
    }

    // stress-strain matrix
    C = gsl_matrix_calloc(stress_strain_size, stress_strain_size);
    // thermal expansion vector
    if (distribution_alpha.defined) {
      et = gsl_vector_calloc(stress_strain_size);
    }  
    
    // if E y nu are scalar variables (i.e. uniform in space) then we compute C once
    // because they do not depend on space and will always be the same for all elements
    if (uniform_properties) {
      wasora_call(fino_break_compute_C(NULL, NULL));
    }
    
  }
  
  // if there are no elemental objects allocated or the volumetric element type
  // is different from the last one we need to re-allocate
  if (J != element->type->nodes) {
    J = element->type->nodes;
    gsl_matrix_free(B);
    B = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);

    gsl_matrix_free(CB);
    CB = gsl_matrix_alloc(stress_strain_size, fino.degrees*J);
    
    if (distribution_alpha.defined) {
      gsl_vector_free(Cet);
      Cet = gsl_vector_alloc(stress_strain_size);
    }  
  }  
  
  // matrix H is the standard one but B is not due to the reduced formulation of the multi-DOF problem
  // i.e the C matrix is 6x6 (or 3x3 in 2D) when it should be 9x9 for the full non-reduced problem
  gsl_matrix_set_zero(B);

  // TODO: use function pointers instead of a big if-elseif
  for (j = 0; j < J; j++) {
    if (fino.problem_kind == problem_kind_full3d) {
      gsl_matrix_set(B, 0, 3*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 3*j+1, gsl_matrix_get(dhdx, j, 1));
      
      gsl_matrix_set(B, 2, 3*j+2, gsl_matrix_get(dhdx, j, 2));
    
      gsl_matrix_set(B, 3, 3*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 3, 3*j+1, gsl_matrix_get(dhdx, j, 0));

      gsl_matrix_set(B, 4, 3*j+1, gsl_matrix_get(dhdx, j, 2));
      gsl_matrix_set(B, 4, 3*j+2, gsl_matrix_get(dhdx, j, 1));

      gsl_matrix_set(B, 5, 3*j+0, gsl_matrix_get(dhdx, j, 2));
      gsl_matrix_set(B, 5, 3*j+2, gsl_matrix_get(dhdx, j, 0));
    
    } else if (fino.problem_kind == problem_kind_axisymmetric) {

      r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
      
      // equation 3.5 AFEM CH.03 sec 3.3.2 pag 3.5
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(dhdx, j, 1));

      if (fino.symmetry_axis == symmetry_axis_y) {
        gsl_matrix_set(B, 2, 2*j+0, h[j]/r_for_axisymmetric);
      } else if (fino.symmetry_axis == symmetry_axis_x) {
        gsl_matrix_set(B, 2, 2*j+1, h[j]/r_for_axisymmetric);
      }
      
      gsl_matrix_set(B, 3, 2*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 3, 2*j+1, gsl_matrix_get(dhdx, j, 0));

    } else  {
      // plane stress and plane strain are the same
      // see equation 14.18 IFEM CH.14 sec 14.4.1 pag 14-11
      gsl_matrix_set(B, 0, 2*j+0, gsl_matrix_get(dhdx, j, 0));
      
      gsl_matrix_set(B, 1, 2*j+1, gsl_matrix_get(dhdx, j, 1));
    
      gsl_matrix_set(B, 2, 2*j+0, gsl_matrix_get(dhdx, j, 1));
      gsl_matrix_set(B, 2, 2*j+1, gsl_matrix_get(dhdx, j, 0));
      
    }
    
    if ((fino.problem_family == problem_family_mechanical) &&
        (distribution_fx.defined != 0 || distribution_fy.defined != 0 || distribution_fz.defined != 0)) {
      // the volumetric force vector
      c = r_for_axisymmetric * element->w[v] * h[j];
      mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
      if (distribution_fx.defined) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+0, c * fino_distribution_evaluate(&distribution_fx, material, element->x[v]));
      }
      if (distribution_fy.defined) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+1, c * fino_distribution_evaluate(&distribution_fy, material, element->x[v]));
      }
      if (distribution_fz.defined && fino.degrees == 3) {
        gsl_vector_add_to_element(fino.bi, fino.degrees*j+2, c * fino_distribution_evaluate(&distribution_fz, material, element->x[v]));
      }
    }
  }
  
  // if E and nu are scalar variables then C is uniform and we already have it
  // but if E or nu are functions or material properties, we need to re-compute C
  if (uniform_properties == 0) {
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    wasora_call(fino_break_compute_C(material, element->x[v]));
  }

  // compute the elemental B'*C*B
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, r_for_axisymmetric, C, B, 0, CB);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v], B, CB, 1.0, fino.Ki);

#ifdef FEM_DUMP  
  printf("Ki(%d,%d) = \n", element->index, v);
  fino_print_gsl_matrix(fino.Ki, stdout);
#endif
  
/*  
  // see if we need to do hourglass control
  if (fino.mesh->integration == integration_reduced && fino.hourglass_epsilon > 0) {  
    gsl_matrix *gamma = NULL;  // lowercase = matrix of plain gamma vectors , one vector per row
    gsl_matrix *Gamma = NULL;  // uppercase = matrix of gamma vectors separated by degree of freedoms
    gsl_matrix_view H;         // hourglass vectors, one vector per row (not to confuse with shape functions)
    gsl_matrix *X = NULL;      // matrix with the coordinates of the nodes
    gsl_matrix *HX = NULL;     // product H*X
    gsl_matrix *BBt = NULL;    // product (B*B') (the elemental matrix is B'*B)
    double *ptr_h = NULL;
    double trace;
    double lambda, mu;
    double eps_tilde;
    
    int n_h = 0;
    int m, alpha;
    
    
    if (element->type->dim == 2 && element->type->nodes == 4) {
      // quad4
      n_h = 1;
      ptr_h = hourglass_2d;
    } else if (element->type->dim == 3 && element->type->nodes == 8) {
      // hex8
      n_h = 4;
      ptr_h = hourglass_3d;
    }
    
    if (n_h != 0) {
      gamma = gsl_matrix_calloc(n_h, element->type->nodes);  
      Gamma = gsl_matrix_calloc(n_h*fino.degrees, element->type->nodes*fino.degrees);  
      H = gsl_matrix_view_array(ptr_h, n_h, element->type->nodes);
      X = gsl_matrix_calloc(element->type->nodes, fino.dimensions);   // either this or the transpose
      HX = gsl_matrix_calloc(n_h, fino.dimensions);
      BBt = gsl_matrix_calloc(fino.dimensions, fino.dimensions);

      // coordinates 
      for (j = 0; j < element->type->nodes; j++) {
        for (m = 0; m < fino.dimensions; m++) {
          gsl_matrix_set(X, j, m, element->node[j]->x[m]);
        }
      }
      
      printf("H %d =\n", element->tag);
      fino_print_gsl_matrix(&H.matrix, stdout);

      printf("X %d =\n", element->tag);
      fino_print_gsl_matrix(X, stdout);

      gsl_matrix_memcpy(gamma, &H.matrix);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &H.matrix, X, 0.0, HX);

      printf("HX %d =\n", element->tag);
      fino_print_gsl_matrix(HX, stdout);

      printf("dhdx %d =\n", element->tag);
      fino_print_gsl_matrix(element->dhdx[v], stdout);
      
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, HX, element->dhdx[v], 1.0, gamma);
      printf("gamma %d =\n", element->tag);
      fino_print_gsl_matrix(gamma, stdout);
      
      for (alpha = 0; alpha < n_h; alpha++) {
        for (m = 0; m < fino.degrees; m++) {
          for (j = 0; j < element->type->nodes; j++) {
            gsl_matrix_set(Gamma, alpha*fino.degrees+m, fino.degrees*j+m, gsl_matrix_get(gamma, alpha, j));
          }
        }
      }
      printf("Gamma %d =\n", element->tag);
      fino_print_gsl_matrix(Gamma, stdout);
  
      // trace of B*B' for normalization
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, element->dhdx[v], element->dhdx[v], 1.0, BBt);
      printf("B*B' %d =\n", element->tag);
      fino_print_gsl_matrix(BBt, stdout);
      
      for (m = 0; m < fino.dimensions; m++) {
        trace += gsl_matrix_get(BBt, m, m);
      }
    
      lambda = E*nu/((1+nu)*(1-2*nu));
      mu = 0.5*E/(1+nu);
      eps_tilde = fino.hourglass_epsilon * (lambda + 2*mu)/3.0 * trace;
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * eps_tilde, Gamma, Gamma, 1.0, fino.Ki);
      

      gsl_matrix_free(Gamma);
    }  
  }
*/  
  // thermal expansion 
  if (distribution_alpha.defined != 0) {
    // note that alpha should be the mean expansion coefficient in the range
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    alphaDT = fino_distribution_evaluate(&distribution_alpha, material, element->x[v]);
    if (alphaDT != 0) {
      alphaDT *= fino_distribution_evaluate(&distribution_T, material, element->x[v])-T0;
      gsl_vector_set(et, 0, alphaDT);
      gsl_vector_set(et, 1, alphaDT);
      if (fino.dimensions == 3) {
        gsl_vector_set(et, 2, alphaDT);
      }  
      gsl_blas_dgemv(CblasTrans, r_for_axisymmetric, C, et, 0, Cet);
      gsl_blas_dgemv(CblasTrans, element->w[v], B, Cet, 1.0, fino.bi);
    }
  }
  
  if (fino.M != NULL) {
    // compute the mass matrix H'*rho*H
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    mesh_compute_H_at_gauss(element, v, fino.degrees, fino.mesh->integration);
    // TODO: check for uniform density
    rho = fino_distribution_evaluate(&distribution_rho, material, element->x[v]);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * rho, element->H[v], element->H[v], 1.0, fino.Mi);
  } 
  
  return WASORA_RUNTIME_OK;
  
}


int fino_break_compute_C(material_t *material, const double *x) {
  
  // TODO: function pointers
  switch (fino.material_type) {
    case material_type_linear_isotropic:
      wasora_call(fino_break_compute_C_isotropic(material, x));
    break; 
    case material_type_linear_orthotropic:
      wasora_call(fino_break_compute_C_orthotropic(material, x)); 
    break;  
  }
  
  return WASORA_RUNTIME_OK;
  
}

int fino_break_compute_C_isotropic(material_t *material, const double *x) {
  
  double E, nu;
  double lambda, mu, lambda2mu;
  
  if ((E = fino_distribution_evaluate(&distribution_E, material, x)) <= 0) {
    wasora_push_error_message("E is not positive (%g)", E);
    return WASORA_RUNTIME_ERROR;
  }

  nu = fino_distribution_evaluate(&distribution_nu, material, x);
  if (nu > 0.499) {
    wasora_push_error_message("nu is greater than 1/2");
    return WASORA_RUNTIME_ERROR;
  } else if (nu < 0) {
    wasora_push_error_message("nu is negative");
    return WASORA_RUNTIME_ERROR;
  }
  
  // using Lame coefficients is more elegant and efficient, but the reference is table 4.3 pag 194 Bathe
  lambda = E*nu/((1+nu)*(1-2*nu));
  mu = 0.5*E/(1+nu);
  lambda2mu = lambda + 2*mu;

  if (fino.problem_kind == problem_kind_full3d) {

    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    gsl_matrix_set(C, 0, 2, lambda);

    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);
    gsl_matrix_set(C, 1, 2, lambda);

    gsl_matrix_set(C, 2, 0, lambda);
    gsl_matrix_set(C, 2, 1, lambda);
    gsl_matrix_set(C, 2, 2, lambda2mu);
  
    gsl_matrix_set(C, 3, 3, mu);
    gsl_matrix_set(C, 4, 4, mu);
    gsl_matrix_set(C, 5, 5, mu);
    
  } else if (fino.problem_kind == problem_kind_plane_stress) {
    
    double c1, c2;
    
    c1 = E/(1-nu*nu);
    c2 = nu * c1;
    gsl_matrix_set(C, 0, 0, c1);
    gsl_matrix_set(C, 0, 1, c2);
    
    gsl_matrix_set(C, 1, 0, c2);
    gsl_matrix_set(C, 1, 1, c1);

    gsl_matrix_set(C, 2, 2, c1*0.5*(1-nu));
    
  } else if (fino.problem_kind == problem_kind_plane_strain) {
    
    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    
    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);

    gsl_matrix_set(C, 2, 2, mu);
    
  } else if (fino.problem_kind == problem_kind_axisymmetric) {
    
    gsl_matrix_set(C, 0, 0, lambda2mu);
    gsl_matrix_set(C, 0, 1, lambda);
    gsl_matrix_set(C, 0, 2, lambda);
    
    gsl_matrix_set(C, 1, 0, lambda);
    gsl_matrix_set(C, 1, 1, lambda2mu);
    gsl_matrix_set(C, 1, 2, lambda);

    gsl_matrix_set(C, 2, 0, lambda);
    gsl_matrix_set(C, 2, 1, lambda);
    gsl_matrix_set(C, 2, 2, lambda2mu);

    gsl_matrix_set(C, 3, 3, mu);
    
  }

  return WASORA_RUNTIME_OK;
}    

int fino_break_compute_C_orthotropic(material_t *material, const double *x) {
  
  double Ex, Ey, Ez;
  double nuxy, nuyz, nuzx;
  double Gxy, Gyz, Gzx;
  
  if ((Ex = fino_distribution_evaluate(&distribution_Ex, material, x)) <= 0) {
    wasora_push_error_message("Ex is not positive (%g)", Ex);
    return WASORA_RUNTIME_ERROR;
  }
  if ((Ey = fino_distribution_evaluate(&distribution_Ey, material, x)) <= 0) {
    wasora_push_error_message("Ey is not positive (%g)", Ey);
    return WASORA_RUNTIME_ERROR;
  }
  if ((Ez = fino_distribution_evaluate(&distribution_Ez, material, x)) <= 0) {
    wasora_push_error_message("Ez is not positive (%g)", Ez);
    return WASORA_RUNTIME_ERROR;
  }

  nuxy = fino_distribution_evaluate(&distribution_nuxy, material, x);
  nuyz = fino_distribution_evaluate(&distribution_nuyz, material, x);
  nuzx = fino_distribution_evaluate(&distribution_nuzx, material, x);
  
  Gxy = fino_distribution_evaluate(&distribution_Gxy, material, x);
  Gyz = fino_distribution_evaluate(&distribution_Gyz, material, x);
  Gzx = fino_distribution_evaluate(&distribution_Gzx, material, x);
  
  if (fino.problem_kind == problem_kind_full3d) {

    gsl_matrix *reducedS = gsl_matrix_alloc(3, 3);  // reduced compliance matrix (only the normal-stress stuff)
    gsl_matrix *reducedC = gsl_matrix_alloc(3, 3);  // reduced stiffness matrix
    
    // fill the compliance matrix
    gsl_matrix_set(reducedS, 0, 0, 1.0/Ex);
    gsl_matrix_set(reducedS, 1, 1, 1.0/Ey);
    gsl_matrix_set(reducedS, 2, 2, 1.0/Ez);

    gsl_matrix_set(reducedS, 0, 1, -nuxy/Ex);
    gsl_matrix_set(reducedS, 1, 0, -nuxy/Ex);
    
    gsl_matrix_set(reducedS, 0, 2, -nuzx/Ez);
    gsl_matrix_set(reducedS, 2, 0, -nuzx/Ez);

    gsl_matrix_set(reducedS, 1, 2, -nuyz/Ey);
    gsl_matrix_set(reducedS, 2, 1, -nuyz/Ey);

    // compute the stiffness by inverting the compliance
    wasora_call(mesh_inverse(reducedS, reducedC));
    
    // fill the full C
    gsl_matrix_set(C, 0, 0, gsl_matrix_get(reducedC, 0, 0));
    gsl_matrix_set(C, 0, 1, gsl_matrix_get(reducedC, 0, 1));
    gsl_matrix_set(C, 0, 2, gsl_matrix_get(reducedC, 0, 2));
    gsl_matrix_set(C, 1, 0, gsl_matrix_get(reducedC, 1, 0));
    gsl_matrix_set(C, 1, 1, gsl_matrix_get(reducedC, 1, 1));
    gsl_matrix_set(C, 1, 2, gsl_matrix_get(reducedC, 1, 2));
    gsl_matrix_set(C, 2, 0, gsl_matrix_get(reducedC, 2, 0));
    gsl_matrix_set(C, 2, 1, gsl_matrix_get(reducedC, 2, 1));
    gsl_matrix_set(C, 2, 2, gsl_matrix_get(reducedC, 2, 2));
    
    gsl_matrix_set(C, 3, 3, Gyz);
    gsl_matrix_set(C, 4, 4, Gzx);
    gsl_matrix_set(C, 5, 5, Gxy);
    
    gsl_matrix_free(reducedC);
    gsl_matrix_free(reducedS);
    

  } else {
    wasora_push_error_message("orthotropic problems are allowed only in 3D");
    return WASORA_RUNTIME_ERROR;
  }

  return WASORA_RUNTIME_OK;
}    


int fino_break_compute_nodal_stresses(element_t *element, int j, double lambda, double mu, double alpha, double *sigmax, double *sigmay, double *sigmaz, double *tauxy, double *tauyz, double *tauzx) {
  
  double dudx = 0;
  double dudy = 0;
  double dudz = 0;

  double dvdx = 0;
  double dvdy = 0;
  double dvdz = 0;
  
  double dwdx = 0;
  double dwdy = 0;
  double dwdz = 0;
  
  double ex = 0;
  double ey = 0;
  double ez = 0;
  double gammaxy = 0;
  double gammayz = 0;
  double gammazx = 0;
  
  double nu = 0;
  double E = 0;
  double DT = 0;
  
  double xi = 0;
  
  // cuter names
  dudx = gsl_matrix_get(element->dphidx_node[j], 0, 0);
  dudy = gsl_matrix_get(element->dphidx_node[j], 0, 1);

  dvdx = gsl_matrix_get(element->dphidx_node[j], 1, 0);
  dvdy = gsl_matrix_get(element->dphidx_node[j], 1, 1);

  if (fino.dimensions == 3) {
    dudz = gsl_matrix_get(element->dphidx_node[j], 0, 2);
    dvdz = gsl_matrix_get(element->dphidx_node[j], 1, 2);

    dwdx = gsl_matrix_get(element->dphidx_node[j], 2, 0);
    dwdy = gsl_matrix_get(element->dphidx_node[j], 2, 1);
    dwdz = gsl_matrix_get(element->dphidx_node[j], 2, 2);
  }

  // strain tensor
  ex = dudx;
  ey = dvdy;

  if (fino.problem_kind == problem_kind_full3d) {
    ez = dwdz;
  } else if (fino.problem_kind == problem_kind_axisymmetric) {
    if (fino.symmetry_axis == symmetry_axis_y) {
      // etheta = u/r
      if (element->node[j]->x[0] > 1e-6) {
        ez = element->node[j]->phi[0]/element->node[j]->x[0];
      }
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      // etheta = v/r
      if (element->node[j]->x[1] > 1e-6) {
        ez = element->node[j]->phi[1]/element->node[j]->x[1];
      }
    }
  }

  gammaxy = dudy + dvdx;
  if (fino.problem_kind == problem_kind_full3d) {
    gammayz = dvdz + dwdy;
    gammazx = dwdx + dudz;
  }

  // so we have derivatives and strains, let's go for the stresses
  // TODO: function pointer
  if (fino.material_type == material_type_linear_isotropic) {
    if (fino.problem_kind != problem_kind_plane_stress) {
    
      // normal stresses
      xi = ex + ey + ez;
      *sigmax = lambda * xi + 2*mu * ex;
      *sigmay = lambda * xi + 2*mu * ey;
      *sigmaz = lambda * xi + 2*mu * ez;  // esta es sigmatheta en axi

      // esfuerzos de corte
      *tauxy =  mu * gammaxy;
      if (fino.dimensions == 3) {
        *tauyz =  mu * gammayz;
        *tauzx =  mu * gammazx;
      } else {
        *tauyz = 0;
        *tauzx = 0;
      }
    } else {
    
      // plane stress es otra milonga
      double c1, c2;
    
      E = mu*(3*lambda + 2*mu)/(lambda+mu);
      nu = lambda / (2*(lambda+mu));
    
      c1 = E/(1-nu*nu);
      c2 = nu * c1;

      *sigmax = c1 * ex + c2 * ey;
      *sigmay = c2 * ex + c1 * ey;
      *tauxy = c1*0.5*(1-nu) * gammaxy;
    
    }  
  } else if (fino.material_type == material_type_linear_orthotropic) {
    
    *sigmax = gsl_matrix_get(C, 0, 0) * ex + gsl_matrix_get(C, 0, 1) * ey + gsl_matrix_get(C, 0, 2) * ez;
    *sigmay = gsl_matrix_get(C, 1, 0) * ex + gsl_matrix_get(C, 1, 1) * ey + gsl_matrix_get(C, 1, 2) * ez;
    *sigmaz = gsl_matrix_get(C, 2, 0) * ex + gsl_matrix_get(C, 2, 1) * ey + gsl_matrix_get(C, 2, 2) * ez;
    *tauxy = gsl_matrix_get(C, 3, 3) * gammaxy;
    *tauyz = gsl_matrix_get(C, 4, 4) * gammayz;
    *tauzx = gsl_matrix_get(C, 5, 5) * gammazx;
    
  }
  
  // subtract the thermal contribution to the normal stresses (see IFEM.Ch30)
  if (alpha != 0) {
    DT = fino_distribution_evaluate(&distribution_T, element->physical_entity->material, element->node[j]->x) - T0;
    E = mu*(3*lambda + 2*mu)/(lambda+mu);
    nu = lambda / (2*(lambda+mu));
    xi = E/(1-2*nu) * alpha * DT;
    
    *sigmax -= xi;
    *sigmay -= xi;
    *sigmaz -= xi;
  }
  
  return WASORA_RUNTIME_OK;
}


int fino_break_compute_stresses(void) {
  
  double sigmax = 0;
  double sigmay = 0;
  double sigmaz = 0;
  double tauxy = 0;
  double tauyz = 0;
  double tauzx = 0;
  double sigma = 0;
  double sigma1 = 0;
  double sigma2 = 0;
  double sigma3 = 0;
  double tresca = 0;
  
  double dsigma_normal = 0;
  double dsigma_shear = 0;
  double delta_sigma = 0;
  
  double displ2 = 0;
  double max_displ2 = 0;
  
  double dudx, dudy, dudz;
  double dvdx, dvdy, dvdz;
  double dwdx, dwdy, dwdz;

  double delta_dudx, delta_dudy, delta_dudz;
  double delta_dvdx, delta_dvdy, delta_dvdz;
  double delta_dwdx, delta_dwdy, delta_dwdz;

  double nu = 0;
  double E = 0;
  double alpha = 0;
  double lambda = 0;
  double mu = 0;
  
  double lambda_max = 0;
  double mu_max = 0;
  
  int j;
  int i, g, m, n;
  int j_global;
  int progress = 0;
  int step = 0; 
  int ascii_progress_chars = 0;

  mesh_t *mesh;
  element_t *element;  
  element_list_item_t *associated_element;
  node_t *node;
  
  // the mesh for "rough" mode is different
  mesh = (fino.rough == 0) ? fino.mesh : fino.mesh_rough;
  
  // number of steps we need to make for progress bar
  if ((step = ceil((double)(mesh->n_elements+mesh->n_nodes)*wasora.nprocs/100.0)) < 1) {
    step = 1;
  }
  
  
  if (fino.sigma->data_value == NULL) {
    // gradient vector
    for (g = 0; g < fino.degrees; g++) {
      for (m = 0; m < fino.dimensions; m++) {
        // derivative of the degree of freedom g with respect to dimension m
        fino_fill_result_function(gradient[g][m]);
        fino_fill_result_function(delta_gradient[g][m]);
      }
    }
    
    // stress tensor
    fino_fill_result_function(sigmax);
    fino_fill_result_function(sigmay);
    fino_fill_result_function(tauxy);
    
    if (fino.dimensions == 3) {
      fino_fill_result_function(tauyz);
      fino_fill_result_function(tauzx);
      fino_fill_result_function(sigmaz);
    }

    // principal stresses
    fino_fill_result_function(sigma1);
    fino_fill_result_function(sigma2);
    fino_fill_result_function(sigma3);

    // von mises
    fino_fill_result_function(sigma);
    fino_fill_result_function(delta_sigma);
    
    // tresca
    fino_fill_result_function(tresca);
  }
  
  // evaluate nu, E and alpha, if they are uniform these values are kept
  if (distribution_nu.variable != NULL ) {
    nu = fino_distribution_evaluate(&distribution_nu, NULL, NULL);
    if (nu > 0.499) {
      wasora_push_error_message("nu is greater than 1/2");
      return WASORA_RUNTIME_ERROR;
    } else if (nu < 0) {
      wasora_push_error_message("nu is negative");
      return WASORA_RUNTIME_ERROR;
    }
  }
  if (distribution_E.variable != NULL) {
    if ((E = fino_distribution_evaluate(&distribution_E, NULL, NULL)) <= 0) {
      wasora_push_error_message("E is not positive (%g)", E);
      return WASORA_RUNTIME_ERROR;
    }
  }
  if (distribution_alpha.variable != NULL) {
    alpha = fino_distribution_evaluate(&distribution_alpha, NULL, NULL);
  }
  
  if (uniform_properties && fino.material_type == material_type_linear_isotropic) {
    // this works for the rest of the routine
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = 0.5*E/(1+nu);
  }  

  // step 1. sweep elements and compute tensors at each node of each element
  for (i = 0; i < mesh->n_elements; i++) {
    if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }

    element = &mesh->element[i];
    if (element->type->dim == fino.dimensions) {

      wasora_call(fino_compute_gradients_at_nodes(mesh, element));

      // if nu, E and/or alpha are not uniform, we need to evalaute them at the nodes
      if (uniform_properties == 0 || distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL) {
        if (element->property_node == NULL) {
          element->property_node = calloc(element->type->nodes, sizeof(double *));
        }  
      }
        
      for (j = 0; j < element->type->nodes; j++) {      
        
        // if nu, E and/or alpha are not uniform, we need to evaluate them at the nodes
        if (fino.material_type == material_type_linear_isotropic &&
            (uniform_properties == 0 ||
             distribution_alpha.function != NULL ||
             distribution_alpha.physical_property != NULL)) {
          
          
          if (uniform_properties == 0) {
            if (distribution_E.variable == NULL) {
              if ((E = fino_distribution_evaluate(&distribution_E, element->physical_entity->material, element->node[j]->x)) <= 0) {
                wasora_push_error_message("E is not positive (%g)", E);
                return WASORA_RUNTIME_ERROR;
              }      
            }
            if (distribution_nu.variable == NULL ) {
              nu = fino_distribution_evaluate(&distribution_nu, element->physical_entity->material, element->node[j]->x);
              if (nu > 0.499) {
                wasora_push_error_message("nu is greater than 1/2");
                return WASORA_RUNTIME_ERROR;
              } else if (nu < 0) {
                wasora_push_error_message("nu is negative");
                return WASORA_RUNTIME_ERROR;
              }
            }
            
            lambda = E*nu/((1+nu)*(1-2*nu));
            mu = 0.5*E/(1+nu);
            
          }
          
          if (distribution_alpha.variable == NULL &&
              (distribution_alpha.function != NULL ||
               distribution_alpha.physical_property != NULL)) {
            alpha = fino_distribution_evaluate(&distribution_alpha, element->physical_entity->material, element->node[j]->x);
          }

          if (element->property_node[j] == NULL) {
            element->property_node[j] = calloc(ISOTROPIC_ELASTIC_PROPERTIES, sizeof(double));
          }  
          element->property_node[j][LAMBDA] = lambda;
          element->property_node[j][MU] = mu;
          element->property_node[j][ALPHA] = alpha;
          
        }
        
        if (fino.rough) {
          // if we are in rough mode, we already have the nodal values
          element->node[j]->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
          element->node[j]->dphidx = element->dphidx_node[j];

          element->node[j]->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));
          fino_break_compute_nodal_stresses(element, j, lambda, mu, alpha, &sigmax, &sigmay, &sigmaz, &tauxy, &tauyz, &tauzx);
          element->node[j]->f[SIGMAX] = sigmax;
          element->node[j]->f[SIGMAY] = sigmay;
          element->node[j]->f[SIGMAZ] = sigmaz;
          element->node[j]->f[TAUXY] = tauxy;
          element->node[j]->f[TAUYZ] = tauyz;
          element->node[j]->f[TAUZX] = tauzx;
        }
      }
    }
  }

  // intermediate step for roughish
  if (fino.roughish) {
    int k;
    double sigmax, sigmay, sigmaz;
    double tauxy, tauyz, tauzx;
    double den;
    node_t *node;
    element_list_item_t *associated_element = NULL;
    element_t *element_smooth;
    element_t *element_rough;
    fino_roughish_avg_t *item, *tmp;
    fino_roughish_avg_t **avg = calloc(fino.mesh->physical_tag_max+1, sizeof(fino_roughish_avg_t *));

    for (j = 0; j < fino.mesh->n_nodes; j++) {
      // para cada nodo de la mall original barro los elementos asociados y hago una linked list
      // para cada physical entity mas o menos asi
      //   nodo original j_global
      //     physical entity 1: elemento i1, nodo local j1
      //                        elemento i2, nodo local j2
      //                        elemento i3, nodo local j3
      //     physical entity 2: elemento i4, nodo local j4
      //                        elemento i5, nodo local j5
      LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element) {
        element_smooth = associated_element->element;
        if (element_smooth != NULL &&
            element_smooth->type->dim == fino.dimensions &&
            element_smooth->physical_entity != NULL) {

          item = calloc(1, sizeof(fino_roughish_avg_t));
          item->smooth_element = element_smooth->tag;
          item->local_node = mesh_compute_local_node_index(element_smooth, j);
          LL_APPEND(avg[element_smooth->physical_entity->tag], item);
        }
      }

      // ahora tengo que buscar todos los ids de los nodo de la malla rough que correspondan
      // al j global de la mall smooth y promediarlos por physical entity
      
      
      for (k = 0; k < fino.mesh->physical_tag_max+1; k++) {
        dudx = 0;
        dudy = 0;
        dvdx = 0;
        dvdy = 0;
        
        sigmax = 0;
        sigmay = 0;
        tauxy = 0;

//        if (fino.dimensions == 3) {
          dudz = 0;
          dvdz = 0;
          
          dwdx = 0;
          dwdy = 0;
          dwdz = 0;
          
          sigmaz = 0;
          tauyz = 0;
          tauzx = 0;
//        }  
        
        den = 0;
        LL_FOREACH(avg[k], item) {
          // la malla rough solamente tiene elementos volumetricos, pero los tiene en orden
          // asi que hay un offset igual al tag del primer elemento rough menos uno
          element_smooth = &fino.mesh->element[item->smooth_element - 1];
          element_rough = &fino.mesh_rough->element[item->smooth_element - fino.mesh_rough->element[0].tag];
          if (element_smooth->tag != element_rough->tag) {
            wasora_push_error_message("numbering mismatch when computing roughish derivatives");
            return WASORA_RUNTIME_ERROR;
          }
          node = element_rough->node[item->local_node];
          
          dudx += gsl_matrix_get(node->dphidx, 0, 0);
          dudy += gsl_matrix_get(node->dphidx, 0, 1);
          dvdx += gsl_matrix_get(node->dphidx, 1, 0);
          dvdy += gsl_matrix_get(node->dphidx, 1, 1);
          
          sigmax += node->f[SIGMAX];
          sigmay += node->f[SIGMAY];
          tauxy += node->f[TAUXY];

          if (fino.dimensions == 3) {
            dudz = gsl_matrix_get(node->dphidx, 0, 2);
            dvdz = gsl_matrix_get(node->dphidx, 1, 2);
          
            dwdx = gsl_matrix_get(node->dphidx, 2, 0);
            dwdy = gsl_matrix_get(node->dphidx, 2, 1);
            dwdz = gsl_matrix_get(node->dphidx, 2, 2);
          
            sigmaz += node->f[SIGMAZ];
            tauyz += node->f[TAUYZ];
            tauzx += node->f[TAUZX];
          }

          den += 1.0;
        }
        
        if (den != 0) {
          dudx /= den;
          dudy /= den;
          dvdx /= den;
          dvdy /= den;
          
          sigmax /= den;
          sigmay /= den;
          tauxy /= den;

          if (fino.dimensions == 3) {
            dudz /= den;
            dvdz /= den;
          
            dwdx /= den;
            dwdy /= den;
            dwdz /= den;
          
            sigmaz /= den;
            tauyz /= den;
            tauzx /= den;
          }
        }  
        
        LL_FOREACH(avg[k], item) {
          element_rough = &fino.mesh_rough->element[item->smooth_element - fino.mesh_rough->element[0].tag];
          node = element_rough->node[item->local_node];
          
          gsl_matrix_set(node->dphidx, 0, 0, dudx);
          gsl_matrix_set(node->dphidx, 0, 1, dudy);
          gsl_matrix_set(node->dphidx, 1, 0, dvdx);
          gsl_matrix_set(node->dphidx, 1, 1, dvdy);
          
          node->f[SIGMAX] = sigmax;
          node->f[SIGMAY] = sigmay;
          node->f[TAUXY] = tauxy;

          if (fino.dimensions == 3) {
            gsl_matrix_set(node->dphidx, 0, 2, dudz);
            gsl_matrix_set(node->dphidx, 1, 2, dvdz);

            gsl_matrix_set(node->dphidx, 2, 0, dwdx);
            gsl_matrix_set(node->dphidx, 2, 1, dwdy);
            gsl_matrix_set(node->dphidx, 2, 2, dwdz);
          
            node->f[SIGMAZ] = sigmaz;
            node->f[TAUYZ] = tauyz;
            node->f[TAUZX] = tauzx;
          }
        }
        
        LL_FOREACH_SAFE(avg[k], item, tmp) {
          LL_DELETE(avg[k], item);
          free(item);
        }  
      }
    }
    free(avg);
  }

  // step 2. sweep nodes of the output mesh (same as input in smooth, rough in rough)
  for (j_global = 0; j_global < mesh->n_nodes; j_global++) {
    if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }

    node = &mesh->node[j_global];
    
    if (fino.rough == 0) {
      double *mean;
      double *current;
      double delta;
      double sum_weight = 0;
      double rel_weight = 0;
      gsl_matrix *m2 = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      int found = 0;
      
      // en iterative si no hacemos esto estamos leakando
      if (node->dphidx == NULL) {
        node->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      } else {
        gsl_matrix_set_zero(node->dphidx);
      }  
      if (node->delta_dphidx == NULL) {
        node->delta_dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      } else {
        gsl_matrix_set_zero(node->delta_dphidx);
      }
      if (node->f == NULL) {
        node->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));
      }
      
      
      n = 0;
      // esto lo hacemos asi para quedarnos con los mayores lambda y mu
      // por si llegamos a tener que calcular dispersion en una interfaz
      lambda_max = 0;
      mu_max = 0;
      LL_FOREACH(mesh->node[j_global].associated_elements, associated_element) {
        element = associated_element->element;
        if (element->dphidx_node != NULL) {
          found = 0;
          for (j = 0; !found && j < element->type->nodes; j++) {
            if (element->node[j]->index_mesh == j_global) {

              n++;
              found = 1;
              sum_weight += element->weight;
              rel_weight = element->weight / sum_weight;
              
              // las derivadas con sus incertezas segun weldford
              for (g = 0; g < fino.degrees; g++) {
                for (m = 0; m < fino.dimensions; m++) {
                  mean = gsl_matrix_ptr(node->dphidx, g, m);
                  current = gsl_matrix_ptr(element->dphidx_node[j], g, m);
                  delta = *current - *mean;
                  *mean += rel_weight * delta;
                  gsl_matrix_add_to_element(m2, g, m, element->weight * delta * ((*current)-(*mean)));
                  gsl_matrix_set(node->delta_dphidx, g, m, sqrt(gsl_matrix_get(m2, g, m)/sum_weight));
                }
              }
              
              if (element->property_node != NULL) {
                lambda = element->property_node[j][LAMBDA];
                mu = element->property_node[j][MU];
                alpha = element->property_node[j][ALPHA];
              }

              fino_break_compute_nodal_stresses(element, j, lambda, mu, alpha, &sigmax, &sigmay, &sigmaz, &tauxy, &tauyz, &tauzx);

              node->f[SIGMAX] += rel_weight*(sigmax - node->f[SIGMAX]);
              node->f[SIGMAY] += rel_weight*(sigmay - node->f[SIGMAY]);
              node->f[SIGMAZ] += rel_weight*(sigmaz - node->f[SIGMAZ]);
              node->f[TAUXY] += rel_weight*(tauxy - node->f[TAUXY]);
              node->f[TAUYZ] += rel_weight*(tauyz - node->f[TAUYZ]);
              node->f[TAUZX] += rel_weight*(tauzx - node->f[TAUZX]);
              
              // necesitamos el peor lambda y el peor mu para calcular el delta de von mises
              if (uniform_properties == 0) {
                if (lambda > lambda_max) {
                  lambda_max = lambda;
                }
                if (mu > mu_max) {
                  mu_max = mu;
                }
              }  
              
            }
          }
        }
      }
      
      gsl_matrix_free(m2);
      
    }
    
    // ya tenemos los promedios ahora, rellenamos las funciones
    // nombres lindos porque los volvemos a usar para calcular la incerteza de von mises
    dudx = gsl_matrix_get(node->dphidx, 0, 0);
    dudy = gsl_matrix_get(node->dphidx, 0, 1);

    dvdx = gsl_matrix_get(node->dphidx, 1, 0);
    dvdy = gsl_matrix_get(node->dphidx, 1, 1);

    if (fino.dimensions == 3) {
      dudz = gsl_matrix_get(node->dphidx, 0, 2);
      dvdz = gsl_matrix_get(node->dphidx, 1, 2);

      dwdx = gsl_matrix_get(node->dphidx, 2, 0);
      dwdy = gsl_matrix_get(node->dphidx, 2, 1);
      dwdz = gsl_matrix_get(node->dphidx, 2, 2);
    } else {
      // con esto sacamos un uninitialized value de valgrind en 2d
      dudz = dvdz = dwdx = dwdy = dwdz = 0;
    }
    
    fino.gradient[0][0]->data_value[j_global] = dudx;
    fino.gradient[0][1]->data_value[j_global] = dudy;
    
    fino.gradient[1][0]->data_value[j_global] = dvdx;
    fino.gradient[1][1]->data_value[j_global] = dvdy;

    if (fino.dimensions == 3) {
      fino.gradient[0][2]->data_value[j_global] = dudz;
      fino.gradient[1][2]->data_value[j_global] = dvdz;

      fino.gradient[2][0]->data_value[j_global] = dwdx;
      fino.gradient[2][1]->data_value[j_global] = dwdy;
      fino.gradient[2][2]->data_value[j_global] = dwdz;
    }
    
    sigmax = node->f[SIGMAX];
    sigmay = node->f[SIGMAY];
    tauxy = node->f[TAUXY];
    
    if (fino.dimensions == 3) {
      sigmaz = node->f[SIGMAZ];
      tauyz = node->f[TAUYZ];
      tauzx = node->f[TAUZX];
    }
    
    fino.sigmax->data_value[j_global] = sigmax;
    fino.sigmay->data_value[j_global] = sigmay;

    fino.tauxy->data_value[j_global] = tauxy;

    if (fino.dimensions == 3) {
      fino.sigmaz->data_value[j_global] = sigmaz;
      
      fino.tauyz->data_value[j_global] = tauyz;
      fino.tauzx->data_value[j_global] = tauzx;
    }

    wasora_call(fino_compute_principal_stress(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx,
                                              &sigma1, &sigma2, &sigma3));
    
    fino.sigma1->data_value[j_global] = sigma1;
    fino.sigma2->data_value[j_global] = sigma2;
    fino.sigma3->data_value[j_global] = sigma3;

    // tresca
    tresca = fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);
    fino.tresca->data_value[j_global] = tresca;

    // von mises
    sigma = fino_compute_vonmises_from_principal(sigma1, sigma2, sigma3);
//    sigma = fino_compute_vonmises_from_strains(lambda, mu, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
    
    // incertezas
    // las incertezas
    if (fino.rough == 0) {
      delta_dudx = gsl_matrix_get(node->delta_dphidx, 0, 0);
      delta_dudy = gsl_matrix_get(node->delta_dphidx, 0, 1);

      delta_dvdx = gsl_matrix_get(node->delta_dphidx, 1, 0);
      delta_dvdy = gsl_matrix_get(node->delta_dphidx, 1, 1);

      if (fino.dimensions == 3) {
        delta_dudz = gsl_matrix_get(node->delta_dphidx, 0, 2);
        delta_dvdz = gsl_matrix_get(node->delta_dphidx, 1, 2);

        delta_dwdx = gsl_matrix_get(node->delta_dphidx, 2, 0);
        delta_dwdy = gsl_matrix_get(node->delta_dphidx, 2, 1);
        delta_dwdz = gsl_matrix_get(node->delta_dphidx, 2, 2);
      }
    
      fino.delta_gradient[0][0]->data_value[j_global] = delta_dudx;
      fino.delta_gradient[0][1]->data_value[j_global] = delta_dudy;
    
      fino.delta_gradient[1][0]->data_value[j_global] = delta_dvdx;
      fino.delta_gradient[1][1]->data_value[j_global] = delta_dvdy;

      if (fino.dimensions == 3) {
        fino.delta_gradient[0][2]->data_value[j_global] = delta_dudz;
        fino.delta_gradient[1][2]->data_value[j_global] = delta_dvdz;

        fino.delta_gradient[2][0]->data_value[j_global] = delta_dwdx;
        fino.delta_gradient[2][1]->data_value[j_global] = delta_dwdy;
        fino.delta_gradient[2][2]->data_value[j_global] = delta_dwdz;
      } else {
        delta_dudz = delta_dvdz = delta_dwdx = delta_dwdy = delta_dwdz = 0;
      }
      
      if (uniform_properties == 0) {
        // esto es para evaluar incertezas en interfaces
        lambda = lambda_max;
        mu = mu_max;
      }  
      
      // derivadas de fino_compute_vonmises_from_strains() calculadas con maxima
      dsigma_normal =  ((-6*(dwdz+dvdy+dudx)*gsl_pow_2(lambda))+2*(dwdz+dvdy+dudx)*lambda
                                       *(3*lambda+4*mu)
                                     -8*(dwdz+dvdy+dudx)*mu*lambda
                                     -4*(dwdz+dvdy)*gsl_pow_2(mu)+8*dudx*gsl_pow_2(mu))
                      /(2*sqrt((-3*gsl_pow_2(dwdz+dvdy+dudx)*gsl_pow_2(lambda))
           +gsl_pow_2(dwdz+dvdy+dudx)*lambda*(3*lambda+4*mu)
           -4*gsl_pow_2(dwdz+dvdy+dudx)*mu*lambda+4*(gsl_pow_2(dwdz)+gsl_pow_2(dvdy)+gsl_pow_2(dudx))*gsl_pow_2(mu)
           -4*(dvdy*dwdz+dudx*dwdz+dudx*dvdy)*gsl_pow_2(mu)
           +3*(gsl_pow_2(dwdy+dvdz)+gsl_pow_2(dwdx+dudz)+gsl_pow_2(dvdx+dudy))*gsl_pow_2(mu)));
      
      dsigma_shear = (3*(dwdy+dvdz)*gsl_pow_2(mu))/sqrt((-3*gsl_pow_2(dwdz+dvdy+dudx)*gsl_pow_2(lambda))
                                  +gsl_pow_2(dwdz+dvdy+dudx)*lambda*(3*lambda+4*mu)
                                  -4*gsl_pow_2(dwdz+dvdy+dudx)*mu*lambda
                                  +4*(gsl_pow_2(dwdz)+gsl_pow_2(dvdy)+gsl_pow_2(dudx))*gsl_pow_2(mu)
                                  -4*(dvdy*dwdz+dudx*dwdz+dudx*dvdy)*gsl_pow_2(mu)
                                  +3*(gsl_pow_2(dwdy+dvdz)+gsl_pow_2(dwdx+dudz)
                                                   +gsl_pow_2(dvdx+dudy))*gsl_pow_2(mu));

      // el libro de taylor de fisica experimental! se me pianta un lagrimoooon
      delta_sigma = sqrt(gsl_pow_2(dsigma_normal*delta_dudx) +
                         gsl_pow_2(dsigma_normal*delta_dvdy) +
                         gsl_pow_2(dsigma_normal*delta_dwdz) +
                         gsl_pow_2(dsigma_shear*delta_dudy) +
                         gsl_pow_2(dsigma_shear*delta_dudz) +
                         gsl_pow_2(dsigma_shear*delta_dvdx) +
                         gsl_pow_2(dsigma_shear*delta_dvdz) +
                         gsl_pow_2(dsigma_shear*delta_dwdx) +
                         gsl_pow_2(dsigma_shear*delta_dwdy));
      fino.delta_sigma->data_value[j_global] = delta_sigma;
    }
        
    if ((fino.sigma->data_value[j_global] = sigma) > wasora_var(fino.vars.sigma_max)) {
      wasora_var(fino.vars.sigma_max) = sigma;
      wasora_var(fino.vars.delta_sigma_max) = delta_sigma;
      
      wasora_var(fino.vars.sigma_max_x) = mesh->node[j_global].x[0];
      wasora_var(fino.vars.sigma_max_y) = mesh->node[j_global].x[1];
      wasora_var(fino.vars.sigma_max_z) = mesh->node[j_global].x[2];
      
      wasora_var(fino.vars.u_at_sigma_max) = mesh->node[j_global].phi[0];
      wasora_var(fino.vars.v_at_sigma_max) = mesh->node[j_global].phi[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_sigma_max) = mesh->node[j_global].phi[2];
      }
    }
    
    displ2 = 0;
    for (g = 0; g < fino.degrees; g++) {
      displ2 += gsl_pow_2(mesh->node[j_global].phi[g]);
    }
    
    // el >= es porque si en un parametrico se pasa por cero tal vez no se actualice displ_max
    if (displ2 >= max_displ2) {
      max_displ2 = displ2;
      wasora_var(fino.vars.displ_max) = sqrt(displ2);
      wasora_var(fino.vars.displ_max_x) = mesh->node[j_global].x[0];
      wasora_var(fino.vars.displ_max_y) = mesh->node[j_global].x[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.displ_max_z) = mesh->node[j_global].x[2];
      }
      
      wasora_var(fino.vars.u_at_displ_max) = mesh->node[j_global].phi[0];
      wasora_var(fino.vars.v_at_displ_max) = mesh->node[j_global].phi[1];
      if (fino.dimensions == 3) {
        wasora_var(fino.vars.w_at_displ_max) = mesh->node[j_global].phi[2];
      }
    }
  }

  if (fino.progress_ascii == PETSC_TRUE) {
    if (wasora.nprocs == 1) {
      while (ascii_progress_chars++ < 100) {
        printf(CHAR_PROGRESS_GRADIENT);
      }
    }
    if (wasora.rank == 0) {
      printf("\n");  
      fflush(stdout);
    }  
  }

  return WASORA_RUNTIME_OK;
}


int fino_break_set_neumann(element_t *element, bc_t *bc) {
  
  double x0, y0, z0;
  double Ix, Iy, Iz;
  double r_for_axisymmetric;
  int i, j, v;
  

  if (bc->type_phys ==  bc_phys_force) {
    if (element->physical_entity->volume == 0) {
      wasora_push_error_message("physical group '%s' has zero volume", element->physical_entity->name);
      return WASORA_RUNTIME_ERROR;
    }
  } else if (bc->type_phys == bc_phys_pressure_normal || bc->type_phys == bc_phys_pressure_real) {
    if (((fino.dimensions - element->type->dim) != 1)) {
      wasora_push_error_message("pressure BCs can only be applied to surfaces");
      return WASORA_RUNTIME_ERROR;
    }
  } else if (bc->type_phys == bc_phys_moment)   {
    
    double xix, xiy, xiz;
    physical_entity_t *physical_entity = element->physical_entity;
    element_t *tmp_element;
  
    if ((fino.dimensions-element->type->dim != 1)) {
      wasora_push_error_message("moment BCs can only be applied to surfaces");
      return WASORA_RUNTIME_ERROR;
    }
    
    // el centro
    if (bc->expr[4].n_tokens != 0) {
      x0 = wasora_evaluate_expression(&bc->expr[3]);
    } else {
      x0 = physical_entity->cog[0];
    }
    if (bc->expr[5].n_tokens != 0) {
      y0 = wasora_evaluate_expression(&bc->expr[4]);
    } else {
      y0 = physical_entity->cog[1];
    }
    if (bc->expr[6].n_tokens != 0) {
      z0 = wasora_evaluate_expression(&bc->expr[5]);
    } else {
      z0 = physical_entity->cog[2];
    }

    // the polar moment of intertia
    Ix = 0;
    Iy = 0;
    Iz = 0;
    for (i = 0; i < physical_entity->n_elements; i++) {
      tmp_element = &fino.mesh->element[physical_entity->element[i]];
      for (v = 0; v < element->type->gauss[fino.mesh->integration].V; v++) {
        mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);

        xix = 0;
        xiy = 0;
        xiz = 0;
        for (j = 0; j < element->type->nodes; j++) {
          xix += element->type->gauss[fino.mesh->integration].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[0]);
          xiy += element->type->gauss[fino.mesh->integration].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[1]);
          xiz += element->type->gauss[fino.mesh->integration].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[2]);
        }

        Ix += element->w[v] * xix;
        Iy += element->w[v] * xiy;
        Iz += element->w[v] * xiz;
      }
    }
  }  
  
  
  for (v = 0; v < element->type->gauss[fino.mesh->integration].V; v++) {
    mesh_compute_integration_weight_at_gauss(element, v, fino.mesh->integration);
    mesh_compute_H_at_gauss(element, v, fino.degrees, fino.mesh->integration);
    mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    mesh_update_coord_vars(element->x[v]);
    r_for_axisymmetric = fino_compute_r_for_axisymmetric(element, v);
    
    if (bc->type_phys == bc_phys_stress) {
      
      gsl_vector_set(fino.Nb, bc->dof, wasora_evaluate_expression(&bc->expr[0]));
      
    } else if (bc->type_phys ==  bc_phys_force) {
      
      gsl_vector_set(fino.Nb, bc->dof, wasora_evaluate_expression(&bc->expr[0])/element->physical_entity->volume);
      
    } else if (bc->type_phys == bc_phys_pressure_normal || bc->type_phys == bc_phys_pressure_real) {
      
      double p;
      
      // la p chica es la proyeccion del vector tension sobre la normal, lo que uno espera en matematica
      p = wasora_evaluate_expression(&bc->expr[0]);
      // la P grande es presion positiva cuando comprime, como lo que uno espera en ingenieria      
      if (bc->type_phys == bc_phys_pressure_real) {
        p = -p;
      }
      
      gsl_vector_set(fino.Nb, 0, wasora_var_value(wasora_mesh.vars.nx) * p);
      gsl_vector_set(fino.Nb, 1, wasora_var_value(wasora_mesh.vars.ny) * p);
      if (fino.dimensions == 3) {
        gsl_vector_set(fino.Nb, 2, wasora_var_value(wasora_mesh.vars.nz) * p);
      }
      
    } else if (bc->type_phys == bc_phys_moment)   {
      
      double dx, dy, dz;
      double d, theta, F;
      
      dx = element->x[v][0] - x0;
      dy = element->x[v][1] - y0;
      dz = element->x[v][2] - z0;
    
      // los tres primeros tienen las componentes Mx My y Mz
      if (bc->expr[0].n_tokens != 0) {
        d = gsl_hypot(dy, dz);
        theta = atan2(dy, dz);
        F = wasora_evaluate_expression(&bc->expr[0]) * 0.5 * d / (Iy+Iz);
        // ty = cos(theta)*dz
        gsl_vector_add_to_element(fino.Nb, 1, -F*cos(theta));
        // tz = -sin(theta)*dy
        gsl_vector_add_to_element(fino.Nb, 2, +F*sin(theta));
      }
    
      if (bc->expr[1].n_tokens != 0) {
        d = gsl_hypot(dx, dz);
        theta = atan2(dx, dz);
        F = wasora_evaluate_expression(&bc->expr[1]) * 0.5 * d / (Ix+Iz);
        gsl_vector_add_to_element(fino.Nb, 0, -F*cos(theta));
        gsl_vector_add_to_element(fino.Nb, 2, +F*sin(theta));
      }
    
      if (bc->expr[2].n_tokens != 0) {
        d = gsl_hypot(dx, dy);
        theta = atan2(dx, dy);
        F = wasora_evaluate_expression(&bc->expr[2]) * 0.5 * d / (Ix+Iy);
        gsl_vector_add_to_element(fino.Nb, 0, -F*cos(theta));
        gsl_vector_add_to_element(fino.Nb, 1, +F*sin(theta));
      }
    
    }  

    gsl_blas_dgemv(CblasTrans, r_for_axisymmetric*element->w[v], element->H[v], fino.Nb, 1.0, fino.bi); 
    
    
  }

  return WASORA_RUNTIME_OK;
}

int fino_compute_principal_stress(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx, double *sigma1, double *sigma2, double *sigma3) {
  
  double I1, I2, I3;
  double c1, c2, c3;
  double phi;
  
  // stress invariants
  // https://en.wikiversity.org/wiki/Principal_stresses
  I1 = sigmax + sigmay + sigmaz;
  I2 = sigmax*sigmay + sigmay*sigmaz + sigmaz*sigmax - gsl_pow_2(tauxy) - gsl_pow_2(tauyz) - gsl_pow_2(tauzx);
  I3 = sigmax*sigmay*sigmaz - sigmax*gsl_pow_2(tauyz) - sigmay*gsl_pow_2(tauzx) - sigmaz*gsl_pow_2(tauxy) + 2*tauxy*tauyz*tauzx;

  // principal stresses
  c1 = sqrt(fabs(gsl_pow_2(I1) - 3*I2));
  phi = 1.0/3.0 * acos((2.0*gsl_pow_3(I1) - 9.0*I1*I2 + 27.0*I3)/(2.0*gsl_pow_3(c1)));
  if (isnan(phi)) {
    phi = 0;
  }

  c2 = I1/3.0;
  c3 = 2.0/3.0 * c1;
  if (sigma1 != NULL) {
    *sigma1 = c2 + c3 * cos(phi);
  }
  if (sigma2 != NULL) {
    *sigma2 = c2 + c3 * cos(phi - 2.0*M_PI/3.0);
  }
  if (sigma3 != NULL) {
    *sigma3 = c2 + c3 * cos(phi - 4.0*M_PI/3.0);
  }
  
  return WASORA_RUNTIME_OK;
  
}

double fino_compute_vonmises_from_principal(double sigma1, double sigma2, double sigma3) {
  
  return sqrt(0.5*(gsl_pow_2(sigma1-sigma2) + gsl_pow_2(sigma2-sigma3) + gsl_pow_2(sigma3-sigma1)));
  
}

double fino_compute_vonmises_from_stress_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {
  
  return sqrt(0.5*(gsl_pow_2(sigmax-sigmay) + gsl_pow_2(sigmay-sigmaz) + gsl_pow_2(sigmaz-sigmax) +
                       6.0 * (gsl_pow_2(tauxy) + gsl_pow_2(tauyz) + gsl_pow_2(tauzx))));
  
}

double fino_compute_vonmises_from_strains(double lambda, double mu,
                                          double dudx, double dudy, double dudz,
                                          double dvdx, double dvdy, double dvdz,
                                          double dwdx, double dwdy, double dwdz) {
  
  return sqrt(lambda*gsl_pow_2(dudx+dvdy+dwdz)*(3*lambda+4*mu)
      + 4*gsl_pow_2(mu)*(gsl_pow_2(dudx)+gsl_pow_2(dvdy)+gsl_pow_2(dwdz))
      - (3*gsl_pow_2(lambda*(dudx+dvdy+dwdz)) + 4*lambda*mu*gsl_pow_2(dudx+dvdy+dwdz)
      + 4*gsl_pow_2(mu)*(dudx*dvdy+dvdy*dwdz+dwdz*dudx))
      + 3*(gsl_pow_2(mu)*(gsl_pow_2(dudy+dvdx)+gsl_pow_2(dvdz+dwdy)+gsl_pow_2(dwdx+dudz))));
  
}


double fino_compute_tresca_from_principal(double sigma1, double sigma2, double sigma3) {

  double S12 = fabs(sigma1-sigma2);
  double S23 = fabs(sigma2-sigma3);
  double S31 = fabs(sigma3-sigma1);

  if (S31 >= S12 && S31 >= S23) {
    return S31;
  } else if (S12 >= S23 && S12 >= S31) {
    return S12;
  } else if (S23 >= S12 && S23 >= S31) {
    return S23;
  }

  return 0;

}

double fino_compute_tresca_from_stress_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {

  double sigma1, sigma2, sigma3;

  wasora_call(fino_compute_principal_stress(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx, &sigma1, &sigma2, &sigma3));
  return fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);

}


int fino_compute_strain_energy(void) {

  PetscScalar e;
  Vec Kphi;
  
  petsc_call(VecDuplicate(fino.phi, &Kphi));
  petsc_call(MatMult(fino.K_nobc, fino.phi, Kphi));
  petsc_call(VecDot(fino.phi, Kphi, &e));
  wasora_var(fino.vars.strain_energy) = 0.5*e;
  
  if (fino.problem_kind == problem_kind_axisymmetric) {
    wasora_var(fino.vars.strain_energy) *= 2*M_PI;
  }
  petsc_call(VecDestroy(&Kphi));

  if (wasora_var_value(wasora_special_var(in_static))) {
    if ((int)(wasora_var(wasora_special_var(step_static))) == 1) {
      *fino.vars.strain_energy->initial_static = *fino.vars.strain_energy->value;
    }
    *fino.vars.strain_energy->initial_transient = *fino.vars.strain_energy->value;
  }
  
  return WASORA_RUNTIME_OK;
}
