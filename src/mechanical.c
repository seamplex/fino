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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#define LAMBDA             0
#define MU                 1
#define ALPHA              2
#define ELASTIC_PROPERTIES 3

#define SIGMAX             0
#define SIGMAY             1
#define SIGMAZ             2
#define TAUXY              3
#define TAUYZ              4
#define TAUZX              5
#define ELASTIC_FUNCTIONS  6


#define fino_fill_result_function(fun_nam) {\
        fino.fun_nam->mesh = fino.rough==0?fino.mesh:fino.mesh_rough; \
        fino.fun_nam->data_argument = fino.fun_nam->mesh->nodes_argument;   \
        fino.fun_nam->data_size = fino.fun_nam->mesh->n_nodes; \
        fino.fun_nam->data_value = calloc(fino.fun_nam->mesh->n_nodes, sizeof(double));}
        
fino_distribution_t distribution_E;     // modulo de young
fino_distribution_t distribution_nu;    // coef de poisson
fino_distribution_t distribution_rho;   // densidad
fino_distribution_t distribution_fx;    // fuerza volumetrica en x
fino_distribution_t distribution_fy;    // fuerza volumetrica en y
fino_distribution_t distribution_fz;    // fuerza volumetrica en z
fino_distribution_t distribution_alpha; // coeficiente de expansion termica
fino_distribution_t distribution_T;     // temperatura

// este es un escalar pero lo ponemos como dist para ver si ya lo inicializamos
fino_distribution_t distribution_T0;    // temperatura de referencia (i.e. sin deformacion)
double T0;  // este es el escalar que usamos para evaluar

int uniform_properties;  // flag global para saber si hay que evaluar cada vez o no



#undef  __func__
#define __func__ "fino_bc_process_mechanical"
int fino_bc_process_mechanical(bc_t *bc, char *name, char *expr) {

  int i;
  bc_t *base_bc = NULL;
  
  PetscFunctionBegin;
  
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

    if (name[0] == 'u') bc->dof = 0;
    if (name[0] == 'v') bc->dof = 1;
    if (name[0] == 'w') bc->dof = 2;

    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "symmetry") == 0 || strcmp(name, "tangential") == 0) {

    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_displacement_symmetry;

  } else if (strcmp(name, "radial") == 0 ||
             (base_bc != NULL && base_bc->type_phys == bc_phys_displacement_radial &&
               (strcmp(name, "x0") == 0 ||
                strcmp(name, "y0") == 0 ||
                strcmp(name, "z0") == 0))) {

    // radial puede tener tres expresiones
    // asi que las alocamos: x0 y0 z0 en la primera de las BCs
    if (base_bc == NULL) {
      bc->type_math = bc_math_dirichlet;
      bc->type_phys = bc_phys_displacement_radial;

      base_bc = bc;
      base_bc->expr = calloc(3, sizeof(expr_t));
    }

    if (strcmp(name, "radial") != 0) {
      i = -1;  // si alguna no aparece es cero (que por default es el baricentro de la entidad)
      if (name[0] == 'x') i = 0;
      if (name[0] == 'y') i = 1;
      if (name[0] == 'z') i = 2;
      if (i == -1) {
        wasora_push_error_message("expecting 'x0', 'y0' or 'z0' instead of '%s'", name);
        return WASORA_PARSER_ERROR;
      }
      wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));          
    }

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
             strcmp(name, "Mz") == 0 ||
             (base_bc != NULL && base_bc->type_phys == bc_phys_moment &&
               (strcmp(name, "x0") == 0 ||
                strcmp(name, "y0") == 0 ||
                strcmp(name, "z0") == 0))) {

    // M necesita seis expresiones
    // asi que las alocamos: Mx My Mz x0 y0 z0 en la primera de las BCs
    if (base_bc == NULL) {
      // solo ponemos el tipo a la base, las otras no hay que procesarlas en bulk
      bc->type_math = bc_math_neumann;
      bc->type_phys = bc_phys_moment;

      base_bc = bc;
      base_bc->expr = calloc(6, sizeof(expr_t));
    }

    i = -1;  // si alguna no aparece es cero (que por default es el baricentro de la entidad)
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

  } else {
    wasora_push_error_message("unknown boundary condition type '%s'", name);
    PetscFunctionReturn(WASORA_PARSER_ERROR);
  }
  
  
  return WASORA_RUNTIME_OK;
}



#undef  __func__
#define __func__ "fino_break_build_element"
int fino_break_build_element(element_t *element, int v) {

  static size_t J;            // cantidad de nodos locales
  // TODO: hacer un campo descripcion en fino_distribution_t para documentar
 
  static size_t stress_strain_size = 0;
  // matrices de la formulacion del problema
  static gsl_matrix *C = NULL;
  static gsl_matrix *B = NULL;
  // vector para calcular las tensiones termicas
  static gsl_vector *et = NULL;

  // matriz intermedia
  static gsl_matrix *CB;
  // vector intermedio
  static gsl_vector *Cet;
  
  double *h;
  gsl_matrix *dhdx;
  
  material_t *material;

  double E, nu;
  double rho;
  double c;
  double alphaDT;
  
  double r_for_axisymmetric = 1.0;
  int j;

  PetscFunctionBegin;
  
  material = (element->physical_entity != NULL)?element->physical_entity->material:NULL;
  
  mesh_compute_integration_weight_at_gauss(element, v);
  mesh_compute_dhdx_at_gauss(element, v);

  dhdx = element->dhdx[v];
  h = element->type->gauss[GAUSS_POINTS_CANONICAL].h[v];
  
  // if the stress-strain matrix C is null then we have to allocate some stuff
  // and resolve the special distributions with material properties
  if (C == NULL) {
    wasora_call(fino_distribution_init(&distribution_E, "E"));
    wasora_call(fino_distribution_init(&distribution_nu, "nu"));
    wasora_call(fino_distribution_init(&distribution_rho, "rho"));
    wasora_call(fino_distribution_init(&distribution_fx, "fx"));
    wasora_call(fino_distribution_init(&distribution_fy, "fy"));
    wasora_call(fino_distribution_init(&distribution_fz, "fz"));
    wasora_call(fino_distribution_init(&distribution_alpha, "alpha"));
    wasora_call(fino_distribution_init(&distribution_T, "T"));
    
    // T0 (the temperature at wich no expansion occurs) might be a function or a constant
    wasora_call(fino_distribution_init(&distribution_T0, "T0"));
    if (distribution_T0.defined) {
      T0 = fino_distribution_evaluate(&distribution_T0, NULL, NULL);
    } else {
      T0 = 0;
    }
    
    if (distribution_E.defined == 0) {
      wasora_push_error_message("cannot find Young modulus 'E'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    } else if (distribution_nu.defined == 0) {
      wasora_push_error_message("cannot find Poisson coefficient 'nu'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    if (fino.math_type == math_type_eigen && distribution_rho.defined == 0) {
      wasora_push_error_message("cannot find density 'rho'");
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
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
    if (distribution_E.variable != NULL && distribution_nu.variable != NULL) {
      uniform_properties = 1;
      if ((E = fino_distribution_evaluate(&distribution_E, material, NULL)) <= 0) {
        wasora_push_error_message("E is not positive (%g)", E);
        return WASORA_RUNTIME_ERROR;
      }

      nu = fino_distribution_evaluate(&distribution_nu, material, NULL);
      if (nu > 0.499) {
        wasora_push_error_message("nu is greater than 1/2");
        return WASORA_RUNTIME_ERROR;
      } else if (nu < 0) {
        wasora_push_error_message("nu is negative");
        return WASORA_RUNTIME_ERROR;
      }
      wasora_call(fino_break_compute_C(C, E, nu));
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
  
  // matrix H is the standard one by B is not due to the reduced formulation of the multi-DOF problem
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
      mesh_compute_x_at_gauss(element, v);
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
  
  // if E and nu are scalar variables C is uniform and we already have it
  // but if E or nu are functions or material properties, we need to re-compute C
  if (uniform_properties == 0) {
    mesh_compute_x_at_gauss(element, v);
    wasora_call(fino_break_compute_C(C,
        fino_distribution_evaluate(&distribution_E,  material, element->x[v]),
        fino_distribution_evaluate(&distribution_nu, material, element->x[v])));
  }

  // compute the elementals Bt*C*B
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, r_for_axisymmetric, C, B, 0, CB);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v], B, CB, 1.0, fino.Ki);

  // thermal expansion 
  if (distribution_alpha.defined != 0) {
    // note that alpha should be the mean expansion coefficient in the range
    mesh_compute_x_at_gauss(element, v);
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
  
  if (fino.has_mass) {
    // calculamos la matriz de masa Ht*rho*H
    mesh_compute_x_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    rho = fino_distribution_evaluate(&distribution_rho, material, element->x[v]);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, element->w[v] * r_for_axisymmetric * rho, element->H[v], element->H[v], 1.0, fino.Mi);
  } 
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_C"
int fino_break_compute_C(gsl_matrix *C, double E, double nu) {
  
  double lambda, mu, lambda2mu;
  
  PetscFunctionBegin;
  
  // esto es mas elegante (y eficiente) pero la referencia posta es tabla 4.3 pag 194 Bathe
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

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}    

#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_nodal_stresses"
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
  
  // nombres lindos
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

  // el tensor de deformaciones
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

  // ya tenemos derivadas y strains, ahora las tensiones
  
  // tensiones normales
  if (fino.problem_kind != problem_kind_plane_stress) {
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
  
  // subtract the thermal contribution to the normal stresses (see IFEM.Ch30)
  if (alpha != 0) {
    DT = fino_distribution_evaluate(&distribution_T, element->physical_entity->material, fino.mesh->node[j].x) - T0;
    E = mu*(3*lambda + 2*mu)/(lambda+mu);
    nu = lambda / (2*(lambda+mu));
    xi = E/(1-2*nu) * alpha * DT;

    *sigmax -= xi;
    *sigmay -= xi;
    *sigmaz -= xi;
  }
  
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_compute_stresses"
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
  
  int v, V;
  int i, j, g, m, n;
  int j_global, j_global_prime;
  int j_local_prime;
  int progress = 0;
  int step = 0; 
  int ascii_progress_chars = 0;

  mesh_t *mesh;
  node_relative_t *parent;
  element_t *element;  
  element_list_item_t *associated_element;
  node_t *node;
  
  PetscFunctionBegin;
  
  mesh = (fino.rough == 0) ? fino.mesh : fino.mesh_rough;
  
  
  // depende de si es gauss o node
  if (fino.gradient_evaluation == gradient_gauss_extrapolated) {
    step = ceil((double)(2*mesh->n_elements+mesh->n_nodes)*wasora.nprocs/100.0);
  } else if (fino.gradient_evaluation == gradient_at_nodes) {
    step = ceil((double)(mesh->n_elements+mesh->n_nodes)*wasora.nprocs/100.0);
  }  
  if (step < 1) {
    step = 1;
  }
  
  
  if (fino.sigma->data_value == NULL) {
    // gradiente
    for (g = 0; g < fino.degrees; g++) {
      for (m = 0; m < fino.dimensions; m++) {
        // derivada del grado de libertad g con respecto a la dimension m
        fino_fill_result_function(gradient[g][m]);
        fino_fill_result_function(delta_gradient[g][m]);
      }
    }
    
    // tensor de tensiones
    fino_fill_result_function(sigmax);
    fino_fill_result_function(sigmay);
    fino_fill_result_function(tauxy);
    
    if (fino.dimensions == 3) {
      fino_fill_result_function(tauyz);
      fino_fill_result_function(tauzx);
      fino_fill_result_function(sigmaz);
    }

    // tensiones principales
    fino_fill_result_function(sigma1);
    fino_fill_result_function(sigma2);
    fino_fill_result_function(sigma3);

    // von mises
    fino_fill_result_function(sigma);
    fino_fill_result_function(delta_sigma);
    
    // tresca
    fino_fill_result_function(tresca);
  }
  
  // evaluamos nu, E y alpha, si son uniformes esto ya nos sirve para siempre
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
  
  if (uniform_properties) {
    // esto ya sirve para toda la cosecha  
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = 0.5*E/(1+nu);
  }  


  // paso 0 (solo si es gauss extrapolate): calculamos las derivadas en los puntos de gauss
  if (fino.gradient_evaluation == gradient_gauss_extrapolated) {
    for (i = 0; i < mesh->n_elements; i++) {
      if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
        printf(CHAR_PROGRESS_GRADIENT);  
        fflush(stdout);
        ascii_progress_chars++;
      }
      
      element = &mesh->element[i];
      if (element->type->dim == fino.dimensions){
        
        V = element->type->gauss[GAUSS_POINTS_CANONICAL].V;
        element->dphidx_gauss = calloc(V, sizeof(gsl_matrix *));
        
        for (v = 0; v < V; v++) {
        
          element->dphidx_gauss[v] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
          mesh_compute_dhdx_at_gauss(element, v);

          // aca habria que hacer una matriz con los phi globales
          // (de j y g, que de paso no depende de v asi que se podria hacer afuera del for de v)
          // y ver como calcular la matriz dphidx como producto de dhdx y esta matriz
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              for (j = 0; j < element->type->nodes; j++) {
                j_global = element->node[j]->index_mesh;
                gsl_matrix_add_to_element(element->dphidx_gauss[v], g, m, gsl_matrix_get(element->dhdx[v], j, m) * mesh->node[j_global].phi[g]);
              }
            }
          }
        }
      }
    }  
  }  
  
  // paso 1. barremos elementos y calculamos los tensores en cada nodo de cada elemento
  for (i = 0; i < mesh->n_elements; i++) {
    if ((fino.progress_ascii == PETSC_TRUE) && (progress++ % step) == 0) {
      printf(CHAR_PROGRESS_GRADIENT);  
      fflush(stdout);
      ascii_progress_chars++;
    }
    
    element = &mesh->element[i];
    if (element->type->dim == fino.dimensions) {

      element->dphidx_node = calloc(element->type->nodes, sizeof(gsl_matrix *));
      V = element->type->gauss[GAUSS_POINTS_CANONICAL].V;

      if (fino.rough == 0) {
        if (fino.gradient_element_weight == gradient_weight_volume) {
          // ok, this looks odd but still I’d rather use C than C++
          element->type->element_volume(element);
          element->weight = element->volume;
        } else if (fino.gradient_element_weight == gradient_weight_quality) {
          mesh_compute_quality(mesh, element);
          element->weight = element->quality;
        } else if (fino.gradient_element_weight == gradient_weight_volume_times_quality) {
          element->type->element_volume(element);
          mesh_compute_quality(mesh, element);
          element->weight = element->volume*GSL_MAX(element->quality, 1);;
        } else if (fino.gradient_element_weight == gradient_weight_flat) {
          element->weight = 1;
        }
      }
      
      // si nu, E y/o alpha no son uniformes, los tenemos que evaluar en los nodos
      if (uniform_properties == 0 || distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL) {
        element->property_node = calloc(element->type->nodes, sizeof(double *));
      }
      
      for (j = 0; j < element->type->nodes; j++) {
      
        element->dphidx_node[j] = gsl_matrix_calloc(fino.degrees, fino.dimensions);
        j_global = element->node[j]->index_mesh;
        
        if (fino.gradient_evaluation == gradient_gauss_extrapolated && j < V && element->type->gauss[GAUSS_POINTS_CANONICAL].extrap != NULL) {
          gsl_vector *inner = gsl_vector_alloc(V);
          gsl_vector *outer = gsl_vector_alloc(V);
          
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              
              for (v = 0; v < V; v++) {
                gsl_vector_set(inner, v, gsl_matrix_get(element->dphidx_gauss[v], g, m));
              }  
                
              gsl_blas_dgemv(CblasNoTrans, 1.0, element->type->gauss[GAUSS_POINTS_CANONICAL].extrap, inner, 0, outer);
              gsl_matrix_set(element->dphidx_node[j], g, m, gsl_vector_get(outer, j));
              
            }
          }
          gsl_vector_free(inner);
          gsl_vector_free(outer);
          
        } else if (element->type->node_parents != NULL && element->type->node_parents[j] != NULL && fino.gradient_highorder_nodes == gradient_average) {
          // promedio de padres
          double den = 0;
          LL_FOREACH(element->type->node_parents[j], parent) {
            den += 1.0;
            for (g = 0; g < fino.degrees; g++) {
              for (m = 0; m < fino.dimensions; m++) {
                gsl_matrix_add_to_element(element->dphidx_node[j], g, m, gsl_matrix_get(element->dphidx_node[parent->index], g, m));
              }
            }  
          }
          
          gsl_matrix_scale(element->dphidx_node[j], 1.0/den);          
          
        } else {
          
          // evaluacion directa en los nodos
          gsl_matrix *dhdx = gsl_matrix_calloc(element->type->nodes, fino.dimensions); // esto esta al vesre
          mesh_compute_dhdx(element, element->type->node_coords[j], NULL, dhdx);
          
          // las nueve derivadas (o menos)
          // TODO: como arriba, aunque hay que pelar ojo si hay menos DOFs
          for (g = 0; g < fino.degrees; g++) {
            for (m = 0; m < fino.dimensions; m++) {
              for (j_local_prime = 0; j_local_prime < element->type->nodes; j_local_prime++) {
                j_global_prime = element->node[j_local_prime]->index_mesh;
                gsl_matrix_add_to_element(element->dphidx_node[j], g, m, gsl_matrix_get(dhdx, j_local_prime, m) * mesh->node[j_global_prime].phi[g]);
              }
            }
          }
          gsl_matrix_free(dhdx);
        }
        

        
        // si nu, E y/o alpha no son uniformes, los tenemos que evaluar en los nodos
        if (uniform_properties == 0 || distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL) {
          
          element->property_node[j] = calloc(ELASTIC_PROPERTIES, sizeof(double));
          wasora_var_value(wasora_mesh.vars.x) = mesh->node[j_global].x[0];
          wasora_var_value(wasora_mesh.vars.y) = mesh->node[j_global].x[1];
          wasora_var_value(wasora_mesh.vars.z) = mesh->node[j_global].x[2];
          
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
          
          if (distribution_alpha.variable == NULL && (distribution_alpha.function != NULL || distribution_alpha.physical_property != NULL)) {
            alpha = fino_distribution_evaluate(&distribution_alpha, element->physical_entity->material, element->node[j]->x);
          }

          element->property_node[j][LAMBDA] = lambda;
          element->property_node[j][MU] = mu;
          element->property_node[j][ALPHA] = alpha;
          
        }
        
        if (fino.rough) {
          // si estamos en rough ya calculamos los valores nodales y ya
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

  // paso intermedio: roughish
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

  // paso 2. barremos nodos de la malla de salida (la misma en smooth, rough en rough)
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
      if (node->dphidx != NULL) {
        gsl_matrix_free(node->dphidx);
      }
      if (node->delta_dphidx != NULL) {
        gsl_matrix_free(node->delta_dphidx);
      }
      if (node->f != NULL) {
        free(node->f);
      }
      node->dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      node->delta_dphidx = gsl_matrix_calloc(fino.degrees, fino.dimensions);
      node->f = calloc(ELASTIC_FUNCTIONS, sizeof(double));
      
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

  PetscFunctionReturn(WASORA_RUNTIME_OK);
}


#undef  __FUNCT__
#define __FUNCT__ "fino_break_set_neumann"
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

    // el momento polar de inercia
    Ix = 0;
    Iy = 0;
    Iz = 0;
    for (i = 0; i < physical_entity->n_elements; i++) {
      tmp_element = &fino.mesh->element[physical_entity->element[i]];
      for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
        mesh_compute_integration_weight_at_gauss(element, v);

        xix = 0;
        xiy = 0;
        xiz = 0;
        for (j = 0; j < element->type->nodes; j++) {
          xix += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[0]);
          xiy += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[1]);
          xiz += element->type->gauss[GAUSS_POINTS_CANONICAL].h[v][j] * gsl_pow_2(tmp_element->node[j]->x[2]);
        }

        Ix += element->w[v] * xix;
        Iy += element->w[v] * xiy;
        Iz += element->w[v] * xiz;
      }
    }
  }  
  
  
  for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
    mesh_compute_integration_weight_at_gauss(element, v);
    mesh_compute_H_at_gauss(element, v, fino.degrees);
    mesh_compute_x_at_gauss(element, v);
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

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_principal_stress"
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

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_vonmises_from_principal"
double fino_compute_vonmises_from_principal(double sigma1, double sigma2, double sigma3) {
  
  return sqrt(0.5*(gsl_pow_2(sigma1-sigma2) + gsl_pow_2(sigma2-sigma3) + gsl_pow_2(sigma3-sigma1)));
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_vonmises_from_stress_tensor"
double fino_compute_vonmises_from_stress_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {
  
  return sqrt(0.5*(gsl_pow_2(sigmax-sigmay) + gsl_pow_2(sigmay-sigmaz) + gsl_pow_2(sigmaz-sigmax) +
                       6.0 * (gsl_pow_2(tauxy) + gsl_pow_2(tauyz) + gsl_pow_2(tauzx))));
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_vonmises_from_strains"
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


#undef  __FUNCT__
#define __FUNCT__ "fino_compute_tresca_from_principal"
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

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_tresca_from_stress_tensor"
double fino_compute_tresca_from_stress_tensor(double sigmax, double sigmay, double sigmaz, double tauxy, double tauyz, double tauzx) {

  double sigma1, sigma2, sigma3;

  wasora_call(fino_compute_principal_stress(sigmax, sigmay, sigmaz, tauxy, tauyz, tauzx, &sigma1, &sigma2, &sigma3));
  return fino_compute_tresca_from_principal(sigma1, sigma2, sigma3);

}


#undef  __FUNCT__
#define __FUNCT__ "fino_compute_strain_energy"
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
