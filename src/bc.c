/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's boundary conditions
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
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifndef _FINO_H
#include "fino.h"
#endif

typedef struct {
  expr_t *expr;
  int dof;
} fino_gsl_function_of_uvw_params_t;


int fino_bc_string2parsed(void) {

  physical_entity_t *physical_entity;
  bc_t *bc;
  char *equal_sign;
  char *name;
  char *expr;
  
  if (fino.mesh == NULL) {
    return WASORA_RUNTIME_OK;
  }
  
  // sweep physical entities and map strings to numerical values
  for (physical_entity = fino.mesh->physical_entities; physical_entity != NULL; physical_entity = physical_entity->hh.next) {
    LL_FOREACH(physical_entity->bcs, bc) {

      // see if there is a "=>" that implies that the BC has a condition
      if ((name = strstr(bc->string, "=>")) != NULL) {
        *name = '\0';
        name += 2;
        wasora_call(wasora_parse_expression(bc->string, &bc->condition));
      } else {
        name = bc->string;
      }
      
      // if there is an equal sign there is an expression, otherwise not
      if ((equal_sign = strchr(name, '=')) != NULL) {
       *equal_sign = '\0';
       expr = equal_sign+1;
      } else {
        expr = NULL;
      } 

      // TODO: function pointers
      if (fino.problem_family == problem_family_mechanical || fino.problem_family == problem_family_modal) {
        fino_bc_process_mechanical(bc, name, expr);
      } else if (fino.problem_family == problem_family_thermal) {
        fino_bc_process_thermal(bc, name, expr, equal_sign);
      }

      // restore the equal sign
      if (equal_sign != NULL) {
        *equal_sign = '=';
      }

    }
  }
  
  PetscFunctionReturn(WASORA_RUNTIME_OK);
}

#undef  __FUNCT__
#define __FUNCT__ "fino_bc_process_mechanical"
int fino_bc_process_mechanical(bc_t *bc, char *name, char *expr) {

  int i;
  bc_t *base_bc = NULL;
  
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

#undef  __FUNCT__
#define __FUNCT__ "fino_bc_process_thermal"
int fino_bc_process_thermal(bc_t *bc, char *name, char *expr, char *equal_sign) {

  int i;
  bc_t *base_bc = NULL;
  bc_t *tmp = NULL;

  if (strcmp(name, "T") == 0) {
    bc->type_math = bc_math_dirichlet;
    bc->type_phys = bc_phys_temperature;
    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if (strcmp(name, "q") == 0 || strcmp(name, "Q") == 0) {
    bc->type_math = bc_math_neumann;
    if (strcmp(name, "Q") == 0) {
      bc->type_phys = bc_phys_heat_total;
    } else {
      bc->type_phys = bc_phys_heat_flux;
    }
    bc->expr = calloc(1, sizeof(expr_t));
    wasora_call(wasora_parse_expression(expr, &bc->expr[0]));

  } else if ((strcmp(name, "h") == 0) ||
             (strcmp(name, "Tref") == 0) ||
             (strcmp(name, "Tinf") == 0)) {

    bc->type_math = bc_math_robin;
    bc->type_phys = bc_phys_convection;

    // conveccion necesita dos expresione
    // las alocamos en la primera de las BCs
    base_bc = bc;
    base_bc->expr = calloc(2, sizeof(expr_t));

    do {
      // volvemos a poner el equal sign, en la primera pasada
      // es para volver a parser, en las siguientes para no romper
      // la ultima se vuelve a arreglar fuera del loop grande
      if (equal_sign != NULL) {
        *equal_sign = '=';
      }
      fino_bc_read_name_expr(bc, &name, &expr, &equal_sign);
      i = -1;
      if (name[0] == 'h') i = 0;
      if (name[0] == 'T') i = 1;
      if (i == -1) {
        wasora_push_error_message("expecting 'h' or 'Tref' instead of '%s'", name);
        return WASORA_PARSER_ERROR;
      }             
      wasora_call(wasora_parse_expression(expr, &base_bc->expr[i]));
      tmp = bc; // esto es para "volver para atras"
    } while ((bc = bc->next) != NULL);

    // ahora bc quedo apuntando a null, tenemos que volver para atras
    // sino el FOREACH de arriba palma
    bc = tmp;

  } else {
    wasora_push_error_message("unknown boundary condition type '%s'", name);
    PetscFunctionReturn(WASORA_PARSER_ERROR);
  }
  
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_bc_read_name_expr"
void fino_bc_read_name_expr(bc_t *bc, char **name, char **expr, char **equal_sign) {

  // si hay signo igual hay expresion, sino no
  *name = bc->string;
  if ((*equal_sign = strchr(bc->string, '=')) != NULL) {
    **equal_sign = '\0';
    *expr = (*equal_sign+1);
  } else {
    *expr = NULL;
  }
  return;
}

int fino_dirichlet_eval(Mat K, Vec b) {
 
  physical_entity_t *physical_entity = NULL;
  physical_entity_t *physical_entity_last = NULL;
  element_list_item_t *associated_element = NULL;
  bc_t *bc = NULL;
  int assembly_needed = 0;

  double n[3] = {0, 0, 0};

  size_t n_bcs = 0;
  size_t current_size = 0;
  
  int j, d;
  int k = 0;
  
  Vec vec_rhs;

  

  if (fino.n_dirichlet_rows != 0) {
    // if we are here then we know more or less the number of BCs we need
    n_bcs = fino.n_dirichlet_rows;
  } else {
    n_bcs = fino.degrees * (fino.last_node - fino.first_node);
  }  
  current_size = n_bcs;

  fino.dirichlet_indexes = calloc(n_bcs, sizeof(PetscInt));
  fino.dirichlet_values = calloc(n_bcs, sizeof(PetscScalar));
  
  for (j = fino.first_node; j < fino.last_node; j++) {

// TODO: arreglar esta logica, los nodos de alto orden terminan con una constante de penalidad diferente que los de primer orden    
//    physical_entity_last = NULL;
    
    LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element) {
      if (associated_element->element != NULL &&
          associated_element->element->type->dim < fino.dimensions &&
          associated_element->element->physical_entity != NULL &&
          associated_element->element->physical_entity != physical_entity_last) {
        LL_FOREACH(associated_element->element->physical_entity->bcs, bc) {
          if (bc->type_math == bc_math_dirichlet) {
//            physical_entity_last = physical_entity; // esto es para no pasar varias veces por lo mismo
            physical_entity = associated_element->element->physical_entity;

            // por que habre querido mirar si n_dirichlet_rows es cero?
//            if (fino.n_dirichlet_rowsk >= (current_size-16)) {
            if (k >= (current_size-16)) {            
              current_size += n_bcs;
              fino.dirichlet_indexes = realloc(fino.dirichlet_indexes, current_size * sizeof(PetscInt));
              fino.dirichlet_values  = realloc(fino.dirichlet_values,  current_size * sizeof(PetscScalar));
            }

            // TODO: see if the outward normal is needed or not
            if ((fino.dimensions - associated_element->element->type->dim) == 1 && associated_element->element->type->dim > 0) {
              wasora_call(mesh_compute_normal(associated_element->element));
              n[0] = wasora_var_value(wasora_mesh.vars.nx);
              n[1] = wasora_var_value(wasora_mesh.vars.ny);
              n[2] = wasora_var_value(wasora_mesh.vars.nz);
            } else {
              n[0] = 0;
              n[1] = 0;
              n[2] = 0;
            }

            // TODO: see if the node coordinates are needed or not
            wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
            wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
            wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
            
            // if there is a condition we evaluate it now
            if (bc->condition.n_tokens == 0 || fabs(wasora_evaluate_expression(&bc->condition)) > 1e-3) {

              // let's see what the user asked for
              // TODO: call functions inside each physics implementation file
              if (bc->type_phys == bc_phys_displacement_fixed) {  

                for (d = 0; d < fino.degrees; d++) {
                  fino.dirichlet_row[k].physical_entity = physical_entity;
                  fino.dirichlet_row[k].dof = d;
                  fino.dirichlet_indexes[k] = fino.mesh->node[j].index_dof[d];
                  fino.dirichlet_values[k] = 0;
                  k++;
                }
                
              } else if (bc->type_phys == bc_phys_displacement ||
                         bc->type_phys == bc_phys_temperature) {

                fino.dirichlet_indexes[k] = fino.mesh->node[j].index_dof[bc->dof];

                if (fino.math_type != math_type_eigen && (strcmp(bc->expr[0].string, "0") != 0)) {
                  fino.dirichlet_values[k] = wasora_evaluate_expression(&bc->expr[0]);
//                  printf("%d bc = %g\n", k, fino.dirichlet_values[k]);
                } else {
                  fino.dirichlet_values[k] = 0;
                }
                
                k++;
                
// temporarily disabled
/*                
              } else if (bc->type_phys == bc_phys_displacement_mimic) {

                if (K != NULL) {
                  // TODO: improve this, check that the master and slave have the same number of nodes
                  gsl_matrix *c;
                  gsl_matrix *P;
                  int l[2];
                  int i, target_index;

                  c = gsl_matrix_calloc(1, 2);
                  P = gsl_matrix_calloc(2, 2);

                  gsl_matrix_set(c, 0, 0, +1);
                  gsl_matrix_set(c, 0, 1, -1);

                  // what to mimic
                  target_index = -1;
                  for (i = 0; i < fino.mesh->n_elements; i++) {
                    if (fino.mesh->element[i].physical_entity != NULL &&
                    fino.mesh->element[i].physical_entity == bc->slave) {
                      target_index = fino.mesh->element[i].node[0]->index_dof[bc->dof];
                      break;
                    }
                  }
                }

                  if (target_index == -1) {
                    wasora_push_error_message("cannot find who to mimic");
                    return WASORA_RUNTIME_ERROR;
                  }

                  l[0] = fino.mesh->node[j].index_dof[bc->dof];
                  l[1] = target_index;

                  if (l[0] != l[1]) {
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    // esto lo necesitamos porque en mimic ponemos cualquier otra estructura diferente a la que ya pusimos antes
                    petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
                    petsc_call(MatSetValues(K, 2, l, 2, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
                    petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
                  }

                  gsl_matrix_free(c);
                  gsl_matrix_free(P);
                  
                }  
*/
              } else if (bc->type_phys == bc_phys_displacement_symmetry) {

                int coordinate_direction = -1;

                // check if the outward normal coincides with a coordinate axis
                for (d = 0; d < 3; d++) {
                  if (fabs(n[d]) > 1-1e-4) {
                    coordinate_direction = d;
                  }
                }

                if (coordinate_direction != -1) {
                  // traditional dirichlet on a single DOF
                  fino.dirichlet_indexes[k] = fino.mesh->node[j].index_dof[coordinate_direction];
                  fino.dirichlet_values[k] = 0;
                  k++;

                } else {
                  
                  // generic multi-freedom
                  if (K != NULL) {
                    
                    int l[3];
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);

                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                      gsl_matrix_set(c, 0, d, n[d]);
                    }

                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));

                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
                    assembly_needed = 1;
                  }

                }  

              } else if (bc->type_phys == bc_phys_displacement_radial) {

                if (K != NULL) {
                  int l[3];
                  double x[3];
                  double eps = 1e-2;

                  for (d = 0; d < 3; d++) {
                    if (bc->expr[d].n_tokens != 0) {
                      x[d] = fino.mesh->node[j].x[d]-wasora_evaluate_expression(&bc->expr[d]);
                    } else {
                      x[d] = fino.mesh->node[j].x[d]-physical_entity->cog[d];
                    }
                  }
  
                  // x-y
                  if (fabs(x[0]) > eps && fabs(x[1]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 0, +x[1]);
                    gsl_matrix_set(c, 0, 1, -x[0]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
  
                  }
  
                  // x-z
                  if (fabs(x[0]) > eps && fabs(x[2]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 0, +x[2]);
                    gsl_matrix_set(c, 0, 2, -x[0]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
                  }
  
                  // y-z
                  if (fabs(x[1]) > eps && fabs(x[2]) > eps) {
  
                    gsl_matrix *c = gsl_matrix_calloc(1, fino.degrees);
                    gsl_matrix *P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                    for (d = 0; d < fino.degrees; d++) {
                      l[d] = fino.mesh->node[j].index_dof[d];
                    }
                    gsl_matrix_set(c, 0, 1, +x[2]);
                    gsl_matrix_set(c, 0, 2, -x[1]);
  
                    wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                    petsc_call(MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES));
  
                    gsl_matrix_free(c);
                    gsl_matrix_free(P);
                  }
                  
                  assembly_needed = 1;
  
                } else if (bc->type_phys == bc_phys_displacement_constrained) {
  
                  gsl_matrix *c;
                  gsl_matrix *P;
                  int l[3];
  
                  gsl_function F;
                  fino_gsl_function_of_uvw_params_t params;
                  double result, abserr;
                  double h = 1e-5;
  
                  params.expr = bc->expr;
  
                  F.function = fino_gsl_function_of_uvw;
                  F.params = &params;
  
                  c = gsl_matrix_calloc(1, fino.degrees);
                  P = gsl_matrix_calloc(fino.degrees, fino.degrees);
  
                  wasora_var_value(fino.vars.U[0]) = 0;
                  wasora_var_value(fino.vars.U[1]) = 0;
                  wasora_var_value(fino.vars.U[2]) = 0;
  
                  for (d = 0; d < fino.degrees; d++) {
                    l[d] = fino.mesh->node[j].index_dof[d];
                    params.dof = d;
                    gsl_deriv_central(&F, 0, h, &result, &abserr);
                    gsl_matrix_set(c, 0, d, -result);
                  }
  
                  wasora_call(gsl_blas_dgemm(CblasTrans, CblasNoTrans, wasora_var(fino.vars.penalty_weight), c, c, 0, P));
                  MatSetValues(K, fino.degrees, l, fino.degrees, l, gsl_matrix_ptr(P, 0, 0), ADD_VALUES);
  
                  // TODO: non-uniform RHS
                  gsl_matrix_free(c);
                  gsl_matrix_free(P);
                  assembly_needed = 1;
                }
              }  
            }  
          }
        }
      }
    }
  }

  // now we know how many rows we need to change
  if (fino.n_dirichlet_rows != k) {
    fino.n_dirichlet_rows = k;
    
    // if k == 0 this like freeing
    fino.dirichlet_indexes = realloc(fino.dirichlet_indexes, fino.n_dirichlet_rows * sizeof(PetscInt));
    fino.dirichlet_values = realloc(fino.dirichlet_values, fino.n_dirichlet_rows * sizeof(PetscScalar));
  }

  if (K != NULL && assembly_needed) {
    // this is needed with multi-freedom BCs
    wasora_call(fino_assembly());
  }  
  
  return WASORA_RUNTIME_OK;
}


// K - stiffness matrix: needs a one in the diagonal and the value in b and keep symmetry
// b - RHS: needs to be updated when modifying K
int fino_dirichlet_set_K(Mat K, Vec b) {
  
  int k;
  PetscScalar diag;
  Vec rhs;
  
  // sometimes there are free nodes with no associated volumes
  // this can trigger zeros on the diagonal and MatZeroRowsColumns complains
  // TODO: change to MatGetDiagonal
  for (k = fino.first_row; k < fino.last_row; k++) {
    petsc_call(MatGetValues(K, 1, &k, 1, &k, &diag));
    if (diag == 0) {
      petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
      petsc_call(MatSetValue(K, k, k, 1.0, INSERT_VALUES));
      petsc_call(MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
      wasora_call(fino_assembly());
    }  
  }
    
  // this vector holds the dirichlet values and is used to re-write
  // the actual rhs vector b in order to keep the symmetry of K
  petsc_call(MatCreateVecs(fino.K, NULL, &rhs));
  petsc_call(VecSetValues(rhs, fino.n_dirichlet_rows, fino.dirichlet_indexes, fino.dirichlet_values, INSERT_VALUES));
  // TODO: scale up the diagonal!
  // see alpha in https://scicomp.stackexchange.com/questions/3298/appropriate-space-for-weak-solutions-to-an-elliptical-pde-with-mixed-inhomogeneo/3300#3300
  petsc_call(MatZeroRowsColumns(K, fino.n_dirichlet_rows, fino.dirichlet_indexes, 1.0, rhs, b));
  petsc_call(VecDestroy(&rhs));
  
  return WASORA_RUNTIME_OK;
}
  
// M - mass matrix: needs a zero in the diagonal and the same symmetry scheme that K
int fino_dirichlet_set_M(Mat M) {

  // the mass matrix is like the stiffness one but with zero instead of one
  petsc_call(MatZeroRowsColumns(M, fino.n_dirichlet_rows, fino.dirichlet_indexes, 0.0, NULL, NULL));

  return WASORA_RUNTIME_OK;
}

// J - Jacobian matrix: needs a one in the diagonal but does not need to keep the symmetry
int fino_dirichlet_set_J(Mat J) {

  // the jacobian is exactly one for the dirichlet values and zero otherwise without keeping symmetry
  petsc_call(MatZeroRows(J, fino.n_dirichlet_rows, fino.dirichlet_indexes, 1.0, NULL, NULL));

  return WASORA_RUNTIME_OK;
}

// JM - Jacobian of the time-derivative matrix (the part of sigma*d(R)/d(phi_dot)): needs a zero in the diagonal but does not need to keep the symmetry
int fino_dirichlet_set_dRdphi_dot(Mat M) {

  // the jacobian is exactly zero for the dirichlet values and zero otherwise without keeping symmetry
  petsc_call(MatZeroRows(M, fino.n_dirichlet_rows, fino.dirichlet_indexes, 0.0, NULL, NULL));

  return WASORA_RUNTIME_OK;
}


// phi - solution: the BC values are set directly in order to be used as a initial condition or guess
int fino_dirichlet_set_phi(Vec phi) {

  // this should be used only to set initial conditions and guesses
  petsc_call(VecSetValues(phi, fino.n_dirichlet_rows, fino.dirichlet_indexes, fino.dirichlet_values, INSERT_VALUES));

  return WASORA_RUNTIME_OK;
}

// r - residual: the BC values are set to the difference between the value and the solution
int fino_dirichlet_set_r(Vec r, Vec phi) {

  int k;
  
  PetscScalar *diff = calloc(fino.n_dirichlet_rows, sizeof(PetscScalar));
  petsc_call(VecGetValues(phi, fino.n_dirichlet_rows, fino.dirichlet_indexes, diff));
  for (k = 0; k < fino.n_dirichlet_rows; k++) {
    diff[k] -= fino.dirichlet_values[k];
  }
  
  petsc_call(VecSetValues(r, fino.n_dirichlet_rows, fino.dirichlet_indexes, diff, INSERT_VALUES));
  free(diff);

  return WASORA_RUNTIME_OK;
}



// esto sirve para calcular derivadas con GSL
double fino_gsl_function_of_uvw(double x, void *params) {

  double y;

  fino_gsl_function_of_uvw_params_t *dummy = (fino_gsl_function_of_uvw_params_t *)params;

  wasora_var_value(fino.vars.U[0]) = 0;
  wasora_var_value(fino.vars.U[1]) = 0;
  wasora_var_value(fino.vars.U[2]) = 0;
  wasora_var_value(fino.vars.U[dummy->dof]) = x;
  
  y = wasora_evaluate_expression(dummy->expr);

  if (gsl_isnan(y) || gsl_isinf(y)) {
    wasora_nan_error();
  }

  return y;

}
