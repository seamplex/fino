/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's boundary conditions
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
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#undef  __FUNCT__
#define __FUNCT__ "fino_read_bcs"
int fino_read_bcs(void) {

  physical_entity_t *physical_entity;
  
  // barremos los physical entities y mapeamos cadenas a valores enteros
  // TODO: elegir si hay que ser lightweight o quejarnos si no hay un mapeo 1 a 1
  LL_FOREACH(wasora_mesh.physical_entities, physical_entity) {
    if (physical_entity->bc_strings != NULL) {
      
      bc_string_based_t *bc;
      char *equal_sign = NULL;
      char *name = NULL;
      char *expr = NULL;
      
      // TODO: poner las BCs de diferentes DOFs en un solo vector asi hacemos una sola cuenta y ya
      // TODO: ver https://scicomp.stackexchange.com/questions/3298/appropriate-space-for-weak-solutions-to-an-elliptical-pde-with-mixed-inhomogeneo/3300#3300      
      LL_FOREACH(physical_entity->bc_strings, bc) {
        
        // si creamos una bc dummy entonces tiene string NULL
        // pero no importa porque ya la procesamos
        if (bc->string == NULL) {
          continue;
        }
        
        // si hay signo igual hay expresion, sino no
        name = bc->string;
        if ((equal_sign = strchr(bc->string, '=')) != NULL) {
          *equal_sign = '\0';
          expr = (equal_sign+1);
        } else {
          expr = NULL;
        }

        if (fino.problem == problem_break) {
          if (strcmp(name, "u") == 0) {
            bc->dof = 0;
            bc->bc_type_int = BC_DIRICHLET;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
          } else if (strcmp(name, "v") == 0) {
            bc->dof = 1;
            bc->bc_type_int = BC_DIRICHLET;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
          } else if (strcmp(name, "w") == 0) {
            bc->dof = 2;
            bc->bc_type_int = BC_DIRICHLET;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
          // fixed es como u=0 v=0 w=0
          } else if (strcmp(name, "fixed") == 0) {
            physical_entity->bc_type_int = BC_DIRICHLET_NULL;

          // se le dan o las componentes del vector t [tx, ty, tz]
          // le ponemos a disposicion las componentes del vector n [nx, ny, nz]
          } else if (strcmp(name, "tx") == 0) {
            bc->dof = 0;
            bc->bc_type_int = BC_NEUMANN;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
          } else if (strcmp(name, "ty") == 0) {
            bc->dof = 1;
            bc->bc_type_int = BC_NEUMANN;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
          } else if (strcmp(name, "tz") == 0) {
            bc->dof = 2;
            bc->bc_type_int = BC_NEUMANN;
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
          } else if (strcmp(name, "p") == 0) {
            // TODO: capaz que conviene hacer esto al evaluar
            bc_string_based_t *bc_v, *bc_w;
            char *buff = malloc(strlen(expr) + 8);
            
            bc->dof = 0;
            bc->bc_type_int = BC_NEUMANN;
            snprintf(buff, strlen(expr) + 7, "-nx*(%s)", expr);
            wasora_call(wasora_parse_expression(buff, &bc->expr));
            
            // agregamos dos mas de prepo para v y w
            bc_v = calloc(1, sizeof(bc_string_based_t));
            bc_v->dof = 1;
            bc_v->bc_type_int = BC_NEUMANN;
            snprintf(buff, strlen(expr) + 7, "-ny*(%s)", expr);
            wasora_call(wasora_parse_expression(buff, &bc_v->expr));
            LL_APPEND(physical_entity->bc_strings, bc_v);
            
            bc_w = calloc(1, sizeof(bc_string_based_t));
            bc_w->dof = 2;
            bc_w->bc_type_int = BC_NEUMANN;
            snprintf(buff, strlen(expr) + 7, "-nz*(%s)", expr);
            wasora_call(wasora_parse_expression(buff, &bc_w->expr));
            LL_APPEND(physical_entity->bc_strings, bc_w);
            
            free(buff);

          } else if (strcmp(name, "implicit") == 0) {
            char *dummy;
            
            physical_entity->bc_type_int = BC_DIRICHLET_ALG;
            bc->bc_type_int = BC_DIRICHLET_ALG;
            
            // el cuento es asi: aca quisieramos que el usuario escriba algo en funcion
            // de x,y,z pero tambien de u,v y w. Pero u,v,w ya son funciones, asi que no
            // se pueden usar como variables
            // mi solucion: definir variables U,V,W y reemplazar u,v,w por U,V,W en
            // esta expresion

            // TODO: ver que haya un separador antes y despues            
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
            
            wasora_call(wasora_parse_expression(expr, &bc->expr));
            
          } else {
            wasora_push_error_message("unknown boundary condition type '%s'", name);
            return WASORA_PARSER_ERROR;
          }
        }
        
        // restauramos el signo igual porque en parametrico en una epoca pasaba de nuevo por aca vamos a volver a pasar por aca
        // ahora ya no pero por si acaso
        if (equal_sign != NULL) {
          *equal_sign = '=';
        }

      }
    }
  }
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_count_bc_expressions"
int fino_count_bc_expressions(expr_t *bc_args) {
  expr_t *bc_arg;
  int n = 0;
  
  LL_FOREACH(bc_args, bc_arg) {
    n++;
  }
  
  return n;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_evaluate_bc_expressions"
int fino_evaluate_bc_expressions(physical_entity_t *physical_entity, node_t *node, int degrees, double multiplier, double *result) {

  int g;
  expr_t *bc_arg = physical_entity->bc_args;

  // las cc dependen de la posicion
  wasora_var(wasora_mesh.vars.x) = node->x[0];
  wasora_var(wasora_mesh.vars.y) = node->x[1];
  wasora_var(wasora_mesh.vars.z) = node->x[2];

  for (g = 0; g < degrees; g++) {
    if (bc_arg == NULL) {
      wasora_push_error_message("boundary condition '%s' needs %d expressions", physical_entity->name, degrees);
      return WASORA_RUNTIME_ERROR;
    }

    result[g] = multiplier * wasora_evaluate_expression(bc_arg);
    bc_arg = bc_arg->next;
  }
  
  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_set_essential_bc"
int fino_set_essential_bc(void) {

  PetscScalar *rhs;
  PetscInt *indexes;
  physical_entity_t *physical_entity;
  element_list_item_t *associated_element;
  
  const PetscInt *cols;
  const PetscScalar *vals;
  size_t current_size = fino.problem_size;
  size_t current_threshold = BC_FACTOR*fino.problem_size - 2*fino.degrees;

  PetscInt ncols;
  Vec vec_rhs;

  int found = 0;
  int non_homogeneous_bc = 0;
  int k = 0;
  int i, j, d;
  
  bc_string_based_t *bc;
  

  double y1, y2;
  double h = 1e-3;
  
  if (wasora_var(fino.vars.dirichlet_diagonal) == 0) {
    wasora_var(fino.vars.dirichlet_diagonal) = 1;
  }
  
  rhs = malloc(BC_FACTOR*fino.problem_size * sizeof(PetscScalar));
  indexes = malloc(BC_FACTOR*fino.problem_size * sizeof(PetscInt));
  fino.dirichlet_row = calloc(BC_FACTOR*fino.problem_size, sizeof(dirichlet_row_t));
   
  for (j = 0; j < fino.mesh->n_nodes; j++) {
    found = 0;
    LL_FOREACH(fino.mesh->node[j].associated_elements, associated_element) {
      if (!found && 
          associated_element->element->type->dim != fino.dimensions &&
          (physical_entity = associated_element->element->physical_entity) != NULL) {
          
        if (k >= current_threshold) {
          current_size += BC_FACTOR*fino.problem_size;
          current_threshold = current_size - 2*fino.degrees;
          indexes = realloc(indexes, current_size * sizeof(PetscInt));
          rhs = realloc(rhs, current_size * sizeof(PetscScalar));
          fino.dirichlet_row = realloc(fino.dirichlet_row, current_size * sizeof(dirichlet_row_t));
        }

        if (physical_entity->bc_type_int == BC_DIRICHLET_NULL) {
          for (d = 0; d < fino.degrees; d++) {
            fino.dirichlet_row[k].physical_entity = physical_entity;
            fino.dirichlet_row[k].dof = d;
            indexes[k] = fino.mesh->node[j].index[d];
            rhs[k] = 0;
            k++;
          }
          found = 1;
          

        } else if (physical_entity->bc_strings != NULL) {

          LL_FOREACH(physical_entity->bc_strings, bc) {

            switch(bc->bc_type_int) {
              case BC_DIRICHLET_ALG:

                fino.dirichlet_row[k].physical_entity = physical_entity;

                wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
                wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
                wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];
                wasora_var_value(fino.vars.U[0]) = 0;
                wasora_var_value(fino.vars.U[1]) = 0;
                wasora_var_value(fino.vars.U[2]) = 0;

                if ((rhs[k] = wasora_evaluate_expression(&bc->expr)) != 0) {
                  non_homogeneous_bc = 1;
                }

                fino.dirichlet_row[k].alg_col = calloc(fino.degrees, sizeof(PetscInt));
                fino.dirichlet_row[k].alg_val = calloc(fino.degrees, sizeof(PetscScalar));

                for (d = 0; d < fino.degrees; d++) {
                  fino.dirichlet_row[k].alg_col[d] = fino.mesh->node[j].index[d];
                  wasora_var_value(fino.vars.U[d]) = +h;
                  y1 = wasora_evaluate_expression(&bc->expr);
                  wasora_var_value(fino.vars.U[d]) = -h;
                  y2 = wasora_evaluate_expression(&bc->expr);
                  wasora_var_value(fino.vars.U[d]) = 0;
                  fino.dirichlet_row[k].alg_val[d] = -(y1-y2)/(2.0*h);
                }

                // el indice es el que tiene el coeficiente mayor (vaya uno a saber por que)
                for (d = 0; d < fino.dimensions; d++) {
                  if (fabs(fino.dirichlet_row[k].alg_val[d]) >= fabs(fino.dirichlet_row[k].alg_val[0]) &&
                      fabs(fino.dirichlet_row[k].alg_val[d]) >= fabs(fino.dirichlet_row[k].alg_val[1]) && 
                      fabs(fino.dirichlet_row[k].alg_val[d]) >= fabs(fino.dirichlet_row[k].alg_val[2])) {
                    indexes[k] = fino.mesh->node[j].index[d];
                    fino.dirichlet_row[k].dof = d;
                  }
                }
                
                k++;
                found = 1;

              break;

              case BC_DIRICHLET:

                if (bc->bc_type_int == BC_DIRICHLET) {
                  fino.dirichlet_row[k].physical_entity = associated_element->element->physical_entity;
                  fino.dirichlet_row[k].dof = bc->dof;

                  indexes[k] = fino.mesh->node[j].index[bc->dof];

                  if (fino.math_type == math_linear && (strcmp(bc->expr.string, "0") != 0)) {
                    wasora_var_value(wasora_mesh.vars.x) = fino.mesh->node[j].x[0];
                    wasora_var_value(wasora_mesh.vars.y) = fino.mesh->node[j].x[1];
                    wasora_var_value(wasora_mesh.vars.z) = fino.mesh->node[j].x[2];

                    if ((rhs[k] = wasora_var(fino.vars.dirichlet_diagonal) * wasora_evaluate_expression(&bc->expr)) != 0) {
                      non_homogeneous_bc = 1;
                    }
                  } else {
                    rhs[k] = 0;
                  }

                  k++;
                }
                found = 1;
                
              break;
            }
          }
        }
      }
    }
  }

  // antes de romper las filas de dirichlet, nos las acordamos para calcular las reacciones  
  // ojo! aca estamos contando varias veces el mismo nodo, porque un nodo pertenece a varios elementos
  // TODO: hacer lo que dijo barry
  for (i = 0; i < fino.n_dirichlet_rows; i++) {
    petsc_call(MatGetRow(fino.A, indexes[i], &ncols, &cols, &vals));
    fino.dirichlet_row[i].ncols = ncols;
    fino.dirichlet_row[i].cols = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscInt *));
    fino.dirichlet_row[i].vals = calloc(fino.dirichlet_row[i].ncols, sizeof(PetscScalar *));
    // por si acaso igualamos en lugar de hacer memcpy
    for (j = 0; j < ncols; j++) {
      fino.dirichlet_row[i].cols[j] = cols[j];
      fino.dirichlet_row[i].vals[j] = vals[j];
    }
    petsc_call(MatRestoreRow(fino.A, indexes[i], &ncols, &cols, &vals));
  }
  
  if (non_homogeneous_bc) {
    if (fino.math_type == math_eigen) {
      wasora_push_error_message("boundary conditions ought to be homogeneous in eigenvalue problems");
      return WASORA_RUNTIME_ERROR;
    }
    petsc_call(MatCreateVecs(fino.A, &vec_rhs, NULL));
    petsc_call(VecSetValues(vec_rhs, k, indexes, rhs, INSERT_VALUES));
    petsc_call(MatZeroRowsColumns(fino.A, k, indexes, wasora_var(fino.vars.dirichlet_diagonal), vec_rhs, fino.b));
    petsc_call(VecDestroy(&vec_rhs));
  } else {
    petsc_call(MatZeroRowsColumns(fino.A, k, indexes, wasora_var(fino.vars.dirichlet_diagonal), PETSC_NULL, PETSC_NULL));
  }
  
  if (fino.math_type == math_eigen) {
    petsc_call(MatZeroRowsColumns(fino.B, k, indexes, 0.0, PETSC_NULL, PETSC_NULL));
  }
  
  // TODO: hacer un array ya listo para hacer un unico MatSetValuesS
  for (i = 0; i < fino.n_dirichlet_rows; i++) {
    if (fino.dirichlet_row[i].alg_val != NULL) {
      for (d = 0; d < fino.degrees; d++) {
        petsc_call(MatSetValues(fino.A, 1, &indexes[i], fino.degrees, fino.dirichlet_row[i].alg_col, fino.dirichlet_row[i].alg_val, INSERT_VALUES));
      }
    }
  }
    
  free(indexes);
  free(rhs);

  wasora_call(fino_assembly());

  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_surface_objects"
// las condiciones de contorno son
// dphi/dn = a * phi + b
// para neumann, a = 0
int fino_build_surface_objects(element_t *element, expr_t *bc_a, expr_t *bc_b) {
  int v, g;
  double w_gauss;
  expr_t *current_bc_a;
  expr_t *current_bc_b;
  gsl_vector *Nb;
  gsl_matrix *Na;
  gsl_matrix *NaH;
  double *cc;

  if (element->type->dim == 0) {

    // si es un punto es un punto y chau
    cc = malloc(fino.degrees * sizeof(double));
    
    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[0], fino.degrees, 1, cc));
    VecSetValues(fino.b, fino.degrees, element->node[0]->index, cc, INSERT_VALUES);
    
    free(cc);
    
  } else if (element->type->dim == 1 && fino.dimensions == 3) {
    // las expresiones estan dadas en "lo que sea" por unidad de longitud, asi que le damos la mitad a cada nodo
    double halflength = 0.5 * gsl_hypot3(element->node[1]->x[0] - element->node[0]->x[0],
                                         element->node[1]->x[1] - element->node[0]->x[1],
                                         element->node[1]->x[2] - element->node[0]->x[2]);
    // si es una linea y el problema es 3 hacemos como si fuesen dos puntos
    cc = malloc(fino.degrees * sizeof(double));

    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[0], fino.degrees, halflength, cc));
    VecSetValues(fino.b, fino.degrees, element->node[0]->index, cc, INSERT_VALUES);

    wasora_call(fino_evaluate_bc_expressions(element->physical_entity, element->node[1], fino.degrees, halflength, cc));
    VecSetValues(fino.b, fino.degrees, element->node[1]->index, cc, INSERT_VALUES);
    
    free(cc);
    
  } else {
    // sino hacemos la cuenta general
  
    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }
  
    Na = gsl_matrix_calloc(fino.degrees, fino.degrees);
    NaH = gsl_matrix_calloc(fino.degrees, fino.n_local_nodes);
    Nb = gsl_vector_calloc(fino.degrees);
  
    gsl_matrix_set_zero(fino.Ai);
    gsl_vector_set_zero(fino.bi);

    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);    
      current_bc_a = bc_a;
      current_bc_b = bc_b;
      for (g = 0; g < fino.degrees; g++) {
        if (current_bc_a != NULL) {
          gsl_matrix_set(Na, g, g, -wasora_evaluate_expression(current_bc_a));
        }
        gsl_vector_set(Nb, g, wasora_evaluate_expression(current_bc_b));

        if (current_bc_a != NULL) {
          current_bc_a = current_bc_a->next;
        }
        if (current_bc_b != NULL) {
          current_bc_b = current_bc_b->next;
        }
      }

      if (bc_a != NULL) {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Na, fino.mesh->fem.H, 0, NaH);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, w_gauss, fino.mesh->fem.H, NaH, 1, fino.Ai);
      }
    
      gsl_blas_dgemv(CblasTrans, -w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
    }
    
    MatSetValues(fino.A, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Ai, 0, 0), ADD_VALUES);
    VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    
    gsl_vector_free(Nb);
    gsl_matrix_free(Na);
    gsl_matrix_free(NaH);
  }
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "fino_add_single_surface_term_to_rhs"
int fino_add_single_surface_term_to_rhs(element_t *element, bc_string_based_t *bc) {
  int v;
  double w_gauss;
  gsl_vector *Nb;
  double cc;
  double a[3], b[3], n[3], surface_center[3], volumetric_neighbor_center[3];
  element_t *volumetric_neighbor;
    
  if (element->type->dim == 0) {

    // si es un punto es un punto y chau
    // las cc dependen de la posicion
    wasora_var(wasora_mesh.vars.x) = element->node[0]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[0]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[0]->x[2];    
    cc = wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[0]->index[bc->dof], &cc, ADD_VALUES);
    
  } else if (element->type->dim == 1 && fino.dimensions == 3) {
    // las expresiones estan dadas en "lo que sea" por unidad de longitud, asi que le damos la mitad a cada nodo
    double halflength = 0.5 * gsl_hypot3(element->node[1]->x[0] - element->node[0]->x[0],
                                         element->node[1]->x[1] - element->node[0]->x[1],
                                         element->node[1]->x[2] - element->node[0]->x[2]);
    
    // si es una linea y el problema es 3 hacemos como si fuesen dos puntos
    wasora_var(wasora_mesh.vars.x) = element->node[0]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[0]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[0]->x[2];    
    cc = halflength * wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[0]->index[bc->dof], &cc, ADD_VALUES);

    wasora_var(wasora_mesh.vars.x) = element->node[1]->x[0];
    wasora_var(wasora_mesh.vars.y) = element->node[1]->x[1];
    wasora_var(wasora_mesh.vars.z) = element->node[1]->x[2];    
    cc = halflength * wasora_evaluate_expression(&bc->expr);
    VecSetValues(fino.b, 1, &element->node[1]->index[bc->dof], &cc, ADD_VALUES);
    
  } else {
    // sino hacemos la cuenta general

    // calculamos el vector normal para las variables nx ny y nx
    mesh_subtract(element->node[0]->x, element->node[1]->x, a);
    mesh_subtract(element->node[0]->x, element->node[2]->x, b);
    mesh_normalized_cross(a, b, n);
    
    // ahora tenemos que ver si la normal que elegimos es efectivamente la outward
    // para eso primero calculamos el centro del elemento de superficie
    wasora_call(mesh_compute_element_barycenter(element, surface_center));

    // y despues el centro del elemento de volumen
    if ((volumetric_neighbor = mesh_find_element_volumetric_neighbor(element)) == NULL) {
      wasora_push_error_message("cannot find any volumetric neighbor for surface element %d", element->id);
      PetscFunctionReturn(WASORA_RUNTIME_ERROR);
    }
    
    volumetric_neighbor = mesh_find_element_volumetric_neighbor(element);
    wasora_call(mesh_compute_element_barycenter(volumetric_neighbor, volumetric_neighbor_center));

    // calculamos el producto entre la normal propuesta y la resta de estos dos vectores
    // si elegimos la otra direccion, la damos tavuel
    if (mesh_subtract_dot(volumetric_neighbor_center, surface_center, n) > 0) {
      n[0] = -n[0];
      n[1] = -n[1];
      n[2] = -n[2];
    }    
        
    wasora_var_value(fino.vars.nx) = n[0];
    wasora_var_value(fino.vars.ny) = n[1];
    wasora_var_value(fino.vars.nz) = n[2];
    
    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }
  
    Nb = gsl_vector_calloc(fino.degrees);
    gsl_vector_set_zero(fino.bi);

    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      w_gauss = mesh_compute_fem_objects_at_gauss(fino.mesh, element, v);    
      gsl_vector_set(Nb, bc->dof, wasora_evaluate_expression(&bc->expr));
      gsl_blas_dgemv(CblasTrans, w_gauss, fino.mesh->fem.H, Nb, 1.0, fino.bi); 
    }
    
    VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    
    gsl_vector_free(Nb);
  }
  
  return WASORA_RUNTIME_OK;
}