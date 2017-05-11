/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's bulk routines
 *
 *  Copyright (C) 2015-2016 jeremy theler
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#define NAME_SIZE 32

#undef  __FUNCT__
#define __FUNCT__ "fino_allocate_elemental_objects"
int fino_allocate_elemental_objects(element_t *element) {

  fino.n_local_nodes = element->type->nodes;
  fino.elemental_size = element->type->nodes * fino.degrees;
  
  // TODO: esta las tendria que alocar mesh
  gsl_matrix_free(fino.mesh->fem.H);
  fino.mesh->fem.H = gsl_matrix_calloc(fino.mesh->degrees_of_freedom, fino.elemental_size);
  
  gsl_matrix_free(fino.mesh->fem.B);
  fino.mesh->fem.B = gsl_matrix_calloc(fino.mesh->degrees_of_freedom * fino.mesh->bulk_dimensions, fino.elemental_size);

  gsl_matrix_free(fino.Ai);
  fino.Ai = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  gsl_matrix_free(fino.Bi);
  fino.Bi = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  gsl_vector_free(fino.bi);
  fino.bi = gsl_vector_calloc(fino.elemental_size);
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_bulk"
int fino_build_bulk(void) {

  int i;
  int step = (fino.mesh->n_elements > 99)?fino.mesh->n_elements/100:1;
  
  for (i = 0; i < fino.mesh->n_elements; i++) {

// ------ progress bar ------------------------------------------    
    if ((i % step) == 0) {
      if (fino.shmem_memory != NULL) {
        getrusage(RUSAGE_SELF, &fino.resource_usage);
        *fino.shmem_memory = (double)(1024.0*fino.resource_usage.ru_maxrss);
      }
      if (fino.shmem_progress_build != NULL) {
        *fino.shmem_progress_build = (double)i/(double)fino.mesh->n_elements;
      }
    }
// --------------------------------------------------------------    

    if (fino.mesh->element[i].type->dim == fino.dimensions) {
      // solo los elementos que tengan la dimension del problema
      // son los que usamos para las matrices elementales
      wasora_call(fino_build_element_volumetric(&fino.mesh->element[i]));
      
    } else if (fino.math_type != math_eigen &&
               fino.mesh->element[i].type->dim < fino.dimensions &&
               fino.mesh->element[i].physical_entity != NULL &&
               (fino.mesh->element[i].physical_entity->bc_type_math == bc_math_neumann ||
                fino.mesh->element[i].physical_entity->bc_type_math == bc_math_robin)) {

      // los otros tienen (o pueden tener) condiciones de contorno de neumann
      // las de dirichlet set ponen en set_essential despues de ensamblar
      wasora_call(fino_build_element_bc(&fino.mesh->element[i]));
    }
  }

  // ver si esto chupa memoria
  wasora_call(fino_assembly());
  
  if (fino.shmem_progress_build != NULL) {
    *fino.shmem_progress_build = 1.0;
  }

  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_element_volumetric"
int fino_build_element_volumetric(element_t *element) {
  int v;           // indice del punto de gauss

  if (element->physical_entity == NULL) {
    // esto (deberia) pasar solo en malla estructuradas
    wasora_push_error_message("volumetric element %d does not have an associated physical entity", element->id);
    return WASORA_RUNTIME_ERROR;
    
  } else {

    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }  
    
    // inicializamos en cero los objetos elementales
    gsl_matrix_set_zero(fino.Ai);
    gsl_matrix_set_zero(fino.Bi);
    gsl_vector_set_zero(fino.bi);

    // TODO: hacer el loop de gauss adentro de cada funcion asi podemos
    // hacer evaluaciones nodo por nodo o lo que sea
    // para cada punto de gauss
    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      // armamos las matrices
      if (fino.problem == problem_shake || fino.problem == problem_break) {
        wasora_call(fino_break_build_element(element, v));
      } else if (fino.problem == problem_bake) {
        wasora_call(fino_build_bake(element, v));
      }
    }

    MatSetValues(fino.A, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Ai, 0, 0), ADD_VALUES);
    if (fino.math_type == math_linear) {
      VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    } else if (fino.math_type == math_eigen)  {
      MatSetValues(fino.B, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Bi, 0, 0), ADD_VALUES);
    }

/*    
    if (fino.dump != NULL) {
      fprintf(fino.dump, "\nelement %d\n", element->id);
      fprintf(fino.dump, "%s\n", fino.matrix_name);
      fino_print_gsl_matrix(fino.Ai, fino.dump);
      fprintf(fino.dump, "%s\n", fino.vector_name);
      fino_print_gsl_vector(fino.bi, fino.dump);
      fprintf(fino.dump, "\n");
    }
 */ 
    
    
  }

  return WASORA_RUNTIME_OK;
}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_element_volumetric"
int fino_build_element_bc(element_t *element) {
  
  double n[3];
  
  if (fino.problem == problem_break) {
    wasora_call(mesh_compute_outward_normal(element, n));
    wasora_var_value(fino.vars.nx) = n[0];
    wasora_var_value(fino.vars.ny) = n[1];
    wasora_var_value(fino.vars.nz) = n[2];
    
    if (element->physical_entity->bc_type_phys == bc_phys_stress) {
      wasora_call(fino_break_add_stress(element));
    } else if (element->physical_entity->bc_type_phys == bc_phys_force) {
      wasora_call(fino_break_add_force(element));
    } else if (element->physical_entity->bc_type_phys == bc_phys_pressure) {
      wasora_call(fino_break_add_pressure(element));
    }
  }
  
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_print_gsl_vector"
int fino_print_gsl_vector(gsl_vector *b, FILE *file) {

  double xi;
  int i;

  for (i = 0; i < b->size; i++) {
    xi = gsl_vector_get(b, i);
    if (xi != 0) {
      fprintf(file, "% .1e ", xi);
    } else {
      fprintf(file, "    0    ");
    }
    fprintf(file, "\n");
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_print_gsl_matrix"
int fino_print_gsl_matrix(gsl_matrix *A, FILE *file) {

  double xi;
  int i, j;

  for (i = 0; i < A->size1; i++) {
    for (j = 0; j < A->size2; j++) {
      xi = gsl_matrix_get(A, i, j);
      if (xi != 0) {
        fprintf(file, "% .1e ", xi);
      } else {
        fprintf(file, "    0    ");
      }
    }
    fprintf(file, "\n");
  }
  
  return WASORA_RUNTIME_OK;

}
