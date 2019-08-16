/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's bulk routines
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
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "fino.h"

#define NAME_SIZE 32

#undef  __FUNCT__
#define __FUNCT__ "fino_allocate_elemental_objects"
int fino_allocate_elemental_objects(element_t *element) {

  fino_free_elemental_objects();
      
  fino.n_local_nodes = element->type->nodes;
  fino.elemental_size = element->type->nodes * fino.degrees;
  
  // TODO: estas las tendria que alocar mesh
  fino.mesh->fem.H = gsl_matrix_calloc(fino.mesh->degrees_of_freedom, fino.elemental_size);
  fino.mesh->fem.B = gsl_matrix_calloc(fino.mesh->degrees_of_freedom * fino.mesh->bulk_dimensions, fino.elemental_size);
  
  fino.Ki = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  fino.Mi = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  fino.bi = gsl_vector_calloc(fino.elemental_size);
  
  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "fino_free_elemental_objects"
int fino_free_elemental_objects(void) {

  if (fino.n_local_nodes != 0 && fino.elemental_size != 0) {
    gsl_matrix_free(fino.mesh->fem.H);
    gsl_matrix_free(fino.mesh->fem.B);
    
    gsl_matrix_free(fino.Ki);
    gsl_matrix_free(fino.Mi);
    gsl_vector_free(fino.bi);
  }
  
  fino.n_local_nodes = 0;
  fino.elemental_size = 0;
  
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_bulk"
int fino_build_bulk(void) {

  bc_t *bc;
  int i;
  int step = (fino.mesh->n_elements > 99)?ceil((double)fino.mesh->n_elements/100.0):1;
  int ascii_progress_chars = 0;
  
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
      if (fino.progress_ascii) {
        printf(CHAR_PROGRESS_BUILD);  
        fflush(stdout);
        ascii_progress_chars++;
      }
    }
// --------------------------------------------------------------    

    if (fino.mesh->element[i].type->dim == fino.dimensions) {
      
      // solo los elementos que tengan la dimension del problema
      // son los que usamos para las matrices elementales
      wasora_call(fino_build_element_volumetric(&fino.mesh->element[i]));
      
    } else if (fino.math_type != math_type_eigen &&
               fino.mesh->element[i].type->dim < fino.dimensions &&
               fino.mesh->element[i].physical_entity != NULL) {
      LL_FOREACH(fino.mesh->element[i].physical_entity->bcs, bc) {
        
        if (bc->type_math == bc_math_neumann || bc->type_math == bc_math_robin) {
          // las de dirichlet set ponen en set_essential despues de ensamblar
          wasora_call(fino_build_element_bc(&fino.mesh->element[i], bc));
        }
        
      }
    }
  }

  // ver si esto chupa memoria
  wasora_call(fino_assembly());
  
  if (fino.shmem_progress_build != NULL) {
    *fino.shmem_progress_build = 1.0;
  }
  if (fino.progress_ascii) {
    while (ascii_progress_chars++ < 100) {
      printf(CHAR_PROGRESS_BUILD);
    }
    printf("\n");  
    fflush(stdout);
  }

  // aca tambien perdemos a C y a et porque son static
  wasora_call(fino_free_elemental_objects());
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "fino_build_element_volumetric"
int fino_build_element_volumetric(element_t *element) {
  int v;           // indice del punto de gauss

  if (element->physical_entity == NULL) {
    // esto (deberia) pasar solo en malla estructuradas
    wasora_push_error_message("volumetric element %d does not have an associated physical entity", element->tag);
    return WASORA_RUNTIME_ERROR;
    
  } else {

    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }  
    
    // inicializamos en cero los objetos elementales
    gsl_matrix_set_zero(fino.Ki);
    gsl_matrix_set_zero(fino.Mi);
    gsl_vector_set_zero(fino.bi);

    // TODO: hacer el loop de gauss adentro de cada funcion asi podemos
    // hacer evaluaciones nodo por nodo o lo que sea
    // para cada punto de gauss
    for (v = 0; v < element->type->gauss[GAUSS_POINTS_CANONICAL].V; v++) {
      // armamos las matrices
      if (fino.problem_family == problem_family_break || fino.problem_family == problem_family_shake) {
        wasora_call(fino_break_build_element(element, v));
      } else if (fino.problem_family == problem_family_bake) {
        wasora_call(fino_build_bake(element, v));
      }
    }

    MatSetValues(fino.K, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Ki, 0, 0), ADD_VALUES);
    if (fino.math_type == math_type_linear) {
      VecSetValues(fino.b, fino.elemental_size, fino.mesh->fem.l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
    }
    if (fino.has_mass)  {
      MatSetValues(fino.M, fino.elemental_size, fino.mesh->fem.l, fino.elemental_size, fino.mesh->fem.l, gsl_matrix_ptr(fino.Mi, 0, 0), ADD_VALUES);
    }

/*    
    printf("\nelement %d\n", element->id);
    fino_print_gsl_matrix(fino.Ki, stdout);
    printf("\n");
//    fino_print_gsl_matrix(fino.Mi, stdout);
    fino_print_gsl_vector(fino.bi, stdout);
    printf("\n");
 */
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
#define __FUNCT__ "fino_build_element_bc"
int fino_build_element_bc(element_t *element, bc_t *bc) {
  
  double n[3];
  
  if (fino.problem_family == problem_family_break) {
    // TODO: poner un flag si se necesita
    wasora_call(mesh_compute_outward_normal(element, n));
    wasora_var_value(fino.vars.nx) = n[0];
    wasora_var_value(fino.vars.ny) = n[1];
    wasora_var_value(fino.vars.nz) = n[2];
    
    // TODO: unificar todos como break_neumann
    if (bc->type_phys == bc_phys_stress) {
      wasora_call(fino_break_set_stress(element, bc));
    } else if (bc->type_phys == bc_phys_force) {
      wasora_call(fino_break_set_force(element, bc));
    } else if (bc->type_phys == bc_phys_pressure_normal ||
               bc->type_phys == bc_phys_pressure_real) {
      wasora_call(fino_break_set_pressure(element, bc));
    } else if (bc->type_phys == bc_phys_moment) {
      wasora_call(fino_break_set_moment(element, bc));
    }
  } else if (fino.problem_family == problem_family_bake) {
    if (bc->type_phys == bc_phys_heat_flux || bc->type_phys == bc_phys_heat_total) {
      if (strcmp(bc->expr[0].string, "0") != 0) { // para no tener que hacer cuentas si es adiabatico
        wasora_call(fino_bake_set_heat_flux(element, bc));
      }
    } else if (bc->type_phys == bc_phys_convection) {
      wasora_call(fino_bake_set_convection(element, bc));
    }
  }
  
/*  
  printf("\nelement %d\n", element->id);
  fino_print_gsl_vector(fino.bi, stdout);
  printf("\n");
 */  
  
  return WASORA_RUNTIME_OK;
  
}

#undef  __FUNCT__
#define __FUNCT__ "fino_compute_r_for_axisymmetric"
double fino_compute_r_for_axisymmetric(void) {

  double r_for_axisymmetric = 1.0;
  
  if (fino.problem_kind == problem_kind_axisymmetric) {
    if (fino.symmetry_axis == symmetry_axis_y) {
      if ((r_for_axisymmetric = gsl_vector_get(fino.mesh->fem.x, 0)) < 0) {
        wasora_push_error_message("axisymmetric problems with respect to y cannot have nodes with x < 0");
        return WASORA_RUNTIME_ERROR;
      }
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      if ((r_for_axisymmetric = gsl_vector_get(fino.mesh->fem.x, 1)) < 0) {
        wasora_push_error_message("axisymmetric problems with respect to x cannot have nodes with y < 0");
        return WASORA_RUNTIME_ERROR;
      }
    }
  }
  
  return r_for_axisymmetric;
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
