/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's bulk routines
 *
 *  Copyright (C) 2015-2020 Seamplex
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
#include <stdio.h>

#ifndef _FINO_H
#include "fino.h"
#endif

#define NAME_SIZE 32

int fino_allocate_elemental_objects(element_t *element) {

  fino_free_elemental_objects();
      
  fino.n_local_nodes = element->type->nodes;
  fino.elemental_size = element->type->nodes * fino.degrees;
  
  fino.Ki = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  fino.Mi = gsl_matrix_calloc(fino.elemental_size, fino.elemental_size);
  fino.bi = gsl_vector_calloc(fino.elemental_size);
  
  return WASORA_RUNTIME_OK;

}


int fino_free_elemental_objects(void) {

  if (fino.n_local_nodes != 0 && fino.elemental_size != 0) {
    gsl_matrix_free(fino.Ki);
    gsl_matrix_free(fino.Mi);
    gsl_vector_free(fino.bi);
  }
  
  fino.n_local_nodes = 0;
  fino.elemental_size = 0;
  
  
  return WASORA_RUNTIME_OK;

}

int fino_build_bulk(void) {

  bc_t *bc;
  int i;
  int step = 0;
  int ascii_progress_chars = 0;
  
  step = ceil((double)fino.mesh->n_elements/100.0);
  if (step < 1) {
    step = 1;
  }
  
  // empty global objects
  if (fino.K != NULL) {
    petsc_call(MatZeroEntries(fino.K));
  }
  if (fino.M != NULL) {
    petsc_call(MatZeroEntries(fino.M));
  }
  if (fino.b != NULL) {
    petsc_call(VecZeroEntries(fino.b));
  }  
  
  for (i = fino.first_element; i < fino.last_element; i++) {

// ------ progress bar ------------------------------------------    
    // only the first time we build
    if (fino.progress_ascii == PETSC_TRUE &&
        fino.first_build == PETSC_TRUE &&
        (i % step) == 0) {
      printf(CHAR_PROGRESS_BUILD);  
      fflush(stdout);
      ascii_progress_chars++;
    }
// --------------------------------------------------------------    

    if (fino.mesh->element[i].type->dim == fino.dimensions) {
       
      // only the elements that have the problem dimension are used to build
      // the elemental matrices of the bulk problem
      wasora_call(fino_build_element_volumetric(&fino.mesh->element[i]));
      
    } else if (fino.mesh->element[i].type->dim < fino.dimensions &&
               fino.mesh->element[i].physical_entity != NULL) {
      
      LL_FOREACH(fino.mesh->element[i].physical_entity->bcs, bc) {
        
        if (bc->type_math == bc_math_neumann || bc->type_math == bc_math_robin) {
          if (bc->condition.n_tokens == 0 || fabs(wasora_evaluate_expression(&bc->condition)) > 1e-3) {
            // only non-dirichlet BCs here, dirichlet is handled in set_essential after the assembly
            wasora_call(fino_build_element_bc(&fino.mesh->element[i], bc));
          }  
        }
      }
    }
  }

  wasora_call(fino_assembly());

  // copy the assembled stiffness matrix before setting dirichlet BCs
  // this is used to compute reactions or resultants at arbitrary locations
  // and to compute the strain energy for non-homogeneous BCs
  if (fino.K_nobc == NULL) {
    petsc_call(MatDestroy(&fino.K_nobc));
    petsc_call(MatDuplicate(fino.K, MAT_DO_NOT_COPY_VALUES, &fino.K_nobc));
    if (fino.b != NULL) {
      petsc_call(VecDuplicate(fino.b, &fino.b_nobc));
    }  
  }
  
  // just in case we need the RHS vector before setting the BCs
  petsc_call(MatCopy(fino.K, fino.K_nobc, SAME_NONZERO_PATTERN));
  if (fino.b != NULL) {
    petsc_call(VecCopy(fino.b, fino.b_nobc));
  }  

  if (fino.progress_ascii == PETSC_TRUE && fino.first_build == PETSC_TRUE) {
    if (wasora.nprocs == 1) {
      while (ascii_progress_chars++ < 100) {
        printf(CHAR_PROGRESS_BUILD);
      }
    }
    if (wasora.rank == 0) {
      printf("\n");  
      fflush(stdout);
    }  
  }
  
  // mark a flag that we already built the matrices
  fino.first_build = PETSC_FALSE;
  
  // C and et are lost here (they are static)
  wasora_call(fino_free_elemental_objects());
  
  return WASORA_RUNTIME_OK;

}

int fino_build_element_volumetric(element_t *element) {
  int V;           // total number of Gauss points
  int v;           // gauss point index 

  if (element->physical_entity == NULL) {
    // this (should) happen only in structured grids
    wasora_push_error_message("volumetric element %d does not have an associated physical group", element->tag);
    return WASORA_RUNTIME_ERROR;
    
  } else {

    V = element->type->gauss[fino.mesh->integration].V;
    
    if (fino.n_local_nodes != element->type->nodes) {
      wasora_call(fino_allocate_elemental_objects(element));
    }
    
    if (element->B == NULL) {
      element->B = calloc(V, sizeof(gsl_matrix *));
    }
    
    // initialize to zero the elemental objects
    gsl_matrix_set_zero(fino.Ki);
    gsl_matrix_set_zero(fino.Mi);
    gsl_vector_set_zero(fino.bi);

    // TODO: hacer el loop de gauss adentro de cada funcion asi podemos
    // hacer evaluaciones nodo por nodo o lo que sea para cada punto de gauss
    for (v = 0; v < V; v++) {
      
      // build elementary matrices
      // TODO: use function pointers
      if (fino.problem_family == problem_family_mechanical || fino.problem_family == problem_family_modal) {
        wasora_call(fino_break_build_element(element, v));
      } else if (fino.problem_family == problem_family_thermal) {
        wasora_call(fino_thermal_build_element(element, v));
      }
    }
    
    // the indices of the DOFs to ensamble
    mesh_compute_l(fino.mesh, element);

    petsc_call(MatSetValues(fino.K, fino.elemental_size, element->l, fino.elemental_size, element->l, gsl_matrix_ptr(fino.Ki, 0, 0), ADD_VALUES));
    if (fino.b != NULL) {
      petsc_call(VecSetValues(fino.b, fino.elemental_size, element->l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES));
    }
    if (fino.M != NULL)  {
      petsc_call(MatSetValues(fino.M, fino.elemental_size, element->l, fino.elemental_size, element->l, gsl_matrix_ptr(fino.Mi, 0, 0), ADD_VALUES));
    }
  }

  return WASORA_RUNTIME_OK;
}

int fino_build_element_bc(element_t *element, bc_t *bc) {

  if (fino.n_local_nodes != element->type->nodes) {
    wasora_call(fino_allocate_elemental_objects(element));
  }
  
  if (fino.Nb == NULL) {
    fino.Nb = gsl_vector_calloc(fino.degrees);
  }
  gsl_vector_set_zero(fino.Nb);
  gsl_vector_set_zero(fino.bi);

  if (fino.problem_family == problem_family_mechanical) {
    
    // TODO: poner un flag si se necesita
    if (element->type->dim == 1 || element->type->dim == 2) {
      wasora_call(mesh_compute_normal(element));
    }  
    
    wasora_call(fino_break_set_neumann(element, bc));
    
  } else if (fino.problem_family == problem_family_thermal) {
    
    if (bc->type_phys == bc_phys_heat_flux || bc->type_phys == bc_phys_heat_total) {
      if (strcmp(bc->expr[0].string, "0") != 0) { // para no tener que hacer cuentas si es adiabatico
        wasora_call(fino_thermal_set_heat_flux(element, bc));
      }
    } else if (bc->type_phys == bc_phys_convection) {
      wasora_call(fino_thermal_set_convection(element, bc));
    }
    
  }
  
  mesh_compute_l(fino.mesh, element);
  VecSetValues(fino.b, fino.elemental_size, element->l, gsl_vector_ptr(fino.bi, 0), ADD_VALUES);
  
  return WASORA_RUNTIME_OK;
  
}

inline double fino_compute_r_for_axisymmetric(element_t *element, int v) {

  double r_for_axisymmetric = 1.0;
  
  if (fino.problem_kind == problem_kind_axisymmetric) {
    if (element->x == NULL || element->x[v] == NULL) {
      mesh_compute_x_at_gauss(element, v, fino.mesh->integration);
    }
  
    if (fino.symmetry_axis == symmetry_axis_y) {
      if ((r_for_axisymmetric = element->x[v][0]) < ZERO) {
        wasora_push_error_message("axisymmetric problems with respect to y cannot have nodes with x <~ 0");
        return WASORA_RUNTIME_ERROR;
      }
    } else if (fino.symmetry_axis == symmetry_axis_x) {
      if ((r_for_axisymmetric = element->x[v][1]) < ZERO) {
        wasora_push_error_message("axisymmetric problems with respect to x cannot have nodes with y <~ 0");
        return WASORA_RUNTIME_ERROR;
      }
    }
  }
  
  return r_for_axisymmetric;
}

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
