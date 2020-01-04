/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's routines for export of post-processing views
 *
 *  Copyright (C) 2020 Seamplex
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

#include "fino.h"


#undef  __FUNCT__
#define __FUNCT__ "fino_instruction_post"
int fino_instruction_post(void *arg) {
  
  mesh_t *mesh;
  element_t *element;
  node_t *node;
  element_list_item_t *element_list;
  
//  int d;
  int i, i_global;
  int j, j_global;
  
  // primera pasada, copiamos lo que podemos, alocamos y contamos la cantidad de nodos
  mesh = calloc(1, sizeof(mesh_t));
  mesh->bulk_dimensions = fino.mesh->bulk_dimensions;
  mesh->n_elements = fino.mesh->n_elements;
  mesh->element = calloc(mesh->n_elements, sizeof(element_t));
  mesh->n_nodes = 0;
  i_global = 0;
  for (i = 0; i < mesh->n_elements; i++) {
    if (fino.mesh->element[i].type->dim == mesh->bulk_dimensions) {
      element = &mesh->element[i_global];

      element->index = i_global;
      element->tag = i+1;
      element->type = fino.mesh->element[i].type;
      element->physical_entity = fino.mesh->element[i].physical_entity;

      mesh->n_nodes += element->type->nodes;
      i_global++;
    }  
  }
  
  mesh->n_elements = i_global;
  mesh->element = realloc(mesh->element, mesh->n_elements*sizeof(element_t));
  
  // segunda pasada, creamos los nodos
  j_global = 0;
  mesh->node = calloc(mesh->n_nodes, sizeof(node_t));
  for (i = 0; i < mesh->n_elements; i++) {
    element = &mesh->element[i];
    element->node = calloc(element->type->nodes, sizeof(node_t));
    
    for (j = 0; j < element->type->nodes; j++) {
      node = &mesh->node[j_global];
      node->tag = j_global+1;
      node->index_mesh = j_global;
      node->x[0] = fino.mesh->element[element->tag-1].node[j]->x[0];
      node->x[1] = fino.mesh->element[element->tag-1].node[j]->x[1];
      node->x[2] = fino.mesh->element[element->tag-1].node[j]->x[2];
      
/*      
      node->phi = calloc(fino.degrees, sizeof(double));
      for (d = 0; d < fino.degrees; d++) {
        node->phi[d] = fino.mesh->element[element->tag-1].node[j]->phi[d];
      }
 */
      node->dphidx = fino.mesh->element[element->tag-1].dphidx_node[j];
      
      element_list = calloc(1, sizeof(element_list_item_t));
      element_list->element = element;
      LL_APPEND(node->associated_elements, element_list);
      
      element->node[j] = node;
      
      j_global++;
    }
  }
  
  mesh_vtk_write_header(stdout);
  mesh_vtk_write_mesh(mesh, 0, stdout);
  
  printf("POINT_DATA %d\n", mesh->n_nodes);
  printf("SCALARS sigmax double\n");
  printf("LOOKUP_TABLE default\n");
  for (j = 0; j < mesh->n_nodes; j++) {
    printf("%.1f\n", 1.0*j);
  }
      
      
  
  return WASORA_RUNTIME_OK;

}