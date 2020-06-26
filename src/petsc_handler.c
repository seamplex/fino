/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's petsc error handling routine
 *
 *  Copyright (C) 2015 jeremy theler
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

PetscErrorCode fino_handler(MPI_Comm comm, int line, const char *func, const char *file, PetscErrorCode n, PetscErrorType p, const char *mess, void *ctx) {
  wasora_push_error_message("PETSc error %d-%d '%s' in %s %s:%d", n, p, mess, file, func, line);
  wasora_runtime_error();
  return WASORA_RUNTIME_OK;
}


