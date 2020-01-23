/*------------ -------------- -------- --- ----- ---   --       -            -
 *  fino's version banner
 *
 *  Copyright (C) 2015--2020 jeremy theler
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
#include <stdlib.h>
#include <signal.h>

#include "fino.h"
#include "version.h"

// global static so the compiler locates these strings in the text section
// so when the plugin_* functions return pointers to the strings there is
// no need to free them afterward
const char finoname[] = "Fino";
const char finodescription[] = "a free finite-element thermo-mechanical solver";
char finoshortversion[128];
char finolongversion[2048];
const char finousage[] = "no commandline options needed";

const char finocopyright[] = "\
Fino is copyright (c) 2015--2020 jeremy theler\n\
and licensed under GNU GPL version 3 or later.\n\
Fino is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.";

const char finohmd5[] = PLUGIN_HEADERMD5;


const char *plugin_name(void) {
  return finoname;
}

const char *plugin_longversion(void) {
  char slepcversion[BUFFER_SIZE];
  char petscversion[BUFFER_SIZE];
  char petscarch[BUFFER_SIZE];
  
  PetscGetVersion(petscversion, BUFFER_SIZE);
  PetscGetArchType(petscarch, BUFFER_SIZE);
#ifdef HAVE_SLEPC
  SlepcGetVersion(slepcversion, BUFFER_SIZE);  
#else
  snprintf(slepcversion, BUFFER_SIZE-1, "none");
#endif  

  sprintf(finolongversion,"\
Last commit date   : %s\n\
Build date         : %s\n\
Build architecture : %s\n\
Builder            : %s@%s\n\
Compiler           : %s\n\
Compiler flags     : %s\n\
PETSc version      : %s\n\
PETSc arch         : %s\n\
SLEPc version      : %s",
   PLUGIN_VCS_DATE,
   COMPILATION_DATE,
   COMPILATION_ARCH,
   COMPILATION_USERNAME,
   COMPILATION_HOSTNAME,
   CCOMPILER_VERSION,
   CCOMPILER_FLAGS,
   petscversion,
   petscarch,
   slepcversion);
  
  return finolongversion;
}

const char *plugin_wasorahmd5(void) {
  return finohmd5;
}
const char *plugin_copyright(void) {
  return finocopyright;
}


const char *plugin_version(void) {
#ifdef PLUGIN_VCS_BRANCH
  sprintf(finoshortversion, "%s%s %s", PLUGIN_VCS_VERSION,
                                       (PLUGIN_VCS_CLEAN==0)?"":"+Δ",
                                       strcmp(PLUGIN_VCS_BRANCH, "master")?PLUGIN_VCS_BRANCH:"");
#else
  sprintf(finoshortversion, "%s", PACKAGE_VERSION);
#endif

  return finoshortversion;
}

const char *plugin_description(void) {
  return finodescription;
}

const char *plugin_usage(void) {
  return finousage;
}
