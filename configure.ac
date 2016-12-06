dnl Process this file with autoconf to produce a configure script.
dnl
dnl This file is free software; as a special exception the author gives
dnl unlimited permission to copy and/or distribute it, with or without
dnl modifications, as long as this notice is preserved.
dnl
dnl This program is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
dnl implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

m4_include([version.m4])
AC_INIT(pluginname,pluginversion,[jeremy@seamplex.com])
AC_CONFIG_SRCDIR(src/pluginname.c)

AC_CONFIG_MACRO_DIR([wasora/m4])

WASORA_PLUGIN_INIT_C
WASORA_OPT_FLAGS_C

AC_CONFIG_FILES([Makefile src/Makefile])

WASORA_PLUGIN_VERSION_H

WASORA_CHECK_BASIC_HEADERS
WASORA_CHECK_BASIC_LIBS
WASORA_CHECK_GSL
WASORA_CHECK_IDA
WASORA_CHECK_READLINE

# fino needs petsc and optionally slepc
WASORA_CHECK_PETSC([required])
WASORA_CHECK_SLEPC([optional])

AC_OUTPUT

AS_BOX([Configuration summary for fino standalone executable])
WASORA_RESUME_LIBS
WASORA_RESUME_PETSC
WASORA_RESUME_SLEPC

# para saber que incluir en el makefile, si PETSc solito o SLEPc
AC_SUBST([PETSC_DIR])
AC_SUBST([PETSC_ARCH])
AC_SUBST([SLEPC_DIR])
AS_IF([test $slepc_found -ne 0],[
  echo "include \$(SLEPC_DIR)/lib/slepc/conf/slepc_variables" > petscslepc.mak
],[
  echo "include \$(PETSC_DIR)/lib/petsc/conf/variables" > petscslepc.mak
])
# listorti


AS_IF([ test $readline -eq 0 ], [
  WASORA_FINAL_WARN([0])
],[
  WASORA_FINAL_WARN([1])
])


