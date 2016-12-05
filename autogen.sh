#!/bin/sh
# 
# Execute this script to generate a configure script
#
# This file is free software: you are free to change and redistribute it.
# There is NO WARRANTY, to the extent permitted by law.
#
plugin=fino

# check for needed tools (cannot put this into an m4 macro
# because we do not even know if we have m4 available)
for i in m4 autoconf libtoolize xargs; do
 if [ -z "`which $i`" ]; then
  echo "error: $i not installed"
  exit 1
 fi
done

# ideally we would use wasora's m4/bootstrap.m4 but we do not
# know where wasora is at this point! Sadly we have to use a certain
# bootstrap.m4 included in the plugin's tree
m4 m4/bootstrap.m4 - << EOF | sh -s $1 || exit 1
plugin=${plugin}
WASORA_CHECK_VCS
PLUGIN_AUTOCLEAN
PLUGIN_FIND_WASORA
PLUGIN_VERSION_VCS
PLUGIN_README_INSTALL
PLUGIN_COPY_M4
EOF


# build standalone's makefile.am
am="src/wasora/Makefile.am"
echo -n "building $am... "
cat << EOF > $am
include ../../petscslepc.mak

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = \$(ACLOCAL_FLAGS)

bin_PROGRAMS = ${plugin}

${plugin}_INCLUDES = \$(WASORA_INCLUDE) \$(SLEPC_CC_INCLUDES) \$(PETSC_CC_INCLUDES) \$(CC_INCLUDES) \$(all_includes)
${plugin}_CFLAGS = -DHARDCODEDPLUGIN
${plugin}_LDADD = ../.libs/lib${plugin}.a \$(STANDALONELIBS) \$(SLEPC_LIB) \$(PETSC_LIB) \$(all_libraries)
${plugin}_LDFLAGS = -rdynamic

${plugin}_SOURCES = `cd src/wasora; find . -maxdepth 2 \( -name "*.c" -o -name "*.h" \) | xargs; cd ../..`
EOF
echo "done"

# build plugins's makefile.am
am="src/Makefile.am"
echo -n "building $am... "
cat << EOF > $am
include ../petscslepc.mak

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = \$(ACLOCAL_FLAGS)

# no se por que no camina como el de arriba
DEFAULT_INCLUDES = \$(WASORA_INCLUDE) \$(SLEPC_CC_INCLUDES) \$(PETSC_CC_INCLUDES) \$(CC_INCLUDES)

lib${plugin}_la_INCLUDES = \$(all_includes)
lib_LTLIBRARIES = lib${plugin}.la
# lib${plugin}_la_LDFLAGS =
lib${plugin}_la_LIBADD = \$(SLEPC_LIB) \$(PETSC_LIB)

include \$(SLEPC_DIR)/lib/slepc/conf/slepc_variables

lib${plugin}_la_SOURCES = ./version.h `cd src; find . -maxdepth 1 \( -name "*.c" -o -name "*.h" \) | xargs; cd ..`

version.\$(OBJEXT): version.h
version.h: Makefile
	./version.sh
EOF
echo "done"

echo "calling autoreconf... "
touch petscslepc.mak
autoreconf -i
echo "done"

