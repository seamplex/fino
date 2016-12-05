include ../petscslepc.mak

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = $(ACLOCAL_FLAGS)

# no se por que no camina como el de arriba
DEFAULT_INCLUDES = $(WASORA_INCLUDE) $(SLEPC_CC_INCLUDES) $(PETSC_CC_INCLUDES) $(CC_INCLUDES)

libfino_la_INCLUDES = $(all_includes)
lib_LTLIBRARIES = libfino.la
# libfino_la_LDFLAGS =
libfino_la_LIBADD = $(SLEPC_LIB) $(PETSC_LIB)

include $(SLEPC_DIR)/lib/slepc/conf/slepc_variables

libfino_la_SOURCES = ./version.h ./fino.c ./bc.c ./mat2sng.h ./eigen_slepc.c ./version.c ./fino.h ./parser.c ./times.c ./bulk.c ./init.c ./bake.c ./distributions.c ./gradient.c ./petschandler.c ./mat2sng.c ./linear_petsc.c ./breakshake.c ./debug.c ./version-vcs.h

version.$(OBJEXT): version.h
version.h: Makefile
	./version.sh
