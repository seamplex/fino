include ../petscslepc.mak

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = $(ACLOCAL_FLAGS)

bin_PROGRAMS = fino

fino_CFLAGS = -I../wasora/src $(SLEPC_CC_INCLUDES) $(PETSC_CC_INCLUDES) $(CC_INCLUDES) $(all_includes) -DHARDCODEDPLUGIN 
fino_LDADD = $(STANDALONELIBS) $(SLEPC_LIB) $(PETSC_LIB) $(all_libraries)
fino_LDFLAGS = -rdynamic

fino_SOURCES = \
./fino.c ./bc.c ./mat2sng.h ./eigen_slepc.c ./version.c ./fino.h ./parser.c ./times.c ./bulk.c ./init.c ./bake.c ./distributions.c ./gradient.c ./petschandler.c ./mat2sng.c ./linear_petsc.c ./breakshake.c ./debug.c ./version-vcs.h\
../wasora/src/history.c ../wasora/src/io.c ../wasora/src/interface.c ../wasora/src/print.c ../wasora/src/realtime.c ../wasora/src/builtindecl.h ../wasora/src/algebra.c ../wasora/src/vector.c ../wasora/src/builtinfunctions.c ../wasora/src/randomline.c ../wasora/src/builtinfunctionals.c ../wasora/src/flow.c ../wasora/src/version.h ../wasora/src/error.c ../wasora/src/builtin.h ../wasora/src/shmem.c ../wasora/src/cleanup.c ../wasora/src/parseaux.c ../wasora/src/wasora.h ../wasora/src/shell.c ../wasora/src/fit.c ../wasora/src/assignment.c ../wasora/src/version.c ../wasora/src/file.c ../wasora/src/parser.c ../wasora/src/init.c ../wasora/src/matrix.c ../wasora/src/multiminf.c ../wasora/src/version-conf.h ../wasora/src/interface.h ../wasora/src/function.c ../wasora/src/thirdparty/utlist.h ../wasora/src/thirdparty/kdtree.h ../wasora/src/thirdparty/kdtree.c ../wasora/src/thirdparty/uthash.h ../wasora/src/builtinvectorfunctions.c ../wasora/src/multiminfdf.c ../wasora/src/alias.c ../wasora/src/m4.c ../wasora/src/getptr.c ../wasora/src/call.c ../wasora/src/debug.c ../wasora/src/parametric.c ../wasora/src/dae.c ../wasora/src/version-vcs.h ../wasora/src/multiroot.c ../wasora/src/minimize.c ../wasora/src/mesh/element.c ../wasora/src/mesh/findmax.c ../wasora/src/mesh/fem.c ../wasora/src/mesh/fillvector.c ../wasora/src/mesh/neighbors.c ../wasora/src/mesh/point.c ../wasora/src/mesh/integrate.c ../wasora/src/mesh/prism.c ../wasora/src/mesh/vtk.c ../wasora/src/mesh/parser.c ../wasora/src/mesh/init.c ../wasora/src/mesh/geom.c ../wasora/src/mesh/cell.c ../wasora/src/mesh/line.c ../wasora/src/mesh/interpolate.c ../wasora/src/mesh/tetrahedron.c ../wasora/src/mesh/hexahedron.c ../wasora/src/mesh/struct.c ../wasora/src/mesh/post.c ../wasora/src/mesh/gmsh.c ../wasora/src/mesh/quad.c ../wasora/src/mesh/mesh.c ../wasora/src/mesh/gauss.c ../wasora/src/mesh/define.c ../wasora/src/mesh/triang.c ../wasora/src/plugin.c ../wasora/src/define.c ../wasora/src/wasora.c

version.$(OBJEXT): version.h
version.h: Makefile
	./version.sh
