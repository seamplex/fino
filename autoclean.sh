#!/bin/sh
# 
# Execute this script to clean the directory and leave it as
# a fresh git repo
#
# This file is free software: you are free to change and redistribute it.
# There is NO WARRANTY, to the extent permitted by law.
#
if test -e Makefile; then
  make clean
fi
rm -f *~ src/*~ src/thirdparty/*~
rm -f src/*.o src/*.lo src/*~
rm -rf src/.deps src/.libs src/.dirstamp src/stamp-h1 src/config.h.in src/config.h
rm -f README INSTALL
rm -f aclocal.m4 configure config.log config.status compile depcomp install-sh missing ltmain.sh config.guess config.sub libtool libtool test-driver
rm -f src/version.h src/version-vcs.h src/version-conf.h version.m4
rm -rf autom4te.cache
rm -rf src/.deps src/.libs src/.dirstamp src/stamp-h1 src/config.h.in src/config.h
rm -f Makefile Makefile.in src/Makefile src/Makefile.in
rm -rf ${plugin}.so ${plugin} src/lib${plugin}.la 
rm -rf src/wasora
cd examples
./clean.sh
cd ..
