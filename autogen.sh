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
for i in git m4 autoconf uname xargs md5sum; do
 if [ -z "`which $i`" ]; then
  echo "error: $i not installed"
  exit 1
 fi
done

echo -n "autogen: cleaning... "
./autoclean.sh
echo "ok"

echo -n "autogen: updating wasora, "
if [ -e ../wasora/.git ]; then
  WASORA_REPO=../wasora/.git
else 
  WASORA_REPO=https://github.com/seamplex/wasora.git
fi

if [ ! -e wasora ]; then
  echo -n "cloning... ";
  git clone ${WASORA_REPO} || exit 1
else
  echo -n "pulling... "
  cd wasora; git checkout -q master; git pull || exit 1; cd ..
fi
echo "ok"

if [ -e sha_wasora ]; then
  echo -n "autogen: checking out... "
  cd wasora; git checkout -q master; git checkout -q `cat ../sha_wasora`; cd ..
  echo "ok"
fi

echo -n "autogen: bootstrapping... "
m4 wasora/m4/bootstrap.m4 - << EOF | sh -s $1 || exit 1
plugin=${plugin}
WASORA_CHECK_VCS
PLUGIN_VERSION_VCS
WASORA_README_INSTALL
EOF
echo "ok"

echo -n "autogen: building makefile.am... "
rm -f petscslepc.mak
touch petscslepc.mak

# check if texinfo is installed
if [ ! -z "`which texi2any`" ]; then
  docsubdir="doc"
  finoinfo="doc/fino.info"
fi

am="Makefile.am"
cat << EOF > ${am}
SUBDIRS = src ${docsubdir}
DIST_SUBDIRS = src wasora doc test

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO \\
                doc/fino.texi doc/fino.svg doc/fino.xml doc/reference.txt ${finoinfo}

dist_man_MANS = doc/fino.1

EXTRA_DIST = examples locateruntest.sh \\
             src/version.sh

TESTS = examples/test-tensile.sh \\
        examples/test-piston.sh \\
        examples/test-wire.sh \\
        examples/test-heated-cylinder.sh \\
        examples/test-cantilever.sh
        
all-local:
	if [ -e src/fino ]; then cp src/fino .; fi
	if [ -e fino ]; then cd examples; ln -sf ../fino; cd;  fi

clean-local:
	rm -f fino src/fino
EOF


am="src/Makefile.am"
echo -n "building ${am}... "

comment=""
if [ "x`uname`" = "xDarwin" ]; then
  comment="#"
fi

cat << EOF > ${am}
include ../petscslepc.mak
${comment} undefine DESTDIR  # this variable is set by petsc somewhere and we need it empty to make install" >> ${am}

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = \$(ACLOCAL_FLAGS)

bin_PROGRAMS = ${plugin}

${plugin}_CFLAGS = -I../wasora/src \$(SLEPC_CC_INCLUDES) \$(PETSC_CC_INCLUDES) \$(CC_INCLUDES) \$(all_includes) -DHARDCODEDPLUGIN 
${plugin}_LDADD = \$(STANDALONELIBS) \$(SLEPC_LIB) \$(PETSC_LIB) \$(all_libraries)
${plugin}_LDFLAGS = -rdynamic

${plugin}_SOURCES = \\
EOF

cd src
find . -maxdepth 1 \( -name "*.c" -o -name "*.h" \) | xargs echo -n >> ../${am}
echo "\\" >> ../${am}
find ../wasora/src -maxdepth 2 \( -name "*.c" -o -name "*.h" \) | xargs echo >> ../${am}
cd ..

# cat << EOF >> ${am}
# 
# version.\$(OBJEXT): version.h
# version.h: Makefile
# 	./version.sh
# EOF
echo "ok"

echo -n "autogen: autoreconfiguring... "
autoreconf -i
echo "ok"

