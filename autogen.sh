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
for i in git m4 autoconf xargs md5sum; do
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
am="src/Makefile.am"
echo -n "building $am... "
cat << EOF > $am
include ../petscslepc.mak
undefine DESTDIR  # this variable is set by petsc somewhere and we need it empty to make install

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = \$(ACLOCAL_FLAGS)

bin_PROGRAMS = ${plugin}

${plugin}_CFLAGS = -I../wasora/src \$(SLEPC_CC_INCLUDES) \$(PETSC_CC_INCLUDES) \$(CC_INCLUDES) \$(all_includes) -DHARDCODEDPLUGIN 
${plugin}_LDADD = \$(STANDALONELIBS) \$(SLEPC_LIB) \$(PETSC_LIB) \$(all_libraries)
${plugin}_LDFLAGS = -rdynamic

${plugin}_SOURCES = \\
EOF

cd src
find . -maxdepth 1 \( -name "*.c" -o -name "*.h" \) | xargs echo -n >> ../$am
echo "\\" >> ../$am
find ../wasora/src -maxdepth 2 \( -name "*.c" -o -name "*.h" \) | xargs echo >> ../$am
cd ..

# cat << EOF >> $am
# 
# version.\$(OBJEXT): version.h
# version.h: Makefile
# 	./version.sh
# EOF
echo "ok"

echo -n "autogen: autoreconfiguring... "
autoreconf -i
echo "ok"

