#!/bin/sh
# set -e

#
# This script is meant for quick & easy install via:
#   'curl https://www.seamplex.com/fino/get.sh | sh'
# or:
#   'wget -O- https://www.seamplex.com/fino/get.sh | sh'

# loosely based on docker's get.sh script at https://get.docker.com/
# jeremy theler

name=fino
petsc_version=3.7.4

# check for needed tools
command_exists() {
  command -v "$@" > /dev/null 2>&1
}


# step 0: needed commands
# if wasora-suite does not exist, create
if test "`whoami`" = "root"; then
  echo "error: do not run me as root"
  exit 1
fi


rm -f get.log
echo -n "0. checking pre-requisistes..." | tee -a get.log; echo >> get.log
for i in m4 autoconf xargs gcc git wget make; do
 if ! command_exists $i; then
  echo
  echo "error: $i not installed"
  echo "install the following packages using your favorite package manager:"
  echo "m4 autoconf gcc git make"
  exit 1
 fi
done
echo "ok!"

# step 1: wasora-suite directory
# if wasora-suite does not exist, create
if test ! -d ./wasora-suite; then
  echo -n "1. creating wasora-suite subdir..." | tee -a get.log; echo >> get.log
  mkdir wasora-suite
  echo "ok!"
else
  echo "1. wasora-suite subdir exists, good!"  | tee -a get.log; echo >> get.log
fi
  
cd wasora-suite

# step 2: clone or update repo
if test ! -d ./${name}; then
  echo -n "2. cloning ${name} repository..."  | tee -a ../get.log; echo >> ../get.log
  git clone https://bitbucket.org/seamplex/${name}.git >> ../get.log; echo >> ../get.log
  cd ${name}
  echo "ok!"
else
  echo -n "2. pulling and updating ${name} repository..."  | tee -a ../get.log; echo >> ../get.log
  cd ${name}
  git pull  >> ../../get.log; echo >> ../../get.log
  echo "ok!"
fi


# step 3: bootstrap repo
echo -n "3. bootstrapping ${name} source tree..."  | tee -a ../../get.log; echo >> ../../get.log
./autogen.sh 2>&1 >> ../../get.log; echo >> ../../get.log
echo "ok!"

# step 4: check for GSL
echo -n "4. checking for GSL..." | tee -a ../../get.log; echo >> ../../get.log
cat << EOF > test.c
#include <gsl/gsl_version.h>

int main(void) {
  return 0;
}
EOF

gcc test.c -lgsl -lgslcblas > /dev/null 2>&1
nodynamicgsl=$?
rm -rf test.c a.out

if test ${nodynamicgsl} -eq 1; then
 echo "no"
 echo
 echo "The development version of GSL is not installed, you may:"
 echo "  a. let me try to install the package for yourself (needs sudo)"
 echo "  b. ask ${name} to try to download  and compile the GSL source"
 echo "  c. abort this script now, manually install package libgsl0-dev and re-run"
 echo -n "What do you want to do? (a/b/c) "
 
 read response
 
 case ${response} in
  "a")
    if command_exists apt-get; then
     sudo apt-get install libgsl0-dev libsundials-serial-dev
     if test $? -ne 0; then
       echo "error: there was some error with sudo, I cannot install GSL for you"
       echo "please manually install package libgsl0-dev and re-run this script"
       exit 1
     fi
    else
     echo "error: I could not find apt-get so I cannot install GSL for you"
     echo "please manually install package libgsl0-dev and re-run this script"
     exit 1
    fi;;
  "b")
    echo "ok, I will ask configure to download and compile it" | tee -a ../../get.log; echo >> ../../get.log
    echo "note that step number 5 may take some time" | tee -a ../../get.log; echo >> ../../get.log
    configflag="--enable-download-gsl"
    ;;
  *)
  echo "please manually install package libgsl0-dev and re-run this script"
  exit;;
 esac
 
else
 configflag=""
 echo "already installed!" | tee -a ../../get.log; echo >> ../../get.log
fi

# step 4.5: check for PETSc/SLEPc
echo -n "4.5. checking for PETSc/SLEPc..." | tee -a ../../get.log
if test -z "$PETSC_DIR"; then
 echo "PETSC_DIR variable not found, so I assume PETSc is not installed. You may:"
 echo "  a. let me try to download and compile the library (may take some time)"
 echo "  b. abort this script now, manually install the library or set the enviornment variables re-run"
 echo -n "What do you want to do? (a/b) "

 read response

 case ${response} in
  "a")
    mkdir -p libs
    cd libs
    if test ! -f petsc-lite-${petsc_version}.tar.gz; then
     if ! command_exists wget; then
      echo "error: wget not installed"
      exit 1
     fi
     wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${petsc_version}.tar.gz
    fi
    tar xvzf petsc-lite-${petsc_version}.tar.gz 
    cd petsc-${petsc_version}/
    export PETSC_DIR=$PWD
    export PETSC_ARCH=arch-linux2-c-opt
    ./configure --download-fblaslapack --with-mpi=0 --with-debugging=0
    make
    make test
    cat >> $HOME/.bashrc << EOF
export PETSC_DIR=${PETSC_DIR}      # added by ${name} get.sh
export PETSC_ARCH=${PETSC_ARCH}    # added by ${name} get.sh
EOF
    cd ../..
    ;;
  *)
   echo "please manually install or fix PETSc and re-run this script"
   exit;;
 esac

elif test ! -d ${PETSC_DIR}; then
 echo "PETSC_DIR is ${PETSC_DIR} but it is not a directory"
 echo "please manually install or fix PETSc re-run this script"
 exit

elif test ! -z "${PETSC_ARCH}"; then
 echo "PETSC_ARCH is empty"
 echo "please manually install or fix PETSc re-run this script"
 exit
fi

# step 5. configure
echo -n "5. configuring ${name}..."  | tee -a ../../get.log; echo >> ../../get.log
./configure ${configflag} >> ../../get.log; echo >> ../../get.log
echo "ok!"

# step 6. compile
echo -n "6. compiling ${name}..."  | tee -a ../../get.log; echo >> ../../get.log
make >> ../../get.log; echo >> ../../get.log
echo "ok!"

# step 7. link to $HOME/bin
echo -n "7. linking $HOME/bin/${name} to the recently compiled binary..."
if test ! -d $HOME/bin; then
  mkdir -p $HOME/bin
fi
rm -rf $HOME/bin/${name}
ln -s $PWD/${name} $HOME/bin/${name}
echo "ok!"

# step 8. add $HOME/bin to the path
if test -z "`echo $PATH | grep $HOME/bin`"; then
  echo -n "8. adding $HOME/bin to the path in ~/.bashrc"
  cat >> $HOME/.bashrc << EOF
export PATH=\$PATH:\$HOME/bin      # added by ${name} get.sh
EOF
  export PATH=\$PATH:\$HOME/bin
  echo "ok!"
else
  echo "8. directory $HOME/bin is already in the PATH envar, good!"
fi

cd ../..

# test if it worked!
echo
if command_exists ${name}; then
  echo "congratulations! ${name} compiled and installed correctly!"
  echo "type ‘${name}’ (without the quotes) if you do not believe me."
else
  echo "something went wrong as ${name} does not run :-("
  echo "share your get.log file with the mailing list for help."
  echo "see <https://groups.google.com/a/seamplex.com/forum/#!forum/wasora>"
fi

