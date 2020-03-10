#!/bin/sh
# this file should not be called directly from the shell,
# it ought to be included from within run.sh for each case

if [[ -x ./fino ]] && [[ ! -h ./fino ]]; then
 finobin=./fino
 testdir="examples/"
elif [[ -x ../fino ]] && [[ ! -d ../fino ]]; then
 finobin="../fino"
 testdir="./"
elif [ ! -z "`which wasora`" ]; then
 if [ -x ./fino.so ]; then
  finobin="wasora -p fino.so"
  testdir="examples/"
 elif [ -x ../fino.so ]; then
  finobin="wasora -p ../fino.so"
  testdir=""
 fi
else
 echo "do not know how to run fino :("
 exit 1
fi

# runs fino
function runfino {
 ${finobin} ${testdir}${1} ${2} ${3} ${4} ${5} ${6}
 outcome=$?
 if [ ${outcome} -ne 0 ]; then
   exit ${outcome}
 fi
}

# calls gnuplot with the provided command if it is installed
function plot {
 if [ "x`which gnuplot`" != "x" ]; then
  if [ "x`uname | cut -c-6`" = "xCYGWIN" ]; then
   if [ "x`ps -e | grep X | wc -l`" = "x0" ]; then
    XWin.exe -multiwindow -clipboard -silent-dup-error > /dev/null &
    sleep 2
   fi
   export DISPLAY=:0.0
  fi
  gnuplot -p -e "$1"
 fi
}

# checks if gmsh is installed
function checkgmsh {
 if [ -z "`which gmsh`" ]; then
  echo "gmsh is not installed, skipping test"
  exit 77
 fi
}

# checks if m4 is installed
function checkm4 {
 if [ -z "`which m4`" ]; then
  echo "m4 is not installed, skipping test"
  exit 77
 fi
}

# calls gmsh (if exists)
function callgmsh {
 if [ ! -z "`which gmsh`" ]; then
  gmsh $1 &
 fi
}


# calls pandoc (or markdown) and shows the result
function callpandoc {
 if [ ! -z "`which pandoc`" ]; then
   pandoc $1.txt -o $1.html -s
   pandoc $1.txt -o $1.pdf
 elif [ ! -z "`which markdown`" ]; then
   markdown $1.txt > $1.html
 fi

 if [ -f $1.html ]; then
   if [ ! -z "`which x-www-browser`" ]; then
     x-www-browser $1.html &
   elif [ ! -z "`which xdg-open`" ]; then
     xdg-open $1.html &
   fi
 fi
}
