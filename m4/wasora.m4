dnl This file is part of wasora and/or one of its plugins
dnl GPL v3+ (c) 2009-2015 jeremy theler
dnl <http://bitbucket.org/wasora/wasora>
dnl

AC_DEFUN([WASORA_PLUGIN_INIT_C],[
AM_INIT_AUTOMAKE
LT_INIT
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CC
AM_PROG_CC_C_O
AC_CANONICAL_HOST
])


AC_DEFUN([WASORA_CHECK_BASIC_HEADERS],[
# checks for header files.
#AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stdlib.h string.h sys/time.h unistd.h sys/select.h stddef.h],[],AC_MSG_ERROR([required header not found.]))
#AC_CHECK_FUNCS([floor ftruncate gettimeofday memset mkdir munmap pow rint sqrt strcasecmp strchr strcspn strdup strspn strerror strpbrk strrchr strstr select clock_gettime strtok_r sem_open],[],AC_MSG_ERROR([required C library function not found.]))
AC_CHECK_FUNCS([strtok_r usleep])
])


AC_DEFUN([WASORA_OPT_FLAGS_C],[
# default is optimized without debugging symbols
AS_IF([test "$CFLAGS" = "-g -O2"], [CFLAGS="-O2"])
])

AC_DEFUN([WASORA_OPT_FLAGS_F],[
# default is optimized without debugging symbols
AS_IF([test "$FCFLAGS" = "-g -O2"], [FCFLAGS="-O2"])
])

AC_DEFUN([WASORA_OPT_FLAGS_F77],[
# default is optimized without debugging symbols
AS_IF([test "$FFLAGS" = "-g -O2"], [FFLAGS="-O2"])
])

AC_DEFUN([WASORA_OPT_FLAGS_CXX],[
# default is optimized without debugging symbols
AS_IF([test "$CXXFLAGS" = "-g -O2"], [CXXFLAGS="-O2"])
])


AC_DEFUN([WASORA_CHECK_BASIC_LIBS],[
# libraries
# si tenemos rt bien, sino no importa
AC_CHECK_LIB([rt],[shm_open])
# idem
AC_CHECK_LIB([pthread],[pthread_create])
#AC_CHECK_HEADERS([sys/mman.h sys/stat.h fcntl.h],[],AC_MSG_ERROR([headers for shm_open not found]))

# parece que dlopen puede estar en la libc
AC_SEARCH_LIBS([dlopen],[dl dld],[])
#AC_CHECK_HEADER([dlfcn.h],[],AC_MSG_ERROR([dlfcn.h not found]))

# lo mismo -lm
AC_SEARCH_LIBS([cos],[m],[],AC_MSG_ERROR([libm not found]))
#AC_CHECK_HEADER([math.h],[],AC_MSG_ERROR([math.h not found]))
])


AC_DEFUN([WASORA_CHECK_GSL],[
AC_CHECK_HEADER([gsl/gsl_vector.h])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_ARG_ENABLE([download-gsl],
  [AS_HELP_STRING([--enable-download-gsl], [try to automatically download GSL @<:@default=no@:>@])],
  [download_gsl=yes],
  [download_gsl=no])

AS_IF([test "x$ac_cv_header_gsl_gsl_vector_h" != "xyes" -o "x$ac_cv_lib_gslcblas_cblas_dgemm" != "xyes" -o "x$ac_cv_lib_gsl_gsl_blas_dgemm" != "xyes"],
 gsldist=gsl-1.16
 gslmirror=http://ftpmirror.gnu.org/gsl/${gsldist}.tar.gz

 AS_IF([test ! -e ${gsldist}.tar.gz],[
   AS_IF([test "x$download_gsl" = "xyes"],[
     AS_IF([test "x`which wget`" != "x"],[
       AC_MSG_NOTICE([downloading ${gslmirror}])
       wget -c http://ftpmirror.gnu.org/gsl/${gsldist}.tar.gz
     ],[
       AC_MSG_ERROR([file ${gsldist}.tar.gz not found and wget not installed])
     ])
   ])
 ])

 AS_IF([test -e ${gsldist}.tar.gz],
  [
   AS_IF([test ! -e ${gsldist}],[
    AC_MSG_NOTICE([uncompressing ${gsldist}.tar.gz])
    tar xzf ${gsldist}.tar.gz
   ])
   AS_IF([test ! -e ${gsldist}/.libs/libgsl.a],[
    AC_MSG_NOTICE([configuring ${gsldist}])
    cd ${gsldist} 
    ./configure --prefix=${prefix} --host=${host}
    AC_MSG_NOTICE([compiling ${gsldist}])
    make
    cd ..
   ],[
    AC_MSG_NOTICE([using already-compiled GSL library ${gsldist}/.libs/libgsl.a])
   ])
   gsl_version="1.16 (embedded)"
   AC_DEFINE([HAVE_GLFIXED_TABLE], [1])
   WASORALIBS="$WASORALIBS ../${gsldist}/.libs/libgsl.a ../${gsldist}/cblas/.libs/libgslcblas.a"
   STANDALONELIBS="$STANDALONELIBS ../../${gsldist}/.libs/libgsl.a ../../${gsldist}/cblas/.libs/libgslcblas.a"
   AC_SUBST([WASORALIBS], [${WASORALIBS}])
   AC_SUBST([STANDALONELIBS], [${STANDALONELIBS}])
   CFLAGS="$CFLAGS -I ../${gsldist} -I ../../${gsldist}"
  ],[
#   AS_BOX([GNU GSL not found])
   AC_MSG_NOTICE([GNU Scientific Library development files not found.
Please install the libgsl-dev (or equivalent) package or download the
${gsldist}.tar.gz file into the current directory and then try to configure again.
You can have me download it automatically by using the --enable-download-gsl
option in the command line.])
   AC_MSG_ERROR([cannot find GNU Scientific Library])
  ]
 )
,
 [
  AS_IF([test -e ${prefix}/bin/gsl-config ],
         [gsl_version=`${prefix}/bin/gsl-config --version`],
        [test "x`which gsl-config`" != "x" ],
         [gsl_version=`gsl-config --version`],
        [gsl_version=unknown])
  AC_CHECK_TYPE([gsl_integration_glfixed_table],
    [AC_DEFINE([HAVE_GLFIXED_TABLE], [1],
      [Define if the GSL library has type integration_glfixed_table (some old versions do not.)])],
    [],
    [#include <gsl/gsl_integration.h>])
  ])
 ])
])



AC_DEFUN([WASORA_CHECK_IDA],[

# the default is read from the macro argument, but the help string
# does not expand variables so it always states that it is "check"
ida_default=m4_default([$1],[check])

AC_ARG_WITH([ida],
  [AS_HELP_STRING([--with-ida],
    [support differential equations @<:@default=check@:>@])],
  [],
  [with_ida=${ida_default}])

AS_IF([test "x$with_ida" != xno],[
   AC_CHECK_HEADERS([sundials/sundials_types.h ida/ida.h], [],
    [AS_IF([test "x$with_ida" != xcheck],
       [AC_MSG_FAILURE([--with-ida was given, but test for ida headers failed])],
       [AC_MSG_WARN([sundials ida headers not found])])
    ])
   AC_CHECK_LIB([sundials_ida], [IDAInit],,
    [AS_IF([test "x$with_ida" != xcheck],
       [AC_MSG_FAILURE([--with-ida was given, but test for ida libray failed])],
       [AC_MSG_WARN([sundials ida library (libsundials-ida) not found])])
    ])
   AC_CHECK_HEADER([nvector/nvector_serial.h], [],
    [AS_IF([test "x$with_ida" != xcheck],
       [AC_MSG_FAILURE([--with-ida was given, but test for sundials nvecserial headers failed])],
       [AC_MSG_WARN([sundials ida headers not found])])
    ])
   AC_CHECK_LIB([sundials_nvecserial], [N_VNew_Serial],,
    [AS_IF([test "x$with_ida" != xcheck],
       [AC_MSG_FAILURE([--with-ida was given, but test for sundials nvecserial libray failed])],
       [AC_MSG_WARN([sundials nvecserial library (libsundials-nvecserial) not found])])
    ])
  ])

# check if we have everything
AS_IF([test "x$ac_cv_lib_sundials_ida_IDAInit" = xyes -a "x$ac_cv_header_sundials_sundials_types_h" = xyes -a "x$ac_cv_lib_sundials_nvecserial_N_VNew_Serial" = xyes -a "x$ac_cv_header_nvector_nvector_serial_h" = xyes ],
  [
   ida=1
   #AS_IF([test "x`which sundials-config`" != x],
           #[ida_include=`sundials-config -lc -mida -ts | grep I | awk '{print [$]1}' | cut -c3-`
            #ida_version=`cat ${ida_include}/sundials/sundials_config.h | grep VERSION | awk '{print [$]3}' | sed s/\"//g`],
         #[ida_version=unknown])
   ida_version=unknown
   AC_DEFINE(HAVE_IDA)
  ],[
   ida=0
  ])
])


AC_DEFUN([WASORA_CHECK_READLINE],[

# the default is read from the macro argument, but the help string
# does not expand variables so it always states that it is "check"
readline_default=m4_default([$1],[check])

AC_ARG_WITH([readline],
  [AS_HELP_STRING([--with-readline],
    [support interactive debug mode @<:@default=check@:>@])],
  [],
  [with_readline=${readline_default}])

AS_IF([test "x$with_readline" != xno],[
   AC_CHECK_HEADER([readline/readline.h], [],
    [AS_IF([test "x$with_readline" != xcheck],
       [AC_MSG_FAILURE([--with-readline was given, but test for readline headers failed])],
       [AC_MSG_WARN([GNU readline headers (libreadline-dev) not found.])])
    ])
   AC_CHECK_LIB([readline], [readline],,
    [AS_IF([test "x$with_readline" != xcheck],
       [AC_MSG_FAILURE([--with-readline was given, but test for readline libray failed])],
       [AC_MSG_WARN([GNU readline library (libreadline) not found.])])
    ])
  ])

# check if we have everything
AS_IF([test "x$ac_cv_lib_readline_readline" = xyes -a "x$ac_cv_header_readline_readline_h" = xyes ],
  [
   readline=1
   AS_IF([test -e /usr/include/readline/readline.h],
           [readline_version=`cat /usr/include/readline/readline.h | grep VERSION_MAJOR | awk '{print [$]3}'`.`cat /usr/include/readline/readline.h | grep VERSION_MINOR | awk '{print [$]3}'`],
         [test -e /usr/local/include/readline/readline.h],
           [readline_version=`cat /usr/local/include/readline/readline.h | grep VERSION_MAJOR | awk '{print [$]3}'`.`cat /usr/local/include/readline/readline.h | grep VERSION_MINOR | awk '{print [$]3}'`],
         [readline_version=unknown])
   AC_DEFINE(HAVE_READLINE)
  ],[
   readline=0
  ])
])


AC_DEFUN([WASORA_CHECK_WASORA_DIR],[

AC_ARG_WITH([wasora],
  [AS_HELP_STRING([--with-wasora=DIR],
    [look for src/wasora.h in the specified directory])])

# buscamos wasora.h
AC_MSG_CHECKING([for wasora directory])
# la opcion --with-wasora sobre-escribe la variable de entorno
AS_IF([test -n "$with_wasora"], [WASORA_DIR=$with_wasora])

# si todavia esta vacia, usamos la embebida
AS_IF([test -n  "${WASORA_DIR}"],
 [AS_IF([test -e "${WASORA_DIR}/src/wasora.h"],
  [AC_MSG_RESULT([yes])
   #AC_CHECK_HEADER([wasora.h],,AC_MSG_ERROR([It seems that the wasora.h header does not work.]))
   WASORA_INCLUDE="${WASORA_DIR}/src -I${WASORA_DIR}/gsl-1.16"
  ],[
   AC_MSG_RESULT([cannot find ${WASORA_DIR}/src/wasora.h])
   AC_MSG_ERROR([
The directory should be set to the wasora directory that contains src, i.e. either
export WASORA_DIR=$HOME/wasora
or
./configure --with-wasora=$HOME/wasora
and not the src subdirectory itself. Please check and configure again.])
 ])],[
  AC_MSG_RESULT([no, proceeding with embedded wasora source])
  WASORA_INCLUDE="./wasora"
 ])
AC_SUBST([WASORA_INCLUDE], [-I${WASORA_INCLUDE}])
])



AC_DEFUN([WASORA_HOST_VERSION_H],[
# nos fabricamos el version-conf.h
VERSIONH=src/$1version-conf.h
AC_MSG_NOTICE([creating $VERSIONH])
cat << EOF > $VERSIONH
#define COMPILATION_ARCH     "${host_os} ${host_cpu}"
#define COMPILER_VERSION     "`$CC --version | head -n1`"
#define COMPILER_CFLAGS      "$CFLAGS"
EOF

dnl ejecutamos version.sh (sin el $1, porque en un plugin no
dnl queremos meternos en src/wasora sino solo en src)
cd src
./version.sh
cd ..
])

AC_DEFUN([WASORA_PLUGIN_VERSION_H],[
# nos fabricamos el version.h del plugin
PLUGINVERSIONH=src/version-conf.h

AC_MSG_NOTICE([creating $PLUGINVERSIONH])

cat << EOF > $PLUGINVERSIONH
#define COMPILATION_ARCH     "${host_os} ${host_cpu}"
#define CCOMPILER_VERSION    "`$CC --version | head -n1`"
#define CCOMPILER_FLAGS      "$CFLAGS"
EOF

cd src
./version.sh
cd ..
])

AC_DEFUN([WASORA_RESUME_LIBS],[
echo "  GSL library (required): yes, version ${gsl_version}"
echo

if [[ $ida -eq 1 ]]; then
  ida_message="yes, version ${ida_version}"
  ida_not=""
else
  ida_message="no"
  ida_not=" NOT"
fi
echo "  IDA library (optional): ${ida_message}"
echo "    differential-algebraic systems will${ida_not} be solved"
echo

if [[ $readline -eq 1 ]]; then
  readline_message="yes, version ${readline_version}"
  readline_not=""
else
  readline_message="no"
  readline_not=" NOT"
fi
echo "  Readline library (opt): ${readline_message}"
echo "    run-time debugging-like capabilities will${readline_not} be provided"
echo
])


AC_DEFUN([WASORA_FINAL_WARN],[
  all_libs=m4_default([$1],[0])
  make_shell=m4_default([$2],[0])

# si no esta gardel con todos los guitarristas, warning
AS_IF([ test ${all_libs} -eq 0],[
  AC_MSG_WARN([there is at least one optional library missing.])
  echo "If this was not the desired result, check config.log for clues."
  echo
])
AS_IF([ test ${make_shell} -eq 0],[
  echo "Now proceed to compile with 'make'"
  echo
],[
  echo "Now proceed to compile with 'make' (or 'make SHELL=/bin/bash' on Debian)"
  echo
])
])
