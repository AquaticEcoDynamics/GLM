#!/bin/bash

if [ "$GLM_CONFIGURED" != "true" ] ; then
  . ./GLM_CONFIG
fi

while [ $# -gt 0 ] ; do
  case $1 in
    --debug)
      export DEBUG=true
      ;;
    --fence)
      export FENCE=true
      ;;
    *)
      ;;
  esac
  shift
done

export OSTYPE=`uname -s`

if [ "$FC" = "" ] ; then
  export FC="ifort"
fi

echo $FC
if [ "$FC" = "ifort" ] ; then
   # for fabm   
   FORTRAN_COMPILER="IFORT"

   if [ `uname -m` = "i686" ] ; then
      CPU="ia32"
   else
      CPU="intel64"
   fi

  if [ -d /opt/intel/bin ] ; then
    . /opt/intel/bin/compilervars.sh $CPU
  fi
  which ifort >& /dev/null
  if [ $? != 0 ] ; then
    echo ifort compiler requested, but not found
    exit 1
  fi
  export PATH="/opt/intel/bin:$PATH"
  export NETCDFHOME=/opt/intel
else
   # for fabm
   # if FC is not ifort assume that it is a variant of gfortran
   FORTRAN_COMPILER="GFORTRAN"

  if [ "$OSTYPE" == "Darwin" ] ; then
     if [ "${HOMEBREW}" = "true" ] ; then
       export NETCDFHOME=/usr/local
     else
       export NETCDFHOME=/opt/local
     fi
  else
     export NETCDFHOME=/usr
  fi
fi

export F77=$FC
export F90=$FC
export F95=$FC

export MPI=OPENMPI

export NETCDFINC=$NETCDFHOME/include
export NETCDFINCL=${NETCDFINC}
export NETCDFLIBDIR=$NETCDFHOME/lib
export NETCDFLIB=${NETCDFLIBDIR}
export NETCDFLIBNAME="-lnetcdff -lnetcdf"

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi
if [ "$PLOTDIR" = "" ] ; then
  export PLOTDIR=../libplot
fi
if [ "$UTILDIR" = "" ] ; then
  export UTILDIR=../libutil
fi

if [ "$FABM" = "true" ] ; then
  if [ "$DEBUG" = "true" ] ; then
    export COMPILATION_MODE=debug
  else
    export COMPILATION_MODE=production
  fi

  if [ ! -d $FABMDIR ] ; then
    echo "FABM directory not found"
    export FABM=false
  else
    which cmake >& /dev/null
    if [ $? != 0 ] ; then
      echo "cmake not found - FABM cannot be built"
      export FABM=false
    fi
  fi
  if [ "$FABM" = "false" ] ; then
    echo build with FABM requested but FABM cannot be built
    exit 1
  fi

  export FABMHOST=glm
  cd ${FABMDIR}
  if [ ! -d build ] ; then
    mkdir build
  fi
  cd build
  export EXTRA_FFLAGS+=-fPIC
  if [ "${USE_DL}" = "true" ] ; then
    cmake ${FABMDIR}/src -DBUILD_SHARED_LIBS=1 || exit 1
  else
    cmake ${FABMDIR}/src || exit 1
  fi
  make || exit 1
fi

if [ "${AED2}" = "true" ] ; then
  cd ${AED2DIR}
  make || exit 1
fi

if [ "$WITH_PLOTS" = "true" ] ; then
  cd ${PLOTDIR}
  make || exit 1
fi

cd ${UTILDIR}
make || exit 1

if [ -f ${CURDIR}/src/glm ] ; then
  /bin/rm ${CURDIR}/src/glm
fi
cd ${CURDIR}
make || exit 1

exit 0
