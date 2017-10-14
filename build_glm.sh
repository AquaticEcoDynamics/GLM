#!/bin/bash

if [ "$GLM_CONFIGURED" != "true" ] ; then
  . ./GLM_CONFIG

  export OSTYPE=`uname -s`

  if [ "$FORTRAN_COMPILER" = "IFORT" ] ; then
    if [ -d /opt/intel/bin ] ; then
      . /opt/intel/bin/compilervars.sh intel64
    fi
    which ifort >& /dev/null
    if [ $? != 0 ] ; then
      echo ifort compiler requested, but not found
      exit 1
    fi
  fi

  if [ "$FORTRAN_COMPILER" = "IFORT" ] ; then
    export PATH="/opt/intel/bin:$PATH"
    export FC=ifort
    export NETCDFHOME=/opt/intel
  else
    export FC=gfortran
    export NETCDFHOME=/usr
    if [ "$OSTYPE" == "Darwin" ] ; then
      export NETCDFHOME=/opt/local
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
else
  export CURDIR=`pwd`
fi

if [ 1 = 0 ] ; then
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

fi

if [ -f ${CURDIR}/src/glm ] ; then
  /bin/rm ${CURDIR}/src/glm
fi
cd ${CURDIR}
make || exit 1

cd ${CURDIR}

if [ "$OSTYPE" = "Linux" ] ; then
  VERSION=`grep GLM_VERSION src/glm.h | cut -f2 -d\"`
  echo glm version $VERSION
  if [ -f /etc/debian_version ] ; then
    VERSDEB=`head -1 debian/changelog | cut -f2 -d\( | cut -f1 -d-`
    echo debian version $VERSDEB
    if [ "$VERSION" != "$VERSDEB" ] ; then
      echo updating debian version
      dch --newversion ${VERSION}-0 "new version ${VERSION}"
    fi

    fakeroot make -f debian/rules binary || exit 1
  fi

  cd ${CURDIR}/win
  ${CURDIR}/vers.sh $VERSION
  cd ${CURDIR}/win-dll
  ${CURDIR}/vers.sh $VERSION

  cd ${CURDIR}
  if [ ! -d bin/ubuntu/$(lsb_release -rs) ] ; then
    mkdir -p bin/ubuntu/$(lsb_release -rs)/
  fi
  mv ../glm*.deb bin/ubuntu/$(lsb_release -rs)/
fi
if [ "$OSTYPE" = "Darwin" ] ; then
  if [ ! -d ${CURDIR}/bin/macos ] ; then
     mkdir -p ${CURDIR}/bin/macos
  fi
  cd ${CURDIR}/macos
  /bin/bash macpkg.sh ${HOMEBREW}
  mv ${CURDIR}/macos/glm_*.zip ${CURDIR}/bin/macos/
fi

exit 0
