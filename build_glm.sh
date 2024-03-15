#!/bin/sh

if [ "$GLM_CONFIGURED" != "true" ] ; then
  . ./GLM_CONFIG
fi

case `uname` in
  "Darwin"|"Linux"|"FreeBSD")
    export OSTYPE=`uname -s`
    ;;
  MINGW*)
    export OSTYPE="Msys"
    ;;
esac

if [ "$OSTYPE" = "FreeBSD" ] ; then
  export FC=flang
  export CC=clang
  export MAKE=gmake
else
  export FC=gfortran
  export CC=gcc
  export MAKE=make
fi

ARGS=""
while [ $# -gt 0 ] ; do
  ARGS="$ARGS $1"
  case $1 in
    --debug)
      export DEBUG=true
      ;;
    --mdebug)
      export MDEBUG=true
      ;;
    --fence)
      export FENCE=true
      ;;
    --fabm)
      export FABM=true
      ;;
    --gfort)
      export FC=gfortran
      ;;
    --ifort)
      export FC=ifort
      ;;
    --flang)
      export FC=flang
      ;;
    --no-gui)
      export WITH_PLOTS=false
      export WITH_XPLOTS=false
      ;;
    *)
      ;;
  esac
  shift
done

if [ "$OSTYPE" = "Darwin" ] ; then
  if [ "$HOMEBREW" = "" ] ; then
    brew -v > /dev/null 2>&1
    if [ $? != 0 ] ; then
      which port > /dev/null 2>&1
      if [ $? != 0 ] ; then
        echo no ports and no brew
      else
        export MACPORTS=true
      fi
    else
      export HOMEBREW=true
    fi
  fi
fi

# see if FC is defined, if not look for gfortran at least v8
if [ "$FC" = "" ] ; then
  gfortran -v > /dev/null 2>&1
  if [ $? != 0 ] ; then
    export FC=ifort
  else
    VERS=`gfortran -dumpversion | cut -d\. -f1`
    if [ $VERS -ge 8 ] ; then
      export FC=gfortran
    else
      gfortran-8 -v > /dev/null 2>&1
      if [ $? != 0 ] ; then
        export FC=ifort
      else
        export gfortran-8
      fi
    fi
  fi
fi

if [ "$FC" = "ifort" ] ; then
  if [ "$OSTYPE" = "Linux" ] ; then
    start_sh="$(ps -p "$$" -o  command= | awk '{print $1}')"
    # ifort config scripts wont work with /bin/sh
    # so we restart using bash
    if [ "$start_sh" = "/bin/sh" ] ;  then
       echo Restart using bash because ifort cant use /bin/sh
       /bin/bash $0 $ARGS
       exit $?
    fi
  fi

  if [ -x /opt/intel/setvars.sh ] ; then
     . /opt/intel/setvars.sh
  elif [ -d /opt/intel/oneapi ] ; then
     . /opt/intel/oneapi/setvars.sh
  elif [ -d /opt/intel/bin ] ; then
     . /opt/intel/bin/compilervars.sh intel64
  fi
  which ifort > /dev/null 2>&1
  if [ $? != 0 ] ; then
     echo ifort compiler requested, but not found
     exit 1
  fi
fi

export F77=$FC
export F90=$FC
export F95=$FC

export MPI=OPENMPI

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
  if [ ! -d $FABMDIR ] ; then
    echo "FABM directory not found"
    export FABM=false
  else
    which cmake > /dev/null 2>&1
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
  export FFLAGS+=-fPIC
  if [ "${USE_DL}" = "true" ] ; then
    cmake ${FABMDIR}/src -DBUILD_SHARED_LIBS=1 || exit 1
  else
    cmake ${FABMDIR}/src || exit 1
  fi
  ${MAKE} || exit 1
fi

if [ "${AED2}" = "true" ] ; then
  cd "${AED2DIR}"
  ${MAKE} || exit 1
  cd ..
  if [ "${AED2PLS}" != "" ] ; then
    if [ -d "${AED2PLS}" ] ; then
      cd "${AED2PLS}"
      ${MAKE} || exit 1
      cd ..
    fi
  fi
fi

if [ "${AED}" = "true" ] ; then
  echo "build libaed-water"
  cd "${CURDIR}/../libaed-water"
  ${MAKE} || exit 1
  DAEDWATDIR=`pwd`
  if [ -d "${CURDIR}/../libaed-benthic" ] ; then
    echo build libaed-benthic
    cd "${CURDIR}/../libaed-benthic"
    ${MAKE} || exit 1
    DAEDBENDIR=`pwd`
  fi
  if [ -d "${CURDIR}/../libaed-demo" ] ; then
    echo build libaed-demo
    cd "${CURDIR}/../libaed-demo"
    ${MAKE} || exit 1
    DAEDDMODIR=`pwd`
  fi
  if [ -d "${CURDIR}/../libaed-riparian" ] ; then
    echo build libaed-riparian
    cd "${CURDIR}/../libaed-riparian"
    ${MAKE} || exit 1
    DAEDRIPDIR=`pwd`
  fi
  if [ -d "${CURDIR}/../libaed-light" ] ; then
    echo build libaed-light
    cd "${CURDIR}/../libaed-light"
    ${MAKE} || exit 1
    DAEDLGTDIR=`pwd`
  fi
  if [ -d "${CURDIR}/../libaed-dev" ] ; then
    if [ -d "${CURDIR}/../phreeqcrm" ] ; then
      PHREEQDIR="${CURDIR}/../phreeqcrm"
    fi
    echo build libaed-dev
    cd "${CURDIR}/../libaed-dev"
    ${MAKE} PHREEQDIR=$PHREEQDIR || exit 1
    DAEDDEVDIR=`pwd`
  fi
fi

if [ -d "${UTILDIR}" ] ; then
  echo "making libutil"
  cd "${UTILDIR}"
  ${MAKE} || exit 1
  cd "${CURDIR}/.."
fi

if [ "$OSTYPE" = "FreeBSD" ] ; then
  echo making flang extras
  cd ancillary/freebsd
  ./fetch.sh
  ${MAKE} || exit 1
elif [ "$OSTYPE" = "Msys" ] ; then
  if [ ! -d ancillary/windows/msys ] ; then
    echo making windows ancillary extras
    cd ancillary/windows/Sources
    ./build_all.sh || exit 1
  fi
fi

if [ "$WITH_PLOTS" = "true" ] ; then
  echo "making libplot"
  cd "${PLOTDIR}"
  ${MAKE} || exit 1
fi

cd "${CURDIR}"
if [ -f obj/aed_external.o ] ; then
  /bin/rm obj/aed_external.o
fi

# Update versions in resource files
VERSION=`grep GLM_VERSION src/glm.h | cut -f2 -d\"`
cd "${CURDIR}/win"
${CURDIR}/vers.sh $VERSION
#cd ${CURDIR}/win-dll
#${CURDIR}/vers.sh $VERSION
cd "${CURDIR}"

${MAKE} AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR || exit 1
if [ "${DAEDDEVDIR}" != "" ] ; then
  if [ -d "${DAEDDEVDIR}" ] ; then
    echo now build plus version
    /bin/rm obj/aed_external.o
    ${MAKE} glm+ AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR AEDRIPDIR=$DAEDRIPDIR AEDLGTDIR=$DAEDLGTDIR AEDDEVDIR=$DAEDDEVDIR PHREEQDIR=$PHREEQDIR || exit 1
  fi
fi

exit 0
