#!/bin/bash

. GLM_CONFIG

cd ${CURDIR}
make distclean

if [ `uname -s` == "Linux" ] ; then
  if [ -f /etc/debian_version ] ; then
    cd ${CURDIR}
    fakeroot make -f debian/rules clean
  fi
fi
if [ `uname -s` == "Darwin" ] ; then
  /bin/rm -rf ${CURDIR}/macos/glm.app
  /bin/rm -rf ${CURDIR}/macos/glm+.app
fi

clean_outputs() {
   cd ${CURDIR}/$1
   LIST1=`find . -name WQ\*.csv`
   LIST2=`find . -name output.nc`
   LIST3=`find . -name lake.csv`
   LIST4=`find . -name outlet_\?\?.csv`
   LIST5=`find . -name overflow.csv`
   LIST6=`find . -name stress_dbg.csv`
   LIST=`echo $LIST1 $LIST2 $LIST3 $LIST4 $LIST5 $LIST6`
   #echo $LIST1
   #echo $LIST2
   #echo $LIST3
   #echo $LIST
   if [ "$LIST" != "" ] ; then
     /bin/rm $LIST
   fi

}

if [ -d "Examples" ] ; then
   clean_outputs "Examples"
fi
