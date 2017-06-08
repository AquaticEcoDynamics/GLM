#!/bin/bash
#
# This script is used to bundle the glm binaries and uncommon library dependancies into an app.
#  The bundle rpaths are modified to find the libraries in the bundle.

if [ "$1" = "true" ] ; then
   BASEDIR=usr
else
   BASEDIR=opt
fi

/bin/rm -rf glm.app

VERSION=`grep GLM_VERSION ../src/glm.h | cut -f2 -d\"`
BUILDDATE=`date +%Y%m%d`

mkdir glm.app
mkdir glm.app/Contents
mkdir glm.app/Contents/Resources
mkdir glm.app/Contents/Resources/English.lproj
mkdir glm.app/Contents/MacOS
cp Info.plist glm.app/Contents
sed -i '' -e "s/VERSION/${VERSION}/" glm.app/Contents/Info.plist
sed -i '' -e "s/BUILDDATE/${BUILDDATE}/" glm.app/Contents/Info.plist

cp PkgInfo    glm.app/Contents
cp ../glm glm.app/Contents/MacOS
cp glm.icns   glm.app/Contents/Resources
cp glm_files.icns  glm.app/Contents/Resources
cp InfoPlist.strings glm.app/Contents/Resources/English.lproj

# find_libs path bin
echo "BASEDIR is ${BASEDIR}" 1>&2
find_libs () {
   echo "*** find_libs \"$1\" \"$2\"" 1>&2
   L2=`otool -L glm.app/Contents/MacOS/$2 | grep \/${BASEDIR}\/$1 | cut -d\  -f1 | grep -o '[^/]*$'`
   LIST=""
   while [ "$L2" != "$LIST" ] ; do
      LIST=$L2
      for i in $LIST ; do
         xx=`find \/${BASEDIR} -name $i 2> /dev/null`
         echo "Looking for \"$i\" ($xx)" 1>&2
         if [ ! -f glm.app/Contents/MacOS/$i ] ; then
            if [ "$xx" = "" ] ; then
               cp /${BASEDIR}/local/lib/$i glm.app/Contents/MacOS
            else
               cp $xx glm.app/Contents/MacOS
            fi
            if [ $? != 0 ] ; then
               echo " ####### Failed to copy $i" 1>&2
            else
               chmod +w glm.app/Contents/MacOS/$i
            fi
         fi
         if [ "$xx" = "" ] ; then
            NLST=`otool -L /${BASEDIR}/$1/lib/$i | grep -v $i | grep \/${BASEDIR}\/$1 | cut -d\  -f1 | grep -o '[^/]*$'`
         else
            NLST=`otool -L $xx | grep -v $i | grep \/${BASEDIR}\/$1 | cut -d\  -f1 | grep -o '[^/]*$'`
         fi
         for j in $NLST ; do
            echo $L2 | grep $j > /dev/null 2>&1
            if [ $? != 0 ] ; then
               echo '**** Adding ' $j 1>&2
               L2+=" $j"
            fi
         done
      done
   done
   echo $LIST
}


LIBS1=`find_libs local glm`
#LIBS2=`find_libs intel glm`

if [ "$FORTRAN_COMPILER" = "IFORT" ] ; then
  PATH2=/opt/intel/lib
  PATH3=
  LIBS2="libifcore.dylib libsvml.dylib libimf.dylib libintlc.dylib"
else
  PATH2=/usr/local/lib
  PATH3=/usr/local/lib/
  LIBS2="libgfortran.3.dylib"
fi

echo "LIBS1 = $LIBS1"
echo "LIBS2 = $LIBS2"

# These general libraries
for i in $LIBS1 ; do
   echo "*** Configuring : $i ***"
   xx=`find /${BASEDIR} -name $i 2> /dev/null`
   #cp /${BASEDIR}/local/lib/$i glm.app/Contents/MacOS
   cp $xx glm.app/Contents/MacOS
   install_name_tool -id $i glm.app/Contents/MacOS/$i
   #install_name_tool -change /${BASEDIR}/local/lib/$i '@executable_path/'$i glm.app/Contents/MacOS/glm
   install_name_tool -change $xx '@executable_path/'$i glm.app/Contents/MacOS/glm
   echo '****' install_name_tool -change $xx '@executable_path/'$i glm.app/Contents/MacOS/glm
   if [ "${BASEDIR}" = "usr" ] ; then
      # This is probably a HOMEBREW setup, so there might be references into the Cellar
      NLST=`otool -L glm.app/Contents/MacOS/$i | grep \/${BASEDIR}\/local/Cellar | cut -d\  -f1`
      for j in $NLST ; do
         k=`echo $j | grep -o '[^/]*$'`
         install_name_tool -change $j '@executable_path/'$k glm.app/Contents/MacOS/$i
      done
   fi
done

# These fortran libraries
for i in $LIBS2 ; do
   cp $PATH2/$i glm.app/Contents/MacOS
# These are redundant since it seems intel fortran dylibs exclude the path from their names
   install_name_tool -id $i glm.app/Contents/MacOS/$i
   install_name_tool -change ${PATH3}$i '@executable_path/'$i glm.app/Contents/MacOS/glm
done


# now update these paths in the libraries as well
export LIBS=`\ls glm.app/Contents/MacOS/ | grep dylib`
echo "======================================================================================"
echo "Now checking libs for dependancies"
echo "LIBS=$LIBS"

for file in $LIBS ; do
  echo "********** $file"

  L2=`otool -L glm.app/Contents/MacOS/$file | grep \/${BASEDIR}\/local | cut -d\  -f1`

  for j in $L2 ; do
    lib=`echo $j | grep -o '[^/]*$'`
    #echo "********** $file : $j ($lib)"
    xx=`find /${BASEDIR} -name $lib 2> /dev/null`
    if [ "$xx" != "" ] ; then
      echo ' YYYY' install_name_tool -change $xx '@executable_path/'$lib glm.app/Contents/MacOS/$file
      install_name_tool -change $xx '@executable_path/'$lib glm.app/Contents/MacOS/$file
    fi
  done

  for j in $LIBS2 ; do
    xx=`otool -L glm.app/Contents/MacOS/$file | grep $j`
    if [ "$xx" != "" ] ; then
      echo ' XXXX' install_name_tool -change $j '@executable_path/'$j glm.app/Contents/MacOS/$file
      install_name_tool -change $j '@executable_path/'$j glm.app/Contents/MacOS/$file
    fi
  done
done

# ln -s glm.app/Contents/MacOS/glm glm
zip -r glm_${VERSION}_macos.zip glm.app # glm

exit 0
