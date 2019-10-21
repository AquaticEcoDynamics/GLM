#!/bin/bash
#
# This script is used to bundle the glm binaries and uncommon library dependancies into an app.
#  The bundle rpaths are modified to find the libraries in the bundle.

MOSLINE=`grep 'SOFTWARE LICENSE AGREEMENT FOR ' '/System/Library/CoreServices/Setup Assistant.app/Contents/Resources/en.lproj/OSXSoftwareLicense.rtf'`
# pre Lion :   MOSNAME=`echo ${MOSLINE} | awk -F 'Mac OS X ' '{print $NF}'  | tr -d '\\' | tr ' ' '_'`
# pre Sierra : MOSNAME=`echo ${MOSLINE} | awk -F 'OS X ' '{print $NF}'  | tr -d '\\' | tr ' ' '_'`
MOSNAME=`echo ${MOSLINE} | awk -F 'macOS ' '{print $NF}'  | tr -d '\\' | tr ' ' '_'`

if [ "$1" = "true" ] ; then
   BASEDIR=usr
else
   BASEDIR=opt
fi
if [ "$2" = "" ] ; then
   PKG=glm
else
   PKG=$2
fi

echo '**************************************************************************************'
echo "Building package for ${PKG}"
echo '**************************************************************************************'

/bin/rm -rf ${PKG}.app

VERSION=`grep GLM_VERSION ../src/glm.h | cut -f2 -d\"`
BUILDDATE=`date +%Y%m%d`
buildapp=`echo ${PKG} | tr [A-Z] [a-z]`
BUILDAPP=`echo ${PKG} | tr [a-z] [A-Z]`

mkdir ${PKG}.app
mkdir ${PKG}.app/Contents
mkdir ${PKG}.app/Contents/Resources
mkdir ${PKG}.app/Contents/Resources/English.lproj
mkdir ${PKG}.app/Contents/MacOS
cp Info.plist ${PKG}.app/Contents
sed -i '' -e "s/VERSION/${VERSION}/" ${PKG}.app/Contents/Info.plist
sed -i '' -e "s/BUILDDATE/${BUILDDATE}/" ${PKG}.app/Contents/Info.plist
sed -i '' -e "s/BUILDAPP/${BUILDAPP}/" ${PKG}.app/Contents/Info.plist
sed -i '' -e "s/buildapp/${buildapp}/" ${PKG}.app/Contents/Info.plist

cp PkgInfo    ${PKG}.app/Contents
cp ../${PKG} ${PKG}.app/Contents/MacOS
cp ${PKG}.icns   ${PKG}.app/Contents/Resources
cp glm_files.icns  ${PKG}.app/Contents/Resources
cp InfoPlist.strings ${PKG}.app/Contents/Resources/English.lproj

if [ "${PKG}" = "glm+" ] ; then
  sed -i '' -e "s/GLM2/XLM2/" ${PKG}.app/Contents/Info.plist
  echo -n "APPLXLM2" > ${PKG}.app/Contents/PkgInfo
  sed -i '' -e "s/au.edu.uwa.science.aquatic.glm/au.edu.uwa.science.aquatic.glm+/" ${PKG}.app/Contents/Info.plist
fi

# find_libs path bin
#echo "BASEDIR is ${BASEDIR}" 1>&2
find_libs () {
   #echo "*** find_libs \"$1\" \"$2\"" 1>&2
   L2=`otool -L ${PKG}.app/Contents/MacOS/$2 | grep \/${BASEDIR}\/$1 | cut -d\  -f1 | grep -o '[^/]*$'`
   LIST=""
   while [ "$L2" != "$LIST" ] ; do
      LIST=$L2
      for i in $LIST ; do
         xx=`find \/${BASEDIR} -name $i 2> /dev/null | tail -n 1`
         #echo "Looking for \"$i\" ($xx)" 1>&2
         if [ ! -f ${PKG}.app/Contents/MacOS/$i ] ; then
            if [ "$xx" = "" ] ; then
               /bin/cp -f /${BASEDIR}/local/lib/$i ${PKG}.app/Contents/MacOS
            else
               /bin/cp -f $xx ${PKG}.app/Contents/MacOS
            fi
            if [ $? != 0 ] ; then
               echo " ####### Failed to copy $i" 1>&2
               exit 1
            else
               chmod +w ${PKG}.app/Contents/MacOS/$i
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


LIBS1=`find_libs local ${PKG}`
#LIBS2=`find_libs intel ${PKG}`

if [ "$FC" = "ifort" ] ; then
  PATH2=/opt/intel
  PATH3=
  LIBS2="libifcore.dylib libsvml.dylib libimf.dylib libintlc.dylib"
  if [ "${PKG}" = "glm+" ] ; then
    LIBS2="${LIBS2} libifport.dylib"
  fi
else
  PATH2=/usr/local/
  PATH3=/usr/local/
  LIBS2="libgfortran.5.dylib"
fi

#echo "LIBS1 = $LIBS1"
#echo "LIBS2 = $LIBS2"

# These general libraries
for i in $LIBS1 ; do
   #echo "*** Configuring : $i ***"
   xx=`find /${BASEDIR} -name $i 2> /dev/null | tail -n 1`
   #cp /${BASEDIR}/local/lib/$i ${PKG}.app/Contents/MacOS
   if [ ! -f ${PKG}.app/Contents/MacOS/$i ] ; then

      /bin/cp -f $xx ${PKG}.app/Contents/MacOS
      if [ $? != 0 ] ; then
         echo " ####### Failed to copy(2) $i" 1>&2
         exit 1
      else
         chmod +w ${PKG}.app/Contents/MacOS/$i
      fi
      install_name_tool -id $i ${PKG}.app/Contents/MacOS/$i
      #install_name_tool -change /${BASEDIR}/local/lib/$i '@executable_path/'$i ${PKG}.app/Contents/MacOS/${PKG}
      install_name_tool -change $xx '@executable_path/'$i ${PKG}.app/Contents/MacOS/${PKG}
      #echo '****' install_name_tool -change $xx '@executable_path/'$i ${PKG}.app/Contents/MacOS/${PKG}
      if [ "${BASEDIR}" = "usr" ] ; then
         # This is probably a HOMEBREW setup, so there might be references into the Cellar
         NLST=`otool -L ${PKG}.app/Contents/MacOS/$i | grep \/${BASEDIR}\/local/Cellar | cut -d\  -f1`
         for j in $NLST ; do
            k=`echo $j | grep -o '[^/]*$'`
            install_name_tool -change $j '@executable_path/'$k ${PKG}.app/Contents/MacOS/$i
         done
      fi
   fi
done

# These fortran libraries
for i in $LIBS2 ; do
   if [ ! -f ${PKG}.app/Contents/MacOS/$i ] ; then
      echo cp $PATH2/$i ${PKG}.app/Contents/MacOS
      src=`find $PATH2 -name $i 2> /dev/null | tail -n 1`
      if [ "$src" != "" ] ; then
         /bin/cp -f $src ${PKG}.app/Contents/MacOS
         if [ $? != 0 ] ; then
            echo " ####### Failed to copy(3) $i" 1>&2
            echo " ####### /bin/cp -f $src ${PKG}.app/Contents/MacOS" 1>&2
            exit 1
         else
            chmod +w ${PKG}.app/Contents/MacOS/$i
         fi
      fi
      # These are redundant for intel fortran since it seems intel fortran dylibs exclude the path from their names
      install_name_tool -id $i ${PKG}.app/Contents/MacOS/$i
      install_name_tool -change ${PATH3}$i '@executable_path/'$i ${PKG}.app/Contents/MacOS/${PKG}
   fi
done


# now update these paths in the libraries as well
export LIBS=`\ls ${PKG}.app/Contents/MacOS/ | grep dylib`
export FILES=`\ls ${PKG}.app/Contents/MacOS/`
echo "======================================================================================"
echo "Now checking libs for dependancies"
echo "LIBS=$LIBS"

for file in $FILES ; do
# echo "********** $file"

  for lib in $LIBS ; do
    xx=`otool -L ${PKG}.app/Contents/MacOS/$file | grep $lib | grep -v ${PKG}.app | cut -d\  -f1`
    if [ "$xx" != "" ] ; then
      install_name_tool -change $xx '@executable_path/'$lib ${PKG}.app/Contents/MacOS/$file
    fi
  done
done

zip -r ${PKG}_${VERSION}_macos_${MOSNAME}.zip ${PKG}.app # ${PKG}

exit 0
