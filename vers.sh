#!/bin/bash

vers=$1
if [ "$vers" = "" ] ; then
  echo cannot change to a non-version
  exit 1
fi


for FILE in ./glm.rc ./glm+.rc ; do
  OFV=`grep FILEVERSION ${FILE} | sed 's/^[ \t]*//' | cut -f2 -d\ `
  OPV=`grep PRODUCTVERSION ${FILE} | sed 's/^[ \t]*//' | cut -f2 -d\ `
  ver2=`echo $vers | sed "s/\./,/g"`

  if [ "$ver2" != "$OFV" ] ; then
    echo sed -i "s/${OFV}/${ver2}/" ${FILE}
    sed -i "s/${OFV}/${ver2}/" ${FILE}
    OFV=`grep FileVersion ${FILE} | sed 's/^[ \t]*//' | cut -f3 -d\ `
    echo sed -i "s/${OFV}/${vers}/" ${FILE}
    sed -i "s/${OFV}/${vers}/" ${FILE}
  else
    echo no change to version number in ${FILE}
  fi
done

exit 0
