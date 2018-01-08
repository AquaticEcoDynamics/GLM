#!/bin/bash

vers=$1
if [ "$vers" = "" ] ; then
  echo cannot change to a non-version
  exit 1
fi

OSTYPE=`uname -s`
if [ "${OSTYPE}" = "Darwin" ] ; then
  EXTN=' ""'
else
  EXTN=''
fi

for FILE in ./glm.rc ./glm+.rc ; do
  OFV=`grep FILEVERSION ${FILE} | sed 's/^[ \t]*//' | cut -f2 -d\ `
  OPV=`grep PRODUCTVERSION ${FILE} | sed 's/^[ \t]*//' | cut -f2 -d\ `
  ver2=`echo $vers | sed "s/\./,/g"`

  if [ "$ver2" != "$OFV" ] ; then
    echo sed -e "s/${OFV}/${ver2}/" -i${EXTN} ${FILE}
    sed -e "s/${OFV}/${ver2}/" -i${EXTN} ${FILE}
    OFV=`grep FileVersion ${FILE} | sed 's/^[ \t]*//' | cut -f3 -d\ `
    echo sed -e "s/${OFV}/${vers}/" -i${EXTN} ${FILE}
    sed -e "s/${OFV}/${vers}/" -i${EXTN} ${FILE}
  else
    echo no change to version number in ${FILE}
  fi
done

exit 0
