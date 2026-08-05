#!/bin/sh

vers=$1
if [ "$vers" = "" ] ; then
  echo cannot change to a non-version
  exit 1
fi

OSTYPE=`uname -s`
if [ "${OSTYPE}" = "Darwin" ] || [ "${OSTYPE}" = "FreeBSD" ] ; then
  EXTN='.x'
else
  EXTN=''
fi

N1=`echo $vers | cut -f1 -d.`
N2=`echo $vers | cut -f2 -d.`
N3=`echo $vers | cut -f3 -d.`
N4=`echo $vers | cut -f4 -d.`

if [ "$N4" = "" ] ; then
  T1=`echo $N3 |  tr 'a' ' '`
  if [ "$T1" = "$N3" ] ; then
    T1=`echo $N3 |  tr 'b' ' '`
    if [ "$T1" = "$N3" ] ; then
      S=' '
    else
      S='b'
    fi
  else
    S='a'
  fi
  if [ "$S" = ' ' ] ; then
    N4=
  else
    N4=`echo $N3 | cut -f2- -d$S`
    N3=`echo $N3 | cut -f1 -d$S`
    N4=',0x'$S$N4
  fi
fi
OPV=$N1\,$N2\,$N3$N4

 echo vers = $vers - OPV  = $OPV

for FILE in ./glm.rc ./glm+.rc ; do
  OFV=`grep FILEVERSION ${FILE} | sed 's/^[ \t]*//' | cut -f2 -d\ `
  echo OFV  = $OFV - OPV  = $OPV

  if [ "$OPV" != "$OFV" ] ; then # new version number
    echo sed -e "s/${OFV}/${OPV}/" -i${EXTN} ${FILE}
         sed -e "s/${OFV}/${OPV}/" -i${EXTN} ${FILE}

    if [ "${OSTYPE}" = "Darwin" ] || [ "${OSTYPE}" = "FreeBSD" ] ; then
      /bin/rm ${FILE}${EXTN}
    fi

    OFV=`grep FileVersion ${FILE} | sed 's/^[ \t]*//' | cut -f3 -d\ | tr -d '\"'`
   echo OFV2 = \'$OFV\' - vers = \'$vers\'

    if [ "$vers" != "$OFV" ] ; then # new version number
      echo sed -e "s/${OFV}/${vers}/" -i${EXTN} ${FILE}
           sed -e "s/${OFV}/${vers}/" -i${EXTN} ${FILE}

      if [ "${OSTYPE}" = "Darwin" ] || [ "${OSTYPE}" = "FreeBSD" ] ; then
        /bin/rm ${FILE}${EXTN}
      fi
    fi
  else
    echo no change to version number in ${FILE}
  fi
done

exit 0
