#!/bin/bash

vers=$1
if [ "$vers" = "" ] ; then
  echo cannot change to a non-version
  exit 1
fi

vs2008 () {
 for f in `find . -name glm.vcproj` ; do
     grep -w Version $f
 done
}

for f in `find . -name glm.vcproj` ; do
  for i in `vs2008 | sort -u | cut -f2 -d\"` ; do
    COUNT=`grep -w Version $f | grep $i | wc -l`
    if [ $COUNT = 4 ] ; then
      echo $i is the version for $f
      if [ $i != $vers ] ; then
        echo sed -i "s/Version=\"${i}\"/Version=\"${vers}\"/" $f
        sed -i "s/Version=\"${i}\"/Version=\"${vers}\"/" $f
      else
        echo no change to version number for $f
      fi
    fi
  done
done

for f in `find . -name glm.vcxproj` ; do
  i=`grep -w Version $f | sort -u | cut -f2 -d\> | cut -f1 -d\<`
  echo $i is the version for $f
  if [ $i != $vers ] ; then
    echo sed -i "s!<Version>${i}</Version>!<Version>${vers}</Version>!" $f
    sed -i "s!<Version>${i}</Version>!<Version>${vers}</Version>!" $f
  else
    echo no change to version number for $f
  fi
done

for f in `find . -name glm-\*.vfproj` ; do
  COUNT=`grep -w Version $f | grep VFLinkerTool | wc -l`
  if [ $COUNT = 4 ] ; then
    i=`grep -w Version $f | grep VFLinkerTool | cut -f4 -d\" | cut -f1 -d\" | sort -u`
    echo $i is the version for $f
    if [ $i != $vers ] ; then
      echo sed -i "s/Version=\"${i}\"/Version=\"${vers}\"/" $f
      sed -i "s/Version=\"${i}\"/Version=\"${vers}\"/" $f
    else
      echo no change to version number for $f
    fi
  fi
done

