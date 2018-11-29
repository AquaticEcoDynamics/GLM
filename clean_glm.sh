#!/bin/bash

make distclean

if [ `uname -s` == "Linux" ] ; then
  if [ -f /etc/debian_version ] ; then
    fakeroot make -f debian/rules clean
  fi
fi
if [ `uname -s` == "Darwin" ] ; then
  /bin/rm -rf macos/glm.app
  /bin/rm -rf macos/glm+.app
fi

