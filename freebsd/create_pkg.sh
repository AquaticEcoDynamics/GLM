#!/bin/sh

# This derived from :
#
# http://lastsummer.de/creating-custom-packages-on-freebsd/

# First, a couple of files will have to be set up in this imaginary package:
export PKGNAME=$1
if [ "${PKGNAME}" != "glm" ] && [ "${PKGNAME}" != "glm+" ] ; then
  echo cannot create package called \"${PKGNAME}\"
  exit 1;
fi

STAGEDIR=/tmp/stage
rm -rf ${STAGEDIR}
mkdir -p ${STAGEDIR}

cat >> ${STAGEDIR}/+PRE_DEINSTALL <<EOF
# careful here, this may clobber your system
echo "Resetting root shell"
pw usermod -n root -s /bin/csh
EOF

cat >> ${STAGEDIR}/+POST_INSTALL <<EOF
# careful here, this may clobber your system
echo "Registering root shell"
pw usermod -n root -s /bin/sh
EOF

# There is +PRE_DEINSTALL, +POST_DEINSTALL, +PRE_INSTALL as well as
# +POST_INSTALL. Generally, you want to run post-install and pre-deinstall
# tasks as it makes sure your files are in place when you do the system
# manipulation, but that's just a rule of thumb. Next, the
# master file +MANIFEST is written:

cat >> ${STAGEDIR}/+MANIFEST <<EOF
name: ${PKGNAME}
version: "3.3.0a5"
origin: science/${PKGNAME}
comment: "automates stuff"
desc: "automates tasks which can also be undone later"
maintainer: casper@ambinet.com.au
www: https://aed.science.uwa.edu.au/
prefix: /
EOF

# This will get you going, but what if the package has dependencies that
# should be enforced? Maybe the package is simply a meta-package that holds
# all essential packages so you can design fast customised system setups?
# If the packages are installed on the build system, then pkg-query(8) may
# be used to generate these dependencies:

echo "deps: {" >> ${STAGEDIR}/+MANIFEST
pkg query "  %n: { version: \"%v\", origin: %o }" portlint >> ${STAGEDIR}/+MANIFEST
pkg query "  %n: { version: \"%v\", origin: %o }" poudriere >> ${STAGEDIR}/+MANIFEST
echo "}" >> ${STAGEDIR}/+MANIFEST

# The benefit of pkg-query(8) is that it will format and output the correct
# information so one doesn't have to triple-check for the correct info
# (trust me, I've been there). Now, shipping files inside the package
# requires another file, namely plist:

mkdir -p ${STAGEDIR}/usr/local/bin
cp ../glm ${STAGEDIR}/usr/local/bin
echo "/usr/local/bin/glm" > ${STAGEDIR}/plist

# The last step uses pkg-create(8) to build a package:

pkg create -m ${STAGEDIR}/ -r ${STAGEDIR}/ -p ${STAGEDIR}/plist -o .
if [ "${PKGNAME}" = "glm+" ] ; then
  cp ../glm+ ${STAGEDIR}/usr/local/bin
  echo "/usr/local/bin/glm+" > ${STAGEDIR}/plist
  pkg create -m ${STAGEDIR}/ -r ${STAGEDIR}/ -p ${STAGEDIR}/plist -o .
fi

# From here on out, you'll be able to use the package archive file with
# pkg-add(8). In the next couple of posts, we'll look deeper into
# pkg(8), e.g. how to create a custom package mirror, how to properly
# shield your build system using chroot(8), etc.
