#!/bin/sh
# $Id: autogen.sh 8435 2012-10-23 17:23:07Z sander $

# barf on errors
set -e

# may be used to force a certain automake-version e.g. 1.7
AMVERS=

# everybody who checks out the CVS wants the maintainer-mode to be enabled
# (should be off for source distributions, this should happen automatically)
DEFAULTCONFOPT="--enable-maintainer-mode"

# check if automake-version was set
if test "x$AMVERS" != x ; then
  echo Warning: explicitly using automake version $AMVERS
  # binaries are called automake-$AMVERS
  AMVERS="-$AMVERS"
fi

## run autotools

echo "--> libtoolize..."
# this script won't rewrite the files if they already exist. This is a
# PITA when you want to upgrade libtool, thus I'm setting --force
libtoolize --force || glibtoolize --force

# prepare everything
echo "--> aclocal..."
aclocal$AMVERS -I m4

# applications should provide a config.h for now
echo "--> autoheader..."
autoheader

echo "--> automake..."
automake$AMVERS --add-missing

echo "--> autoconf..."
autoconf
