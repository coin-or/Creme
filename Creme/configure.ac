## $Id$
#
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti (from a template by Andreas Waechter)

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Crème],[0.0trunk],[http://projects.coin-or.org/Creme])

AC_COPYRIGHT(
[Copyright 2011-2021 Pietro Belotti.
All Rights Reserved.
This file is part of the open source package Creme which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Base/main.c)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Creme)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Creme)

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

AC_COIN_HAS_PKGCONFIG

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHECK_LIBM(CremeLib)

#AC_SUBST(CREMELIBS)
#CREMELIBS=

########################################################################
##          Create Makefiles and other stuff                          ##
########################################################################

AC_CONFIG_FILES([Makefile
		 src/Base/Makefile
		 test/Makefile
		 test/grandTest.sh
		 creme.pc creme-uninstalled.pc])

AC_CONFIG_HEADER([inc/config_creme.h])

AC_COIN_FINALIZE
