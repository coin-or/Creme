## $Id$
#
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti (from a template by Andreas Waechter)

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Crème],[0.0trunk],[http://projects.coin-or.org/Creme/newticket])

AC_COPYRIGHT([
Copyright 2011 Pietro Belotti.
All Rights Reserved.
This file is part of the open source package Creme which is distributed
under the Eclipse Public License.])

AC_REVISION(0.9)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

AC_COIN_CHECK_GNU_BZLIB

#############################################################################
#              Create the libtool script for the subdirectories             #
#############################################################################

AC_COIN_CREATE_LIBTOOL

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_PACKAGEDIR(Creme)

AC_COIN_HAS_PKGCONFIG

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile doxydoc/doxygen.conf])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
