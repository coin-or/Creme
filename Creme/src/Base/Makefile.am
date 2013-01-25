## $Id$
#
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti (from a template by Andreas Waechter)

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libCreme.la
bin_PROGRAMS = creme

libCreme_la_SOURCES = 	\
	chooseblock.c	cmdline.c init.c	isfeas.c			linopt.c \
	locsrch.c 	lpio.c	  misc.c	move.c 			readbz.c	rtr.c	writelp.c

creme_SOURCES = main.c

creme_LDADD = libCreme.la

creme_DEPENDENCIES = libCreme.la

libCreme_la_LDFLAGS = $(LT_LDFLAGS) -lm -lbz2

includecremedir = $(includedir)/coin
includecreme_HEADERS = chooseblock.h cmdline.h 	init.h isfeas.h 	linopt.h 	locsrch.h \
	lpio.h misc.h move.h rtr.h 	sparse.h

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ../../creme.pc

DISTCLEANFILES =
