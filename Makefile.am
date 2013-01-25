## $Id$
#
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti (from a template by Andreas Waechter)

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = doxydoc/doxygen.conf

# Subdirectories ########################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = Creme/src/Base
#$(subdirs)

# Extra targets #########################################################

test: all
	cd Creme && $(MAKE) test

unitTest: test

#tests: all
#	for dir in $(subdirs); do \
#	  if test -r $$dir/test/Makefile; then \
#	    (cd $$dir; $(MAKE) test) \
#	  fi; \
#	done

#unitTests: tests

doxydoc:
	cd $(srcdir) && doxygen doxydoc/doxygen.conf

.PHONY: test unitTest doxydoc 
#.PHONY: test unitTest tests unitTests doxydoc cuter

# Maintainer ############################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = coin_subdirs.txt

include BuildTools/Makemain.inc
