## $Id$
#
# This file is distributed under the Eclipse Public License.
#
# Author: Pietro Belotti (from a template by Andreas Waechter)

AUTOMAKE_OPTIONS = foreign

SUBDIRS = src/Base

#if ALWAYS_FALSE
  SUBDIRS += test
#endif

#EXTRA_DIST =

EXTRA_DIST = test/grandTest.sh.in

.PHONY: test unitTest

# Stuff to clean up from "make test"

DISTCLEANFILES =

tmpdir = tmpdir

libcreme = libcreme.a

test: unitTest

unitTest: all
	cd test && $(MAKE) test

install-exec-local: install-doc

uninstall-local: uninstall-doc

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = creme.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Creme

# setup addlibs files
# the addlibs_cpp file is created by running pkg-config on the installed creme.pc, if pkg-config available,
#   or by using variables setup in configure
# the addlibs_c file is created by adding @CXXLIBS@ to the content of the addlibs_cpp file
# the addlibs_f file is creating by taking the content of the addlibs_cpp file and removing the FLIBS
#   since pkg-config may have reordered the flags from FLIBS, we take them out one by one
install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"

if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs creme > $(addlibsdir)/creme_addlibs_cpp.txt
else

endif
	addlibs=`cat $(addlibsdir)/creme_addlibs_cpp.txt` ; \
	echo "$$addlibs @CXXLIBS@" > $(addlibsdir)/creme_addlibs_c.txt ; \
	for i in @FLIBS@ coin_dummy ; do \
		addlibs=`echo -n " $$addlibs " | sed -e "s! $$i ! !g"` ; \
	done ; \
	echo "$$addlibs @CXXLIBS@" > $(addlibsdir)/creme_addlibs_f.txt

uninstall-hook:
	rm -f $(addlibsdir)/creme_addlibs_{cpp,c,f}.txt

include BuildTools/Makemain.inc
