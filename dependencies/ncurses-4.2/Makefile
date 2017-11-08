# Generated automatically from Makefile.in by configure.
# $Id: Makefile.in,v 1.15 1998/02/11 12:13:39 tom Exp $
##############################################################################
# Copyright (c) 1998 Free Software Foundation, Inc.                          #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, distribute    #
# with modifications, sublicense, and/or sell copies of the Software, and to #
# permit persons to whom the Software is furnished to do so, subject to the  #
# following conditions:                                                      #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE ABOVE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
# Except as contained in this notice, the name(s) of the above copyright     #
# holders shall not be used in advertising or otherwise to promote the sale, #
# use or other dealings in this Software without prior written               #
# authorization.                                                             #
##############################################################################
#
# Author: Thomas E. Dickey <dickey@clark.net> 1996,1997
#
# Master Makefile for ncurses library.

SHELL = /bin/sh

INSTALL_PREFIX=
CF_MFLAGS = $(MFLAGS) INSTALL_PREFIX="$(INSTALL_PREFIX)"



NCURSES_MAJOR	= 4
NCURSES_MINOR	= 2
NCURSES_PATCH	= 980228

prefix		= /usr/local
exec_prefix	= ${prefix}

bindir		= ${exec_prefix}/bin
datadir		= ${prefix}/share
includedir	= ${prefix}/include
libdir		= ${exec_prefix}/lib
mandir		= ${prefix}/man

DIRS_TO_MAKE = lib objects obj_g

all ::	$(DIRS_TO_MAKE)

$(DIRS_TO_MAKE) :
	mkdir $@

preinstall :
	@ echo ''
	@ echo '** Configuration summary for NCURSES $(NCURSES_MAJOR).$(NCURSES_MINOR) $(NCURSES_PATCH):'
	@ echo ''
	@ echo '      bin directory: '$(bindir)
	@ echo '      lib directory: '$(libdir)
	@ echo '  include directory: '$(includedir)
	@ echo '      man directory: '$(mandir)
	@ echo ' terminfo directory: '$(datadir)/terminfo
	@ echo ''
	@ if test "$(includedir)" != "$(prefix)/include" ; then \
		echo '** Include-directory is not in a standard location'; fi

# Put the common rules here so that we can easily construct the list of
# directories to visit.
all \
clean \
distclean \
mostlyclean \
realclean \
sources \
uninstall \
install ::
	cd man && $(MAKE) $(CF_MFLAGS) $@
	cd include && $(MAKE) $(CF_MFLAGS) $@
	cd ncurses && $(MAKE) $(CF_MFLAGS) $@
	cd progs && $(MAKE) $(CF_MFLAGS) $@
	cd panel && $(MAKE) $(CF_MFLAGS) $@
	cd menu && $(MAKE) $(CF_MFLAGS) $@
	cd form && $(MAKE) $(CF_MFLAGS) $@
	cd misc && $(MAKE) $(CF_MFLAGS) $@
	cd test && $(MAKE) $(CF_MFLAGS) $@

install.libs \
uninstall.libs \
install.includes \
uninstall.includes ::
	cd include && $(MAKE) $(CF_MFLAGS) $@

lint \
lintlib \
install.libs \
uninstall.libs \
install.ncurses \
uninstall.ncurses ::
	cd ncurses && $(MAKE) $(CF_MFLAGS) $@

lint \
lintlib \
install.libs \
uninstall.libs \
install.progs \
uninstall.progs ::
	cd progs && $(MAKE) $(CF_MFLAGS) $@

install.includes \
uninstall.includes \
lint \
lintlib \
install.libs \
uninstall.libs \
install.panel \
uninstall.panel ::
	cd panel && $(MAKE) $(CF_MFLAGS) $@

install.includes \
uninstall.includes \
lint \
lintlib \
install.libs \
uninstall.libs \
install.menu \
uninstall.menu ::
	cd menu && $(MAKE) $(CF_MFLAGS) $@

install.includes \
uninstall.includes \
lint \
lintlib \
install.libs \
uninstall.libs \
install.form \
uninstall.form ::
	cd form && $(MAKE) $(CF_MFLAGS) $@

lint \
lintlib \
install.libs \
uninstall.libs \
install.test \
uninstall.test ::
	cd test && $(MAKE) $(CF_MFLAGS) $@

install.data ::
	cd misc && $(MAKE) $(CF_MFLAGS) $@

install.man ::
	cd man && $(MAKE) $(CF_MFLAGS) $@

distclean ::
	rm -f config.cache config.log config.status Makefile include/ncurses_cfg.h
	rm -f headers.sh headers.sed
	rm -rf $(DIRS_TO_MAKE)
