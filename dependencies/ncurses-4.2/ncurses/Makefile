# Generated automatically from Makefile.in by configure.
# $Id: Makefile.in,v 1.46 1998/02/11 12:13:57 tom Exp $
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
# Makefile for ncurses source code.
#
# This makes the following:
#	programs
#	includes
#	libraries (normal/debug/profile/shared)
#
# The variable 'srcdir' refers to the source-distribution, and can be set with
# the configure script by "--srcdir=DIR".
#
# The rules are organized to produce the libraries for the configured models,
# and the programs with the configured default model.

# turn off _all_ suffix rules; we'll generate our own
.SUFFIXES:

SHELL		= /bin/sh
THIS		= Makefile

x		= 

CF_MFLAGS 	= $(MFLAGS)


MODEL		= normal
INSTALL_PREFIX	= 
srcdir		= .
prefix		= /usr/local
exec_prefix	= ${prefix}
libdir		= ${exec_prefix}/lib
includedir	= ${prefix}/include
datadir		= ${prefix}/share

ticdir		= $(datadir)/terminfo

INSTALL		= /usr/bin/install -c
INSTALL_LIB	= /usr/bin/install -c -m 644
INSTALL_DATA	= ${INSTALL} -m 644

AR		= ar
AR_OPTS		= rv
AWK		= mawk
LD		= ld
LN_S		= ln -s

CC		= gcc
CPP		= gcc -E
CFLAGS		= -O2 

INCDIR		= $(srcdir)/../include
CPPFLAGS	= -I../ncurses -I$(srcdir)  -DNDEBUG -I. -I../include -I$(includedir) \
		  -DHAVE_CONFIG_H -DTERMINFO=\"$(ticdir)\"

CCFLAGS		= $(CPPFLAGS) $(CFLAGS)

HOSTCC		= gcc
HOSTCCFLAGS	= -O2  $(CPPFLAGS)
HOSTLDFLAGS	=  

CFLAGS_NORMAL	= $(CCFLAGS)
CFLAGS_DEBUG	= $(CCFLAGS) -g -DTRACE
CFLAGS_PROFILE	= $(CCFLAGS) -pg
CFLAGS_SHARED	= $(CCFLAGS) unknown

CFLAGS_DEFAULT	= $(CFLAGS_NORMAL)

LINK		= $(CC)
LDFLAGS		=   

SHLIB_LIST	= 
MK_SHARED_LIB	= echo unknown
REL_VERSION	= 4.2
ABI_VERSION	= 4

RANLIB		= ranlib

LIBRARIES	=  ../lib/libncurses.a ../lib/libncurses_g.a

LINT		= 
LINT_OPTS	= 
LINT_LIBS	= -lncurses 

FALLBACK_LIST	= 

AUTO_SRC = \
	../include/nomacros.h \
	./comp_captab.c \
	./expanded.c \
	./fallback.c \
	keys.tries \
	./lib_keyname.c \
	./lib_gen.c \
	./codes.c \
	./names.c \
	./unctrl.c

TEST_DEPS	= ../lib/libncurses.a
TEST_ARGS	= -L../lib -lncurses
TEST_LDFLAGS	=  $(TEST_ARGS)    

TEST_PROGS = \
	captoinfo$x \
	hardscroll$x \
	hashmap$x \
	lib_mvcur$x

################################################################################
all ::		$(AUTO_SRC) ../lib $(LIBRARIES)

sources:	$(AUTO_SRC)

$(INSTALL_PREFIX)$(libdir) :
	$(srcdir)/../mkinstalldirs $@

../lib : ; mkdir $@

./fallback.c : $(srcdir)/MKfallback.sh
	sh $(srcdir)/MKfallback.sh $(FALLBACK_LIST) >$@

./lib_gen.c : $(srcdir)/MKlib_gen.sh ../include/curses.h
	sh $(srcdir)/MKlib_gen.sh "$(CPP)" "$(AWK)" <../include/curses.h >$@

../include/nomacros.h : $(srcdir)/MKlib_gen.sh ../include/curses.h
	sh $(srcdir)/MKlib_gen.sh "$(CPP)" "$(AWK)" <../include/curses.h | \
		fgrep undef >$@

keys.tries: make_keys$x $(srcdir)/keys.list 
	./make_keys $(srcdir)/keys.list > keys.tries

make_keys$x : \
		$(srcdir)/make_keys.c \
		./names.c
	$(HOSTCC) -o $@ $(HOSTCCFLAGS) $(srcdir)/make_keys.c $(HOSTLDFLAGS)

make_hash$x : \
		$(srcdir)/comp_hash.c \
		../include/hashsize.h
	$(HOSTCC) -o $@ $(HOSTCCFLAGS) -DMAIN_PROGRAM $(srcdir)/comp_hash.c $(HOSTLDFLAGS)

./expanded.c : $(srcdir)/MKexpanded.sh
	sh $(srcdir)/MKexpanded.sh "$(CPP)" $(CPPFLAGS)

./comp_captab.c: \
		make_hash$x \
		../include/hashsize.h \
		$(srcdir)/MKcaptab.awk
	sh $(srcdir)/MKcaptab.awk $(AWK) $(srcdir)/../include/Caps > $@

./lib_keyname.c: $(srcdir)/keys.list $(srcdir)/MKkeyname.awk
	$(AWK) -f $(srcdir)/MKkeyname.awk $(srcdir)/keys.list > $@

./names.c ./codes.c: $(srcdir)/MKnames.awk
	$(AWK) -f $(srcdir)/MKnames.awk $(srcdir)/../include/Caps
	cat namehdr boolnames boolfnames numnames numfnames strnames strfnames nameftr >./names.c
	cat namehdr boolcodes numcodes strcodes codeftr >./codes.c
	-rm -f namehdr nameftr codeftr boolnames boolfnames boolcodes numnames numfnames numcodes strnames strfnames strcodes

./unctrl.c: $(srcdir)/MKunctrl.awk
	echo | $(AWK) -f $(srcdir)/MKunctrl.awk >$@

tags:
	ctags *.[ch]

TAGS:
	etags *.[ch]

mostlyclean ::
	-rm -f core tags TAGS *~ *.ln *.atac trace
	-rm -f $(TEST_PROGS)

clean :: mostlyclean
	-rm -f $(AUTO_SRC)
	-rm -f make_keys
	-rm -f make_hash

distclean :: clean
	-rm -f Makefile

realclean :: distclean

# These rules are used to allow "make -n" to work on a clean directory-tree
../include/hashsize.h \
../include/parametrized.h \
../include/term.h :
	cd ../include; $(MAKE) $(CF_MFLAGS)

# These rules build test-programs for the modules that have test-drivers
test_progs : $(TEST_PROGS)

captoinfo$x : $(srcdir)/captoinfo.c $(TEST_DEPS)
	 $(CC) -o $@ $(CCFLAGS) -DMAIN $(srcdir)/captoinfo.c $(TEST_LDFLAGS)

hardscroll$x : $(srcdir)/hardscroll.c $(TEST_DEPS)
	 $(CC) -o $@ $(CCFLAGS) -DSCROLLDEBUG $(srcdir)/hardscroll.c $(TEST_LDFLAGS)

hashmap$x : $(srcdir)/hashmap.c $(TEST_DEPS)
	 $(CC) -o $@ $(CCFLAGS) -DHASHDEBUG $(srcdir)/hashmap.c $(TEST_LDFLAGS)

lib_mvcur$x : $(srcdir)/lib_mvcur.c $(TEST_DEPS) \
		../objects/dump_entry.o
	 $(CC) -o $@ $(CCFLAGS) -DNCURSES_TEST -I$(srcdir)/../progs $(srcdir)/lib_mvcur.c ../objects/dump_entry.o $(TEST_LDFLAGS)

../objects/dump_entry.o:
	cd ../progs && $(MAKE) ../objects/dump_entry.o:

###############################################################################
# The remainder of this file is automatically generated during configuration
###############################################################################

# generated by mk-0th.awk

C_SRC = \
	$(srcdir)/hardscroll.c \
	$(srcdir)/hashmap.c \
	$(srcdir)/lib_acs.c \
	$(srcdir)/lib_addch.c \
	$(srcdir)/lib_addstr.c \
	$(srcdir)/lib_baudrate.c \
	$(srcdir)/lib_beep.c \
	$(srcdir)/lib_bkgd.c \
	$(srcdir)/lib_box.c \
	$(srcdir)/lib_chgat.c \
	$(srcdir)/lib_clear.c \
	$(srcdir)/lib_clearok.c \
	$(srcdir)/lib_clrbot.c \
	$(srcdir)/lib_clreol.c \
	$(srcdir)/lib_color.c \
	$(srcdir)/lib_data.c \
	$(srcdir)/lib_delch.c \
	$(srcdir)/lib_delwin.c \
	$(srcdir)/lib_doupdate.c \
	$(srcdir)/lib_echo.c \
	$(srcdir)/lib_endwin.c \
	$(srcdir)/lib_erase.c \
	$(srcdir)/lib_flash.c \
	./lib_gen.c \
	$(srcdir)/lib_getch.c \
	$(srcdir)/lib_getstr.c \
	$(srcdir)/lib_hline.c \
	$(srcdir)/lib_immedok.c \
	$(srcdir)/lib_inchstr.c \
	$(srcdir)/lib_initscr.c \
	$(srcdir)/lib_insch.c \
	$(srcdir)/lib_insdel.c \
	$(srcdir)/lib_insstr.c \
	$(srcdir)/lib_instr.c \
	$(srcdir)/lib_isendwin.c \
	$(srcdir)/lib_leaveok.c \
	$(srcdir)/lib_longname.c \
	$(srcdir)/lib_mouse.c \
	$(srcdir)/lib_move.c \
	$(srcdir)/lib_mvcur.c \
	$(srcdir)/lib_mvwin.c \
	$(srcdir)/lib_napms.c \
	$(srcdir)/lib_newterm.c \
	$(srcdir)/lib_newwin.c \
	$(srcdir)/lib_nl.c \
	$(srcdir)/lib_overlay.c \
	$(srcdir)/lib_pad.c \
	$(srcdir)/lib_print.c \
	$(srcdir)/lib_printw.c \
	$(srcdir)/lib_redrawln.c \
	$(srcdir)/lib_refresh.c \
	$(srcdir)/lib_restart.c \
	$(srcdir)/lib_scanw.c \
	$(srcdir)/lib_screen.c \
	$(srcdir)/lib_scroll.c \
	$(srcdir)/lib_scrollok.c \
	$(srcdir)/lib_scrreg.c \
	$(srcdir)/lib_set_term.c \
	$(srcdir)/lib_slk.c \
	$(srcdir)/lib_slkatrof.c \
	$(srcdir)/lib_slkatron.c \
	$(srcdir)/lib_slkatrset.c \
	$(srcdir)/lib_slkattr.c \
	$(srcdir)/lib_slkclear.c \
	$(srcdir)/lib_slkinit.c \
	$(srcdir)/lib_slklab.c \
	$(srcdir)/lib_slkrefr.c \
	$(srcdir)/lib_slkset.c \
	$(srcdir)/lib_slktouch.c \
	$(srcdir)/lib_touch.c \
	$(srcdir)/lib_traceatr.c \
	$(srcdir)/lib_tracedmp.c \
	$(srcdir)/lib_tracemse.c \
	$(srcdir)/lib_tstp.c \
	$(srcdir)/lib_ungetch.c \
	$(srcdir)/lib_vidattr.c \
	$(srcdir)/lib_vline.c \
	$(srcdir)/lib_wattroff.c \
	$(srcdir)/lib_wattron.c \
	$(srcdir)/lib_window.c \
	$(srcdir)/nc_panel.c \
	$(srcdir)/safe_sprintf.c \
	$(srcdir)/write_entry.c \
	$(srcdir)/sigaction.c \
	$(srcdir)/vsscanf.c \
	$(srcdir)/define_key.c \
	./expanded.c \
	$(srcdir)/keyok.c \
	$(srcdir)/lib_adabind.c \
	$(srcdir)/lib_dft_fgbg.c \
	$(srcdir)/lib_freeall.c \
	$(srcdir)/resizeterm.c \
	$(srcdir)/softscroll.c \
	$(srcdir)/wresize.c \
	$(srcdir)/alloc_entry.c \
	$(srcdir)/captoinfo.c \
	./codes.c \
	./comp_captab.c \
	$(srcdir)/comp_error.c \
	$(srcdir)/comp_expand.c \
	$(srcdir)/comp_hash.c \
	$(srcdir)/comp_parse.c \
	$(srcdir)/comp_scan.c \
	./fallback.c \
	$(srcdir)/lib_cur_term.c \
	$(srcdir)/lib_kernel.c \
	./lib_keyname.c \
	$(srcdir)/lib_options.c \
	$(srcdir)/lib_raw.c \
	$(srcdir)/lib_termcap.c \
	$(srcdir)/lib_ti.c \
	$(srcdir)/lib_tparm.c \
	$(srcdir)/lib_tputs.c \
	$(srcdir)/lib_trace.c \
	$(srcdir)/lib_tracechr.c \
	$(srcdir)/lib_twait.c \
	$(srcdir)/name_match.c \
	./names.c \
	$(srcdir)/parse_entry.c \
	$(srcdir)/read_entry.c \
	$(srcdir)/read_termcap.c \
	$(srcdir)/lib_setup.c \
	$(srcdir)/trace_buf.c \
	$(srcdir)/tries.c \
	./unctrl.c

# Producing llib-lncurses is time-consuming, so there's no direct-dependency for
# it in the lintlib rule.  We'll only remove in the cleanest setup.
clean ::
	rm -f llib-lncurses.*

realclean ::
	rm -f llib-lncurses

llib-lncurses : $(C_SRC)
	cproto -a -l -DLINT $(CPPFLAGS) $(C_SRC) >$@

lintlib :
	$(srcdir)/../misc/makellib ncurses $(CPPFLAGS)
lint :
	$(LINT) $(LINT_OPTS) $(CPPFLAGS) $(C_SRC) $(LINT_LIBS)

# generated by mk-1st.awk

NORMAL_OBJS = \
	../objects/hardscroll.o \
	../objects/hashmap.o \
	../objects/lib_acs.o \
	../objects/lib_addch.o \
	../objects/lib_addstr.o \
	../objects/lib_baudrate.o \
	../objects/lib_beep.o \
	../objects/lib_bkgd.o \
	../objects/lib_box.o \
	../objects/lib_chgat.o \
	../objects/lib_clear.o \
	../objects/lib_clearok.o \
	../objects/lib_clrbot.o \
	../objects/lib_clreol.o \
	../objects/lib_color.o \
	../objects/lib_data.o \
	../objects/lib_delch.o \
	../objects/lib_delwin.o \
	../objects/lib_doupdate.o \
	../objects/lib_echo.o \
	../objects/lib_endwin.o \
	../objects/lib_erase.o \
	../objects/lib_flash.o \
	../objects/lib_gen.o \
	../objects/lib_getch.o \
	../objects/lib_getstr.o \
	../objects/lib_hline.o \
	../objects/lib_immedok.o \
	../objects/lib_inchstr.o \
	../objects/lib_initscr.o \
	../objects/lib_insch.o \
	../objects/lib_insdel.o \
	../objects/lib_insstr.o \
	../objects/lib_instr.o \
	../objects/lib_isendwin.o \
	../objects/lib_leaveok.o \
	../objects/lib_longname.o \
	../objects/lib_mouse.o \
	../objects/lib_move.o \
	../objects/lib_mvcur.o \
	../objects/lib_mvwin.o \
	../objects/lib_napms.o \
	../objects/lib_newterm.o \
	../objects/lib_newwin.o \
	../objects/lib_nl.o \
	../objects/lib_overlay.o \
	../objects/lib_pad.o \
	../objects/lib_print.o \
	../objects/lib_printw.o \
	../objects/lib_redrawln.o \
	../objects/lib_refresh.o \
	../objects/lib_restart.o \
	../objects/lib_scanw.o \
	../objects/lib_screen.o \
	../objects/lib_scroll.o \
	../objects/lib_scrollok.o \
	../objects/lib_scrreg.o \
	../objects/lib_set_term.o \
	../objects/lib_slk.o \
	../objects/lib_slkatrof.o \
	../objects/lib_slkatron.o \
	../objects/lib_slkatrset.o \
	../objects/lib_slkattr.o \
	../objects/lib_slkclear.o \
	../objects/lib_slkinit.o \
	../objects/lib_slklab.o \
	../objects/lib_slkrefr.o \
	../objects/lib_slkset.o \
	../objects/lib_slktouch.o \
	../objects/lib_touch.o \
	../objects/lib_traceatr.o \
	../objects/lib_tracedmp.o \
	../objects/lib_tracemse.o \
	../objects/lib_tstp.o \
	../objects/lib_ungetch.o \
	../objects/lib_vidattr.o \
	../objects/lib_vline.o \
	../objects/lib_wattroff.o \
	../objects/lib_wattron.o \
	../objects/lib_window.o \
	../objects/nc_panel.o \
	../objects/safe_sprintf.o \
	../objects/write_entry.o \
	../objects/sigaction.o \
	../objects/vsscanf.o \
	../objects/define_key.o \
	../objects/expanded.o \
	../objects/keyok.o \
	../objects/lib_adabind.o \
	../objects/lib_dft_fgbg.o \
	../objects/lib_freeall.o \
	../objects/resizeterm.o \
	../objects/softscroll.o \
	../objects/wresize.o \
	../objects/alloc_entry.o \
	../objects/captoinfo.o \
	../objects/codes.o \
	../objects/comp_captab.o \
	../objects/comp_error.o \
	../objects/comp_expand.o \
	../objects/comp_hash.o \
	../objects/comp_parse.o \
	../objects/comp_scan.o \
	../objects/fallback.o \
	../objects/lib_cur_term.o \
	../objects/lib_kernel.o \
	../objects/lib_keyname.o \
	../objects/lib_options.o \
	../objects/lib_raw.o \
	../objects/lib_termcap.o \
	../objects/lib_ti.o \
	../objects/lib_tparm.o \
	../objects/lib_tputs.o \
	../objects/lib_trace.o \
	../objects/lib_tracechr.o \
	../objects/lib_twait.o \
	../objects/name_match.o \
	../objects/names.o \
	../objects/parse_entry.o \
	../objects/read_entry.o \
	../objects/read_termcap.o \
	../objects/lib_setup.o \
	../objects/trace_buf.o \
	../objects/tries.o \
	../objects/unctrl.o

$(NORMAL_OBJS) : ../include/ncurses_cfg.h ./curses.priv.h

../lib/libncurses.a : $(NORMAL_OBJS)
	$(AR) $(AR_OPTS) $@ $?
	$(RANLIB) $@

install \
install.libs \
install.ncurses :: $(INSTALL_PREFIX)$(libdir) ../lib/libncurses.a
	@echo installing ../lib/libncurses.a as $(INSTALL_PREFIX)$(libdir)/libncurses.a 
	$(INSTALL_DATA) ../lib/libncurses.a $(INSTALL_PREFIX)$(libdir)/libncurses.a 
	@echo linking libcurses.a to libncurses.a 
	-@rm -f $(INSTALL_PREFIX)$(libdir)/libcurses.a 
	(cd $(INSTALL_PREFIX)$(libdir) && $(LN_S) libncurses.a libcurses.a)
	$(RANLIB) $(INSTALL_PREFIX)$(libdir)/libncurses.a

uninstall \
uninstall.libs \
uninstall.ncurses ::
	@echo uninstalling $(INSTALL_PREFIX)$(libdir)/libncurses.a 
	-@rm -f $(INSTALL_PREFIX)$(libdir)/libncurses.a
	@echo linking libcurses.a to libncurses.a 
	-@rm -f $(INSTALL_PREFIX)$(libdir)/libcurses.a 

clean ::
	rm -f ../lib/libncurses.a

mostlyclean::
	rm -f $(NORMAL_OBJS)

clean ::
	rm -f $(NORMAL_OBJS)

# generated by mk-2nd.awk


../objects/hardscroll.o :	$(srcdir)/hardscroll.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/hardscroll.c

../objects/hashmap.o :	$(srcdir)/hashmap.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/hashmap.c

../objects/lib_acs.o :	$(srcdir)/lib_acs.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_acs.c

../objects/lib_addch.o :	$(srcdir)/lib_addch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_addch.c

../objects/lib_addstr.o :	$(srcdir)/lib_addstr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_addstr.c

../objects/lib_baudrate.o :	$(srcdir)/lib_baudrate.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_baudrate.c

../objects/lib_beep.o :	$(srcdir)/lib_beep.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_beep.c

../objects/lib_bkgd.o :	$(srcdir)/lib_bkgd.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_bkgd.c

../objects/lib_box.o :	$(srcdir)/lib_box.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_box.c

../objects/lib_chgat.o :	$(srcdir)/lib_chgat.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_chgat.c

../objects/lib_clear.o :	$(srcdir)/lib_clear.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_clear.c

../objects/lib_clearok.o :	$(srcdir)/lib_clearok.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_clearok.c

../objects/lib_clrbot.o :	$(srcdir)/lib_clrbot.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_clrbot.c

../objects/lib_clreol.o :	$(srcdir)/lib_clreol.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_clreol.c

../objects/lib_color.o :	$(srcdir)/lib_color.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_color.c

../objects/lib_data.o :	$(srcdir)/lib_data.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_data.c

../objects/lib_delch.o :	$(srcdir)/lib_delch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_delch.c

../objects/lib_delwin.o :	$(srcdir)/lib_delwin.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_delwin.c

../objects/lib_doupdate.o :	$(srcdir)/lib_doupdate.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_doupdate.c

../objects/lib_echo.o :	$(srcdir)/lib_echo.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_echo.c

../objects/lib_endwin.o :	$(srcdir)/lib_endwin.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_endwin.c

../objects/lib_erase.o :	$(srcdir)/lib_erase.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_erase.c

../objects/lib_flash.o :	$(srcdir)/lib_flash.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_flash.c

../objects/lib_gen.o :	./lib_gen.c \
			../include/curses.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_gen.c

../objects/lib_getch.o :	$(srcdir)/lib_getch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_getch.c

../objects/lib_getstr.o :	$(srcdir)/lib_getstr.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_getstr.c

../objects/lib_hline.o :	$(srcdir)/lib_hline.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_hline.c

../objects/lib_immedok.o :	$(srcdir)/lib_immedok.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_immedok.c

../objects/lib_inchstr.o :	$(srcdir)/lib_inchstr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_inchstr.c

../objects/lib_initscr.o :	$(srcdir)/lib_initscr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_initscr.c

../objects/lib_insch.o :	$(srcdir)/lib_insch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_insch.c

../objects/lib_insdel.o :	$(srcdir)/lib_insdel.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_insdel.c

../objects/lib_insstr.o :	$(srcdir)/lib_insstr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_insstr.c

../objects/lib_instr.o :	$(srcdir)/lib_instr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_instr.c

../objects/lib_isendwin.o :	$(srcdir)/lib_isendwin.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_isendwin.c

../objects/lib_leaveok.o :	$(srcdir)/lib_leaveok.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_leaveok.c

../objects/lib_longname.o :	$(srcdir)/lib_longname.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_longname.c

../objects/lib_mouse.o :	$(srcdir)/lib_mouse.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_mouse.c

../objects/lib_move.o :	$(srcdir)/lib_move.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_move.c

../objects/lib_mvcur.o :	$(srcdir)/lib_mvcur.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_mvcur.c

../objects/lib_mvwin.o :	$(srcdir)/lib_mvwin.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_mvwin.c

../objects/lib_napms.o :	$(srcdir)/lib_napms.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_napms.c

../objects/lib_newterm.o :	$(srcdir)/lib_newterm.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_newterm.c

../objects/lib_newwin.o :	$(srcdir)/lib_newwin.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_newwin.c

../objects/lib_nl.o :	$(srcdir)/lib_nl.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_nl.c

../objects/lib_overlay.o :	$(srcdir)/lib_overlay.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_overlay.c

../objects/lib_pad.o :	$(srcdir)/lib_pad.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_pad.c

../objects/lib_print.o :	$(srcdir)/lib_print.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_print.c

../objects/lib_printw.o :	$(srcdir)/lib_printw.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_printw.c

../objects/lib_redrawln.o :	$(srcdir)/lib_redrawln.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_redrawln.c

../objects/lib_refresh.o :	$(srcdir)/lib_refresh.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_refresh.c

../objects/lib_restart.o :	$(srcdir)/lib_restart.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_restart.c

../objects/lib_scanw.o :	$(srcdir)/lib_scanw.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_scanw.c

../objects/lib_screen.o :	$(srcdir)/lib_screen.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_screen.c

../objects/lib_scroll.o :	$(srcdir)/lib_scroll.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_scroll.c

../objects/lib_scrollok.o :	$(srcdir)/lib_scrollok.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_scrollok.c

../objects/lib_scrreg.o :	$(srcdir)/lib_scrreg.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_scrreg.c

../objects/lib_set_term.o :	$(srcdir)/lib_set_term.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_set_term.c

../objects/lib_slk.o :	$(srcdir)/lib_slk.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slk.c

../objects/lib_slkatrof.o :	$(srcdir)/lib_slkatrof.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkatrof.c

../objects/lib_slkatron.o :	$(srcdir)/lib_slkatron.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkatron.c

../objects/lib_slkatrset.o :	$(srcdir)/lib_slkatrset.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkatrset.c

../objects/lib_slkattr.o :	$(srcdir)/lib_slkattr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkattr.c

../objects/lib_slkclear.o :	$(srcdir)/lib_slkclear.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkclear.c

../objects/lib_slkinit.o :	$(srcdir)/lib_slkinit.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkinit.c

../objects/lib_slklab.o :	$(srcdir)/lib_slklab.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slklab.c

../objects/lib_slkrefr.o :	$(srcdir)/lib_slkrefr.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkrefr.c

../objects/lib_slkset.o :	$(srcdir)/lib_slkset.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slkset.c

../objects/lib_slktouch.o :	$(srcdir)/lib_slktouch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_slktouch.c

../objects/lib_touch.o :	$(srcdir)/lib_touch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_touch.c

../objects/lib_traceatr.o :	$(srcdir)/lib_traceatr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_traceatr.c

../objects/lib_tracedmp.o :	$(srcdir)/lib_tracedmp.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tracedmp.c

../objects/lib_tracemse.o :	$(srcdir)/lib_tracemse.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tracemse.c

../objects/lib_tstp.o :	$(srcdir)/lib_tstp.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tstp.c

../objects/lib_ungetch.o :	$(srcdir)/lib_ungetch.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_ungetch.c

../objects/lib_vidattr.o :	$(srcdir)/lib_vidattr.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_vidattr.c

../objects/lib_vline.o :	$(srcdir)/lib_vline.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_vline.c

../objects/lib_wattroff.o :	$(srcdir)/lib_wattroff.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_wattroff.c

../objects/lib_wattron.o :	$(srcdir)/lib_wattron.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_wattron.c

../objects/lib_window.o :	$(srcdir)/lib_window.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_window.c

../objects/nc_panel.o :	$(srcdir)/nc_panel.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/nc_panel.c

../objects/safe_sprintf.o :	$(srcdir)/safe_sprintf.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/safe_sprintf.c

../objects/write_entry.o :	$(srcdir)/write_entry.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/write_entry.c

../objects/sigaction.o :	$(srcdir)/sigaction.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/sigaction.c

../objects/vsscanf.o :	$(srcdir)/vsscanf.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/vsscanf.c

../objects/define_key.o :	$(srcdir)/define_key.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/define_key.c

../objects/expanded.o :	./expanded.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/expanded.c

../objects/keyok.o :	$(srcdir)/keyok.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/keyok.c

../objects/lib_adabind.o :	$(srcdir)/lib_adabind.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_adabind.c

../objects/lib_dft_fgbg.o :	$(srcdir)/lib_dft_fgbg.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_dft_fgbg.c

../objects/lib_freeall.o :	$(srcdir)/lib_freeall.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_freeall.c

../objects/resizeterm.o :	$(srcdir)/resizeterm.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/resizeterm.c

../objects/softscroll.o :	$(srcdir)/softscroll.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/softscroll.c

../objects/wresize.o :	$(srcdir)/wresize.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/wresize.c

../objects/alloc_entry.o :	$(srcdir)/alloc_entry.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/alloc_entry.c

../objects/captoinfo.o :	$(srcdir)/captoinfo.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/captoinfo.c

../objects/codes.o :	./codes.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/codes.c

../objects/comp_captab.o :	./comp_captab.c \
			../include/term.h \
			../include/hashsize.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_captab.c

../objects/comp_error.o :	$(srcdir)/comp_error.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_error.c

../objects/comp_expand.o :	$(srcdir)/comp_expand.c \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_expand.c

../objects/comp_hash.o :	$(srcdir)/comp_hash.c \
			../include/term.h \
			../include/hashsize.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_hash.c

../objects/comp_parse.o :	$(srcdir)/comp_parse.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_parse.c

../objects/comp_scan.o :	$(srcdir)/comp_scan.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/comp_scan.c

../objects/fallback.o :	./fallback.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/fallback.c

../objects/lib_cur_term.o :	$(srcdir)/lib_cur_term.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_cur_term.c

../objects/lib_kernel.o :	$(srcdir)/lib_kernel.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_kernel.c

../objects/lib_keyname.o :	./lib_keyname.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_keyname.c

../objects/lib_options.o :	$(srcdir)/lib_options.c \
			../include/term.h \
			keys.tries
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_options.c

../objects/lib_raw.o :	$(srcdir)/lib_raw.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_raw.c

../objects/lib_termcap.o :	$(srcdir)/lib_termcap.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_termcap.c

../objects/lib_ti.o :	$(srcdir)/lib_ti.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_ti.c

../objects/lib_tparm.o :	$(srcdir)/lib_tparm.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tparm.c

../objects/lib_tputs.o :	$(srcdir)/lib_tputs.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tputs.c

../objects/lib_trace.o :	$(srcdir)/lib_trace.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_trace.c

../objects/lib_tracechr.o :	$(srcdir)/lib_tracechr.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_tracechr.c

../objects/lib_twait.o :	$(srcdir)/lib_twait.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_twait.c

../objects/name_match.o :	$(srcdir)/name_match.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/name_match.c

../objects/names.o :	./names.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/names.c

../objects/parse_entry.o :	$(srcdir)/parse_entry.c \
			../include/term.h \
			../include/parametrized.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/parse_entry.c

../objects/read_entry.o :	$(srcdir)/read_entry.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/read_entry.c

../objects/read_termcap.o :	$(srcdir)/read_termcap.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/read_termcap.c

../objects/lib_setup.o :	$(srcdir)/lib_setup.c \
			../include/term.h
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/lib_setup.c

../objects/trace_buf.o :	$(srcdir)/trace_buf.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/trace_buf.c

../objects/tries.o :	$(srcdir)/tries.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/tries.c

../objects/unctrl.o :	./unctrl.c
	cd ../objects; $(CC) $(CFLAGS_NORMAL) -c ../ncurses/unctrl.c


# generated by mk-1st.awk

DEBUG_OBJS = \
	../obj_g/hardscroll.o \
	../obj_g/hashmap.o \
	../obj_g/lib_acs.o \
	../obj_g/lib_addch.o \
	../obj_g/lib_addstr.o \
	../obj_g/lib_baudrate.o \
	../obj_g/lib_beep.o \
	../obj_g/lib_bkgd.o \
	../obj_g/lib_box.o \
	../obj_g/lib_chgat.o \
	../obj_g/lib_clear.o \
	../obj_g/lib_clearok.o \
	../obj_g/lib_clrbot.o \
	../obj_g/lib_clreol.o \
	../obj_g/lib_color.o \
	../obj_g/lib_data.o \
	../obj_g/lib_delch.o \
	../obj_g/lib_delwin.o \
	../obj_g/lib_doupdate.o \
	../obj_g/lib_echo.o \
	../obj_g/lib_endwin.o \
	../obj_g/lib_erase.o \
	../obj_g/lib_flash.o \
	../obj_g/lib_gen.o \
	../obj_g/lib_getch.o \
	../obj_g/lib_getstr.o \
	../obj_g/lib_hline.o \
	../obj_g/lib_immedok.o \
	../obj_g/lib_inchstr.o \
	../obj_g/lib_initscr.o \
	../obj_g/lib_insch.o \
	../obj_g/lib_insdel.o \
	../obj_g/lib_insstr.o \
	../obj_g/lib_instr.o \
	../obj_g/lib_isendwin.o \
	../obj_g/lib_leaveok.o \
	../obj_g/lib_longname.o \
	../obj_g/lib_mouse.o \
	../obj_g/lib_move.o \
	../obj_g/lib_mvcur.o \
	../obj_g/lib_mvwin.o \
	../obj_g/lib_napms.o \
	../obj_g/lib_newterm.o \
	../obj_g/lib_newwin.o \
	../obj_g/lib_nl.o \
	../obj_g/lib_overlay.o \
	../obj_g/lib_pad.o \
	../obj_g/lib_print.o \
	../obj_g/lib_printw.o \
	../obj_g/lib_redrawln.o \
	../obj_g/lib_refresh.o \
	../obj_g/lib_restart.o \
	../obj_g/lib_scanw.o \
	../obj_g/lib_screen.o \
	../obj_g/lib_scroll.o \
	../obj_g/lib_scrollok.o \
	../obj_g/lib_scrreg.o \
	../obj_g/lib_set_term.o \
	../obj_g/lib_slk.o \
	../obj_g/lib_slkatrof.o \
	../obj_g/lib_slkatron.o \
	../obj_g/lib_slkatrset.o \
	../obj_g/lib_slkattr.o \
	../obj_g/lib_slkclear.o \
	../obj_g/lib_slkinit.o \
	../obj_g/lib_slklab.o \
	../obj_g/lib_slkrefr.o \
	../obj_g/lib_slkset.o \
	../obj_g/lib_slktouch.o \
	../obj_g/lib_touch.o \
	../obj_g/lib_traceatr.o \
	../obj_g/lib_tracedmp.o \
	../obj_g/lib_tracemse.o \
	../obj_g/lib_tstp.o \
	../obj_g/lib_ungetch.o \
	../obj_g/lib_vidattr.o \
	../obj_g/lib_vline.o \
	../obj_g/lib_wattroff.o \
	../obj_g/lib_wattron.o \
	../obj_g/lib_window.o \
	../obj_g/nc_panel.o \
	../obj_g/safe_sprintf.o \
	../obj_g/write_entry.o \
	../obj_g/sigaction.o \
	../obj_g/vsscanf.o \
	../obj_g/define_key.o \
	../obj_g/expanded.o \
	../obj_g/keyok.o \
	../obj_g/lib_adabind.o \
	../obj_g/lib_dft_fgbg.o \
	../obj_g/lib_freeall.o \
	../obj_g/resizeterm.o \
	../obj_g/softscroll.o \
	../obj_g/wresize.o \
	../obj_g/alloc_entry.o \
	../obj_g/captoinfo.o \
	../obj_g/codes.o \
	../obj_g/comp_captab.o \
	../obj_g/comp_error.o \
	../obj_g/comp_expand.o \
	../obj_g/comp_hash.o \
	../obj_g/comp_parse.o \
	../obj_g/comp_scan.o \
	../obj_g/fallback.o \
	../obj_g/lib_cur_term.o \
	../obj_g/lib_kernel.o \
	../obj_g/lib_keyname.o \
	../obj_g/lib_options.o \
	../obj_g/lib_raw.o \
	../obj_g/lib_termcap.o \
	../obj_g/lib_ti.o \
	../obj_g/lib_tparm.o \
	../obj_g/lib_tputs.o \
	../obj_g/lib_trace.o \
	../obj_g/lib_tracechr.o \
	../obj_g/lib_twait.o \
	../obj_g/name_match.o \
	../obj_g/names.o \
	../obj_g/parse_entry.o \
	../obj_g/read_entry.o \
	../obj_g/read_termcap.o \
	../obj_g/lib_setup.o \
	../obj_g/trace_buf.o \
	../obj_g/tries.o \
	../obj_g/unctrl.o

$(DEBUG_OBJS) : ../include/ncurses_cfg.h ./curses.priv.h

../lib/libncurses_g.a : $(DEBUG_OBJS)
	$(AR) $(AR_OPTS) $@ $?
	$(RANLIB) $@

install \
install.libs \
install.ncurses :: $(INSTALL_PREFIX)$(libdir) ../lib/libncurses_g.a
	@echo installing ../lib/libncurses_g.a as $(INSTALL_PREFIX)$(libdir)/libncurses_g.a 
	$(INSTALL_DATA) ../lib/libncurses_g.a $(INSTALL_PREFIX)$(libdir)/libncurses_g.a 
	$(RANLIB) $(INSTALL_PREFIX)$(libdir)/libncurses_g.a

uninstall \
uninstall.libs \
uninstall.ncurses ::
	@echo uninstalling $(INSTALL_PREFIX)$(libdir)/libncurses_g.a 
	-@rm -f $(INSTALL_PREFIX)$(libdir)/libncurses_g.a

clean ::
	rm -f ../lib/libncurses_g.a

mostlyclean::
	rm -f $(DEBUG_OBJS)

clean ::
	rm -f $(DEBUG_OBJS)

# generated by mk-2nd.awk


../obj_g/hardscroll.o :	$(srcdir)/hardscroll.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/hardscroll.c

../obj_g/hashmap.o :	$(srcdir)/hashmap.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/hashmap.c

../obj_g/lib_acs.o :	$(srcdir)/lib_acs.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_acs.c

../obj_g/lib_addch.o :	$(srcdir)/lib_addch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_addch.c

../obj_g/lib_addstr.o :	$(srcdir)/lib_addstr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_addstr.c

../obj_g/lib_baudrate.o :	$(srcdir)/lib_baudrate.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_baudrate.c

../obj_g/lib_beep.o :	$(srcdir)/lib_beep.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_beep.c

../obj_g/lib_bkgd.o :	$(srcdir)/lib_bkgd.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_bkgd.c

../obj_g/lib_box.o :	$(srcdir)/lib_box.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_box.c

../obj_g/lib_chgat.o :	$(srcdir)/lib_chgat.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_chgat.c

../obj_g/lib_clear.o :	$(srcdir)/lib_clear.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_clear.c

../obj_g/lib_clearok.o :	$(srcdir)/lib_clearok.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_clearok.c

../obj_g/lib_clrbot.o :	$(srcdir)/lib_clrbot.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_clrbot.c

../obj_g/lib_clreol.o :	$(srcdir)/lib_clreol.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_clreol.c

../obj_g/lib_color.o :	$(srcdir)/lib_color.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_color.c

../obj_g/lib_data.o :	$(srcdir)/lib_data.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_data.c

../obj_g/lib_delch.o :	$(srcdir)/lib_delch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_delch.c

../obj_g/lib_delwin.o :	$(srcdir)/lib_delwin.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_delwin.c

../obj_g/lib_doupdate.o :	$(srcdir)/lib_doupdate.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_doupdate.c

../obj_g/lib_echo.o :	$(srcdir)/lib_echo.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_echo.c

../obj_g/lib_endwin.o :	$(srcdir)/lib_endwin.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_endwin.c

../obj_g/lib_erase.o :	$(srcdir)/lib_erase.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_erase.c

../obj_g/lib_flash.o :	$(srcdir)/lib_flash.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_flash.c

../obj_g/lib_gen.o :	./lib_gen.c \
			../include/curses.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_gen.c

../obj_g/lib_getch.o :	$(srcdir)/lib_getch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_getch.c

../obj_g/lib_getstr.o :	$(srcdir)/lib_getstr.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_getstr.c

../obj_g/lib_hline.o :	$(srcdir)/lib_hline.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_hline.c

../obj_g/lib_immedok.o :	$(srcdir)/lib_immedok.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_immedok.c

../obj_g/lib_inchstr.o :	$(srcdir)/lib_inchstr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_inchstr.c

../obj_g/lib_initscr.o :	$(srcdir)/lib_initscr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_initscr.c

../obj_g/lib_insch.o :	$(srcdir)/lib_insch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_insch.c

../obj_g/lib_insdel.o :	$(srcdir)/lib_insdel.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_insdel.c

../obj_g/lib_insstr.o :	$(srcdir)/lib_insstr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_insstr.c

../obj_g/lib_instr.o :	$(srcdir)/lib_instr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_instr.c

../obj_g/lib_isendwin.o :	$(srcdir)/lib_isendwin.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_isendwin.c

../obj_g/lib_leaveok.o :	$(srcdir)/lib_leaveok.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_leaveok.c

../obj_g/lib_longname.o :	$(srcdir)/lib_longname.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_longname.c

../obj_g/lib_mouse.o :	$(srcdir)/lib_mouse.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_mouse.c

../obj_g/lib_move.o :	$(srcdir)/lib_move.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_move.c

../obj_g/lib_mvcur.o :	$(srcdir)/lib_mvcur.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_mvcur.c

../obj_g/lib_mvwin.o :	$(srcdir)/lib_mvwin.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_mvwin.c

../obj_g/lib_napms.o :	$(srcdir)/lib_napms.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_napms.c

../obj_g/lib_newterm.o :	$(srcdir)/lib_newterm.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_newterm.c

../obj_g/lib_newwin.o :	$(srcdir)/lib_newwin.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_newwin.c

../obj_g/lib_nl.o :	$(srcdir)/lib_nl.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_nl.c

../obj_g/lib_overlay.o :	$(srcdir)/lib_overlay.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_overlay.c

../obj_g/lib_pad.o :	$(srcdir)/lib_pad.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_pad.c

../obj_g/lib_print.o :	$(srcdir)/lib_print.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_print.c

../obj_g/lib_printw.o :	$(srcdir)/lib_printw.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_printw.c

../obj_g/lib_redrawln.o :	$(srcdir)/lib_redrawln.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_redrawln.c

../obj_g/lib_refresh.o :	$(srcdir)/lib_refresh.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_refresh.c

../obj_g/lib_restart.o :	$(srcdir)/lib_restart.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_restart.c

../obj_g/lib_scanw.o :	$(srcdir)/lib_scanw.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_scanw.c

../obj_g/lib_screen.o :	$(srcdir)/lib_screen.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_screen.c

../obj_g/lib_scroll.o :	$(srcdir)/lib_scroll.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_scroll.c

../obj_g/lib_scrollok.o :	$(srcdir)/lib_scrollok.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_scrollok.c

../obj_g/lib_scrreg.o :	$(srcdir)/lib_scrreg.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_scrreg.c

../obj_g/lib_set_term.o :	$(srcdir)/lib_set_term.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_set_term.c

../obj_g/lib_slk.o :	$(srcdir)/lib_slk.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slk.c

../obj_g/lib_slkatrof.o :	$(srcdir)/lib_slkatrof.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkatrof.c

../obj_g/lib_slkatron.o :	$(srcdir)/lib_slkatron.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkatron.c

../obj_g/lib_slkatrset.o :	$(srcdir)/lib_slkatrset.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkatrset.c

../obj_g/lib_slkattr.o :	$(srcdir)/lib_slkattr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkattr.c

../obj_g/lib_slkclear.o :	$(srcdir)/lib_slkclear.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkclear.c

../obj_g/lib_slkinit.o :	$(srcdir)/lib_slkinit.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkinit.c

../obj_g/lib_slklab.o :	$(srcdir)/lib_slklab.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slklab.c

../obj_g/lib_slkrefr.o :	$(srcdir)/lib_slkrefr.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkrefr.c

../obj_g/lib_slkset.o :	$(srcdir)/lib_slkset.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slkset.c

../obj_g/lib_slktouch.o :	$(srcdir)/lib_slktouch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_slktouch.c

../obj_g/lib_touch.o :	$(srcdir)/lib_touch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_touch.c

../obj_g/lib_traceatr.o :	$(srcdir)/lib_traceatr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_traceatr.c

../obj_g/lib_tracedmp.o :	$(srcdir)/lib_tracedmp.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tracedmp.c

../obj_g/lib_tracemse.o :	$(srcdir)/lib_tracemse.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tracemse.c

../obj_g/lib_tstp.o :	$(srcdir)/lib_tstp.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tstp.c

../obj_g/lib_ungetch.o :	$(srcdir)/lib_ungetch.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_ungetch.c

../obj_g/lib_vidattr.o :	$(srcdir)/lib_vidattr.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_vidattr.c

../obj_g/lib_vline.o :	$(srcdir)/lib_vline.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_vline.c

../obj_g/lib_wattroff.o :	$(srcdir)/lib_wattroff.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_wattroff.c

../obj_g/lib_wattron.o :	$(srcdir)/lib_wattron.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_wattron.c

../obj_g/lib_window.o :	$(srcdir)/lib_window.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_window.c

../obj_g/nc_panel.o :	$(srcdir)/nc_panel.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/nc_panel.c

../obj_g/safe_sprintf.o :	$(srcdir)/safe_sprintf.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/safe_sprintf.c

../obj_g/write_entry.o :	$(srcdir)/write_entry.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/write_entry.c

../obj_g/sigaction.o :	$(srcdir)/sigaction.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/sigaction.c

../obj_g/vsscanf.o :	$(srcdir)/vsscanf.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/vsscanf.c

../obj_g/define_key.o :	$(srcdir)/define_key.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/define_key.c

../obj_g/expanded.o :	./expanded.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/expanded.c

../obj_g/keyok.o :	$(srcdir)/keyok.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/keyok.c

../obj_g/lib_adabind.o :	$(srcdir)/lib_adabind.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_adabind.c

../obj_g/lib_dft_fgbg.o :	$(srcdir)/lib_dft_fgbg.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_dft_fgbg.c

../obj_g/lib_freeall.o :	$(srcdir)/lib_freeall.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_freeall.c

../obj_g/resizeterm.o :	$(srcdir)/resizeterm.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/resizeterm.c

../obj_g/softscroll.o :	$(srcdir)/softscroll.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/softscroll.c

../obj_g/wresize.o :	$(srcdir)/wresize.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/wresize.c

../obj_g/alloc_entry.o :	$(srcdir)/alloc_entry.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/alloc_entry.c

../obj_g/captoinfo.o :	$(srcdir)/captoinfo.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/captoinfo.c

../obj_g/codes.o :	./codes.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/codes.c

../obj_g/comp_captab.o :	./comp_captab.c \
			../include/term.h \
			../include/hashsize.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_captab.c

../obj_g/comp_error.o :	$(srcdir)/comp_error.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_error.c

../obj_g/comp_expand.o :	$(srcdir)/comp_expand.c \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_expand.c

../obj_g/comp_hash.o :	$(srcdir)/comp_hash.c \
			../include/term.h \
			../include/hashsize.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_hash.c

../obj_g/comp_parse.o :	$(srcdir)/comp_parse.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_parse.c

../obj_g/comp_scan.o :	$(srcdir)/comp_scan.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/comp_scan.c

../obj_g/fallback.o :	./fallback.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/fallback.c

../obj_g/lib_cur_term.o :	$(srcdir)/lib_cur_term.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_cur_term.c

../obj_g/lib_kernel.o :	$(srcdir)/lib_kernel.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_kernel.c

../obj_g/lib_keyname.o :	./lib_keyname.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_keyname.c

../obj_g/lib_options.o :	$(srcdir)/lib_options.c \
			../include/term.h \
			keys.tries
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_options.c

../obj_g/lib_raw.o :	$(srcdir)/lib_raw.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_raw.c

../obj_g/lib_termcap.o :	$(srcdir)/lib_termcap.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_termcap.c

../obj_g/lib_ti.o :	$(srcdir)/lib_ti.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_ti.c

../obj_g/lib_tparm.o :	$(srcdir)/lib_tparm.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tparm.c

../obj_g/lib_tputs.o :	$(srcdir)/lib_tputs.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tputs.c

../obj_g/lib_trace.o :	$(srcdir)/lib_trace.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_trace.c

../obj_g/lib_tracechr.o :	$(srcdir)/lib_tracechr.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_tracechr.c

../obj_g/lib_twait.o :	$(srcdir)/lib_twait.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_twait.c

../obj_g/name_match.o :	$(srcdir)/name_match.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/name_match.c

../obj_g/names.o :	./names.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/names.c

../obj_g/parse_entry.o :	$(srcdir)/parse_entry.c \
			../include/term.h \
			../include/parametrized.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/parse_entry.c

../obj_g/read_entry.o :	$(srcdir)/read_entry.c \
			../include/term.h \
			$(INCDIR)/tic.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/read_entry.c

../obj_g/read_termcap.o :	$(srcdir)/read_termcap.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/read_termcap.c

../obj_g/lib_setup.o :	$(srcdir)/lib_setup.c \
			../include/term.h
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/lib_setup.c

../obj_g/trace_buf.o :	$(srcdir)/trace_buf.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/trace_buf.c

../obj_g/tries.o :	$(srcdir)/tries.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/tries.c

../obj_g/unctrl.o :	./unctrl.c
	cd ../obj_g; $(CC) $(CFLAGS_DEBUG) -c ../ncurses/unctrl.c

