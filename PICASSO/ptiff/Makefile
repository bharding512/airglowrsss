#
# "$Id: Makefile.in,v 1.35 2003/09/13 23:00:35 easysw Exp $"
#
#   Makefile for flphoto.
#
#   Copyright 2002-2003 by Michael Sweet
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#

#
# Programs...
#

CC		=	gcc
CP		=	/bin/cp
CXX		=	g++
FLTKCONFIG	=	/usr/local/bin/fltk-config
MKDIR		=	/bin/mkdir -p
MSGFMT		=	/usr/bin/msgfmt
NROFF		=	/usr/bin/nroff
RANLIB		=	ranlib
RM		=	/bin/rm -f
SHELL		=	/bin/sh


#
# Windows crapola...
#
# Note: I do not test flPhoto on Windows and offer no assistance in
#       getting it to work on that platform.
#

EXEEXT		=	


#
# Program options...
#
# OPTIM   defines the common compiler optimization/debugging options.
#

CFLAGS		=	 -I/usr/local/include -I/usr/X11R6/include   $(OPTIM)
CXXFLAGS	=	 -I/usr/local/include -I/usr/X11R6/include   $(OPTIM)
LDFLAGS		=	 $(OPTIM)
LIBS		=	-L/usr/local/lib -L/usr/X11R6/lib -lfltk_images -lpng -lz -ljpeg -lfltk -lm -lXext -lX11 -lsupc++
OPTIM		=	-Wall -Wno-return-type -O -fno-rtti -fno-exceptions -DHAVE_CONFIG_H


#
# Manpage extensions...
#

CAT1EXT		=	1
MAN1EXT		=	1
CAT5EXT		=	5
MAN5EXT		=	5


#
# Directories...
#

bindir		=	${exec_prefix}/bin
datadir		=	/usr/local/share
docdir		=	/usr/local/share/doc/flphoto
exec_prefix	=	/usr/local
mandir		=	${prefix}/man
prefix		=	/usr/local
libdir		=	${exec_prefix}/lib
srcdir		=	.
FLPHOTO_LOCALE	=	/usr/local/share/locale

DESTDIR		=



#
# Rules...
#

.SILENT:

.SUFFIXES:	.c .cxx .h .mo .o .po
.c.o:
	echo Compiling $<...
	$(CC) $(CFLAGS) -I. -I$(srcdir) -c $<
.cxx.o:
	echo Compiling $<...
	$(CXX) $(CXXFLAGS) -I. -I$(srcdir) -c $<
.po:	espmsg$(EXEEXT)
	echo Generating $@...
	./espmsg $< compile


#
# Targets...
#

TRANSLATIONS =	po/de po/en_CA po/en_GB po/es po/fr po/it po/nl po/pt
TARGETS	=	ptiff$(EXEEXT)
CXXFILES =	\
		Fl_Image_Display.cxx \
		Fl_TIFF_Image.cxx \
		ptiffUI.cxx \
		util.cxx \
		main.cxx

OBJS	=	$(CFILES:.c=.o) $(CXXFILES:.cxx=.o)


#
# Make all targets...
#

all:	$(TARGETS)


#
# Clean all object files...
#

clean:
	echo Cleaning out object files...
	$(RM) $(OBJS)
	$(RM) $(TARGETS)
	$(RM) $(TRANSLATIONS)


#
# Make dependencies...
#

depend:
	makedepend -Y -f Dependencies *.c *.cxx

#
# flphoto
#

ptiff$(EXEEXT):	$(OBJS)
	echo Linking ptiff$(EXEEXT)...
	$(CXX) $(CXXFLAGS) -o ptiff$(EXEEXT) $(OBJS) $(LIBS)
	$(FLTKCONFIG) --post tiff



