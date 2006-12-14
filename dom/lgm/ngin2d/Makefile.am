# $Id$

if BUILD2D
  LIB_2D = libngin2d.la
else
  LIB_2D =
endif

noinst_LTLIBRARIES = $(LIB_2D)

# prefix to use instead of 'yy'
PARSERPREFIX = ng2d

# generate header file, rename prefix
AM_YFLAGS = -d -p $(PARSERPREFIX)
AM_LFLAGS = -P$(PARSERPREFIX)
LEX_OUTPUT_ROOT = lex.$(PARSERPREFIX)

BUILT_SOURCES = ngin-yacc.h

libngin2d_la_SOURCES = ng.c ngin-yacc.y ngin-lex.l ngin-lex.h ng2d.h
libngin2d_la_CPPFLAGS = $(UG_2DFLAGS) -I$(top_srcdir)/dom/lgm
