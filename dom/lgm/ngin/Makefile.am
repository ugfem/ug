if BUILD3D
  LIB_3D = libngin.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_3D)

# prefix to use instead of 'yy'
PARSERPREFIX = ng

# generate header file, rename prefix
AM_YFLAGS = -d -p $(PARSERPREFIX)
AM_LFLAGS = -P$(PARSERPREFIX)
LEX_OUTPUT_ROOT = lex.$(PARSERPREFIX)

BUILT_SOURCES = ngin-yacc.hh

libngin_la_SOURCES = ngin-yacc.yy ngin-lex.ll ngin-lex.hh ng.cc ng.h
libngin_la_CPPFLAGS = $(UG_3DFLAGS) -I$(top_srcdir)/dom/lgm
