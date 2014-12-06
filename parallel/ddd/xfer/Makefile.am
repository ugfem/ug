SUBDIRS =

if PARALLEL
  PLIB = libxfer2.la libxfer3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = sll.h xfer.h

EXTRA_DIST = sll.ct sll.ht

SOURCES = cmds.c \
          cmdmsg.c \
          cplmsg.c \
          ctrl.c \
          pack.c \
          supp.c \
          unpack.c \
          xfer.c

libxfer2_la_SOURCES = $(SOURCES)
libxfer2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS) -DXferMemFromHeap

libxfer3_la_SOURCES = $(SOURCES)
libxfer3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS) -DXferMemFromHeap
