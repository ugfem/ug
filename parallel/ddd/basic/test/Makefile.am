# -*- tab-width: 4; indent-tabs-mode: nil -*-

TESTPROGS = testbtree

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

# programs to run when "make check" is used
TESTS = $(TESTPROGS)

# define the programs
testbtree_SOURCES = testbtree.cc
testbtree_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

#include $(top_srcdir)/am/global-rules
