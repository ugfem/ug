#****************************************************************************#
#*																			*#
#* File:	  graph.make													*#
#*																			*#
#* Purpose:   build ug graphics 											*#
#*																			*#
#* Author:	  Peter Bastian 												*#
#*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*#
#*			  Universitaet Heidelberg										*#
#*			  Im Neuenheimer Feld 368										*#
#*			  6900 Heidelberg												*#
#*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*#
#*																			*#
#* History:   24.05.92 begin, ug version 2.0								*#
#*			  13.08.92 update for ug 2.0.2									*#
#*			  05 Sep 92 update for ug 2.1(.0)								*#
#*			  05 Sep 94 update for ug 2.3									*#
#*			  01 Jan 95 update for ug 3.0									*#
#*																			*#
#* Remarks:   MPW 3.2 version												*#
#*																			*#
#****************************************************************************#

# list of source files
OBJECTS = wpm.c.o wop.c.o graph.c.o plotproc.c.o initgraph.c.o

# local C compiler flags
LCFLAGS = -i "::low" -i "::dev" -i "::gm"

# the main rule
all Ä {OBJECTS} graph.make
	Lib -o "::lib:libgraph{LIBSUFFIX}.a" {OBJECTS}

# compile all source files
wpm.c.o Ä  wpm.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} wpm.c

# compile all source files
wop.c.o Ä  wop.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} wop.c

# compile all source files
graph.c.o Ä  graph.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} graph.c

# compile all source files
plotproc.c.o Ä	plotproc.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} plotproc.c

# compile all source files
initgraph.c.o Ä  initgraph.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} initgraph.c

# clean up
clean Ä
	Set Exit 0
	delete -i Å.o
