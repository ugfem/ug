#****************************************************************************#
#*																			*#
#* File:	  ug.make														*#
#*																			*#
#* Purpose:   top level makefile to build the ug library					*#
#*			  for the Macintosh under MPW									*#
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

# set DIM to 2 or 3 
DIM = 2

# set the following switches to 'T' if you need degrees of freedom in the
# corresponding geometric object or to 'F' else.
NODEDATA = T
EDGEDATA = F
SIDEDATA = F
ELEMDATA = F

# These are the command line defines passed to ug source files
UGDEFINES = -d _{DIM} -d Node{NODEDATA} -d Edge{EDGEDATA} -d Side{SIDEDATA} -d Elem{ELEMDATA}

# library name suffix
LIBSUFFIX = {DIM}{NODEDATA}{EDGEDATA}{SIDEDATA}{ELEMDATA}

# dimension dependent compilation
version = {DIM}Dversion

#	C compiler options to use
COPTS = -d __MPW32__ -m -opt speed -model far -warnings on -r -sym on -mc68020 -mc68881 -elems881

# list of modules to build
MODULES = LOW DEV GM NUMERICS GRAPH UI

# list of source files
OBJECTS = initug.c.o

# local C compiler flags
LCFLAGS = -i ":low" -i ":dev" -i ":gm" -i ":numerics" -i ":graph" -i ":ui"

# the main rule
ug Ä {MODULES} {OBJECTS}
	#    Macintosh version of ug 3.0 compiled successfully

# compile all source files
initug.c.o Ä  initug.c
	C {COPTS} {LCFLAGS} initug.c
	Lib -o ":lib:libinit{LIBSUFFIX}.a" {OBJECTS}

# build module low
LOW Ä
	directory low
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f low.make > makeout
	makeout
	delete makeout
	directory ::

DEV Ä
	directory dev
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f dev.make > makeout
	makeout
	delete makeout
	directory ::

GM Ä
	directory gm
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f gm.make {version} > makeout
	makeout
	delete makeout
	directory ::

NUMERICS Ä
	directory numerics
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f numerics.make {version} > makeout
	makeout
	delete makeout
	directory ::

GRAPH Ä
	directory graph 
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f graph.make > makeout
	makeout
	delete makeout
	directory ::

UI Ä
	directory ui
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f ui.make > makeout
	makeout
	delete makeout
	directory ::

NUMERICS Ä
	directory numerics
	make -d UGDEFINES="{UGDEFINES}" -d LIBSUFFIX="{LIBSUFFIX}" -d COPTS="{COPTS}" -f numerics.make > makeout
	makeout
	delete makeout
	directory ::

CleanAll Ä
	Set Exit 0
	delete -i ":lib:libug{LIBSUFFIX}.a"
	delete -i ":lib:libdev.a"
	delete -i Å.makeout
	delete -i :low:Å.o
	delete -i :dev:Å.o
	delete -i :dev:meta:Å.o
	delete -i :dev:mif:Å.o
	delete -i :gm:Å.o
	delete -i :graph:Å.o
	delete -i :ui:Å.o
	delete -i :numerics:Å.o
	echo "Macintosh version of ug 3.0 cleaned up"
