#****************************************************************************#
#*																			*#
#* File:	  gm.make														*#
#*																			*#
#* Purpose:   build the grid manager module 								*#
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

# local C compiler flags
LCFLAGS = -i "::low" -i "::dev"

# object files for both dimensions
OBJECTS = algebra.c.o enrol.c.o evm.c.o ugio.c.o ugm.c.o cw.c.o initgm.c.o elements.c.o

# objects for 2D version
OBJECTS2 = shapes2d.c.o ugm2d.c.o ugrefine2d.c.o

# objects for 3D version
OBJECTS3 = shapes3d.c.o ugm3d.c.o simplex.c.o ugrefine3d.c.o

# trick to switch source files depending on dimension
DOBJECTS = OBJECTS{DIM}

# the main rules
2Dversion Ä {OBJECTS} {OBJECTS2} gm.make
	Lib -o "::lib:libgm{LIBSUFFIX}.a" {OBJECTS} {OBJECTS2}

3Dversion Ä {OBJECTS} {OBJECTS3} gm.make
	Lib -o "::lib:libgm{LIBSUFFIX}.a" {OBJECTS} {OBJECTS3}

# compile all source files
initgm.c.o Ä  initgm.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} initgm.c
elements.c.o Ä  elements.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} elements.c
cw.c.o Ä  cw.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} cw.c
algebra.c.o Ä  algebra.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} algebra.c
enrol.c.o Ä  enrol.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} enrol.c
evm.c.o Ä  evm.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} evm.c
ugio.c.o Ä	ugio.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugio.c
ugm.c.o Ä  ugm.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugm.c
shapes2d.c.o Ä	shapes2d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} shapes2d.c
ugm2d.c.o Ä  ugm2d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugm2d.c
ugrefine2d.c.o Ä  ugrefine2d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugrefine2d.c
shapes3d.c.o Ä	shapes3d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} shapes3d.c
ugm3d.c.o Ä  ugm3d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugm3d.c
simplex.c.o Ä  simplex.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} simplex.c
ugrefine3d.c.o Ä  ugrefine3d.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugrefine3d.c

# clean up
clean Ä
	Set Exit 0
	delete -i Å.o
