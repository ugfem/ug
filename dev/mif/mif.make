#****************************************************************************#
#*																			*#
#* File:	  mif.make														*#
#*																			*#
#* Purpose:   build Macintosh graphical user interface						*#
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
OBJECTS = MacMain.c.o MacShell.c.o MacGraph.c.o MacSurface.c.o

# local C compiler flags
LCFLAGS = -i ":::low" -i "::"

# the main rule
all Ä {OBJECTS} mif.make ug.rsrc
	Lib -o ":::lib:libmif.a" {OBJECTS}

# compile all source files
MacShell.c.o Ä	MacShell.c
	 C {COPTS} {LCFLAGS} MacShell.c

MacMain.c.o Ä  MacMain.c
	 C {COPTS} {LCFLAGS} MacMain.c

MacGraph.c.o Ä	MacGraph.c
	 C {COPTS} {LCFLAGS} MacGraph.c

MacSurface.c.o Ä  MacSurface.c
	 C {COPTS} {LCFLAGS} MacSurface.c

ug.rsrc Ä ugrsrc.txt
	Rez -o ug.rsrc ugrsrc.txt

