#****************************************************************************#
#*																			*#
#* File:	  low.make														*#
#*																			*#
#* Purpose:   make file for low modules library 							*#
#*																			*#
#* Author:	  Peter Bastian 												*#
#*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*#
#*			  Universitaet Heidelberg										*#
#*			  Im Neuenheimer Feld 368										*#
#*			  6900 Heidelberg												*#
#*			  internet: ug@ica3.uni-stuttgart.de					*#
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

#	all object files for ug library
OBJECTS = ugenv.c.o heaps.c.o fifo.c.o misc.c.o initlow.c.o defaults.c.o fileopen.c.o ugstruct.c.o

#	main rule
all Ä low.make {OBJECTS}
	Lib -o "::lib:liblow{LIBSUFFIX}.a" {OBJECTS}


# compile all source files

initlow.c.o Ä  initlow.c
	 C {COPTS} {UGDEFINES} initlow.c

misc.c.o Ä	misc.c
	 C {COPTS} {UGDEFINES} misc.c

defaults.c.o Ä	defaults.c
	 C {COPTS} {UGDEFINES} defaults.c

ugenv.c.o Ä  ugenv.c
	 C {COPTS} {UGDEFINES} ugenv.c

heaps.c.o Ä  heaps.c
	 C {COPTS} {UGDEFINES} heaps.c

fifo.c.o Ä	fifo.c
	 C {COPTS} {UGDEFINES} fifo.c

fileopen.c.o Ä	fileopen.c
	 C {COPTS} {UGDEFINES} fileopen.c

ugstruct.c.o Ä	ugstruct.c
	 C {COPTS} {UGDEFINES} ugstruct.c
	
# clean up
clean Ä
	delete -i Å.o
