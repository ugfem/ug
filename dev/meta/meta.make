#****************************************************************************#
#*																			*#
#* File:	  meta.make 													*#
#*																			*#
#* Purpose:   build metafile device 										*#
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

# list of source files
OBJECTS = metafile.c.o

# local C compiler flags
LCFLAGS = -i ":::low" -i "::"

# the main rule
all Ä {OBJECTS} meta.make
	Lib -o ":::lib:libmeta.a" {OBJECTS}

# compile all source files
metafile.c.o Ä	metafile.c
	 C {COPTS} {LCFLAGS} metafile.c


