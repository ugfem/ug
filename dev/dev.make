#****************************************************************************#
#*																			*#
#* File:	  dev.make														*#
#*																			*#
#* Purpose:   build devices 												*#
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

# list of devices to build
DEVICES = META MIF

# list of source files
OBJECTS = devices.c.o

# local C compiler flags
LCFLAGS = -i "::low"

# the main rule
all Ä {OBJECTS} {DEVICES} dev.make
	Lib -o "::lib:libdev.a" {OBJECTS}

# compile all source files
devices.c.o Ä  devices.c
	 C {COPTS} {LCFLAGS} devices.c

# build devices
META Ä
	directory meta
	make -d COPTS="{COPTS}" -f meta.make > makeout
	makeout
	delete makeout
	directory ::

MIF Ä
	directory mif
	make -d COPTS="{COPTS}" -f mif.make > makeout
	makeout
	delete makeout
	directory ::

# clean up
clean Ä
	Set Exit 0
	delete -i "::lib:libdev.a"
	delete -i Å.o
	delete -i :meta:Å.o
	delete -i :mif:Å.o
