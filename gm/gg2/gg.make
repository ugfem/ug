#****************************************************************************#
#*																			*#
#* File:	  gg.make													*#
#*																			*#
#* Purpose:   build the numerics module     								*#
#*																			*#
#* Author:	  Dirk Feuchter 												*#
#*																			*#
#* History:   23.10.95 begin, ug version 3.1								*#
#*																			*#
#* Remarks:   MPW 3.2 version												*#
#*																			*#
#****************************************************************************#

# local C compiler flags
LCFLAGS = -i "::low" -i "::dev" -i "::gm" -i "::numerics" -i "::ui"

# object files for both dimensions
OBJECTS = basics.c.o  initnumerics.c.o mgsolver.c.o nlsolver.c.o num.c.o smoother.c.o transgrid.c.o ugblas.c.o ugiter.c.o tnlsolve.c.o

# objects for 2D version
OBJECTS2 = fv.c.o

# objects for 3D version
OBJECTS3 = 

# trick to switch source files depending on dimension
DOBJECTS = OBJECTS{DIM}

# the main rules
2Dversion Ÿ {OBJECTS} {OBJECTS2} numerics.make
	Lib -o "::lib:libnumerics{LIBSUFFIX}.a" {OBJECTS} {OBJECTS2}

3Dversion Ÿ {OBJECTS} {OBJECTS3} numerics.make
	Lib -o "::lib:libnumerics{LIBSUFFIX}.a" {OBJECTS} {OBJECTS3}

# compile all source files
basics.c.o Ÿ  basics.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} basics.c

initnumerics.c.o Ÿ  initnumerics.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} initnumerics.c

mgsolver.c.o Ÿ  mgsolver.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} mgsolver.c

nlsolver.c.o Ÿ  nlsolver.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} nlsolver.c

num.c.o Ÿ  num.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} num.c

smoother.c.o Ÿ  smoother.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} smoother.c

transgrid.c.o Ÿ  transgrid.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} transgrid.c

ugblas.c.o Ÿ  ugblas.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugblas.c

ugiter.c.o Ÿ  ugiter.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} ugiter.c

fv.c.o Ÿ  fv.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} fv.c

tnlsolve.c.o Ÿ tnlsolve.c
	 C {COPTS} {LCFLAGS} {UGDEFINES} tnlsolve.c

# clean up
clean Ÿ
	Set Exit 0
	delete -i ‰.o
