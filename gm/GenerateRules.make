#	File:		GenerateRules.make
#	Target: 	GenerateRules
#	Sources:	GenerateRules.c simplex.c
#	Created:	Dienstag, 14. September 1993 14:44:07 Uhr


OBJECTS = GenerateRules.c.o simplex.c.o



GenerateRules ÄÄ GenerateRules.make {OBJECTS}
	Link -d -c 'MPS ' -t MPST -sym on -mf ¶
		{OBJECTS} ¶
		"{CLibraries}"Clib881.o ¶
		#"{CLibraries}"CSANELib881.o ¶
		#"{CLibraries}"Math881.o ¶
		#"{CLibraries}"Complex881.o ¶
		"{CLibraries}"StdClib.o ¶
		"{Libraries}"Stubs.o ¶
		"{Libraries}"Runtime.o ¶
		"{Libraries}"Interface.o ¶
		#"{Libraries}"ToolLibs.o ¶
		-o GenerateRules
GenerateRules.c.o Ä GenerateRules.make GenerateRules.c
	 C -d __MPW32__ -r -sym on -mc68020 -mc68881 -elems881 GenerateRules.c
#	 C -d TEST_GR -d __MPW32__ -r -sym on -mc68020 -mc68881 -elems881 GenerateRules.c
simplex.c.o Ä GenerateRules.make simplex.c
	 C -d __MPW32__ -r -sym on -mc68020 -mc68881 -elems881 simplex.c
