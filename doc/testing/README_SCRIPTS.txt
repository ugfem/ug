Description of the scripts to use when testing UG with dart (client usage)

1. ug_dart_co (tcsh script)

purpose: Cvs checkout of UG into the right directory

usage: ug_dart_co <model>
     
       <model> can be either "Client" or "Server". The default is set to 
       "Client".



2. ug_dart_conf.pl (perl script)

purpose: Creating the "DartConfiguration.tcl" file which is needed to
         run the test cycle.

usage: ug_dart_conf.pl -bd <build_dir_suffix> <model> <ugconf_arg 1 ... ugconfarg n>
     
       If you want to specify another build directory than ".../UG" which is 
       default you have to enter -bd followed by the subdirectory. 
       Example: The build directory should be ".../UG/cd/appl". Then the parameter
       -bg cd/appl is suitable. 			
       <model> may be either "Client" or "Server". The default is set to 
       "Client". 
       <ugconf_arg i> is an argument that should be passed to 
       ugconf. It's important that the the first of these arguments has to be
       the architecture specification. For further information on ugconf 
       see the UG manual. 



3. ug_dart_tests.pl (perl script)

purpose: Creating the needed "DartTestfile.txt" files. That doesn't include
         creating of "DartTestfile.txt" files in which tests are specified.

usage: ug_dart_tests.pl <build_dir_suffix>
     
       <build_dir_suffix>: see the description of "ug_dart_conf.pl". 		 	       



4. ug_dart_client_test (tcsh script)

purpose: Running the whole configure/build/test/submit cycle. Therefor it 
         calls several Tcl scripts provided by Dart. 

usage: ug_dart_client_test <test_type> <coverage> <build_dir_suffix> 

       <test_type> may be "Nightly", "Continuous" or "Experimental".
       <coverage> may be either "COV" if coverage information should be
       evaluated.
       <build_dir_suffix>: see the description of "ug_dart_conf.pl". 		 	       



5. ug_dart.pl (perl script)

purpose: It calls all the scripts above to run a whole cycle.

usage: ug_dart.pl -bd <build_dir_suffix> <test_type> <mode> <ugconf_arg 1 ... ugconfarg n>

       -bg <build_dir_suffix>: see the description of "ug_dart_conf.pl"
       <test_type> may be "Nightly", "Continuous" or "Experimental".	 	
       <mode> can be either "Client" or "Server". The default is set to
       "Client". <ugconf_arg i> is an argument that should be passed to
       ugconf. It's important that the the first of these arguments has to be
       the architecture specification. For further information on ugconf
       see the UG manual



6. ugpart (tcsh script)

purpose: building UG - libraries and some executables

usage: ugpart

       This script is called automatically by Dart. 

