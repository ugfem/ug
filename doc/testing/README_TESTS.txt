Defining and running tests.

1. Defining tests

(i)  There is a file in ".../UG" called "DartTestFiles" which includes the paths
     to the libraries or executables to test. The syntax for an entry is:

     <Build Directory Identifier>:<Build Directory>

     Example:
    
     UGLib:ug
     CD:cd/appl

     There is one entry per line. The build directory is a relative path. The
     Build Directory Identifier has to be unique!

(ii) In the build directory given in (i) there has to be a file called
     "TestCases" which includes the tests that should be run. The syntax for
     an entry is:
        
     <Test Identifier>::<Test Application>::<ECFS>
     
     <ECFS> = <Build Directory Identifier>,<CFS>
     <CFS> = <ARCH>-<MODEL>-<DIM>-<DOM_MODULE>-<DEBUG_MODE>-<OPTIM_MODE>

     The Extended Configuration String (ECFS) contains the Configuration
     Strings (CFS) for the library or application itself and for every library
     needed to build the library or application to test. To specify which CFS
     belongs to a certain library or application one has to add the Build
     Directory Identifier as well. The CFS consists of parameters that ugconf 
     needs to set up a configuration for the build. For further information on
     the parameters please see the UG documentation. The order of the
     parameters given above has to be minded!
 
     Example:	

     UGlib1::::UGlib,PC-SEQ-2-LGM_DOMAIN-DEBUG-OPTIM
     UGlib2::::UGlib,PCGCOV-SEQ-3-LGM_DOMAIN-DEBUG-OPTIM




2. Running tests
   
   There are three steps to be done.

   (i) Open the Dashboard

       On the Dart server run "ug_dart_server_start -t <Test Type>".


   (ii) Run the Dart build/test/submit cycle

	On the Dart client run "ug_dart.pl". For further informations about the
	parameters to pass to "ug_dart.pl" see the file README_SCRIPTS.txt. 

   (iii) Close the Dashboard and generate the HTML pages

	 On the Dart server run "ug_dart_server_stop -t <Test Type>".
	 
   
   Remark: The Test Type may be Experimental, Nightly or Continuous. We will
	   only need Experimental or Nightly. You can run an experimental
	   build whenever you want to test something and you don't want to
	   wait to the next day until the nightly builds are run.
	   When you invoke ug_dart_server_start and ug_dart_server_stop setting 
	   the Test Type to "Experimental" the results of this build will
	   appear on their own "experimental dashboard" and not on the nightly
	   dashboard!     
  
   
   