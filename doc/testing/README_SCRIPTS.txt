Description of the scripts to use when testing UG with dart


1. Client side 

1.1 ug_dart_conf.pl (perl script)

usage:	ug_dart_conf.pl [-b build_dir_suffix] [-m mode] [-i test id] 
[-a architecture] [-c ugconf options]

-b: Specifying another build directory than .../UG by entering the
    sub directory, e.g. cd/appl 
-i: Test identifier
-d: Build directory identifier
-a: Architecture, e.g. PC
-c: Options which ugconf would accept except the architecture

purpose: Creating the "DartConfiguration.tcl" file which includes the 
	 configuration for each test



1.2 ug_dart_tests.pl (perl script)

usage:	ug_dart.pl [-b build directory suffix] [-s source directory suffix] 
                   [-p test application] [-i test id] 

-b: Specifying another build directory than .../UG by entering the
    sub directory, e.g. cd/appl
-s: Source directory
-p: Specifying the application to test
-i: Test identifier

purpose: Creating all "DartTestFile.txt" files. These files specify the 
	 test application and its parameters.



1.3 ug_dart_client_test (tcsh script)

usage: ug_dart_client_test [-t test type] [-b build directory]
		           [-cov coverage] [-test testing]

-t: May be Nightly, Continuous or Experimental
-b: the whole path to the library or application to build
-cov: May be yes (coverage on), or no (coverage off)
-test: May be yes (testing on), or no (testing off)		    

purpose: Running the whole configure/build/test/submit cycle. Therefor it 
         calls several Tcl scripts provided by Dart. 



1.4 ug_dart.pl (perl script)

usage: ug_dart.pl [-c configuration string] [-e test extent] [-t test type]
		   [-d build directory identifier]
-d: Build directory identifier: It is defined in the file ".../Ug/DartTestFiles"
    and specifies the directory in which a library or an application should be
    build		   
-c: Configuration String which is necessary if the test extent
    is either 'test' or 'lib'. It has to be given in the form
   -c <DirId1>:<TestId1,...,TestIdn> ... <DirIdn>:<TestId1,...,TestIdn>,
    e.g. -c 'Appl1:tutor3,tutor5 Appl2:cd36 (test extent = test)
-e: Test extent.May be
    - manually: execute a given test
    - test    : execute either all tests in the given directory
                or only several test cases
    - lib     : building libraries with either all configurations
                or only with several configurations
    - test_all: execute all tests in all directories
    - all     : build the libraries with all configurations
                execute all tests
-t: May be Nightly, Continuous or Experimental
    Default: Nightly
-h: Print help message   

purpose: running all scripts which are necessary for a complete dart
         build-test-submit cycle 



1.5 ug_dart_cron (tcsh script)

usage: ug_dart_cron

purpose: This script is called by the cron job running on the dart client.
	 It calls ug_dart.pl -e all -t Nightly.



2. Server side

2.1 ug_dart_server_start (tcsh script)

usage: ug_dart_server_start [-t test type]

-t: May be Nightly, Continuous or Experimental
    Default: Nightly

purpose: open a dashboard on the Dart server



2.2 ug_dart_server_stop (tcsh script)

usage: ug_dart_server_stop [-t test type]

-t: May be Nightly, Continuous or Experimental
    Default: Nightly

purpose: close the dashboard (generate html code with the XSLT engine)



2.3 ug_dart_copy_html (tcsh script)

usage: ug_dart_copy_html

purpose: copies over the dashboards to "~/public_html" 

