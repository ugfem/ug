Running "Experimental" builds

There are four steps to be done.

1. Open a new dashboard on the server

	ug_dart_server_start

2. Prepare the build

	Create a file called "DartTestfile.txt" in the directory
	where the application to be tested will reside. In this file
	you have to specify the test with the command

		ADD_TEST( test_name make test ). 
	Then add two targets to the according Makefile:
	a) build: building the needed libraries.
	b) test: building and running the apllication. 
	

3. Run the build-test-submit cycle on the client

	ug_dart.pl -bg <build_dir_suffix> Experimental <ugconf_params>

4. Close the dashboard on the server

	ug_dart_server_stop

For further information on using the script "ug_dart.pl" in step 2. please
see the README_SCRIPTS.txt file.

Example: Testing of tutor2d which will be build in the directory ".../UG/tutor/appl"

1. ug_dart_server_start

2. Write 

	ADD_TEST( tutor make test )

   into the file ".../UG/tutor/appl/DartTestfile.txt".
   Add the two targets test and build to the file ".../UG/tutor/appl/Makefile":

3. ug_dart.pl -bg tutor/appl Experimental PC 

4. ug_dart_server_stop


			 
	

