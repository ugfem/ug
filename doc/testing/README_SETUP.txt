Set up Dart for UG (client usage) 

1. Download Dart

Download Dart by using cvs

 cvs -d :pserver:anonymous@public.kitware.com:/Dart/cvsroot login
 (respond with password dart) 

Follow this command by checking out the source code:
 
 cvs -d :pserver:anonymous@public.kitware.com:/Dart/cvsroot co Dart




2. Installing Dart

You have to install Dart by copying the downloaded directory "Dart" anywhere
to your system. Do NOT delete the subdirectory "Source"!




3. Setting some environment variables in shell resource file

There are two environment variables that have to be set in .tcshrc, .cshrc, ..., 
one for Dart and one for UG.

(i)  DART_HOME has to be set to the root of your Dart installation (../Dart).
(ii) UGROOT has to be set to "../UG/ug"

Further more the environmental path has to be extended that it contains 
the directory "../UG/ug/bin". This is also described in the UG manual.




4. Setting up Dart for UG to use it with our Dart server

First of all you need an installation of Tcl and ssh on your system because
the configure/build/test cycle is driven by Tcl scripts and the data is
copied to the Dart server via scp.
The only thing left to do is defining tests. There are several tests 
predefined, extend by new lines if you want to add your own tests. For further information
see the README_TESTS.txt file.




5. Running the configure/build/test/submit cycle

There are scripts in ../UG/ug/bin which you can use to run a cycle.
For further information on these scripts see the README_SCRIPTS.txt file.




6. Setting up a Dart server for UG 

Please see the README_ADVANCED.txt file.


 
