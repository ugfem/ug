Set up the Dart server

1. Download and install Dart and set some environment variables 

See the README_SETUP.txt file (1.-3.).




2. Reqirements

a)
A web server capable of running a cgi-bin script has to be installed.
The perl script ".../Dart/Source/Server/Dart.pl" has to be copied onto the web
servers cgi-bin directory. The name of the script should match the
"TriggerSite" variable (see 3.a)). There are two variables in that file to
modify. "dropLocation" must match the "DropLocation" variable which is set in
"DartConfiguration.tcl.proto" (see 3.a)) and "destination" is where the XML
reports submitted by the client should be moved. "destination" should be set
to ".../UG/Testing/HTML/TestingResults/Sites" .

b)
An installation of Java.

c)
An installation of Tcl.




3. Modify the "DartConfiguration.tcl.proto" file 

You have to checkout the "DartConfiguration.tcl.proto" file from the UG cvs
repository and modify it so that it fits to your installation. Then write it
back to the repository. The modifications that have to be done are described 
below.

a) Modify the paragraph "Submission information"

# Submission information
DropSite:  
DropLocation: 
DropSiteUser: 
DropSitePassword:
DropSiteMode:
DropMethod: 
TriggerSite: 

Explanation:

-"DropSite" has to match the Dart server location.

-"DropLocation" is the directory on the Dart server where the clients leave
 their build and test results.

-"DropSiteUser" is the username which the clients use to connect to the server.

-"DropSitePassword" is password of the "DropSiteUser".

-"DropSiteMode" is an optional ftp mode and may be "active" or "passive". This
 option may only be used when "DropMethod is set to "ftp". 

-"DropMethod" is the method which is used by clients to submit their build and
 test results to the Dart server.
  
-"TriggerSite" has to match the directory where the "Dart.pl" file resides. See
 below.
 For example, if Dart.pl resides in "/home/ugtest/public-html/cgi-bin/" the
 "TriggerSite" variable has to match "http://hal.iwr.uni-heidelberg.de/~ugtest/
 cgi-bin/Dart.pl".

b) Modify the paragraph "Commands for the build/test/submit cycle"

# Commands for the build/test/submit cycle
ConfigureCommand: ---
CMakeCommand: 
MakeCommand: make build
CVSCommand: 
TclshCommand: 
JavaCommand: 
ScpCommand: 
PurifyCommand: PURIFYCOMMAND-NOTFOUND
ValgrindCommand:
ValgrindCommandOptions:
# Compression commands
GunzipCommand: 
CompressionCommand: 
CompressionType: 

Comment: In the file "DartConfiguration.tcl.proto" a few variables have the
value "---". This value musten'd be changed! It will be replaced by the perl
script "ug_dart_conf.pl".
When an application mentioned above isn't installed on your system you have to 
fill in "...COMMAND-NOTFOUND" (see the PurifyCommand above).
The "MakeCommand" should be "ugpart" except you want to use your own make
script.




4. Prepare the directories needed on the server

Checkout the UG source from the cvs repository. This has to be done only once!




5. Configure Dart

Simply run "ug_dart_conf.pl" to create the "DartConfiguration.tcl" file
using "DartConfiguration.tcl.proto".




6. Setting up cron jobs to run Dart tests automatically

You need a cron job which calls the script "ug_dart_server_start" to open a new
Dashboard. For a Nightly Build you may do this at 00.00. Further you also need
to run "ug_dart_server_stop" to end this Dashboard and build HTML pages out of
the XSL files which were copied onto the server by the clients.
Note: The cron jobs need to be run as the same user as the web server.
