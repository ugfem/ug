Running "Experimental" builds

There are three steps to be done.

1. Open a new Dashboard

On the Dart server run "ug_dart_server_start_ex".


2. Run the Dart build/test/submit cycle

On the Dart client run 

"ug_dart_ex.pl Client <arch> <ugconf_param1 ...  ugconf_paramN>"

 where <arch> is the architecture (e.g. PC,PCI,...) and 
<ugconf_param1 ... ugconf_paramN> are all the parameters that may be passed to
ugconf.


3. Close the Dashboard and generate the HTML pages

On the Dart server run "ug_dart_server_stop_ex".
