#!/usr/bin/perl
#
# ug_dart.pl
# 
# purpose: runs all necessary scripts for either the server (mode=Server)
#          or the client (mode=Client)

# print help message 
if($ARGV[0] eq "-help" || $ARGV[0] eq "--help" || $ARGV[0] eq "-h")
{
	print "usage:	ug_dart.pl -bg build_dir_suffix test_type modeparameter\n";
        print "                    architecture parameter1 ... parametern\n";
	print "	 -bg build_dir_suffix: specifying another build directory than .../UG by entering the\n";
	print "        	              sub directory, e.g. cd/appl\n"; 
	print "  	 test_type: May be Nightly, Continuous or Experimental\n";
	print "         modeparameter: calling this script on the server: modeparameter = Server\n";
	print "                        calling this script on the client: modeparameter = Client\n";
	print "         parameter1 ... parametern: same parameters as ugconf\n";
	print "\n";
	print "purpose: running all scripts which are necessary for a complete dart\n";
	print "         build/test/submit cycle\n";
	die "\n";
}
# set up a variable to count the given parameters
my $param_count = 0;

# specify the build directory
my $build_dir = ''; # default
if($ARGV[0] eq "-bd")
{
	$param_count++;
	$build_dir = $ARGV[$param_count];
	$param_count++;
}

# specify the test type (Nightly, Experimantal, Continuous)
my $test_type = "Nightly"; # default
if($ARGV[$param_count] eq "Experimental" || $ARGV[$param_count] eq "Nightly" || $ARGV[$param_count] eq "Continuous")
{
	$test_type = splice(@ARGV,$param_count,1);
}

# specify the mode (Client, Server)
my $mode = "Client"; # default
if($ARGV[$param_count] eq "Client" || $ARGV[$param_count] eq "Server")
{
	$mode = $ARGV[$param_count];
	$param_count++;
}

# build the ug_dart_conf command
my $command = "ug_dart_conf.pl";
foreach $param (@ARGV)
{
	$command=join(' ',$command,$param);
}

# specify the architecture
my $arch = $ARGV[$param_count];

# specify the build parameter (COV, NOCOV)
my $build_param = "NOCOV";
if($arch =~ /GCOV/)
{
	$build_param = "COV";
}

# run the build- test-submit-cycle either client or open/close a Dashboard
# on the server
 
if($mode eq "Client")
{
       # checkout only for the first build/test/submit cycle			
        #system("ug_dart_co $mode");
	system("$command");
	system("ug_dart_tests.pl $build_dir");
	system("ug_dart_client_test $test_type $build_param $build_dir");
}
elsif($mode eq "Server")
{
       # checkout only for the first build/test/submit cycle 
        #system("ug_dart_co $mode");  
       # generate the build tree only for the first build/test/submit cycle
        #system("ug_dart_tree");
	system("$command");
	system("ug_dart_server_start $test_type $build_dir");
}	
