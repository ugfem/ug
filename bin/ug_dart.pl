#!/usr/bin/perl
#
# ug_dart.pl
# 
# purpose: runs all necessary scripts for either the server (mode=Server)
#          or the client (mode=Client)

if($ARGV[0] eq "-help" || $ARGV[0] eq "-h")
{
	print "usage:   ug_dart.pl modeparameter parameter1 ... parametern\n";
	print "         modeparameter: calling this script on the server: modeparameter = Server\n";
	print "                        calling this script on the client: modeparameter = Client\n";
	print "         parameter1 ... parametern: same parameters as ugconf\n";
	print "\n";
	print "purpose: running all scripts which are necessary for a complete dart\n";
	print "         build/test/submit cycle\n";
	die "\n";
}

my $test_type = "Nightly";
if($ARGV[0] eq "Experimental" || $ARGV[0] eq "Nightly" || $ARGV[0] eq "Continuous")
{
	$test_type = splice(@ARGV,0,1);
}

my $mode = "Client";
if($ARGV[0] eq "Client" || $ARGV[0] eq "Server")
{
	$mode = $ARGV[0];
}
my $command = "ug_dart_conf.pl";
foreach $param (@ARGV)
{
	$command=join(' ',$command,$param);
}
if($mode eq "Client")
{
       # checkout only for the first build/test/submit cycle			
        #system("ug_dart_co $mode");
	system("$command");
	system("ug_dart_tests.pl");
	system("ug_dart_client_test $test_type");
}
elsif($mode eq "Server")
{
       # checkout only for the first build/test/submit cycle 
        #system("ug_dart_co $mode");  
       # generate the build tree only for the first build/test/submit cycle
        #system("ug_dart_tree");
	system("$command");
	system("ug_dart_server_start $test_type");
}	
