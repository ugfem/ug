package ug;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw(ug);
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();

##############################################
## source 
##############################################

use IPC::Open2;
BEGIN
{
	my $debug=0;
	sub debug
	{
		$debug=$_[0];
	}
	sub submit
	{
		if ($_[0]=~/quit/) {die "ERROR: quit command is blocked, use 'end'\n";}
	    print IN $_[0];
		if ($debug) { print "UG-CMD: $_[0]"} ;
	}
}
sub out
{
	my ($line,$ret);
	$ret="";
	while($line=<OUT>)
	{
		if ($line=~/^EOO$/) {last;}
		if ($_[0]) {print $line;}
		$ret.=$line;
	}
	return $ret;
}
sub set
{
	print IN "set $_[0]\n";
    return split /[=\s]+/,out(0);
}
sub ug
{
	my ($i,$cmd,$print,$command,$ui,$stat);

	# basic check
	if (@_<=0) 
	{
		die "ERROR: provide ug command\n";
	} 

	# detect internal print
	$command=$_[0]; $print=0;
	if ($command=~/^\s*print\s+/) 
	{
		$print=1; 
		$command=~s/\s*^print\s+//g;
		if ($command eq "") {die "ERROR: print must come with ug command\n";}
	}

	# check for I/O channels
	if ($command ne "start")
	{
		1==stat IN and 1==stat OUT or die "ERROR: IN/OUT channel missing\n";
		1==stat IN or die "ERROR: IN channel missing\n";
		1==stat OUT or die "ERROR: OUT channel missing\n";
	}

	SWITCH:
	{
		# debug
		if ($command eq "debug")
		{
			if (@_!=2 || ($_[1]!=0 && $_[1]!=1))
            {
                die "ERROR: provide [0|1] with 'debug'\n";
            }
			debug $_[1];
		}

		# command 'end'
		if ($command eq "end")
		{
			if(@_!=1)
        	{   
        	    die "ERROR: don't provide any option with 'end'\n";
        	} 
        	print IN "quit\n";
        	close(IN);
        	close(OUT);
			return;
		}

		# command 'set'
		if ($command eq "set")
		{
			if(@_!=2)
            {   
                die "ERROR: provide one option with 'set'\n";
            }
			print IN "set $_[1]\n";
			return split /[=\s]+/,out($print);
		}

		# command 'start'
		if ($command eq "start")
		{
			if(@_<2 || @_>3 || (@_==3 && $_[2] ne "x") || $_[1]=~/-ui/)
        	{   
        	    die 'ERROR: usage: ug("start", "<program>" [,"x"]);'."\n";
        	} 
			
			if (@_==2)	{$ui="-ui cn";}
			if (@_==3)	{$ui="-ui c";}
        	open2(*OUT,*IN,"$_[1] $ui -perl");
			return out($print);
		}
	
		# command 'ug' 
		if ($command eq "ug")
        {
			if(@_!=2)
            {
                die "ERROR: command 'ug' must come with one argument\n\n";
            }
			submit "$_[1]\n";
			return out($print);
		}

		# std command
		if (@_%2==1) {$cmd="$command "; $i=1;}
		else {$cmd="$command $_[1] "; $i=2;}
		for (;$i<@_;$i+=2)
		{
			$cmd.='$'."$_[$i] $_[$i+1] ";
		}
		$cmd.="\n";
		submit $cmd;
		return out($print);
	}
}

1;




