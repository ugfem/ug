package gnuplot;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw(gnuplot);
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();

##############################################
## source 
##############################################

use IO::File;
use POSIX qw(tmpnam);

BEGIN
{
	my $end='';
	END
	{
		eval $end;
	}
	sub getname
	{
		my $name;
		do { $name=tmpnam(); } until IO::File->new($name,O_RDWR|O_CREAT|O_EXCL);
		$end.="unlink('$name');";
		return $name;
	}
}

sub gnuplot
{
	my ($command,$pid,$fh);

	# basic check
	if (@_<=0) 
	{
		die "ERROR: provide gnuplot command\n";
	} 
	$command=$_[0];

	# check for I channel
	if ($command ne "start")
	{
		1==stat IN or die "ERROR: IN channel missing\n";
	}

	SWITCH:
	{
		# command 'start'
		if ($command eq "start")
		{
			if(@_!=1)
        	{   
        	    die "ERROR: provide no option with 'start'\n";
        	} 
        	$pid=open(IN,"| gnuplot");
			$fh=select(IN);
			$|=1;
			select($fh);
			return $pid;
		}

		# command 'end'
        if ($command eq "end")
        {
            if(@_!=1)
            {
                die "ERROR: provide no option with 'end'\n";
            }
            close(IN);
			return;
        }
	
		# std command
        if (@_==1)
        {
			print IN "$_[0]\n";
            return;
		}
		if (@_>=2)
		{
			my ($i,$j,$s,$t,$name);
			$t=$_[0]; 
			for ($j=1; $j<@_; $j++)
			{
				$s=$_[$j]; 
				$name=getname;
				open(TMP,">$name"); 
				for $i (sort {$a<=>$b} keys %$s) { print TMP "$i $$s{$i}\n"; } 
				close(TMP);
				$t=~s/%/"$name"/;
			}
			print IN "$t\n";
            return;
		}

	}
}

1;




