package ug;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw(ug float tmpfile time_start time_stop time_set time_eval);
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();

##############################################
## source 
##############################################

use IPC::Open2;
use IO::Handle;
use IO::File;
use POSIX qw(tmpnam);


BEGIN
{
	my $debug=0;
	my $end='';
	my $TimeHiRes;
	my %time;

	sub module
	{
		my $i;
		for ($i=0; $i<@INC; $i++) 
		{
			-e "$INC[$i]/$_[0]" and return 1;
		}
		return 0;
	}

	# time
	$TimeHiRes=0; if (module('Time/HiRes.pm')) { require Time::HiRes; $TimeHiRes=1; };
	sub gettimeofday
	{
		if (!$TimeHiRes) { return -1; }
		return Time::HiRes::gettimeofday();
	}
	sub time_start
	{
		my $name=shift;
		
		if (!(defined $time{"$name"})) { $time{"total $name"}=$time{"diff $name"}=0; } 
		$time{$name}=gettimeofday();
		$time{"running $name"}=1;
	}
	sub time_stop
	{
		my $dt;
		my $name=shift;

		defined $time{$name} or return;
		$dt=gettimeofday()-$time{$name};
		$time{"total $name"}+=$dt;
		$time{"diff $name"}+=$dt;
		$time{"running $name"}=0;
	}
	sub time_set
	{
		my $name=shift;
		my $dt=shift;

		if (!(defined $time{"$name"})) { $time{"total $name"}=$time{"diff $name"}=$dt; }
		else { $time{"total $name"}+=$dt; $time{"diff $name"}=$dt; }
		$time{$name}=-1;
		$time{"running $name"}=0;
	}
	sub time_eval
	{
		my $diff;
		my $name=shift;

		$TimeHiRes or return (-1,-1);
		defined $time{$name} or die "time_eval: using evaluation on undefined entry '$name'\n";
		if ($time{"running $name"})
		{
			time_stop($name);
			time_start($name);
		}
		$diff=$time{"diff $name"}; $time{"diff $name"}=0;
		return ($time{"total $name"},$diff);
	}

	STDOUT->autoflush(1);
	sub debug
	{
		$debug=$_[0];
		if ($debug)
		{
			open(DEBUG,">debug.scr");
			DEBUG->autoflush(1);
		}
		elsif (-e 'debug.scr')
		{
			`rm debug.scr`;
		}
	}
	sub submit
	{
		if ($_[0]=~/quit/) {die "ERROR: quit command is blocked, use 'end'\n";}
	    print IN $_[0];
		if ($debug) 
		{ 
			$_[0]=~/(.*)/;
			print DEBUG "$1;\n"; 
		}
	}
	sub tmpfile
	{
        my $name;
        do { $name=tmpnam(); } until IO::File->new($name,O_RDWR|O_CREAT|O_EXCL);
        $end.="unlink('$name');";
        return $name;
	}
	END
	{
		eval $end;
		close(DEBUG);
	}
	sub out
	{
		my ($line,$ret,$error);
		$ret=""; $error=0;
		while($line=<OUT>)
		{
			if ($line=~/ERROR/) {$error=1;}
			if ($line=~/^EOO$/) {last;}
			if ($_[0] || $error) {print $line;}
			$ret.=$line;
		}
		if ($error) {die "ug aborted due to ERROR\n";}
		return $ret;
	}
}

sub usage_error
{
	my $cmd=shift;
	my $usage=shift;
	die "ERROR in usage of '$cmd':\nusage: '$usage'\n";
}

sub set
{
	print IN "set $_[0]\n";
    return split /[=\s]+/,out(0);
}

sub ug
{
	my (@in,$i,$cmd,$print,$command,$ui,$stat,%argv,$dummy);

	# cancel trailing white spaces
	@in=@_;
	for ($i=0; $i<@in; $i++) { $in[$i]=~s/\s*$//g; }

	# basic check
	if (@in<=0) 
	{
		die "ERROR: provide ug command\n";
	} 
	if ($in[0] eq '') { return; }
	if ($in[0]=~/^\s*\#/) { return; }

	# detect internal print
	$command=$in[0]; 
	$print=0;
	if ($command=~/^\s*print\s+/) 
	{
		$print=1; 
		$command=~s/\s*^print\s+//g;
		if ($command eq "") {die "ERROR: print must come with ug command\n";}
	}

	# check for I/O channels
	if ($command ne "start" && $command ne "running")
	{
		1==stat IN and 1==stat OUT or die "ERROR in '$command': IN/OUT channel missing\n";
		1==stat IN or die "ERROR in '$command': IN channel missing\n";
		1==stat OUT or die "ERROR in '$command': OUT channel missing\n";
	}

	# scan arguments (used for some commands)
	@in%2==1 or die "ERROR: odd number of arguments provided with command '@in'\n";
	($dummy,%argv)=@in; 
	SWITCH:
	{
		# running
        if ($command eq "running")
        {
			if (1==stat IN and 1==stat OUT) { return 1; }
			return 0;
        }

		# debug
		if ($command eq "debug")
		{
			debug $argv{'d'};
			return;
		}

		# command 'end'
		if ($command eq "end")
		{
			if(@in!=1)
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
			@in==3 or usage_error('set','set "v"=><name>');
			print IN "set $argv{'v'}\n";
			return split /[=\s]+/,out($print);
		}

		# command 'start'
		if ($command eq "start")
		{
			my ($exec,$appl,$e,$r,$amp,$model,$a_mode,$m_mode,$p_mode,$r_mode,$pre,@name);
			if(@in!=7)
        	{   
        	    die 'ERROR: usage: ug "start", "p"=>"program", "x"=>[0|1], "n"=><# of procs>;'."\n";
        	} 
			-e $argv{'p'} or die "ERROR: program '$argv{'p'}' does not exist\n";
			$argv{'n'}>0 or die "ERROR: nb of processors out of range\n"; 
			$argv{'x'}==0 || $argv{'x'}==1 or die "ERROR: wrong specification of 'x'-option\n";
            @name=split /\//,$argv{'p'}; $appl=$name[@name-1];
			if ($argv{'x'}==1) { $ui="-ui c"; } else { $ui="-ui cn"; }

			# determine appl-mode
			$model=`strings $argv{'p'} | grep 'Model:'`;
			$a_mode='s'; if ($model=~/parallel/) { $a_mode='p'; }

			# determine module_mode
			$m_mode='s'; if (module('ugp.pm')) { require ugp; $m_mode='p'; $ui="-ui cn"; } 

			# determin procs-mode
			$p_mode='s'; if ($argv{'n'}>1) { $p_mode='p'; }

			$amp=$a_mode.$m_mode.$p_mode; $r_mode=0;
			SWITCH:
			{
				if ($amp eq 'sss')
				{
					# classic sequential
					$pre ="################# start ################\n";
					$pre.="application: $appl\n";
					$pre.="mode: running sequential\n";
					$pre.="########################################\n";
					$exec="$argv{'p'} $ui -perl";
					last SWITCH;
				}
				if ($amp eq 'ssp')
				{
					die "ERROR: cannot run sequential code on $argv{'n'} processors\n";
					last SWITCH;
				}
				if ($amp eq 'sps')
				{
					# run sequential code on parallel machine: one processor
					$pre ="################# start ################\n";
					$pre.="application: $appl\n";
					$pre.="mode: running sequential code on one\n";
					$pre.="      processor of parallel machine\n";
					$pre.="########################################\n";
					($e,$exec)=ugp::run('n'=>$argv{'n'},'p'=>"$argv{'p'} $ui -perl",'b'=>0,'r'=>0);
					last SWITCH;
				}
				if ($amp eq 'spp')
				{
					die "ERROR: cannot run sequential code on $argv{'n'} processors\n";
					last SWITCH;
				}
				if ($amp eq 'pss')
				{
					die "ERROR: cannot run classic parallel on one processor: module 'ugp.pm' missing\n";
					last SWITCH;
				}
				if ($amp eq 'psp')
				{
					die "ERROR: cannot run classic parallel: module 'ugp.pm' missing\n";
					last SWITCH;
				}
				if ($amp eq 'pps')
				{
					# parallel code on one processor
					$pre ="################# start ################\n";
					$pre.="application: $appl\n";
					$pre.="mode: running parallel code on one\n";
					$pre.="      processor of parallel machine\n";
					$pre.="########################################\n";
					($e,$exec)=ugp::run('n'=>$argv{'n'},'p'=>"$argv{'p'} $ui -perl",'b'=>0,'r'=>0);
					last SWITCH;
				}
				if ($amp eq 'ppp')
				{
					# classic parallel
					$pre ="################# start ################\n";
					$pre.="application: $appl\n";
					$pre.="mode: running parallel on $argv{'n'} procs\n";
					$pre.="########################################\n";
					($e,$exec)=ugp::run('n'=>$argv{'n'},'p'=>"$argv{'p'} $ui -perl",'b'=>0,'r'=>0);
					last SWITCH;
				}
			}

			# run ug 
        	open2(*OUT,*IN,$exec);
			IN->autoflush(1); OUT->autoflush(1);
			if ($print) { $pre.=out(0); } else { out(0); }
			return $pre;
		}
	
		# command 'ug' 
		if ($command eq "ug")
        {
			if(@in!=2)
            {
                die "ERROR: command 'ug' must come with one argument\n\n";
            }
			submit "$in[1]\n";
			return out($print);
		}

		# std command
		if (@in%2==1) {$cmd="$command "; $i=1;}
		else {$cmd="$command $in[1] "; $i=2;}
		for (;$i<@in;$i+=2)
		{
			$cmd.='$'."$in[$i] $in[$i+1] ";
		}
		$cmd.="\n";
		submit $cmd;
		return out($print);
	}
}
sub float
{
	my $real='[+-]?\d+\.?\d*[eE]?[+-]?\d+|[+-]?\d*\.?\d+[eE]?[+-]?\d+|[+-]?\d+';
	my (@list,$f,$s,$in);

	if (@_==1) { @list=grep /$real/,split /($real)/,$_[0]; }
	elsif (@_==2) 
	{
		$in=' '.$_[1];
		($f,$s)=split /$_[0]/,$in,2;
		@list=grep /$real/,split /($real)/,$s;
	}
	else
	{
		$in=' '.$_[2];
        ($f,$s)=split /$_[0]/,$in,2;
        @list=grep /$real/,split /($real)/,$s;
		if ($_[1]>=@list || $_[1]<0) { return undef;}
		for ($s=0; $s<$_[1]; $s++) { shift @list; } 
    }
	return wantarray ? @list : $list[0];
}

1;




