#!/usr/bin/perl
#
# ug_dart_conf.pl
#
# purpose: build the "DartConfiguration.tcl" file which includes the 
#          configuration for the dart test cycle


##############################################################################
# defining several subroutines
##############################################################################


sub config
{
	my ($ConfParams,$BuildParams) = ("","");
	my @control_params;
	my $arch_type = splice(@_,0,1);
	for(my $i=0;$i<=13;$i++)
	{
		$control_params[$i]=0;
	}
	foreach $param (@_)
	{
    		if($param eq "SEQ" || $param eq "MPI" || $param eq "PVM" ||
		   $param eq "NX" || $param eq "NXLIB" || $param eq "SHMEM" ||
		   $param eq "SHMEMT3D" || $param eq "PARIX")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "SEQ")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}
			$control_params[0]=1;
		}
		elsif($param eq "2" || $param eq "3")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "2")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[1]=1;
		}
		elsif($param eq "GRAPE" || $param eq "NOGRAPE")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOGRAPE")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[2]=1;
		}
		elsif($param eq "COVISE" || $param eq "NOCOVISE")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOCOVISE")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[3]=1;
		}
		elsif($param eq "PV3" || $param eq "NOPV3")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOPV3")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[4]=1;
		}
		elsif($param eq "NETGEN" || $param eq "NONETGEN")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NONETGEN")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[5]=1;
		}
		elsif($param eq "RIF" || $param eq "NORIF")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NORIF")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[6]=1;
		}
		elsif($param eq "SIF" || $param eq "XIF" || $param eq "MIF" ||
		       $param eq "RIF" || $param eq "XRIF" || $param eq "NORIF")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "XIF")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[7]=1;
		}
		elsif($param eq "GUI" || $param eq "NOGUI")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOGUI")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[8]=1;
		}
		elsif($param eq "STD_DOMAIN" || $param eq "LGM_DOMAIN")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "STD_DOMAIN")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[9]=1;
		}
		elsif($param eq "DEBUG" || $param eq "NODEBUG")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "DEBUG")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[10]=1;
		}
		elsif($param eq "OPTIM" || $param eq "NOOPTIM")
		{
			if($arch_type =~ /GCOV/)
			{
				$ConfParams=join(' ',$ConfParams,"NOOPTIM");
				print "\n\nOptimimization can't be used with code coverage!\n\n";
			}
			else
			{
				$ConfParams=join(' ',$ConfParams,$param);
			}
			unless($param eq "NOOPTIM")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[11]=1;
		}
		elsif($param eq "CHACO" || $param eq "NOCHACO")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOCHACO")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[12]=1;
		}
		elsif($param eq "CAD" || $param eq "NOCAD")
		{
			$ConfParams=join(' ',$ConfParams,$param);
			unless($param eq "NOCAD")
			{
				$BuildParams=join('-',$BuildParams,$param);
			}	
			$control_params[13]=1;
		}
	}
	for(my $i=0;$i<=13;$i++)
	{
		if($control_params[$i]==0)
		{
			if($i==0)
			{
				$ConfParams=join(' ',$ConfParams,'SEQ');
			}				
			elsif($i==1)
			{
				$ConfParams=join(' ',$ConfParams,'2');
			}				
			elsif($i==2)
			{
				$ConfParams=join(' ',$ConfParams,'NOGRAPE');
			}				
			elsif($i==3)
			{
				$ConfParams=join(' ',$ConfParams,'NOCOVISE');
			}				
			elsif($i==4)
			{
				$ConfParams=join(' ',$ConfParams,'NOPV3');
			}				
			elsif($i==5)
			{
				$ConfParams=join(' ',$ConfParams,'NONETGEN');
			}				
			elsif($i==6)
			{
				$ConfParams=join(' ',$ConfParams,'NORIF');
			}				
			elsif($i==7)
			{
				$ConfParams=join(' ',$ConfParams,'XIF');
			}				
			elsif($i==8)
			{
				$ConfParams=join(' ',$ConfParams,'NOGUI');
			}				
			elsif($i==9)
			{
				$ConfParams=join(' ',$ConfParams,'STD_DOMAIN');
			}				
			elsif($i==10)
			{
				$ConfParams=join(' ',$ConfParams,'DEBUG');
			}				
			elsif($i==11)
			{
				if($arch_type =~ /GCOV/)
				{
					$ConfParams=join(' ',$ConfParams,'NOOPTIM');
				}
				else
				{
					$ConfParams=join(' ',$ConfParams,'OPTIM');
				}
			}				
			elsif($i==12)
			{
				$ConfParams=join(' ',$ConfParams,'NOCHACO');
			}				
			elsif($i==13)
			{
				$ConfParams=join(' ',$ConfParams,'NOCAD');
			}				
		}
	}
	if($BuildParams eq "")
	{
		$BuildParams="-default";
	}
	my @output;
	push(@output,$ConfParams);
	push(@output,$BuildParams);
	return(@output);
}

# creating the "DartRoot" dart variable for "DartConfiguration.tcl"
sub dartroot
{
	return(join(' ','DartRoot:',$_[0]));
}

# creating the "SourceDirectory" dart variable for "DartConfiguration.tcl"
sub sourcedir
{
	return(join('','SourceDirectory: ',$_[0],'/Source/',$_[1],'/UG'));
}

# creating the "BuildDirectory" dart variable for "DartConfiguration.tcl"
sub builddir
{
	return(join('','BuildDirectory: ',$_[0],'/Source/',$_[1],'/UG'));
}

# creating the "Site" dart variable for "DartConfiguration.tcl"
sub site
{
	my $site;
	
	system("uname -n > ./temp.txt");
	open TEMP, "./temp.txt";
	while (<TEMP>)
	{
    		$site=substr $_, 0, (index $_, "\n");
	}
	system("rm ./temp.txt");
	return(join(' ','Site:',$site));
}

# creating the "BuildName" dart variable for "DartConfiguration.tcl"
sub buildname
{
	my ($string,$kernel_name,$os_release,$compiler);
        $string = '[A-Za-z]\w*';
        # setting the compiler
        open ARCHFILE, (join('',$_[0],'/arch/',$_[1],'/mk.arch'))
         or die "\nArchitecture doesn't exist!\n";
        while (<ARCHFILE>)
        {
        	if ( $_ =~ /ARCH_CC\s*=\s*($string)/ )
                {
			$compiler = $1;
		}
        }
	system("uname -s > ./temp.txt");
	open TEMP, "./temp.txt";
	while (<TEMP>)
	{
    		$kernel_name=substr $_, 0, (index $_, "\n");
	}
	system("uname -r > ./temp.txt");
	open TEMP, "./temp.txt";
	while (<TEMP>)
	{
    		$os_release=substr $_, 0, (index $_, "\n");
	}
	system("rm ./temp.txt");
	return(join('','BuildName: ',$kernel_name,'-',$os_release,'-',$compiler,'-',$_[1],,$_[2]));
}

# creating the "ConfigureCommand" dart variable for "DartConfiguration.tcl"
sub confcomm
{
	my $ConfComm=join('','ConfigureCommand: ',$_[0],'/bin/ugconf',' ',$_[1],$_[2]);
	return($ConfComm);
}

##############################################################################
# main 
##############################################################################

# setting the "dartroot" variable
my $dartroot = $ENV {"DART_HOME"};

# detect whether this script runs on the server or the client(default: client)
my $mode = "Client";
if($#ARGV > 0 && ($ARGV[0] eq "Client" || $ARGV[0] eq "Server"))
{
	$mode = splice(@ARGV,0,1);
}

# setting the "ugroot" variable
my $ugroot = join('',$dartroot,'/Source/',$mode,'/UG/ug');

# setting the arch variable
my $arch = splice(@ARGV,0,1);

# generating the parameters for ugconf and for the BuildName variable
my @config_output=config($arch,@ARGV);
my $ugconf_params=$config_output[0];
my $build_params=$config_output[1];

# open the "DartConfiguration.tcl.proto" file for reading
open IN, (join('',$dartroot,'/Source/',$mode,'/UG/DartConfiguration.tcl.proto'));

# open the "DartConfiguration.tcl" file for writing
open OUT, (join('','>',$dartroot,'/Source/',$mode,'/UG/DartConfiguration.tcl'));

# copy the content of IN to OUT and write the dart variables into OUT
while (<IN>)
{
    if ( $_=~ /DartRoot: ---/)
    {
	print OUT dartroot($dartroot),"\n";
    }
    elsif ( $_=~ /SourceDirectory: ---/)
    {
	print OUT sourcedir($dartroot,$mode),"\n";
    }
    elsif ( $_=~ /BuildDirectory: ---/)
    {
	print OUT builddir($dartroot,$mode),"\n";
    }
    elsif ( $_=~ /BuildName: ---/)
    {
	print OUT buildname($ugroot,$arch,$build_params),"\n";
    }
    elsif ( $_=~ /ConfigureCommand: ---/)
    {
        print OUT confcomm($ugroot,$arch,$ugconf_params),"\n";
    }
    elsif ( $_=~ /Site: ---/)
    {
        print OUT site(),"\n";
    }
    else
    {
        print OUT $_;
    }
}
