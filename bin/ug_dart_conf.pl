#!/usr/bin/perl
#
# ug_dart_conf.pl
#
# purpose: build the "DartConfiguration.tcl" file which includes the 
#          configuration for the dart test cycle

###############################################################################
# modules
###############################################################################
use Getopt::Std;

###############################################################################
# defining several subroutines
###############################################################################

# creating the "DartRoot" dart variable 
sub dartroot # argument: dartroot directory
{
	return(join(' ', 'DartRoot:', $_[0]));
}

# creating the "SourceDirectory" dart variable
sub sourcedir # argument: test directory
{
    return(join('', 'SourceDirectory: ', $_[0]));
}

# creating the "BuildDirectory" dart variable
sub builddir #argument: test directory
{
    return(join('', 'BuildDirectory: ', $_[0]));
}

# creating the "Site" dart variable 
sub site
{
	my $site;
	system("uname -n > ./temp.txt");
	open TEMP, "./temp.txt";
	while (<TEMP>)
	{
    		$site=substr($_, 0, (index $_, "\n"));
	}
	system("rm ./temp.txt");
	return(join(' ','Site:',$site));
}

# creating the "BuildName" dart variable 
sub buildname # arguments: test identifier, ugroot directory, architecture  
{
	my ($string,$kernel_name,$os_release,$compiler);
  $string = '[-+A-Za-z]\w*';
  # setting the compiler
  open ARCHFILE, (join('', $_[1], '/arch/', $_[2], '/mk.arch'))
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
    		$kernel_name=substr($_, 0, (index $_, "\n"));
	}
	system("uname -r > ./temp.txt");
	open TEMP, "./temp.txt";
	while (<TEMP>)
	{
    		$os_release=substr($_, 0, (index $_, "\n"));
	}
	system("rm ./temp.txt");
	return(join('', 'BuildName: ', $_[3], '-', $_[0], '-', $kernel_name, '-', $os_release, '-', $compiler));
}

# creating the "ConfigureCommand" dart variable 
sub confcomm # arguments: architecture, ugconf parameter
{
    my $ConfComm=join('', 'ConfigureCommand: ugconf -v', ' ', $_[0], ' ', $_[1]);
    return($ConfComm);
}

# creating the "MakeCommand" dart variable
sub makecomm
{
    my $makecomm = "MakeCommand: make build";
    return $makecomm;
} 

# print the help message
sub help
{
    print "usage:	ug_dart_conf.pl [-b build_dir_suffix] [-m mode] [-i test id]\n"; 
    print "[-a architecture] [-c ugconf options]\n";
    print "-b: Specifying another build directory than .../UG by entering the\n";
    print "    sub directory, e.g. cd/appl\n"; 
    print "-i: Test identifier\n";
    print "-d: Build directory identifier\n";
    print "-a: Architecture, e.g. PC\n";
    print "-c: Options which ugconf would accept except the architecture\n";
    print "\n";
    print "purpose: running all scripts which are necessary for a complete dart\n";
    print "         build/test/submit cycle\n";
}

###############################################################################
# main
###############################################################################

# read the command line parameters
%option = ();
getopts("b:i:a:c:d:h", \%option);

# print help message if the parameter -h is given
if($option{h})
{
    help();
    exit 0;
}

# exit if no architecture is given
unless($option{a})
{
    die "No architecture specified!";
}

# set up some variables

# dartroot
my $dartroot = $ENV {"DART_HOME"};

# ugroot
my $ugroot = join('', $dartroot, '/Source/Client/UG/ug');

# write DartConfiguration.tcl
# open the "DartConfiguration.tcl.proto" file for reading
open IN, (join('', $dartroot, '/Source/Client/UG/DartConfiguration.tcl.proto'));

# open the "DartConfiguration.tcl" file for writing
open OUT, (join('', '>', $option{b}, '/DartConfiguration.tcl'));

# copy the content of IN to OUT and write the dart variables into OUT
while (<IN>)
{
    if ( $_ =~ /DartRoot: ---/)
    {
        print OUT dartroot($dartroot),"\n";
    }
    elsif ( $_ =~ /SourceDirectory: ---/)
    {
        print OUT sourcedir($option{b}),"\n";
    }
    elsif ( $_ =~ /BuildDirectory: ---/)
    {
        print OUT builddir($option{b}),"\n";
    }
    elsif ( $_ =~ /BuildName: ---/)
    {
        print OUT buildname($option{i},$ugroot,$option{a},$option{d}),"\n";
    }
    elsif ( $_ =~ /ConfigureCommand: ---/)
    {
        print OUT confcomm($option{a},$option{c}),"\n";
    }
    elsif ( $_ =~ /Site: ---/)
    {
        print OUT site(),"\n";
    }
    elsif ( $_ =~ /MakeCommand: ---/)
    {
        print OUT makecomm($option{i}),"\n";
    } 
    else
    {
        print OUT $_;
    }
}





