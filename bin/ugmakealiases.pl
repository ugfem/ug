#!/usr/bin/perl -w

# ignore these directories
%exact = ( "CVS" => 1 );
$pattern = "";

for ($i=0; $i<@ARGV; $i++) {
	
	if ( $ARGV[$i] eq "-help") {
		print "purpose: make aliases file 'ugaliases' for  easy changing\n";
		print "    into ug's directories; e.g. instead of typing on\n";
		print "    your shell command line:\n";
		print "      cd \$UGROOT/gm/gg3\n";
		print "    you can use:\n";
		print "      cdgg3\n";
		print "     (note that no space is used)\n";
		print "     To activate the 'ugaliases' file in your shell\n";
		print "     type on your shell command line:\n";
		print "         source \$UGROOT/bin/ugaliases\n";
		print "     you can also put this line into your shell\n";
		print "     resource file (e.g. ~/.tcshrc) to use the\n";
		print "     shortcutting permanently.\n";
		print_usage();
		exit (0);
	}
	elsif ( $ARGV[$i] eq "-e" ) {
		$exact{"$ARGV[$i+1]"} = 1;
	}
	elsif ( $ARGV[$i] eq "-p" ) {
		if ($pattern eq "") { $pattern = "$ARGV[$i+1]"; }
		else { $pattern = "$pattern|$ARGV[$i+1]"; }
	}
}
	
$ENV{'UGROOT'} or die "environment variable 'UGROOT' not defined\n";
$ugroot = $ENV{'UGROOT'};

# determine shell
$ENV{'SHELL'} or die "couldn't determine what shell you are using.\n";
$shell = $ENV{'SHELL'};

# alias syntax is either csh or sh
if    ( $shell =~ "tcsh" ) { $shell = "csh"; }
elsif ( $shell =~ "csh"  ) { $shell = "csh"; }
elsif ( $shell =~ "bash" ) { $shell = "sh";  }
else  { $shell = "sh"; }  # default is sh


$ugaliasesfile = $ugroot . "/bin/ugaliases";

system("rm -f $ugaliasesfile");

# global hash of aliases
%alias = ();

chdir $ugroot or die "Couldn't change dir to $ugroot";

# create aliases for applications
app_aliases();

# create aliases for ug
ug_aliases();

write_aliasfile();

print "New aliases file ($ugaliasesfile) created.\n";

exit 0;


######### subroutines ############################################

sub app_aliases {

	my $UGdir = $ugroot ."/..";

	opendir UGDIR, $UGdir or die "Couldn't open dir to $ugroot/..";

	my @UGfiles = readdir UGDIR;
	closedir UGDIR;

	foreach $file (@UGfiles) {
		
		# don't work in ug,  .,  .. and specified directories 
		if ( $file eq "ug" || exists($exact{$file}) ||
			 $file =~ /$pattern/ || substr($file,0,1) eq "." ) { 
			next;
		}
		
		$fullpath = "../" . $file;
		if ( -d $fullpath ) {
			my $v;
			($v = $fullpath ) =~ s/[^\w\d]//g;
			if ( !exists($alias{$v}) ) {
				$alias{$v} = "cd \$UGROOT/../" . $file;
			}
			else {
				print "    name clash: $v is either \"cd $fullpath\" or \"$alias{$v}\"\n";
			}
			work_on_directory($file,"..");
		}
	}
}



sub work_on_directory {

	# $path is the path from $ugroot
	# $dir is the directory name
	local($dir,$path) = @_;

	my $pathto;
	if ( $path eq "" )  { $pathto = $dir; }
	else                { $pathto = $path . "/" . $dir; }

    opendir DIR, $pathto or die "Couldn't open dir $pathto";

    local @files = readdir DIR;
    closedir DIR;

	my $file;
    foreach $file (@files) {
	    
		my $thefile = $pathto . "/" . $file;

        # don't use directories starting with . or those specified with -n
        if ( exists($exact{$file}) || $file =~ /$pattern/ || substr($file,0,1) eq "." ) {
			next;
        }

        if ( -d $thefile ) {
			my $v;
			($v = $thefile ) =~ s/[^\w\d]//g;
            if ( !exists($alias{$v}) ) {
				$alias{$v} = "cd \$UGROOT/$thefile";
			}
			else {
				print "    name clash: $v is either \"cd $thefile\" or \"$alias{$v}\"\n";
			}
            work_on_directory($file,$pathto);
        }
    }
}


sub ug_aliases { 
	work_on_directory(".","");
}

sub write_aliasfile {

	open(UGALIAS,"> $ugaliasesfile");
	my $key;
	
	if ( $shell eq "sh") {
		foreach $key (sort keys %alias)  {
			print UGALIAS "alias cd$key=\"$alias{$key}\"\n";
		}
	}
	elsif ( $shell eq "csh" ) {
		foreach $key (sort keys %alias)  {
			print UGALIAS "alias cd$key \"$alias{$key}\"\n";
		}
	}
	close UGALIAS;
}

sub print_usage {
	print "\nusage: ugmakaliases [-help] [-e name] [-p pattern]\n";
	print "-help: prints this help infromation\n";
	print "-e name: Don't work on directories calleded name.\n";
	print "-p pattern: Don't work on directories whose name contains pattern.\n";

	print "Example:\n    ugmakealiases -e test -e dummy -p .bak\n";
	print "In this example, directories called 'test' or 'dummy' are ignored,\n";
	print "as are direktories containing the string '.bak'\n";
}


