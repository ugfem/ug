#!/usr/bin/perl -w

#################################################################################################
# Replacement for ugmake shell script.
# Tries to be fast by using parallel make.
#
# Author:  Volker Reichenberger
#          volker.reichenberger@iwr.uni-heidelberg.de
#
# History: September 27, 2000 begin
#################################################################################################



$debug = 0;

@ugmods = ( "ug", "dev", "meta", "xif", "gm", "sif", "mif",
"graphics", "low", "np", "ui", "gg2", "gg3", "dom", "std", "lgm", "gen", "netgen",
	"diff2d", "diff2da", "cd", "cda", "fem", "fema", "ns", "ns2d", "ns3d", "simple",
	"scalar", "tools"
);

# if '-help' appears on the command line, we show the help, and do nothing else
foreach $arg (@ARGV)	{
	if ($arg eq "-help")	{
		print_usage();
		exit(0);
	}
}

#################################################################################################
# determine ug location
$ENV{'UGROOT'} or die "[ugmake.pl] environment variable 'UGROOT' not defined\n";
$ugroot = $ENV{UGROOT};
if ($debug) { print "UGROOT is $ugroot\n"; }

$ugconf=$ugroot.'/ug.conf';


#################################################################################################
# determine architecture, this defines $arch (e.g. PC, MACOSX)
determine_architecture();

# determine dimension, this defines $dim (2 or 3)
determine_dimension();

# determine parallel, this defines $parallel (0 or 1)
determine_parallel();

# do we want CHACO?
$chaco = 0;
if ($parallel) { determine_chaco(); }




#################################################################################################
# get appropriate make, store it in $make
$mkarch=$ugroot.'/arch/'.$arch.'/mk.arch';
if ($debug) {print "Reading architecture info from $mkarch\n";}
open(MKARCH,$mkarch);
$foundmake = $foundpmake = 0;
while(<MKARCH>)
{
	$line = $_;
	if ( $line =~ /^[ \t]*(ARCH_MAKE)([ \t]*=[ \t]*)([^\n\r#]*)[ \t]*$/ )	{
		$make = $3;
		$foundmake = $1;
	}
	if ( $line =~ /^[ \t]*(ARCH_PMAKE)([ \t]*=[ \t]*)([^\n\r#]*)[ \t]*$/ )	{
		$pmake = $3;
		$foundpmake = 1;
	}
}
close(MKARCH);

$standardmake = "make";
if ($foundpmake) { $standardmake = $make; $make = $pmake; }
else { if (!$foundmake) { print "Couldn't find ARCH_MAKE in $mkarch, using 'make' instead."; } }

if ($debug) {print "'make' is $make\n";}


#################################################################################################
# check options

$globaloptions = "";
@modules = ();

print "*****\n";
print "@ARGV\n";
for ($i=0; $i<@ARGV; $i++) { print "$ARGV[$i]\n"; } 
print "*****\n";

for ($i=0; $i<@ARGV; $i++) {
	$mfound=0;
	foreach $ugmod (@ugmods)	{
		# if the argument is a ug module name, append it to @modules
		if ($ARGV[$i] eq $ugmod) { 
			@modules = ( @modules, $ARGV[$i] );
			$mfound=1;
		}
	}
	if (!$found) { $i++; last; }
	$found = 0;
}

# the rest of the command line is considered as an additional
# make option

$makeclean = 0;	# the $makeclean variable is used below.  Some ug files are
				# explicitly compiled (initug.c and initnp.c), because they 
				# are sometimes forgotten. (TODO: why? has probably to do with Makefile)
				# If $makeclean is on, we can avoid compiling them.
for ( ; $i<@ARGV; $i++)	{
    $globaloptions = $globaloptions." ".$ARGV[$i];
	if ($ARGV[$i] eq "clean") { $makeclean=1; }
}

if ( $debug )	{
	print "Building modules ", commify_series(@modules),"\n";
	print "Global options '$globaloptions'\n";
}

#################################################################################################
# say what we are doing
print "[ugmake.pl] Building for $arch $model, dimension $dim\n";



#################################################################################################
# if no argument was given look for Makefile
if (@ARGV==0) {
	@makefiles = <Makefile*>; 
	if (@makefiles>1) {
		print "Makefile not found or not unique (", commify_series(@makefiles),")\n";
		exit(0);
	}
	print "$make -f @makefiles\n";
	system("$make -f @makefiles")  == 0
		or die "[ugmake.pl] Couldn't make in $ENV{PWD}\n";
	exit(0);
}

#################################################################################################
# In these directories we can savely issue parallel make

%simple = (
		"meta" => "dev/meta",
		"nif" => "dev/nif",
		"ppm" => "dev/ppm",
		"ps" => "dev/ps",
		"rif" => "dev/rif",
		"sif" => "dev/sif",
		"xif" => "dev/xif",
		"gen" => "dom/gen",
		"ngin" => "dom/lgm/ngin",
		"ngin2d" => "dom/lgm/ngin2d",
		"std" => "dom/std",
		"gg2" => "gm/gg2",
		"gg3" => "gm/gg3",
		"netgen" => "gm/gg3/netgen",
		"covise" => "graphics/covise",
		"grape" => "graphics/grape",
		"uggraph" => "graphics/uggraph",
		"low" => "low",
		"algebra" => "np/algebra",
		"famglib" => "np/famglib",
		"field" => "np/field",
		"procs" => "np/procs",
		"udm" => "np/udm",
		"assign" => "parallel/chaco/assign",
		"bpmatch" => "parallel/chaco/bpmatch",
		"coarsen" => "parallel/chaco/coarsen",
		"connect" => "parallel/chaco/connect",
		"eigen" => "parallel/chaco/eigen",
		"graph" => "parallel/chaco/graph",
		"inertial" => "parallel/chaco/inertial",
		"input" => "parallel/chaco/input",
		"klspiff" => "parallel/chaco/klspiff",
		"main" => "parallel/chaco/main",
		"misc" => "parallel/chaco/misc",
		"optimize" => "parallel/chaco/optimize",
		"symmlq" => "parallel/chaco/symmlq",
		"util" => "parallel/chaco/util",
		"dddif" => "parallel/dddif",
		"dddobj" => "parallel/dddobj",
		"ppif" => "parallel/ppif",
		"tools" => "tools",
		"ui" => "ui"
);

### Use normal make here
%standard = (
	"ns" => "../ns/pclib",
	"ns3d" => "../ns/appl3d",
	"ns2d" => "../ns/appl2d",
	"diff2d" => "../diff2d/pclib",
	"diff2da" => "../diff2d/appl",
	"cd" => "../cd/pclib",
	"cda" => "../cd/appl",
	"fem" => "../fem/pclib",
	"fema" => "../fem/appl",
	"tools" => "tools",
	"simple" => "../simple"
);

if ($debug && 0) {
	print "directories without Makefile dependencies into other directories\n";
	while (($key,$value) = each(%simple)) {
		print "modlue name $key, directory $value\n";
	}
}


# check for include directory
if ( ! -e "$ugroot/include" ) {
	print "[ugmake.pl] making include file links\n";
	system("ugmakelinks") == 0 or die "call to ugmakelinks failed\n";
}

#################################################################################################
# process the arguments
foreach $module (@modules) {
	# for some cases, we just call the old ugmake
	if ($module eq "links") {
		system("ugmakelinks") == 0 or die "Call to ugmakelinks failed\n";
	}
	# the default case is to make the module
	else {
		make_module($module);
	}
}

# and that's it
exit(0);

#################################################################################################
#################################################################################################


#################################################################################################
###  subroutines     ############################################################################
#################################################################################################

sub make_module {
	my ($module) = @_;
	print "[ugmake.pl] making ug module $module \n";
	
	# for these modules we can just call make in the directory
	if (exists($simple{$module}))	{
		print "cd $ugroot/$simple{$module}\n";
		@makef = <$ugroot/$simple{$module}/Makefile*>;
		print "$make -f @makef $globaloptions\n";
		system("cd $ugroot/$simple{$module}; $make -f @makef $globaloptions") == 0
			or die "[ugmake.pl] Couldn't make $module\n";
	}
	elsif (exists($standard{$module}))	{
		print "cd $ugroot/$standard{$module}\n";
		@makef = <$ugroot/$standard{$module}/Makefile*>;
		@makef == 1 or die "More than one Makefile found in $ugroot/$standard{$module}\n";
		print "$standardmake -f @makef $globaloptions\n";
		system("cd $ugroot/$standard{$module}; $standardmake -f @makef $globaloptions") == 0
			or die "[ugmake.pl] Couldn't make $module\n";
	}
	# for these modules we need special care
	else {
		if ($module eq "ug")	{
			make_dev();
			make_dom();
			make_gm();
			make_graphics();
			make_module("low");
			make_np();
			make_module("ui");
			if ($parallel) {
				make_parallel();
			}
			if (!$makeclean) {
				system("cd $ugroot; $make initug.o") == 0
					or die "[ugmake.pl] Couldn't make initug.o\n";
			}

		}
		elsif ($module eq "gm")	{
			make_gm();
		}
		elsif ($module eq "np")	{
			make_np();
		}
		elsif ($module eq "parallel")	{
			make_parallel();
		}
		elsif ($module eq "graphics")	{
			make_graphics();
		}
		else {
			print "Unknown ug module $module\n";
		}
	}
	print "[ugmake.pl] done making ug module $module \n";
}

sub make_in_directory {
	my ($dir) = @_;
	print "[ugmakepl] cd $ugroot/$dir; $make $globaloptions\n";
	system("cd $ugroot/$dir; $make $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make in $ugroot/$dir\n";
}

sub make_dev {
	make_module("meta");
	make_module("ppm");
	make_module("ps");
	print "[ugmakepl] cd $ugroot/dev; $make -f Makefile.dev $globaloptions\n";
	system("cd $ugroot/dev; $make -f Makefile.dev $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make dev\n";
}

# TODO: Use parallel make
sub make_dom {
	print "[ugmakepl] cd $ugroot/dom; $make -f Makefile.dom $globaloptions\n";
	system("cd $ugroot/dom; $make -f Makefile.dom $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make dom\n";
}

sub make_gm {
	make_module("gg$dim");
	print "[ugmakepl] cd $ugroot/gm; $make -f Makefile.gm $globaloptions\n";
	system("cd $ugroot/gm; $make -f Makefile.gm $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make gm\n";

}

sub make_graphics {
	make_module("uggraph");
	print "[ugmakepl] cd $ugroot/graphics; $make -f Makefile.graphics $globaloptions\n";
	system("cd $ugroot/graphics; $make -f Makefile.graphics $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make graphics\n";
}

sub make_np {
	make_module("algebra");
	make_module("field");
	make_module("procs");
	make_module("udm");
	print "[ugmakepl] cd $ugroot/np/amglib; $make -f Makefile.amglib $globaloptions\n";
	system("cd $ugroot/np/amglib; $make -f Makefile.amglib $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make amglib\n";
	if (!$makeclean)	{
		print "[ugmakepl] cd $ugroot/np; $make -f Makefile.np initnp.o $globaloptions\n";
		system("cd $ugroot/np; $make -f Makefile.np initnp.o $globaloptions") == 0
			or die "[ugmake.pl] Couldn't make initnp.o\n";
	}
	print "[ugmakepl] cd $ugroot/np; $make -f Makefile.np $globaloptions\n";
	system("cd $ugroot/np; $make -f Makefile.np $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make np\n";
}

sub make_parallel	{
	make_module("dddif");
	make_module("dddobj");
	make_module("ppif");
	if ($chaco) { make_chaco(); }
	print "[ugmakepl] cd $ugroot/parallel; $make -f Makefile.parallel $globaloptions\n";
	system("cd $ugroot/np/parallel; $make -f Makefile.parallel $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make parallel\n";
}

sub make_chaco {
	make_module("assign");
	make_module("assign");
	make_module("bpmatch");
	make_module("coarsen");
	make_module("connect");
	make_module("eigen");
	make_module("graph");
	make_module("inertial");
	make_module("input");
	make_module("klspiff");
	make_module("main");
	make_module("misc");
	make_module("optimize");
	make_module("symmlq");
	make_module("util");
	print "[ugmakepl] cd $ugroot/parallel/chaco; $make -f Makefile.chaco $globaloptions\n";
	system("cd $ugroot/np/parallel/chaco; $make -f Makefile. chaco $globaloptions") == 0
		or die "[ugmake.pl] Couldn't make chaco \n";
}

#################################################################################################

sub determine_architecture	{
    if ($debug) {print "Reading ARCH variable from $ugconf\n";}
    open(UGCONF,$ugconf);
    $arch = "NONE";
    while(<UGCONF>)
    {
        $line = $_;
        if ( $line =~ /^[ \t]*(ARCH)([ \t]*=[ \t]*)(\w*)[ \t]*$/ )	{
            $arch = $3;
        }
    }
    close(UGCONF);
    if ($debug) {print "ARCH is $arch\n";}
}

sub determine_dimension	{
    if ($debug) {print "Reading DIM variable from $ugconf\n";}
    open(UGCONF,$ugconf);
    $dim = 0;
    while(<UGCONF>)
    {
        $line = $_;
        if ( $line =~ /^[ \t]*(DIM)([ \t]*=[ \t]*)([23])[ \t]*$/ )	{
            $dim = $3;
        }
    }
    close(UGCONF);
    if ($debug) {print "DIM is $dim\n";}
}

sub determine_parallel	{
    if ($debug) {print "Check if we need parallel support (looking in $ugconf)\n";}
    open(UGCONF,$ugconf);
    $model = "SEQ"; # default
    while(<UGCONF>)
    {
        $line = $_;
        if ( $line =~ /^[ \t]*(MODEL)([ \t]*=[ \t]*)([23])[ \t]*$/ )	{
            $model = $3;
        }
    }
    if ( $model eq "SEQ" ) { $parallel=0; }
    else { $parallel = 1; }
    close(UGCONF);
    if ($debug) {
        if ($parallel) { print "parallel support is on\n"; }
        else { print "parallel support is off\n"; }
    }
}

sub determine_chaco	{
    if ($debug) {print "Check if we need CHACO (looking in $ugconf)\n";}
    open(UGCONF,$ugconf);
    $chaco = 0;
    while(<UGCONF>)
    {
        $line = $_;
        if ( $line =~ /^[ \t]*(CHACO)([ \t]*=[ \t]*)([23])[ \t]*$/ )	{
            $tmp = $3;
        }
    }
    if ( $tmp eq "OFF" ) { $chaco=0; }
    else { $chaco = 1; }
    close(UGCONF);
    if ($debug) {
        if ($chaco) { print "CHACO is onn\n"; }
        else { print "CHACO is off\n"; }
    }
} 

sub print_usage {
	print "Usage: $0 [ug modules] [make arguments]\n";
	print "       $0 will call the appropriate actions to build\n";
	print "       the specified ug module. If called without arguments, $0\n";
	print "       will make the ug module in the current directory.\n";
	print "       After the list of valid ug modules, a list of arguments can be\n";
	print "       specified, which will be passsed to 'make'.\n";
	print "       \n";
	print "       Valid ug modules are:\n       ";
	for ($i=0; $i<@ugmods; $i++)	{
		print "$ugmods[$i] ";
		if ($i%10==9) { print "\n       "; }
	}
	print "\n";
}



#################################################################################################
# commify_series: Used for separating array entries by commas
# from Perl Cookbook, pp 93

sub commify_series {
	(@_ == 0) ? ''                :
	(@_ == 1) ? $_[0]             :
	(@_ == 2) ? join(" and ",@_)  :
	           join(", ", @_[0 .. ($#_-1)], "and $_[-1]");
}

exit(0);
