#!/usr/bin/perl
# ug_dart_tests.pl
#
# findind "DartTestfile.txt" files which declare tests
# and create all other needed "DartTestfile.txt" files

###############################################################################
# modules
###############################################################################
use File::Find ();
use Getopt::Std;

###############################################################################
# defining several subroutines
###############################################################################
sub help
{
    print "usage:	ug_dart.pl [-b build_dir_suffix] [-t test_type] [-m mode]\n";
    print "[-p test application] [-o test application options] [-i test id]\n"; 
    print "[-a architecture] [-c ugconf options] [-e test extent]\n";
    print "-b: Specifying another build directory than .../UG by entering the\n";
    print "    sub directory, e.g. cd/appl\n";
    print "-s: Source directory\n"; 
    print "-p: Specifying the application to test\n";
    print "-i: Test identifier\n";
    print "\n";
    print "purpose: running all scripts which are necessary for a complete dart\n";
    print "         build/test/submit cycle\n";
}

###############################################################################
# read the command line parameters
###############################################################################
%option = ();
getopts("b:s:p:i:h", \%option);

###############################################################################
# print help message if the parameter -h is given
###############################################################################
if($option{h})
{
    help();
    exit 0;
}

###############################################################################
# set default values
###############################################################################
# test application
unless($option{p})
{
    die "No test application given!\n";
}
# test identifier
unless($option{i})
{
    $option{i} = $option{p};
}

###############################################################################
# set up some variables
###############################################################################
# Dart directory
my $dartdir = $ENV {"DART_HOME"};
# start directory 
my $startdir = $option{b};
# test directory 
my $testdir = $option{s}; 

###############################################################################
# generate the DartTestfile.txt file that defines the tests
###############################################################################
unless($option{p} eq "()")
{
    # open the file
    open(ROOT_TESTFILE, join('','>',$option{b},'/DartTestfile.txt'));
    #write to the file
    print ROOT_TESTFILE join('','ADD_TEST(',$option{i},' ',$option{p},')');
    close(ROOT_TESTFILE);
}
###############################################################################
# create the other DartTestfile.txt files only if it's necessary
###############################################################################
unless($option{s} eq $option{b})
{ 
# saving the length of "startdir"
my $i = 0;
foreach $byte (split //,$startdir)
{
    $i++;
}
# finding all subdirectories of "$startdir" and saving them in "@dirs"
sub find(&@) { &File::Find::find }
*name = *File::Find::name;
find { push(@dirs, $name) if -d } $startdir;
# finding all directories whith a "DartTestfile.txt" file in it and saving them
# in "@temptestdirs"
foreach $dir (@dirs)
{
    opendir(DIR,$dir);
    while(defined($file = readdir(DIR)))
    {
        if($file eq "DartTestfile.txt")
        {
            push(@temptestdirs, substr($dir,$i+1));
        }
    }
    closedir(DIR);
}

# saving only the whole path to the "DartTestfile.txt" files defined in "@testdirs" 
my $check = 0;
foreach $dir1 (@temptestdirs)
{
    foreach $dir2 (@temptestdirs)
    {
        if($dir2 =~/$dir1/)
        {
            $check++;
        }
    }
    if($check < 2)
    {
        push(@testdirs, $dir1);
    }
    $check = 0;
}

# creating the input of the "DartTestfile.txt" file which resides in the 
# build directory and creating the other needed files
my $testfile_input = "SUBDIRS(";
foreach $dir (@testdirs)
{
    my @part_dirs = split '/', $dir;
    my $j = 0;
    foreach $part_dir (@part_dirs)
    {
        $j++;
        if($j == 1)
        {
            $sub_dir = $part_dir;
        }
        else
        {
            $sub_dir = join('/',$sub_dir,$part_dir);
        }
        $write_file = join('/',$startdir,$sub_dir,'DartTestfile.txt');
        system("touch $write_file");
        $testfile_input = join(' ',$testfile_input,$sub_dir);
    }
}
$testfile_input = join(' ',$testfile_input,')');

# creating the "DartTestfile.txt" file which resides in the build directory
open(TESTFILE, join('','>',$startdir,'DartTestfile.txt'));
print TESTFILE $testfile_input;
}












