#!/usr/bin/perl
# ug_dart_tests.pl
#
# findind "DartTestfile.txt" files which declare tests
# and create all other needed "DartTestfile.txt" files

# use the following modules
use File::Find ();

# reading and saving the Dart root directory 
my $dartdir = $ENV {"DART_HOME"};

# defining the start directory 
my $startdir = join('',$dartdir,'/Source/Client/UG'); # default
if($#ARGV > 0)
{
    $startdir = join('/',$startdir,$ARGV[0]);
}

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

# saving only the whole path to the "DartTestfile.txt" files in which tests are
# defined in "@testdirs" 
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
my $root_testfile_input = "SUBDIRS(";
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
        $root_testfile_input = join(' ',$root_testfile_input,$sub_dir);
    }
}
$root_testfile_input = join(' ',$root_testfile_input,')');

# creating the "DartTestfile.txt" file which resides in the build directory
open(ROOT_TESTFILE, join('','>',$startdir,'/DartTestfile.txt'));
print ROOT_TESTFILE $root_testfile_input;
