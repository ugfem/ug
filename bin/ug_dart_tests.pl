#!/usr/bin/perl
# ug_dart_tests.pl
#
# findind "DartTestfile.txt" files which declare tests
# and create all other neede "DartTestfile.txt" files

# use the following modules
use File::Find ();
use File::Copy;

# reading and saving the Dart install directory 
my $dartdir = $ENV {"DART_HOME"};

# defining start directory and saving its "length"
my $i = 0;
$startdir = join('',$dartdir,'/Source/Client/UG');
foreach $byte (split //,$startdir)
{
	$i++;
}

#*****************************************************************************
# searching for existing "DartTestfile.txt" files and saving the directories
# in which they reside in the variable "@testdirs" 
#*****************************************************************************

# searching for "DartTestfile.txt" files and writing the directories
# in which they reside to "@temptestdirs"
sub find(&@) { &File::Find::find }
*name = *File::Find::name;
find { push(@dirs, $name) if -d } $startdir;
foreach $dir (@dirs)
{
	opendir(DIR,$dir);
	while(defined($file = readdir(DIR)))
	{
		if($file eq "DartTestfile.txt")
		{
			push(@temptestdirs,$dir);
		}
	}
	closedir(DIR);
}

# prevent from writing some directories multiple times to "@testdirs"
# and writing directories to "@testdirs"
my ($l,$savedir,$tempdir) = (0,"","");
foreach $dir (@temptestdirs)
{
	$l++;
	if($l != 1)
	{
        	if($dir =~ /$savedir/)
		{
			$tempdir = $dir;
		}
		else
		{
			push(@testdirs,$tempdir);
			$tempdir = $dir;
		}			
	}
	else
	{
		$tempdir = $dir;
	}
	$savedir = $dir;
}
if(l == 1)
{
	push(@testdirs,$savedir);
}
else
{
	push(@testdirs,$savedir);
}

#*****************************************************************************
# extracting the subdirectories from every directory in "@testdirs" and 
# creating the directories which are necessary
#*****************************************************************************

foreach $dir (@testdirs)
{
	my ($j,$k,@I) = (0,0,());
	my $partdir = substr($dir,$i);
	print $partdir,"\n";
	# searching for backslashs in the directory name "$dir"
	foreach $byte (split //,$partdir)
	{
	        $j++;
		if($byte eq "/")
		{
			push(@I,$j);
		}
	}
	# extracting the subdirectories between the backslashs
	foreach $index (@I)
	{
	        $k++;
		if($k != 1)
		{
			$count = $index - $offset - 1;
			push(@partdirs,substr($partdir,$offset,$count));
		}
		$offset = $index;
	}
	$count = $j - $offset;
	push(@partdirs,substr($partdir,$offset,$count)); 
	$l = 0;
	$tempdir2 = $startdir;
	# building the several subdirectories
	foreach $dir (@partdirs)
	{
#		print $dir,"\n\n";
		$l++;
		if($l == 1)
		{
			$tempdir1 = $dir;
		}
		else
		{
			$tempdir1 = join('/',$tempdir1,$dir);
		}	
		$tempdir2 = join('/',$tempdir2,$dir);
		push(@subdirs,$tempdir1);
		push(@tempdirs,$tempdir2);
	}
}

#*****************************************************************************
# creating the needed "DartTestfile.txt" files and the SUB_DIRS command which
# is needed in "DartTestfile.txt" file which resides in the BuildDirectory
#*****************************************************************************

foreach $dir (@tempdirs)
{
      	$testfile = join('/',$dir,'DartTestfile.txt');
	system("touch $testfile");
}
$SUB_DIRS = "SUBDIRS(";
foreach $dir (@subdirs)
{
	$SUB_DIRS = join(' ',$SUB_DIRS,$dir);
}
$SUB_DIRS = join(' ',$SUB_DIRS,')');
$testfile = join('/',$startdir,'DartTestfile.txt');
system("touch $testfile");
open(TEST,join('','>',$testfile));
print TEST $SUB_DIRS
