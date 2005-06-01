#!/usr/bin/perl

###############################################################################
# load modules
###############################################################################
use File::Find;
use Time::localtime;
use Time::Local;

###############################################################################
# global variables
###############################################################################
# half a year and a quarter year in seconds
my ($HY, $QY)=(15811200, 7905600);
# current time in DMYHMS format
my $tm=localtime;
# current time in epoch seconds
my $T=timelocal($tm->sec, $tm->min, $tm->hour, $tm->mday, $tm->mon,$tm->year);
my @DirsToDelete = ();

###############################################################################
# subroutines
###############################################################################
# delete outdated directories
sub DeleteDir
{
    print "\nDeleting: ", $_[0];
    system("rm -rf $_[0]");
}

# get outdated directories
sub GetDirsToDelete
{
    # writetime in epoch seconds
    my $WT = (stat($_))[9];
    # writetime in DMYHMS format
    my $tm_w = localtime($WT);
	    
    # check if the directory is outdated
    if(-d $_)
    {
	if(($T > $WT + $HY) && !($tm_w->wday == 1 && $tm_w->mday <=7))
	{
	    push(@DirsToDelete, $File::Find::name);
	}
	elsif(($T > $WT + $QY) && ($T < $WT + $HY) && !($tm_w->wday == 1))
	{
	    push(@DirsToDelete, $File::Find::name);
	}
    }
}

###############################################################################
# main 
###############################################################################

my @BaseDir = ();

push(@BaseDir, "/home/ugtest/public_html/TestingResults/Sites/speedo");
push(@BaseDir, "/home/ugtest/public_html/TestingResults/Dashboard");

find(\&GetDirsToDelete, @BaseDir);

foreach $Dir (@DirsToDelete)
{
    DeleteDir($Dir);
}
