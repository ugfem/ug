#!/usr/bin/perl
#
# ug_dart.pl
# 
# purpose: runs all necessary scripts on the client 

###############################################################################
# modules
###############################################################################
use File::Find ();
use Getopt::Std;

###############################################################################
# defining subroutines
###############################################################################

### test in the given directories either all test cases or only the given ones
sub Test # arguments: configuration string, test type, directory id 
{ 
    my ($LibTestDirs, $LibDirIds, $TestDirs, $DirIds);
    my (@AllLibTestCases, @AllTestCases);
    my (@TestCaseFamilies, @TestOrder);

    # get the library test directories
    ($LibTestDirs, $LibDirIds) = read_testdirs("Lib");
    
    # get the application test directories
    ($TestDirs, $DirIds) = read_testdirs("Appl");
    
    # get the library test cases and save them 
    for(my $i = 0; $i < @$LibTestDirs; $i++)
    {
        my $DirLibTestCases = read_case($LibTestDirs->[$i], $LibDirIds->[$i]);
        for(my $j = 0; $j < @$DirLibTestCases; $j++)
        {
            push(@AllLibTestCases, $DirLibTestCases->[$j]);
        }
    }
    
    # get the application test cases and save them 
    for(my $i = 0; $i < @$TestDirs; $i++)
    {
	if($DirIds->[$i] eq $_[2])
	{
	    my @DirTestCases = @{ read_case($TestDirs->[$i], $DirIds->[$i])};
	    if($_[1] eq "")
	    {
		for(my $j = 0; $j <= $#DirTestCases; $j++)
		{
		    push(@AllTestCases, $DirTestCases[$j]);
		}
	    }
	    else
	    {
		for(my $j = 0; $j <= $#DirTestCases; $j++)
		{
		    if($DirTestCases->{CFS} eq $_[1])
		    {
			push(@AllTestCases, $DirTestCases[$j]);
		    }
		}
	    }
	}
    }
    
    # divide all application test cases into 'families' which members need 
    # the same libraries
    my $FamilyCount = 0;
    while($#AllTestCases > -1)
    {
        my (@TestCaseFamily, @TestCasePosition, @BackupTestCases) = ((), (), ());
        my $CurrentTestCase = splice(@AllTestCases, 0, 1);
        
        push(@TestCaseFamily, $CurrentTestCase);
        for(my $i = 0; $i <= $#AllTestCases; $i++)
        {
            if(subset($CurrentTestCase, $AllTestCases[$i]))
            {
                push(@TestCasePosition, $i);
            }
        }
        foreach $Position (@TestCasePosition)
        {
            push(@TestCaseFamily, $AllTestCases[$Position]);
        }
        @BackupTestCases = @AllTestCases;
        @AllTestCases = ();
        for(my $i = 0; $i <= $#BackupTestCases; $i++)
        {
            my $Keep = "yes";
            foreach $Position (@TestCasePosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
            }
            if($Keep eq "yes")
            {
                push(@AllTestCases, $BackupTestCases[$i]);
            }
        }
        $TestCaseFamilies[$FamilyCount] = [ @TestCaseFamily ];
        $FamilyCount++;
    }
    
    # get the reference test case for every family 
    for(my $i = 0; $i <= $#TestCaseFamilies; $i++)
    {
        my @RefDirId = @{ $TestCaseFamilies[$i][0]->{"DirId"} };
        @RefPosition[$i] = 0;
        for(my $j = 1; $j <= $#{ $TestCaseFamilies[$i] }; $j++)
        {
            @DirId = @{ $TestCaseFamilies[$i][$j]->{"DirId"} };
            if($#DirId > $#RefDirId)
            {
                @RefDirId = @DirId;
                @RefPosition[$i] =  $j;
            }
        }
    }

    # get the library test cases needed by every reference test case and
    # build '@TestOrder' which contains all test cases to be executed in the
    # right order
    for(my $i = 0; $i <= $#TestCaseFamilies; $i++)
    {
        my @RefDirId = @{ $TestCaseFamilies[$i][$RefPosition[$i]]->{"DirId"} };
        my @RefCFS = @{ $TestCaseFamilies[$i][$RefPosition[$i]]->{"CFS"} };
        my (@RefLibPosition, @DepLibPosition) = ((), ());
        
        for(my $j = 0; $j <=  $#AllLibTestCases; $j++)
        {
            my $LibTestCase = $AllLibTestCases[$j];
            for(my $k = 0; $k <= $#RefDirId; $k++)
            {
                if($RefDirId[$k] eq $LibTestCase->{"MainDirId"} && $RefCFS[$k] eq $LibTestCase->{"MainCFS"})
                {
                    push(@RefLibPosition, $j);
                }
            }
        }
        foreach $LibPosition (@RefLibPosition)
        {
            if($AllLibTestCases[$LibPosition]->{"DirId"})
            {
                my @DepLibDirIds = @{ $AllLibTestCases[$LibPosition]->{"DirId"} };
                my @DepLibCFSs = @{ $AllLibTestCases[$LibPosition]->{"CFS"} };
                for(my $i = 0; $i <=  $#AllLibTestCases; $i++)
                {
                    for(my $j = 0; $j <= $#DepLibDirIds; $j++)
                    {
                        if($DepLibDirIds[$j] eq $AllLibTestCases[$i]->{"MainDirId"} && $DepLibCFSs[$j] eq $AllLibTestCases[$i]->{"MainCFS"})
                        {
                            push(@DepLibPosition, $i);
                        }
                    }
                }    
                
            }
        }
        foreach $Position (@DepLibPosition)
        {
            push(@TestOrder, $AllLibTestCases[$Position]);
        }
        foreach $Position (@RefLibPosition)
        {
            push(@TestOrder, $AllLibTestCases[$Position]);
        }
        my @BackupLibTestCases = ();
        @BackupLibTestCases = @AllLibTestCases;
        @AllLibTestCases = ();
        for(my $i = 0; $i <= $#BackupLibTestCases; $i++)
        {
            my $Keep = "yes";
            foreach $Position (@DepLibPosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
                
            }
            foreach $Position (@RefLibPosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
            }
            if($Keep eq "yes")
            {
                push(@AllLibTestCases, $BackupLibTestCases[$i]);
            }
        }
        push(@TestOrder, @{ $TestCaseFamilies[$i] } );
    }

    if($_[1] eq "all")
    {
	push(@TestOrder, @AllLibTestCases);
    }

    # execute the library builds and the application tests in the right order
    for(my $k = 0; $k <= $#TestOrder; $k++)
    {
        dart_test($_[0], $TestOrder[$k]);
    }
}

### build the libraries either with all configurations or only with
### the given configurations
sub Lib # arguments: configuration string, test type, directory id
{
    my ($LibTestDirs, $LibDirIds);
    my (@AllLibTestCases, @TestOrder);

    # get the library test directories
    ($LibTestDirs, $LibDirIds) = read_testdirs("Lib");
    
    # get the library test cases and save them 
    for(my $i = 0; $i < @$LibTestDirs; $i++)
    {
	if($LibDirIds->[$i] eq $_[2])
	{
	    my @DirLibTestCases = @{ read_case($LibTestDirs->[$i], $LibDirIds->[$i])};
	    if($_[1] eq "")
	    {
		for(my $j = 0; $j <= $#DirLibTestCases; $j++)
		{
		    push(@AllLibTestCases, $DirLibTestCases[$j]);
		}
	    }
	    else
	    {
		for(my $j = 0; $j <= $#DirLibTestCases; $j++)
		{
		    if($DirLibTestCases->{CFS} eq $_[1])
		    {
			push(@AllLibTestCases, $DirLibTestCases[$j]);
		    }
		}
	    }
	}
    }
    for(my $i = 0; $i <= $#AllLibTestCases; $i++)
    {
	my @DepLibPosition = ();
	if($AllLibTestCases[$i]->{"DirId"})
	{
	    my @DepLibDirIds = @{ $AllLibTestCases[$i]->{"DirId"} };
	    my @DepLibCFSs = @{ $AllLibTestCases[$i]->{"CFS"} };
	    for(my $j = 0; $j <=  $#AllLibTestCases; $j++)
	    {
		for(my $k = 0; $k <= $#DepLibDirIds; $k++)
		{
		    if($DepLibDirIds[$k] eq $AllLibTestCases[$j]->{"MainDirId"} && $DepLibCFSs[$k] eq $AllLibTestCases[$j]->{"MainCFS"})
		    {
			push(@DepLibPosition, $j);
		    }
		}
	    }    
	}
	foreach $LibPosition (@DepLibPosition)
	{
	    push(@TestOrder, $AllLibTestCases[$LibPosition]);
	}
	push(@TestOrder, $AllLibTestCases[$i]);
    }
   
    # execute the library builds in the right order
    for(my $k = 0; $k <= $#TestOrder; $k++)
    {
        dart_test($_[0], $TestOrder[$k]);
    }
}

### execute all possible builds of the libraries
sub Liball
{
    my ($LibTestDirs, $LibDirIds);
    my (@AllLibTestCases, @TestOrder);

    # get the library test directories
    ($LibTestDirs, $LibDirIds) = read_testdirs("Lib");
    
    # get the library test cases and save them 
    for(my $i = 0; $i < @$LibTestDirs; $i++)
    {
        my $DirLibTestCases = read_case($LibTestDirs->[$i], $LibDirIds->[$i]);
        for(my $j = 0; $j < @$DirLibTestCases; $j++)
        {
            push(@AllLibTestCases, $DirLibTestCases->[$j]);
        }
    }
    
    for(my $i = 0; $i <= $#AllLibTestCases; $i++)
    {
	my @DepLibPosition = ();
	if($AllLibTestCases[$i]->{"DirId"})
	{
	    my @DepLibDirIds = @{ $AllLibTestCases[$i]->{"DirId"} };
	    my @DepLibCFSs = @{ $AllLibTestCases[$i]->{"CFS"} };
	    for(my $j = 0; $j <=  $#AllLibTestCases; $j++)
	    {
		for(my $k = 0; $k <= $#DepLibDirIds; $k++)
		{
		    if($DepLibDirIds[$k] eq $AllLibTestCases[$j]->{"MainDirId"} && $DepLibCFSs[$k] eq $AllLibTestCases[$j]->{"MainCFS"})
		    {
			push(@DepLibPosition, $j);
		    }
		}
	    }    
	}
	foreach $LibPosition (@DepLibPosition)
	{
	    push(@TestOrder, $AllLibTestCases[$LibPosition]);
	}
	push(@TestOrder, $AllLibTestCases[$i]);
    }
   
    # execute the library builds in the right order
    for(my $k = 0; $k <= $#TestOrder; $k++)
    {
        dart_test($_[0], $TestOrder[$k]);
    }
}

### execute all possible builds of the libraries and all possible tests
sub AllOrTestall # argument: test type
{
    my ($LibTestDirs, $LibDirIds, $TestDirs, $DirIds);
    my (@AllLibTestCases, @AllTestCases);
    my (@TestCaseFamilies, @TestOrder);

    # get the library test directories
    ($LibTestDirs, $LibDirIds) = read_testdirs("Lib");
    
    # get the application test directories
    ($TestDirs, $DirIds) = read_testdirs("Appl");
    
    # get the library test cases and save them 
    for(my $i = 0; $i < @$LibTestDirs; $i++)
    {
        my $DirLibTestCases = read_case($LibTestDirs->[$i], $LibDirIds->[$i]);
        for(my $j = 0; $j < @$DirLibTestCases; $j++)
        {
            push(@AllLibTestCases, $DirLibTestCases->[$j]);
        }
    }
    
    # get the application test cases and save them 
    for(my $i = 0; $i < @$TestDirs; $i++)
    {
        my $DirTestCases = read_case($TestDirs->[$i], $DirIds->[$i]);
        for(my $j = 0; $j < @$DirTestCases; $j++)
        {
            push(@AllTestCases, $DirTestCases->[$j]);
        }
    }
    
    # divide all application test cases into 'families' which members need 
    # the same libraries
    my $FamilyCount = 0;
    while($#AllTestCases > -1)
    {
        my (@TestCaseFamily, @TestCasePosition, @BackupTestCases) = ((), (), ());
        my $CurrentTestCase = splice(@AllTestCases, 0, 1);
        
        push(@TestCaseFamily, $CurrentTestCase);
        for(my $i = 0; $i <= $#AllTestCases; $i++)
        {
            if(subset($CurrentTestCase, $AllTestCases[$i]))
            {
                push(@TestCasePosition, $i);
            }
        }
        foreach $Position (@TestCasePosition)
        {
            push(@TestCaseFamily, $AllTestCases[$Position]);
        }
        @BackupTestCases = @AllTestCases;
        @AllTestCases = ();
        for(my $i = 0; $i <= $#BackupTestCases; $i++)
        {
            my $Keep = "yes";
            foreach $Position (@TestCasePosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
            }
            if($Keep eq "yes")
            {
                push(@AllTestCases, $BackupTestCases[$i]);
            }
        }
        $TestCaseFamilies[$FamilyCount] = [ @TestCaseFamily ];
        $FamilyCount++;
    }
    
    # get the reference test case for every family 
    for(my $i = 0; $i <= $#TestCaseFamilies; $i++)
    {
        my @RefDirId = @{ $TestCaseFamilies[$i][0]->{"DirId"} };
        @RefPosition[$i] = 0;
        for(my $j = 1; $j <= $#{ $TestCaseFamilies[$i] }; $j++)
        {
            @DirId = @{ $TestCaseFamilies[$i][$j]->{"DirId"} };
            if($#DirId > $#RefDirId)
            {
                @RefDirId = @DirId;
                @RefPosition[$i] =  $j;
            }
        }
    }

    # get the library test cases needed by every reference test case and
    # build '@TestOrder' which contains all test cases to be executed in the
    # right order
    for(my $i = 0; $i <= $#TestCaseFamilies; $i++)
    {
        my @RefDirId = @{ $TestCaseFamilies[$i][$RefPosition[$i]]->{"DirId"} };
        my @RefCFS = @{ $TestCaseFamilies[$i][$RefPosition[$i]]->{"CFS"} };
        my (@RefLibPosition, @DepLibPosition) = ((), ());
        
        for(my $j = 0; $j <=  $#AllLibTestCases; $j++)
        {
            my $LibTestCase = $AllLibTestCases[$j];
            for(my $k = 0; $k <= $#RefDirId; $k++)
            {
                if($RefDirId[$k] eq $LibTestCase->{"MainDirId"} && $RefCFS[$k] eq $LibTestCase->{"MainCFS"})
                {
                    push(@RefLibPosition, $j);
                }
            }
        }
        foreach $LibPosition (@RefLibPosition)
        {
            if($AllLibTestCases[$LibPosition]->{"DirId"})
            {
                my @DepLibDirIds = @{ $AllLibTestCases[$LibPosition]->{"DirId"} };
                my @DepLibCFSs = @{ $AllLibTestCases[$LibPosition]->{"CFS"} };
                for(my $i = 0; $i <=  $#AllLibTestCases; $i++)
                {
                    for(my $j = 0; $j <= $#DepLibDirIds; $j++)
                    {
                        if($DepLibDirIds[$j] eq $AllLibTestCases[$i]->{"MainDirId"} && $DepLibCFSs[$j] eq $AllLibTestCases[$i]->{"MainCFS"})
                        {
                            push(@DepLibPosition, $i);
                        }
                    }
                }    
                
            }
        }
        foreach $Position (@DepLibPosition)
        {
            push(@TestOrder, $AllLibTestCases[$Position]);
        }
        foreach $Position (@RefLibPosition)
        {
            push(@TestOrder, $AllLibTestCases[$Position]);
        }
        my @BackupLibTestCases = ();
        @BackupLibTestCases = @AllLibTestCases;
        @AllLibTestCases = ();
        for(my $i = 0; $i <= $#BackupLibTestCases; $i++)
        {
            my $Keep = "yes";
            foreach $Position (@DepLibPosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
                
            }
            foreach $Position (@RefLibPosition)
            {
                if($i == $Position)
                {
                    $Keep = "no";
                }
            }
            if($Keep eq "yes")
            {
                push(@AllLibTestCases, $BackupLibTestCases[$i]);
            }
        }
        push(@TestOrder, @{ $TestCaseFamilies[$i] } );
    }

    if($_[1] eq "all")
    {
	push(@TestOrder, @AllLibTestCases);
    }

    # execute the library builds and the application tests in the right order
    for(my $k = 0; $k <= $#TestOrder; $k++)
    {
        dart_test($_[0], $TestOrder[$k]);
    }
}

### get the test directories
sub read_testdirs
{
    my @Line;
    my @TestDirs;
    my @DirIds;
    
    # save the path to the DartTestFiles file 
    my $TestRootDir = join('',$ENV{"DART_HOME"},'/Source/Client/UG');
    
    # open the DartTestFiles file to get directories where TestCases files reside
    open(TESTFILES, join('',$TestRootDir,'/DartTestFiles'));
    while(<TESTFILES>)
    {
        # cut off the "end-of-line" (\n)
        chomp;
        
        # split the line
        @Line = split(/:/, $_);
        if($_[0] eq "Lib")
        {
            if(uc($Line[0]) =~ /LIB/)
            {
                # save the whole directory in @testdirs
                push(@TestDirs, join('/',$TestRootDir,$Line[1]));
                push(@DirIds, $Line[0]);
            }
        }
        else
        {
            unless(uc($Line[0]) =~ /LIB/)
            {
                # save the whole directory in @testdirs
                push(@TestDirs, join('/',$TestRootDir,$Line[1]));
                push(@DirIds, $Line[0]);
            }
        }
    }
    close(TESTFILES);
    return \(@TestDirs, @DirIds);    
}

### get the test cases
sub read_case # argument: test directory, directory id
{
    my @TestCases;
         
    # open the testcase file and read the entries
    open(TESTCASE, join('',$_[0],'/TestCases'));
    while(<TESTCASE>)
    {
        my (@TestCase, @ECFS, %TestCase) = ((), (), ());
        my (@Arch, @DirIds, @ECFS, @ECFSParts, @CFS, @CFSParts) = ((), (), (), (), (), ());
        my ($CFS, $MainCFS, $MainDirId, $MainArch) = ("", "", "", "");
        my ($WithDependencies, $WithMain) = ("no", "no");

        chomp;
        @TestCase = split(/::/, $_);
        @ECFS = split(/:/, $TestCase[2]);
        foreach $ECFSPart (@ECFS)
        {
            @ECFSParts = split(',', $ECFSPart);
            unless(uc($ECFSParts[0]) eq uc($_[1]))
            {
                push(@DirIds, $ECFSParts[0]);
                @CFSParts = split(/-/, $ECFSParts[1]);
                push(@Arch, splice(@CFSParts, 0, 1));
                $CFS = "";
                foreach $CFSPart (@CFSParts)
                {
                    $CFS = join('', $CFS, $CFSPart, ' ');
                }
                push(@CFS, $CFS);
                $WithDependencies = "yes";
            }
            else
            {
                $MainDirId = $ECFSParts[0];
                @CFSParts = split(/-/, $ECFSParts[1]);
                $MainArch = splice(@CFSParts, 0, 1);
                $MainCFS = "";
                foreach $CFSPart (@CFSParts)
                {
                    $MainCFS = join('', $MainCFS, $CFSPart, ' ');
                }
                $WithMain = "yes";
            }
        }
        $TestCase{"TestDir"}   = $_[0];
        $TestCase{"TestId"}    = $TestCase[0];       
        $TestCase{"TestApp"}   = $TestCase[1];
        if($WithMain eq "yes")
        {
            $TestCase{"MainArch"}  = $MainArch;
            $TestCase{"MainDirId"} = $MainDirId;
            $TestCase{"MainCFS"}   = $MainCFS;
        }
        else
        {
            $TestCase{"MainArch"}  = $Arch[0];
            $TestCase{"MainDirId"} = $_[1];
            $TestCase{"MainCFS"}   = "";
        }
        if($WithDependencies eq "yes")
        {
            $TestCase{"Arch"}      = \@Arch;
            $TestCase{"DirId"}     = \@DirIds;
            $TestCase{"CFS"}       = \@CFS;
        }
        push(@TestCases, \%TestCase);
    }
    close(TESTCASE);
    return \@TestCases;
}

### execute the necessary commands to perform the build-test-submit cycle
sub dart_test # arguments: test type, testcase(hash)
{
    my $TestType;
    my $TestCase;
    my ($TestDir, $TestId, $Arch, $TestApp, $DirId, $CFS);
    my $WithCoverage;
    my $WithTest;
  
    # get the test type
    $TestType = $_[0];
    
    # save the testcase as hash
    $TestCase = $_[1];
    
    # save the hash values in variables to pass them to the other scripts
    $TestDir = $TestCase->{"TestDir"};
    $TestId = $TestCase->{"TestId"};
    $TestApp = $TestCase->{"TestApp"};
    $Arch = $TestCase->{"MainArch"};
    $DirId = $TestCase->{"MainDirId"};
    $CFS = $TestCase->{"MainCFS"};
    
    # coverage: yes or no
    if($Arch =~ /GCOV/)
    {
        $WithCoverage = "yes";
    }
    else 
    {
        $WithCoverage = "no";
    }
    # build-test-submit cycle
    system("ug_dart_conf.pl -b $TestDir -i $TestId -a $Arch -d $DirId -c '$CFS'");
    
    # execute ug_dart_tests.pl and call Test.tcl only if a test application is given
    if($TestApp eq "")
    {
        $WithTest = "no";
        system("ug_dart_client_test $TestType $TestDir $WithCoverage $WithTest");
    }
    else
    {
        $WithTest = "yes";
        system("ug_dart_tests.pl -b $TestDir -s $TestDir -p '$TestApp' -i $TestId");
        system("ug_dart_client_test $TestType $TestDir $WithCoverage $WithTest");
    }
}

### check whether the configuration of two test cases matches or not
sub subset
{
    my $match = -1;
    my @DirId1 = @{ $_[0]->{"DirId"} };
    my @DirId2 = @{ $_[1]->{"DirId"} };
    my @CFS1 = @{ $_[0]->{"CFS"} };
    my @CFS2 = @{ $_[1]->{"CFS"} };
    
    # count the matching directory ids and CFS
    for(my $i = 0; $i <= $#DirId1; $i++)
    {
        for(my $j = 0; $j <= $#DirId2; $j++)
        {
            if($DirId1[$i] eq $DirId2[$j] && $CFS1[$i] eq $CFS2[$j])
            {
                $match++;
            }
        }
    }
    
    # if $match = $#DirId1($#DirId2) then test case configuration 1(2) is a subset of test
    # case configuration 2(1)
    if($match == $#DirId1 || $match == $#DirId2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

### print help message 
sub help
{
    print "usage:	ug_dart.pl [-c configuration string] [-e test extent] [-t test type]\n\n";
    print "-c: Configuration String which is necessary if the test extent\n";
    print "    is either 'test' or 'lib'. It has to be given in the form\n";
    print "    -c <DirId1>:<TestId1,...,TestIdn> ... <DirIdn>:<TestId1,...,TestIdn>,\n";
    print "    e.g. -c 'Appl1:tutor3,tutor5 Appl2:cd36 (est extent = test)\n";
    print "-e: Test extent.May be\n";
    print "      - manually: execute a given test\n";
    print "      - test    : execute either all tests in the given directory\n";
    print "                  or only several test cases\n";
    print "      - lib     : building libraries with either all configurations\n";
    print "                  or only with several configurations\n";
    print "      - test_all: execute all tests in all directories\n";
    print "      - all     : build the libraries with all configurations\n";
    print "                  execute all tests\n";
    print "-t: May be Nightly, Continuous or Experimental\n";
    print "    Default: Nightly\n";
    print "-h: Print this help message\n";    
    print "\n";
    print "purpose: running all scripts which are necessary for a complete dart\n";
    print "         build-test-submit cycle\n";
}

###############################################################################
# main
###############################################################################

# read the command line parameters
my %option = ();
getopts("c:d:e:t:h", \%option);

# print help message if the parameter -h is given
if($option{h})
{
    help();
    exit 0;
}

# set test type "Nightly" if no test type is given
unless($option{t})
{
    $option{t} = "Nightly";
}

# exit if no configuration string is given but it is needed
if($option{e} eq "test" || $option{e} eq "lib")
{
    unless($option{d})
    {
        die "A build directory identifier is recommended!\n";
    }
    unless($option{c})
    {
	$option{c} = "";
    }
}

# distinguish between the different entries of the test_extent variable
if($option{e} eq "test")
{
    Test($option{t}, $option{c}, $option{d});
}
elsif($option{e} eq "lib")
{
    Lib($option{t}, $option{c}, $option{d});
}
elsif($option{e} eq "lib_all")
{
    Liball($option{t});
}
elsif($option{e} eq "all" || $option{e} eq "test_all")
{
    AllOrTestall($option{t}, $option{e});
}
else 
{
    die "The test extent you entered isn't defined!\n";
}

###############################################################################
