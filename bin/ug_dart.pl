#!/usr/bin/perl
#
# ug_dart.pl
# 
# purpose: runs all necessary scripts to perform a dart cycle on the client 

###############################################################################
# modules
###############################################################################
use File::Find ();
use Getopt::Std;

###############################################################################
# global variables
###############################################################################
my $TestRootdir = join('', $ENV{"DART_HOME"}, '/Source/Client/UG');

###############################################################################
# subroutines
###############################################################################

#### test one specific application
###############################################################################
sub Test_App
{
#### parameters
    my ($Type, $TestOrder) = @_;
    
####
    for(my $i = 0; $i < @$TestOrder; $i++)
    {
	RunTest($option{t}, $TestOrder->[$i]);
    }    
}

#### test all applications
###############################################################################
sub Test_AppAll
{
#### parameters
    my ($Type, $TestOrder) = @_;
    
####
    for(my $i = 0; $i < @$TestOrder; $i++)
    {
	RunTest($option{t}, $TestOrder->[$i]);
    }    
}

#### test one specific library
###############################################################################
sub Test_Lib
{
#### parameters
    my ($Type, $TestOrder) = @_;
    
####
    for(my $i = 0; $i < @$TestOrder; $i++)
    {
	RunTest($option{t}, $TestOrder->[$i]);
    }    
}

#### test all libraries
###############################################################################
sub Test_LibAll
{
#### parameters
    my ($Type, $TestOrder) = @_;
    
####
    for(my $i = 0; $i < @$TestOrder; $i++)
    {
	RunTest($option{t}, $TestOrder->[$i]);
    }    
}

#### test all applications and all libraries
###############################################################################
sub Test_All
{
#### parameters
    my ($Type, $TestOrder) = @_;
    
####
    for(my $i = 0; $i <= $#{ $TestOrder }; $i++)
    {
	RunTest($Type, $TestOrder->[$i]);
    }    
}

#### sort the test cases to perform an efficient test cycle
###############################################################################
sub SetupTestOrder
{
#### local variables
    my ($DirId, $TestDir, $TestcaseFile, $LibOrApp) = ReadTestfile();
    my (@LibTest, @AppTest) = ();
    my $NotDefinedCount = 0;
    my @TestOrder = ();
    	
#### collect all test cases and save them either in @LibTest or @AppTest
    for(my $i = 0; $i < @$DirId; $i++)
    {
	my $Test = ReadCasefile($DirId->[$i], $TestDir->[$i], 
				$TestcaseFile->[$i]);
	
	if($LibOrApp->[$i] eq "LIB")
	{
	    for(my $j = 0; $j < @$Test; $j++)
	    {
		push(@LibTest, $Test->[$j]);
	    }
	}
	elsif($LibOrApp->[$i] eq "APP")
	{
	    for(my $j = 0; $j < @$Test; $j++)
	    {
		push(@AppTest, $Test->[$j]);
	    }
	}
    }

#### create needed but not defined library test cases
    for(my $i = 0; $i <= $#AppTest; $i++)
    {
	for(my $j = 0; $j <= $#{ $AppTest[$i]->{"DirId"} }; $j++)
	{
	    my $Match = 0;
	    for(my $k = 0; $k <= $#LibTest; $k++)
	    {
		if($AppTest[$i]->{"DirId"}->[$j] eq $LibTest[$k]->{"SelfDirId"}
		   && $AppTest[$i]->{"Arch"} eq $LibTest[$k]->{"Arch"} 
		   && $AppTest[$i]->{"CFS"}->[$j] eq $LibTest[$k]->{"SelfCFS"})
		{
		    $Match = 1;
		}
	    }
	    unless($Match)
	    {
		my %NotDefinedCase = ();
		my (@DirId, @CFS) = ();
		
		$NotDefinedCount ++;
		for(my $k = 0; $k <= $#{ $DirId }; $k++)
		{
		    if($DirId->[$k] eq $AppTest[$i]->{"DirId"}->[$j])
		    {
			$NotDefinedCase{"Dir"} = $TestDir->[$k];
		    }
		}
		$NotDefinedCase{"Id"} = "ND.$NotDefinedCount";
		$NotDefinedCase{"App"} = "";
		$NotDefinedCase{"Arch"} = $AppTest[$i]->{"Arch"};
		$NotDefinedCase{"SelfDirId"} = $AppTest[$i]->{"DirId"}->[$j];
		$NotDefinedCase{"SelfCFS"} = $AppTest[$i]->{"CFS"}->[$j] ;
		$NotDefinedCase{"DirId"} = [ @DirId ];
		$NotDefinedCase{"CFS"} = [ @CFS ];
		push(@LibTest, { %NotDefinedCase });
	    }
	}
    }

#### divide the test cases into classes
    my $LibClass = SetupTestClasses([ @LibTest ]);
    my $AppClass = SetupTestClasses([ @AppTest ]);

#### get a reference test case per class
    my $LibClassReference = GetClassReference($LibClass);
    my $AppClassReference = GetClassReference($AppClass);
    
#### get the dependencies for every class
    my $LibClassDep = GetClassDep($LibClassReference, [ @LibTest ]);
    my $AppClassDep = GetClassDep($AppClassReference, [ @LibTest ]);

#### get the needed library classes for every class
    my $LibClassDepIndex = GetClassDepIndex($LibClassDep, $LibClassReference);
    my $AppClassDepIndex = GetClassDepIndex($AppClassDep, $LibClassReference);
    
#### find classes with 'overlapping' dependencies

#### build the test order
    my @Match = ();
    for(my $i = 0; $i <= $#{ $AppClassDepIndex }; $i++)
    {
	for(my $k = 0; $k <= $#{ $AppClassDepIndex->[$i] }; $k++)
	{
	    for(my $j = 0; $j <= $#{ $LibClass }; $j++)
	    {
		if($j == $AppClassDepIndex->[$i][$k])
		{
		    $Match[$j] = 1;
		    for(my $l = 0; $l <= $#{ $LibClass->[$j] }; $l++)
		    {
			push(@TestOrder, $LibClass->[$j][$l]);
		    }
		    for(my $l = 0; $l <= $#{ $AppClass->[$i] }; $l++)
		    {
			push(@TestOrder, $AppClass->[$i][$l]);
		    }
		}
	    }
	}
    }
    for(my $i = 0; $i <= $#Match; $i++)
    {
	unless($Match[$i])
	{
	    for(my $j = 0; $j <= $#{ $LibClass->[$i] }; $j++)
	    {
		push(@TestOrder, $LibClass->[$i][$j]);
	    }
	}
    }
    return \@TestOrder;   
}

#### read the test file entries in the "DartTestFile" file
###############################################################################
sub ReadTestfile 
{
#### local variables
    my (@DirId, @TestDir, @TestcaseFile, @LibOrApp) = ((), (), (), ());
    
#### get the test directories, directory identifiers and the 'TestCase.xxx'
#### files
    open(TESTFILE, join('', '< ', $TestRootdir, '/DartTestFiles'));
  LINE: while(<TESTFILE>)
  {
      if($_ =~ /^\s*(.*)\s*:\s*(.*)\s*:\s*(.*)\s*/)
      {
	  my ($DirId, $TestDir, $TestcaseFile) = ($1, $2, $3);
	  
	  push(@DirId, $DirId);
	  push(@TestDir, join('/',$TestRootdir,$TestDir));
	  push(@TestcaseFile, $TestcaseFile);
	  if($DirId =~ /[lL][iI][bB]/)
	  {
	      push(@LibOrApp, "LIB");
	      
	  }
	  else
	  {
	      push(@LibOrApp, "APP");
	  }
      }
      else
      {
	  next LINE;
      }
  }
    close(TESTFILE);
    return \(@DirId, @TestDir, @TestcaseFile, @LibOrApp);    
}

#### read the test case entries in the "TestCases.xxx" file
###############################################################################
sub ReadCasefile
{
#### parameters
    my ($DirId, $TestDir, $TestcaseFile) = @_;
 
#### local variables
    my @Test;

#### get the test case information
    open(TESTCASE, join('', '< ', $TestDir, ,'/', $TestcaseFile));
  LINE: while(<TESTCASE>)
  {
      my %Test = ();
      my (@DirId, @CFS, @CFSParams) = ((), (), ());
      my $CFS = "";
      my $SelfConfig = 0;
      
      $Test{"Dir"} = $TestDir;
      chomp;
      if($_ =~ /^\s*(.*)\s*::\s*(.*)\s*::\s*(.*)\s*$/)
      {
	  my $eCFS = $3;
	  my @eCFS = split(/:/, $eCFS);
	 
	  $Test{"Id"} = $1;
	  $Test{"App"} = $2;
	  for(my $i = 0; $i <= $#eCFS; $i++)
	  {
	      if($eCFS[$i] =~ /^\s*(.*)\s*,\s*(.*)\s*$/)
	      {
		  @CFSParams = split(/-/, $2);
		  if($i == 0)
		  {
		      $Test{"Arch"} = shift(@CFSParams);
		  }
		  else
		  {
		      my $Dummy = shift(@CFSParams);
		  }
		  for(my $j = 0; $j <= $#CFSParams; $j++)
		  {
		      $CFS = join('', $CFS, $CFSParams[$j], ' ');
		  }
		  if($DirId eq $1)
		  {
		      $Test{"SelfDirId"} = $1;
		      $Test{"SelfCFS"} = $CFS;
		      $CFS = "";
		      $SelfConfig = 1;
		  }
		  else
		  {
		      push(@DirId, $1);
		      push(@CFS, $CFS);
		      $CFS = "";
		  }
	      }
	  }
	  $Test{"DirId"} = [ @DirId ];
	  $Test{"CFS"} = [ @CFS ];
      }
      else
      {
	  next LINE;
      }
      if($SelfConfig == 1)
      {
	  push(@Test, { %Test });  
      }
      else
      {
	  print "Error: The test case has got no own configuration string!";
      }
  }
    return \@Test;
}

#### divide all test cases into classes
###############################################################################
sub SetupTestClasses
{
#### parameters
    my $Test = $_[0];

#### local variables
    my ($ClassCount, $EmptyEntryCount) = 0;
    my @Class = ();
   
#### divide all test cases into classes which members have the same library
#### dependencies
    while($EmptyEntryCount < $#{ $Test } + 1)
    {
	my $CaseCount = 0;
      TESTCASE: for(my $i = 0; $i <= $#{ $Test }; $i++)
      {
	  unless($Test->[$i] eq "")
	  {
	      $EmptyEntryCount++;
	      $Class[$ClassCount][$CaseCount] = $Test->[$i];
	      $Test->[$i] = "";
	      last TESTCASE;
	  }
      }
      TESTCASE: for(my $i = 0; $i <= $#{ $Test }; $i++)
      {
	  if($Test->[$i] eq "")
	  {
	      next TESTCASE;
	  }
	  elsif((CheckSubset($Class[$ClassCount][0], $Test->[$i])) ||
		(CheckSubset($Test->[$i], $Class[$ClassCount][0])))
	  {
	      $CaseCount++;
	      $EmptyEntryCount++;
	      $Class[$ClassCount][$CaseCount] = $Test->[$i];
	      $Test->[$i] = "";
	  }
      }
	$ClassCount++;
    }
    return \@Class;
}

#### check whether the configuration of test case 1 is a subset of the one of 
#### the second test case
###############################################################################
sub CheckSubset
{
#### parameters
    my ($Testcase1, $Testcase2) = @_;

#### local variables
    my $MatchCount = -1;

#### check
    for(my $i = 0; $i <= $#{ $Testcase1->{"CFS"} }; $i++)
    {
	for(my $j = 0; $j <= $#{ $Testcase2->{"CFS"} }; $j++)
	{
	    if(($Testcase1->{"DirId"}->[$i] eq $Testcase2->{"DirId"}->[$j]) && 
	       ($Testcase1->{"Arch"} eq $Testcase2->{"Arch"}) &&
	       ($Testcase1->{"CFS"}->[$i] eq $Testcase2->{"CFS"}->[$j]))
	    {
		$MatchCount++;
	    }
	}
    }
    if(($MatchCount == $#{ $Testcase1->{"DirId"} }) && 
       !($#{ $Testcase1->{"DirId"} } == -1))
    {
	return 1;
    }
    else
    {
	return 0;
    }
}

#### get a reference test case per class
###############################################################################
sub GetClassReference
{
#### parameters
    my $Class = $_[0];

#### local variables
    my @ClassReference = ();

#### get the reference case with aid of "CheckSubset" 
    for(my $i = 0; $i <= $#{ $Class }; $i++)
    {
	my $ClassReference = $Class->[$i][0];

	for(my $j = 1; $j <= $#{ $Class->[$i] }; $j++)
	{
	    if(CheckSubset($ClassReference, $Class->[$i][$j]))
	    {
		$ClassReference = $Class->[$i][$j];
	    }
	}
	push(@ClassReference, $ClassReference);
    }
    return \@ClassReference;
}

#### get the library dependencies for every class
###############################################################################
sub GetClassDep
{
#### parameters
    my ($ClassReference, $LibTest) = @_;

#### local variables
    my @ClassDep = ();

    for(my $i = 0; $i <= $#{ $ClassReference }; $i++)
    {
      DIRID: for(my $j = 0; $j <= $#{ $ClassReference->[$i]->{"DirId"} }; $j++)
      {
	  my $CaseCount = 0;
	  
	  for(my $k = 0; $k <= $#{ $LibTest }; $k++)
	  {
	      if(($ClassReference->[$i]->{"DirId"}->[$j] eq 
		  $LibTest->[$k]->{"SelfDirId"}) && 
		 ($ClassReference->[$i]->{"Arch"} eq 
		  $LibTest->[$k]->{"Arch"}) && 
		 ($ClassReference->[$i]->{"CFS"}->[$j] eq 
		  $LibTest->[$k]->{"SelfCFS"}))
	      {
		  $ClassDep[$i][$CaseCount] = $LibTest->[$k];
		  $CaseCount++;
		  next DIRID;
	      }
	  } 
	  
      }
    } 
    return \@ClassDep;
}

#### get the indices of the neeeded libraries per class
###############################################################################
sub GetClassDepIndex
{
#### parameters
    my ($ClassDep, $LibClassReference) = @_;

#### local variables
    my $DepCount = 0;
    my @ClassDepIndex = ();

   
    for(my $i = 0; $i <= $#{ $ClassDep }; $i++)
    {
	for(my $j = 0; $j <= $#{ $ClassDep->[$i] }; $j++)
	{
	    for(my $k = 0; $k <= $#{ $LibClassReference }; $k++)
	    {
		if($ClassDep->[$i][$j]->{"Id"} eq 
		   $LibClassReference->[$k]->{"Id"})
		{
		    $ClassDepIndex[$i][$DepCount] = $k;
		    $DepCount++;
		}
	    }
	}
	$DepCount = 0;    
    }
    return \@ClassDepIndex;
}


#### execute the necessary commands to perform the build-test-submit cycle
###############################################################################
sub RunTest # arguments: test type, test case
{
#### parameters
    my ($Type, $Test) = @_;

#### local variables
    my $Coverage = "no";
# save the hash values in variables to pass them to the other scripts
    foreach $Key ( keys %{ $Test } )
    {
	${ $Key } = $Test->{$Key};
    }
    if($Arch =~ /GCOV/)
    {
        $Coverage = "yes";
    }

#### build-test-submit cycle
    system("ug_dart_conf.pl -b $Dir -i $Id -a $Arch -d $SelfDirId -c '$SelfCFS'");
# execute ug_dart_tests.pl and call Test.tcl only if a test application is 
# given
    if($App eq "")
    {
    	system("ug_dart_client_test $Type $Dir $Coverage no");
    }
    else
    {
       	system("ug_dart_tests.pl -b $Dir -s $Dir -p '$App' -i $Id");
	system("ug_dart_client_test $Type $Dir $Coverage yes");
    }
}

#### print help message
###############################################################################
sub Help
{
    print "\nusage:	ug_dart.pl [-c configuration string] [-e test extent]";
    print " [-t test type]\n\n";
    print "-e: Test extent.May be\n";
    print "      - manually: execute a given test\n";
    print "      - test    : execute either all tests in the given ";
    print "directory\n";
    print "                  or only several test cases\n";
    print "      - lib     : building libraries with either all ";
    print "configurations\n";
    print "                  or only with several configurations\n";
    print "      - test_all: execute all tests in all directories\n";
    print "      - all     : build the libraries with all configurations\n";
    print "                  execute all tests\n\n";
    print "-c: Configurations to test. This is necessary if the test extent\n";
    print "    is either 'test' or 'lib'. It has to be given in the form\n";
    print "    -c <DirId1>:<TestId1,...,TestIdn> ... <DirIdn>:<TestId1,...,";
    print "TestIdn>,\n";
    print "    e.g. -c 'TUTOR:tutor3,tutor5 CD:cd36 (test extent = test)\n\n";
    print "-t: May be Nightly, Continuous or Experimental\n";
    print "    Default: Nightly\n\n";
    print "-h: Print this help message\n";    
    print "\n\n";
    print "purpose: running all scripts which are necessary for a complete ";
    print "dart\n";
    print "         build-test-submit cycle\n\n";
}

###############################################################################
# main
###############################################################################

#### read the command line parameters
my %option = ();
getopts("c:e:t:h", \%option);

#### print help message if the parameter -h is given
if($option{h})
{
    Help();
    exit 0;
}

#### set test type "Nightly" if no test type is given
unless($option{t})
{
    $option{t} = "Nightly";
}

#### exit if no configuration string is given but it is needif($option{e} eq 
#### "test" || $option{e} eq "lib")
{
    unless($option{c})
    {
	$option{c} = "";
    }
}

#### distinguish between the different entries of the test_extent variable
if($option{e} eq "app")
{
    Test_App($option{t}, $option{c});
}
elsif($option{e} eq "app_all")
{
    Test_AppAll($option{t}, SetupTestOrder());
}
elsif($option{e} eq "lib")
{
    Test_Lib($option{t}, $option{c});
}
elsif($option{e} eq "lib_all")
{
    Test_LibAll($option{t});
}
elsif($option{e} eq "all")
{
    Test_All($option{t}, SetupTestOrder());
}
else 
{
    die "The test extent you entered isn't defined!\n";
}

###############################################################################
