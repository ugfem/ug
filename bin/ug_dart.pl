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
### manually entered test is run
sub man # arguments: 
{

}
### test in the given directories either all test cases or only the given ones
sub test # arguments: test directory, test identifier, test type 
{       
    # read the ugroot environment variable and save it
    my $lib_source_dir = $ENV{"UGROOT"};
    # get the test directories
    my @testdirs;
    @testdir_parts = split ' ', $_[0];
    foreach $dir (@testdir_parts)
    {
        push(@testdirs, join('','/home/ugtest/Dart/Source/Client/UG/',$dir));
    }
    # get the test identifiers
    @test_identifiers = split ' ', $_[1];
    # get all test cases for application tests which match the given test identifiers 
    foreach $dir (@testdirs)
    {
        my @testcases = read_case($dir);
        my $hash_elem_count = splice(@testcases, -1);
        my $n = splice(@testcases, -1);
        for(my $i = 1;$i <= $n;$i++)
        {
            %testcase = splice(@testcases, 0, $hash_elem_count);
            foreach $testidentifier (@testidentifiers)
            {
                if($testidentifier eq $testcase{"test_id"})
                {
                    push(@array_of_testcase_hashes, %testcase);    
                }
            }
        }
    }
    # get all test cases for building the libraries
    my @libtestcases = read_case($lib_source_dir);
    my $hash_elem_count = splice(@libtestcases, -1);
    my $n = splice(@libtestcases, -1);
    for(my $i = 1;$i <= $n;$i++)
    {
        %libtestcase = splice(@libtestcases, 0, $hash_elem_count);
        push(@array_of_libtestcase_hashes, %testcase); 
    }
    # compare the test cases for application tests and for building the libraries
    foreach $libtestcase (@array_of_libtestcase_hashes)
    {
        %libtestcase = $libtestcase;
        my $lib_build = 0;
        foreach $testcase (@array_of_testcase_hashes)
        {
            %testcase = $testcase;
            if($testcase{"config_string"} eq $libtestcase{"config_string"})
            {
                if($lib_build == 0)
                {
                    dart_test($_[2], %libtestcase);
                    $lib_build = 1;
                }
                dart_test($_[2], %testcase);
            }
        }
    }
    
}
### build the libraries either with all configurations or only with
### the given configurations
sub lib # arguments: test identifier, test type
{
    if($_[0] eq "")
    {
        # read the ugroot environment variable and save it
        my $lib_source_dir = $ENV{"UGROOT"};
        # get the test cases for building the library
        my @libtestcases = read_case($lib_source_dir);
        my $hash_elem_count = splice(@libtestcases, -1);
        my $n = splice(@libtestcases, -1);
        for(my $i = 1;$i <= $n;$i++)
        {
            %libtestcase = splice(@libtestcases, 0, $hash_elem_count);
            # build the libraries for the configuration
            dart_test($_[1], %libtestcase);
        }
    }
    else
    {
        # get the test identifiers for the tests to be executed
        my @test_identifiers = split ' ', $_[0];
        # read the ugroot environment variable and save it
        my $lib_source_dir = $ENV{"UGROOT"};
        # get all the test cases for building the library given in the TestCases file
        my @libtestcases = read_case($lib_source_dir);
        my $lib_hash_elem_count = splice(@libtestcases, -1);
        my $lib_n = splice(@libtestcases, -1);
        for(my $lib_i = 1;$lib_i <= $lib_n;$lib_i++)
        {
            %libtestcase = splice(@libtestcases, 0, $lib_hash_elem_count);
            # build the libraries which match the given test identifiers
            foreach $test_identifier (@test_identifiers)
            {
                if($test_identifier eq $libtestcase{"test_id"})
                {
                    dart_test($_[1], %libtestcase);
                }
            }
        }
    }
}
### execute all tests (all test cases in all directories)
sub test_all # argument: test type 
{
    # read the ugroot environment variable and save it
    my $lib_source_dir = $ENV{"UGROOT"};
    # get the test directories
    my @testdirs = read_testdirs();
    my @array_of_testcase_hashes;
    # get all test cases for application tests 
    foreach $dir (@testdirs)
    {
        my @testcases = read_case($dir);
        my $hash_elem_count = splice(@testcases, -1);
        my $n = splice(@testcases, -1);
        for(my $i = 1;$i <= $n;$i++)
        {
            %testcase = splice(@testcases, 0, $hash_elem_count);
            push(@array_of_testcase_hashes, %testcase); 
        }
    }
    # get all test cases for building the libraries
    my @libtestcases = read_case($lib_source_dir);
    my $hash_elem_count = splice(@libtestcases, -1);
    my $n = splice(@libtestcases, -1);
    for(my $i = 1;$i <= $n;$i++)
    {
        %libtestcase = splice(@libtestcases, 0, $hash_elem_count);
        push(@array_of_libtestcase_hashes, %testcase); 
    }
    # compare the test cases for application tests and for building the libraries
    foreach $libtestcase (@array_of_libtestcase_hashes)
    {
        %libtestcase = $libtestcase;
        my $lib_build = 0;
        foreach $testcase (@array_of_testcase_hashes)
        {
            %testcase = $testcase;
            if($testcase{"config_string"} eq $libtestcase{"config_string"})
            {
                if($lib_build == 0)
                {
                    dart_test($_[0], %libtestcase);
                    $lib_build = 1;
                }
                dart_test($_[0], %testcase);
            }
        }
    }
}
### execute all possible builds of the libraries and all possible tests
sub all # argument: test type
{
    # read the ugroot environment variable and save it
    my $lib_source_dir = $ENV{"UGROOT"};
    # get the test directories
    my @testdirs = read_testdirs();
    # get the test cases for building the library
    my @libtestcases = read_case($lib_source_dir);
    my $lib_hash_elem_count = splice(@libtestcases, -1);
    my $m = splice(@libtestcases, -1);
    # build all
    for(my $k = 1;$k <= $m;$k++)
    {
        %libtestcase = splice(@libtestcases, 0, $lib_hash_elem_count);
        # build the libraries for the configuration saved in %libtestcase
        dart_test($_[0], %libtestcase);
        # execute the tests whose configuration matches the lib configuration
        foreach $dir (@testdirs)
        {
            # get all test cases
            my @testcases = read_case($dir);
            my $hash_elem_count = splice(@testcases, -1);
            my $n = splice(@testcases, -1);
            for(my $i = 1;$i <= $n;$i++)
            {
                # extract one test case
                %testcase = splice(@testcases, 0, $hash_elem_count);
                # execute the test 
                if($testcase{"config_string"} eq $libtestcase{"config_string"})
                {
                    dart_test($_[0], %testcase);
                }
            }
        }
    }
}
### get the test directories
sub read_testdirs
{
    # read the dartroot environment variable and build the testroot directory
    my $testroot_dir = join('',$ENV{"DART_HOME"},'/Source/Client/UG');
    # open the DartTestFiles file to get directories where TestCases files reside
    open(TESTFILES, join('',$testroot_dir,'/DartTestFiles'));
    while(<TESTFILES>)
    {
        # cut off the "end-of-line" (\n)
        chomp;
        # save the whole directory in @testdirs
        push(@testdirs, join('/',$testroot_dir,$_));
    }
    close(TESTFILES);
    return @testdirs;    
}
### get the test cases
sub read_case # argument: TestCases file entry
{
    my @testcases;
    my $i = 0;
    # open the testcase file and read the entries
    open(TESTCASE, join('',$_[0],'/TestCases'));
    while(<TESTCASE>)
    {
        # cut off the "end-of-line" (\n)
        chomp;
        # divide the test case into test identifier, test application
        # and configuration string
        my @testcase_array = split ':', $_;
        # save the test case in %testcase
        my %testcase =();
        $testcase{"test_id"} = $testcase_array[0];
        # cut off the "(" at the beginning and the ")" at the end
        my @test_app = split(//, $testcase_array[1]); 
        $testcase{"test_app"} = substr($testcase_array[1],1,$#test_app -1);
        $testcase{"config_string"} = $testcase_array[2];
        # divide the configuration string into architecture and other
        # ugconf parameter
        my @ugconf_params = split '-', $testcase_array[2];
        $testcase{"arch"} = splice(@ugconf_params, 0, 1);
        $testcase{"conf_params"} = "";
        # save the test directory
        $testcase{"test_dir"} = $_[0];
        # build the configure parameter list
        foreach $param (@ugconf_params)
        {
            $testcase{"conf_params"} = join('',$testcase{"conf_params"},$param,' ');
        }
        # save the test case in @testcases
        push(@testcases, %testcase);
        $i++;
    }
    close(TESTCASE);
    $hash_elem_count = ($#testcases + 1)/$i;
    push(@testcases, $i);
    push(@testcases, $hash_elem_count);
    return @testcases;
}
### execute the necessary commands to perform the build-test-submit cycle
sub dart_test # arguments: test type, testcase(hash)
{
    # get the test type
    my $test_type = splice(@_,0,1);
    # save the testcase as hash
    my %case = @_;
    # coverage: yes or no
    if($case{"arch"} =~ /GCOV/)
    {
        my $coverage_param = "COV";
    }
    else 
    {
        my $coverage_param = "NOCOV";
    }
    # save the hash values in variables in order to have the possibility
    # to pass them to the other scripts
    my $test_dir = $case{"test_dir"};
    my $test_id = $case{"test_id"};
    my $arch = $case{"arch"};
    my $test_app = $case{"test_app"};
    my $conf_params = $case{"conf_params"};
    # build-test-submit cycle
    system("ug_dart_conf.pl -b $test_dir -i $test_id -a $arch -c '$conf_params'");
    # execute ug_dart_tests.pl and call Test.tcl only if a test application is given
    if($test_app eq "")
    {
        $with_test = "no";
        system("ug_dart_client_test $test_type $test_dir $coverage_param $with_test");
    }
    else
    {
        $with_test = "yes";
        system("ug_dart_tests.pl -b $test_dir -s $test_dir -p '$test_app' -i $test_id");
        system("ug_dart_client_test $test_type $test_dir $coverage_param '$with_test'");
    }
}
### print help message 
sub help
{
    print "usage:	ug_dart.pl [-b build_dir_suffix] [-t test_type] [-m mode]\n";
    print "[-p test application] [-o test application options] [-i test id]\n"; 
    print "[-a architecture] [-c ugconf options] [-e test extent]\n";
    print "-e: Test extent.May be\n";
    print "      - manually: execute a given test\n";
    print "      - test    : execute either all tests in the given directory\n";
    print "                  or only several test cases\n";
    print "      - lib     : building libraries with either all configurations\n";
    print "                  or only with several configurations\n";
    print "      - test_all: execute all tests in all directories\n";
    print "      - all     : build the libraries with all configurations\n";
    print "                  execute all tests\n";
    print "-b: Specifying another build directory than .../UG by entering the\n";
    print "    sub directory, e.g. cd/appl\n";
    print "-s: Specifying another source directory than .../UG by entering\n";
    print "    the sub directory, e.g. cd/appl\n";  
    print "-t: May be Nightly, Continuous or Experimental\n";
    print "-i: Test identifier\n";
    print "-p: Specifying the application to test with optional parameters\n";
    print "    and scripts\n"; 
    print "-c: Options which ugconf would accept\n";
    print "\n";
    print "purpose: running all scripts which are necessary for a complete dart\n";
    print "         build-test-submit cycle\n";
}

###############################################################################
# main
###############################################################################
# read the command line parameters
my %option = ();
getopts("b:s:t:p:e:i:a:c:h", \%option);
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
# distinguish between the different entries of the test_extent variable
if($option{e} eq "manually")
{
    # exit if no test directory, test identifier,
    # architecture or ugconf parameters is given
    unless($option{b})
    {
        die "No test directory given!\n";
    }
    unless($option{i})
    {
        die "No test identifier given!\n";
    }
    unless($option{a})
    {
        die "No architecture given!\n";
    }
    unless($option{c})
    {
        die "No configuration string given!\n";
    }
    # if no source directory is given it is the same as the test directory
    unless($option{s})
    {
        $option{s} = $option{b}
    }
    man($option{b}, $option{s}, $option{t}, $option{p}, $option{i}, $option{a}, $option{c});
}
elsif($option{e} eq "test")
{
    # exit if no test directory is given
    unless($option{b})
    {
        die "No test directory given!\n";
    }
    # set test identifier to "" if it is not given
    unless($option{i})
    {
        $option{i} = "";
    }
    test($option{b}, $option{i}, $option{t});
}
elsif($option{e} eq "lib")
{
    # set test identifier to "" if it is not given
    unless($option{i})
    {
        $option{i} = "";
    }
    lib($option{i}, $option{t});
}
elsif($option{e} eq "test_all")
{
    test_all($option{t});
}
elsif($option{e} eq "all")
{
    all($option{t});
}
else 
{
    die, "The test extent you entered isn't defined!\n";
}

###############################################################################
