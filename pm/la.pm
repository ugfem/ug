package la; use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw(la_mv la_inv la_p la_spec);
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();

##############################################
## source 
##############################################
use IO::File;
use POSIX qw(tmpnam);

sub tmpfile
{
    my $name;
    do { $name=tmpnam(); } until IO::File->new($name,O_RDWR|O_CREAT|O_EXCL);
    return $name;
}

sub float
{
    my $real='[+-]?\d+\.?\d*[eE]?[+-]?\d+|[+-]?\d*\.?\d+[eE]?[+-]?\d+|[+-]?\d+';
    my (@list,$f,$s);

    if (@_==1) { @list=grep /$real/,split /($real)/,$_[0]; }
    else
    {
        ($f,$s)=split /$_[0]/,$_[1],2;
        @list=grep /$real/,split /($real)/,$s;
    }
    return wantarray ? @list : $list[0];
}

sub la_mv
{
	my ($m,$v,$n,$i,$j,$out);

	$m=shift;
	$v=shift;
    $n=$#{$m}+1;
	for ($i=0; $i<$n; $i++)
    {
		$out->[$i]=0;
        for ($j=0; $j<$n; $j++)
        {
			$out->[$i]+=$m->[$i][$j]*$v->[$j];
		}
	}
	return $out;
}

sub la_inv
{
	my ($i,$j,$k,$n,$m,$d,$f,$out,$in,$v);

	$in=shift; 
	$n=$#{$in}+1;
	for ($i=0; $i<$n; $i++)
    {
		for ($j=0; $j<$n; $j++)
		{
			$m->[$i][$j]=$in->[$i][$j];
			if ($i==$j) {$out->[$i][$j]=1;}
			else {$out->[$i][$j]=0;}
		}
	}
	for ($i=0; $i<$n-1; $i++)
	{
		$d=$m->[$i][$i] or die "cannot invert matrix\n";
		for ($j=$i+1; $j<=$n-1; $j++)
		{
			$f=$m->[$j][$i]/$d;
			$m->[$j][$i]=$f;
			for ($k=$i+1; $k<=$n-1; $k++)
			{
				$m->[$j][$k]-=$f*$m->[$i][$k];
			}
	 	}
	}
	for ($v=0; $v<$n; $v++)
	{
		for ($i=1; $i<$n; $i++)
		{
			for ($j=0; $j<$i; $j++)
			{
				$out->[$i][$v]-=$m->[$i][$j]*$out->[$j][$v];
			}
		}
		for ($i=$n-1; $i>=0; $i--)
		{
			for ($j=$i+1; $j<=$n-1; $j++)
			{
				$out->[$i][$v]-=$m->[$i][$j]*$out->[$j][$v];
			}
			$m->[$i][$i]!=0 or die "cannot invert matrix\n";
			$out->[$i][$v]/=$m->[$i][$i];
		}
	}
	return $out;
}

sub la_p
{
    my ($i,$j,$n,$m,$ret);

	$ret='';
    $m=shift; 
    $n=$#{$m}+1; 
	if (ref $m->[0])
	{
    	for ($i=0; $i<$n; $i++)
    	{
			for ($j=0; $j<$n; $j++)
			{
				$ret.=sprintf ("%15e ",$m->[$i][$j]);
			}
			$ret.="\n";
		}
	}
	else
	{
        for ($i=0; $i<$n; $i++)
        {
			$ret.=sprintf ("%15e\n",$m->[$i]);
		}
	}
	$ret.="\n";
}

sub la_spec
{
	my ($i,$j,$m,$s,$n,$name,@r,@f,@real,@imag);

	$m=shift;
	$n=$#{$m}+1;
	$s='a={';
	for ($i=0; $i<$n; $i++)
	{
		for ($j=0; $j<$n; $j++)
		{
			$s.="$m->[$i][$j]";
			if ($j<$n-1) {$s.=','}
			elsif ($i<$n-1 && $j==$n-1) {$s.=";\n"}
		}
	}
	$s.="};\nspec(a)\n"; $name=tmpfile; 
	open(FOO,">$name"); print FOO $s; close(FOO);
	$s=`cat $name | scilab -nw`; `rm $name`;
	($s,$s)=split /ans/,$s; 
	$s=~s/\s+/ /g; $s=~s/! !/!/g; $s=~s/i//g; $s=~s/\+\s*//g; $s=~s/\-\s*/-/g;
	@f=split /!/,$s; shift @f;
	for ($i=0; $i<@f-1; $i++)
	{
		@r=float $f[$i];
		if (@r==1) {$real[$i]=$r[0];$imag[$i]=0;}
		elsif (@r==2) {$real[$i]=$r[0];$imag[$i]=$r[1];}
		else {die "cannot evaluate spectrum from '$s'\n";}
	}
	return(\@real,\@imag);
}


1;












