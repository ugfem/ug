package la; use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw(la_tr la_add la_mv la_mm la_inv la_p la_spec la_norm la_scale la_spec_z);
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

sub la_mm
{
	my ($m1,$m2,$n,$i,$j,$k,$out);

	$m1=shift;
	$m2=shift;
    $n=$#{$m1}+1;
	for ($i=0; $i<$n; $i++)
    {
        for ($j=0; $j<$n; $j++)
        {
			$out->[$i][$j]=0;
        	for ($k=0; $k<$n; $k++)
        	{
				$out->[$i][$j]+=$m1->[$i][$k]*$m2->[$k][$j];
			}
		}
	}
	return $out;
}

sub la_tr
{
    my ($m,$n,$i,$j,$out);

    $m=shift;
    $n=$#{$m}+1;
    for ($i=0; $i<$n; $i++)
    {
        for ($j=0; $j<$n; $j++)
        {
            $out->[$i][$j]=$m->[$j][$i];
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
				$ret.=sprintf ("%16.8e ",$m->[$i][$j]);
			}
			$ret.="\n";
		}
	}
	else
	{
        for ($i=0; $i<$n; $i++)
        {
			$ret.=sprintf ("%16.8e\n",$m->[$i]);
		}
	}
	$ret.="\n";
}

sub la_add
{
    my ($i,$j,$n,$m,$mm,$sum,$ret);

	$ret='';
    $m=shift; 
    $mm=shift; 
    $n=$#{$m}+1; 
	if (ref $m->[0])
	{
    	for ($i=0; $i<$n; $i++)
    	{
			for ($j=0; $j<$n; $j++)
			{
				$sum->[$i][$j]=$m->[$i][$j]+$mm->[$i][$j];
			}
		}
	}
	else
	{
        for ($i=0; $i<$n; $i++)
        {
			$sum->[$i]=$m->[$i]+$mm->[$i];
		}
	}
	return $sum;
}

sub la_spec
{
	my ($s1,$s2,$s3,$i,$j,$m,$s,$n,$name,@r,@f,@real,@imag,$b);

	$m=shift;
	$b=shift; 
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
	$s.="};\n";
	if (defined $b)
	{
		$s.='b={';
		for ($i=0; $i<$n; $i++)
		{
			for ($j=0; $j<$n; $j++)
			{
				$s.="$b->[$i][$j]";
				if ($j<$n-1) {$s.=','}
				elsif ($i<$n-1 && $j==$n-1) {$s.=";\n"}
			}
		}
		$s.="};\n";
	}
	if (!defined $b)
	{
		$s.="spec(a)\n"; 
	}
	else
	{
		$s.="[x,y]=gspec(a,b);\nv=x./y;\nreal(v)\nimag(v)\n";
	}
	$name=tmpfile; open(FOO,">$name"); print FOO $s; close(FOO); 
	$s=`cat $name | scilab -nw`; `rm $name`;
	($s1,$s2,$s3)=split /ans/,$s;
	$s2=~s/\-\s+/-/g; $s2=~s/\+\s+/+/g;
	$s3=~s/\-\s+/-/g; $s3=~s/\+\s+/+/g;
	@real=float $s2; @imag=float $s3;
	return(\@real,\@imag);
}

sub la_spec_z
{
	my ($fac,$evr,$evi,$s1,$s2,$ii,$jj,$s3,$s4,$i,$j,$k,$m,$s,$n,$name,@r,@f,@ri,@real,@imag,$b);

	$m=shift;
	$b=shift; 
	$n=$#{$m}+1;
	$s='';
	for ($i=0; $i<$n; $i++)
	{
		for ($j=0; $j<$n; $j++)
		{
			$ii=$i+1; $jj=$j+1;
			$s.="a($ii,$jj)=$m->[$i][$j];\n";
		}
	}
	if (defined $b)
	{
		for ($i=0; $i<$n; $i++)
		{
			for ($j=0; $j<$n; $j++)
			{
				$ii=$i+1; $jj=$j+1;
				$s.="b($ii,$jj)=$b->[$i][$j];\n";
			}
		}
	}
	if (!defined $b)
	{
		$s.="v=spec(a);\nreal(v)\nimag(v)\n"; 
	}
	else
	{
		-e '/tmp/Sci_SR' and `rm /tmp/Sci_SR`;
		-e '/tmp/Sci_SI' and `rm /tmp/Sci_SI`;
		-e '/tmp/Sci_ER' and `rm /tmp/Sci_ER`;
		-e '/tmp/Sci_EI' and `rm /tmp/Sci_EI`;
		$s.="[x,y,z]=gspec(a,b); v=x./y; 
		     write('/tmp/Sci_SR',real(v),'(e32.16)'); 
			 write('/tmp/Sci_SI',imag(v),'(e32.16)'); 
			 write('/tmp/Sci_ER',real(z),'($n(e32.16))'); 
			 write('/tmp/Sci_EI',imag(z),'($n(e32.16))');";
	}
	open(FOO,">/tmp/Sci_script"); print FOO $s; close(FOO); 
	`cat /tmp/Sci_script | scilab -nw`;  `rm /tmp/Sci_script`;
	$s=`cat /tmp/Sci_SR`; @real=float $s;  `rm /tmp/Sci_SR`;
	$s=`cat /tmp/Sci_SI`; @imag=float $s;  `rm /tmp/Sci_SI`;
	$s=`cat /tmp/Sci_ER`; @ri=float $s;  `rm /tmp/Sci_ER`;
	for ($i=$k=0; $i<$n; $i++) { for ($j=0; $j<$n; $j++) { $evr->[$i][$j]=$ri[$k++]; } }
	$s=`cat /tmp/Sci_EI`; @ri=float $s;  `rm /tmp/Sci_EI`;
	for ($i=$k=0; $i<$n; $i++) { for ($j=0; $j<$n; $j++) { $evi->[$i][$j]=$ri[$k++]; } }
	return(\@real,\@imag,$evr,$evi);
}

sub la_norm
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
	$s.="};\nsqrt(max(spec(a*a')))\n"; $name=tmpfile; 
	open(FOO,">$name"); print FOO $s; close(FOO);
	$s=`cat $name | scilab -nw`; `rm $name`;
	return float "ans",$s;
}

sub la_scale
{
	my ($r,$i,$m,$n,$s);

	$s=shift; $m=shift;
	$n=$#{$m}+1;
	for ($i=0; $i<$n; $i++)
    {
    	$r->[$i]=$s*$m->[$i];
    }
	return $r;
}

1;












