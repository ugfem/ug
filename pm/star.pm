=head1 NAME

d3f - Distributed Density Driven Flow

=head1 DESCRIPTION

This module provides high-level commands of the software-package
d3f.

=head2 Functions

=cut

package star;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use Exporter;
$VERSION = 1.0;
@ISA = qw(Exporter);

@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();

##############################################
## source 
##############################################

use Term::ANSIColor qw(:constants);
use Math::Complex;

sub float
{
    my $real='[+-]?\d+\.?\d*[eE]?[+-]?\d+|[+-]?\d*\.?\d+[eE]?[+-]?\d+|[+-]?\d+';
    my (@list,$f,$s,$in);

    if (@_==1) { @list=grep /$real/,split /($real)/,$_[0]; }
    elsif (@_==2)
    {
        $in=' '.$_[1];
        ($f,$s)=split /$_[0]/,$in,2;
        @list=grep /$real/,split /($real)/,$s;
    }
    else
    {
        $in=' '.$_[2];
        ($f,$s)=split /$_[0]/,$in,2;
        @list=grep /$real/,split /($real)/,$s;
        if ($_[1]>=@list || $_[1]<0) { return undef;}
        for ($s=0; $s<$_[1]; $s++) { shift @list; }
    }
    return wantarray ? @list : $list[0];
}

sub dim 
{
    @_==1 or die "ERROR in dim: provide 1 argument\n";
    my ($key,$dim,@a);
    $dim=-1;
    for $key (keys %{$_[0]})
    {
        if ($dim==-1) { @a=float $key; $dim=@a; }
        else { @a=float $key; if ($dim!=@a) { die "ERROR in s_dim: key-dimension mismatch\n"; } }
    }
    return $dim;
}

sub multiply
{
    @_==2 or die "ERROR in star_multiply: provide 2 arguments\n";
    my $a=shift; my $b=shift;
    my $dim=dim $a; if ($dim!=dim $b) { die "ERROR in s_multiply: dimensions of the arguments mismatch\n"; }

    SWITCH:
    {
        if ($dim==1)
        {
            my ($k1,$k2,$k);
            my %c=();
            for $k1 (keys %$a)
            {
                for $k2 (keys %$b)
                {
                    $k=$k1+$k2;
                    if (defined $c{$k}) { $c{$k}+=$$a{$k1}*$$b{$k2}; }
                    else { $c{$k}=$$a{$k1}*$$b{$k2}; }
                }
            }
            return \%c;
        }
        if ($dim==2)
        {
            my ($k1,$k2,$k,@a1,@a2,@a);
            my %c=();
            for $k1 (keys %$a)
            {
                @a1=float $k1;
                for $k2 (keys %$b)
                {
                    @a2=float $k2;
                    $a[0]=$a1[0]+$a2[0]; $a[1]=$a1[1]+$a2[1]; $k="$a[0] $a[1]";
                    if (defined $c{$k}) { $c{$k}+=$$a{$k1}*$$b{$k2}; }
                    else { $c{$k}=$$a{$k1}*$$b{$k2}; }
                }
            }
            return \%c;
        }
        die "ERROR in star_print: dimension $dim not implemented yet\n";
    }
}

sub subtract
{
    @_==2 or die "ERROR in subtract: provide 2 arguments\n";
    my $a=shift; my $b=shift;
    my $dim=dim $a; if ($dim!=dim $b) { die "ERROR in subtract: dimensions of the arguments mismatch\n"; }
    my $k;
    my %c=();
    for $k (keys %$a)
    {
        if (defined $c{$k}) { $c{$k}+=$$a{$k}; }
        else { $c{$k}=$$a{$k}; }
    }
    for $k (keys %$b)
    {
        if (defined $c{$k}) { $c{$k}-=$$b{$k}; }
        else { $c{$k}=-$$b{$k}; }
    }
    return \%c;
}

sub transpose
{
    @_==1 or die "ERROR in star_transpose: provide 1 argument\n";
    my $a=shift; my $dim=dim $a;
    my %c=();
    SWITCH:
    {
        if ($dim==1)
        {
            my ($k,$l);
            my %c=();
            for $k (keys %$a)
            {
                $l=-$k;
                $c{$l}=$$a{$k};
            }
            return \%c;
        }
        if ($dim==2)
        {
            my ($k,$l,@b);
            my %c=();
            for $k (keys %$a)
            {
                @b=float $k;
                $b[0]*=-1; $b[1]*=-1;
                $l="$b[0] $b[1]";
                $c{$l}=$$a{$k};
            }
            return \%c;
        }
        die "ERROR in star_transpose: dimension $dim not implemented yet\n";
    }
}

sub scale
{
	@_==2 or die "ERROR in star_scale: provide 2 arguments\n";
	my $a=shift; my $dim=dim $a;
	my $scale=shift;
	my %r=();my $k;
	for $k (keys %$a)
	{
		$r{$k}=$scale*$$a{$k};
	}
	return \%r;
}

sub restrict
{
    @_<1 and die "ERROR in star_restict: provide 1+dim arguments\n";
    my $a=shift; my $dim=dim $a;
    if (@_!=$dim) { die "ERROR in star_restict: provide 1+dim arguments\n"; }
    my ($copy,$i,$k,@a);
    my %c=();
    for $k (keys %$a)
    {
        @a=float $k;
        $copy=1;
        for ($i=0; $i<$dim; $i++) { if (abs($a[$i])>$_[$i]) {$copy=0; } }
        if ($copy) { $c{$k}=$$a{$k}; }
    }
    return \%c;
}

sub norm
{
    @_==1 or die "ERROR in norm: provide 1 argument\n";
    my $a=shift; dim $a;
    my $norm=0;
    my $k;
    for $k (keys %$a)
    {
        if (defined $$a{$k}) { $norm+=$$a{$k}*$$a{$k}; }
    }
    return sqrt($norm);
}

sub print
{
    @_==1 or die "ERROR in star_print: provide 1 argument\n";
    my $a=shift;
    my $dim=dim $a;
    SWITCH:
    {
        if ($dim==1)
        {
            my ($key,$i);
            my $max=0; for $key (keys %$a) { if ($max<$key || $max<-$key) { $max=abs($key); } }
            print "| ";
            for ($i=-$max; $i<=$max; $i++)
            {
                if ($i!=0) { print RED; }
                if (defined $$a{$i}) { printf "%+e ",$$a{$i}; }
                else { printf "%+e ",0; }
                if ($i!=0) { print RESET; }
            }
            print "|\n";
            last SWITCH;
        }
        if ($dim==2)
        {
            my (@ar,$key,$i,$j,$max0,$max1);
            $max0=$max1=0;
            for $key (keys %$a) { @ar=float $key; if ($max0<$ar[0] || $max0<-$ar[0]) { $max0=abs($ar[0]); } if ($max1<$ar[1] || $max1<-$ar[1]) { $max1=abs($ar[1]); } }
            print "| ";
            for ($j=$max1; $j>=-$max1; $j--)
            {
            	for ($i=-$max0; $i<=$max0; $i++)
                {
                    if ($i!=0 || $j!=0) { print RED; }
                    $key="$i $j";
                    if (defined $$a{$key}) { printf "%+e ",$$a{$key}; }
                    else { printf "%+e ",0; }
                    if ($i!=0 || $j!=0) { print RESET; }
                }
                print "|\n";
                if ($j>-$max1) { print "| "; }
            }
            last SWITCH;
        }
        die "ERROR in s_print: dimension $dim not implemented yet\n";
    }
}

sub f
{
    @_==1 or die "ERROR in star_f: provide 1 argument\n";
    my $a=shift;
	my ($key,$f,$eta,$c,$sign_eta,$sign_c,$v,$i,$nu,$sign_nu);
    my $dim=dim $a;
    SWITCH:
    {
        if ($dim==1)
        {
			$f='sub {my $eta=shift;return(';
			for $key (keys %$a) 
			{ 
				$c=$$a{$key}; if ($c<0) { $sign_c=-1; } elsif ($c==0) { next; } else { $sign_c=1; } $c=abs($c);
				$eta=float $key; $eta*=3.1415926535; if ($eta<0) { $sign_eta=-1; } elsif ($eta==0) {$sign_eta=0; } else { $sign_eta=1; } $eta=abs($eta);
				if ($sign_eta<0)
				{
					if ($sign_c<0) 	{ $f.="-$c*exp(-$eta*".'$eta*i)'; }
					else 			{ $f.="+$c*exp(-$eta*".'$eta*i)'; }
				}
				elsif ($sign_eta==0)
				{
					if ($sign_c<0) 	{ $f.="-$c"; }
					else 			{ $f.="+$c"; }
				}
				else
				{
					if ($sign_c<0) 	{ $f.="-$c*exp($eta*".'$eta*i)'; }
					else 			{ $f.="+$c*exp($eta*".'$eta*i)'; }
				}
			}
		}
        if ($dim==2)
        {
			$f='sub {my $eta=shift;my $xi=shift;return(';
			for $key (keys %$a) 
			{ 
				$c=$$a{$key}; if ($c<0) { $sign_c=-1; } elsif ($c==0) { next; } else { $sign_c=1; } $c=abs($c);
				($eta,$nu)=float $key; 
				$eta*=3.1415926535; if ($eta<0) { $sign_eta=-1; } elsif ($eta==0) {$sign_eta=0; } else { $sign_eta=1; } $eta=abs($eta);
				$nu*=3.1415926535; if ($nu<0) { $sign_nu=-1; } elsif ($nu==0) {$sign_nu=0; } else { $sign_nu=1; } $nu=abs($nu);
				if ($sign_eta<0)
				{
					if ($sign_c<0) 	{ $f.="-$c*exp(-$eta*".'$eta*i)'; }
					else 			{ $f.="+$c*exp(-$eta*".'$eta*i)'; }
				}
				elsif ($sign_eta==0)
				{
					if ($sign_c<0) 	{ $f.="-$c"; }
					else 			{ $f.="+$c"; }
				}
				else
				{
					if ($sign_c<0) 	{ $f.="-$c*exp($eta*".'$eta*i)'; }
					else 			{ $f.="+$c*exp($eta*".'$eta*i)'; }
				}
				if ($sign_nu<0)
				{
					$f.="*exp(-$nu*".'$nu*i)'; 
				}
				elsif ($sign_nu>0)
				{
					$f.="*exp($nu*".'$nu*i)'; 
				}
			}
		}
	}
	$f.=');}'; $f=eval($f);
	return $f;
}

1;










