package tree;
use strict;
use vars qw (@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use vars qw ($DIM @DEBUG $HUGE $TIGHT_BOXES @TMP1);
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(BBT_NewBBox BBT_NewTree BuildTree NewNode BBoxForBBoxes BBT_TreePointDistance);
@EXPORT_OK = qw();
%EXPORT_TAGS = qw();
@DEBUG=(0,0,0);
$TIGHT_BOXES = 1;

#######################################################################
# BBoxPointDistance2
#######################################################################
sub BBoxPointDistance2
{
	my($bb, $x, $min, $max)=@_;

if($TIGHT_BOXES){
	my($dmax,$dmin,$dll_i, $dur_i,$d,$m,$M,$i,$j);

	# set pointers
	$dmin = $min;
	$dmax = $max;

	for($i=0; $i<$DIM; $i++){
		$TMP1[$i] = 0.0;
	}
	$$dmin = 0.0;
	for($i=0; $i<$DIM; $i++){
		$d = $x->[$i] - $bb->{ll}->[$i];
		$dll_i = $d * $d;
		$d = $x->[$i] - $bb->{ur}->[$i];
		$dur_i = $d * $d;
		if($DEBUG[2]){
			printf("BBoxPointDistance2: d=%f dll_i=%f dur_i=%f\n",$d,$dll_i,$dur_i); }
		if($x->[$i] < $bb->{ll}->[$i]){
			$$dmin += $dll_i;}
		else{
			if($x->[$i] > $bb->{ur}->[$i]){
				$$dmin += $dur_i;
			}else{ $$dmin += 0.0; }
		}
		$m = min($dll_i, $dur_i);
		$M = max($dll_i, $dur_i);
		for($j=0; $j<$DIM; $j++){
			if($j == $i){
				$TMP1[$j] += $m; 
			}else{
				$TMP1[$j] += $M;
			}
		}
	}
	$$dmax = $TMP1[0];
	for($i=0; $i<$DIM; $i++){
		$$dmax = min($$dmax, $TMP1[$i]);
	}

	
}else{
	my ($i,$dmax,$dmin,$dll_i,$dur_i,$d);

	if($DEBUG[2]){
		printf("\nBBoxPointDistance2: x[0]=%f x[1]=%f bb->ll[0]=%f\n",$x->[0],$x->[1],$bb->{ll}->[0]);}

	# set pointers
	$dmin = $min;
	$dmax = $max;

	$$dmax = $$dmin =0.0;
	for($i=0; $i<$DIM; $i++){
		$d = $x->[$i] - $bb->{ll}->[$i];
		$dll_i = $d * $d;
		$d = $x->[$i] - $bb->{ur}->[$i];
		$dur_i = $d * $d;
		if($DEBUG[2]){
			printf("BBoxPointDistance2: d=%f dll_i=%f dur_i=%f\n",$d,$dll_i,$dur_i); }
		if($x->[$i] < $bb->{ll}->[$i]){
			$$dmin += $dll_i;}
		else{
			if($x->[$i] > $bb->{ur}->[$i]){
				$$dmin += $dur_i;
			}else{ $$dmin += 0.0; }
		}
		$$dmax += max($dll_i, $dur_i);
	}
	if($DEBUG[2]){
		printf("BBoxPointDistance2: dmin=%f dmax=%f\n",$$dmin,$$dmax);}
}
}
#######################################################################
# MinMaxBBoxPointDist2 - 
#######################################################################
sub MinMaxBBoxPointDist2 
{
	my ($node,$x,$minmax)=@_;
	my ($lmin,$lmax,$rmin,$rmax);

	if($node->{left} eq "NULL"){
		# AH: TODO
		#if($node->{right} eq "NULL"){ die "node->{right} eq NULL$!\n"; }
		return $minmax;
	}
	BBoxPointDistance2($node->{left}->{bb}, $x, \$lmin, \$lmax);
	if($lmax < $minmax){ $minmax = $lmax; }
	BBoxPointDistance2($node->{right}->{bb}, $x, \$rmin, \$rmax);
	if($rmax < $minmax){ $minmax = $rmax; }
	if($lmin < $rmin){
		if($lmin < $minmax){
			$minmax = MinMaxBBoxPointDist2($node->{left}, $x, $minmax);
			if($rmin < $minmax){
				$minmax = MinMaxBBoxPointDist2($node->{right}, $x, $minmax);
			}
		}
	}else{
		if($rmin < $minmax){
			$minmax = MinMaxBBoxPointDist2($node->{right}, $x, $minmax);
			if($lmin < $minmax){
				$minmax = MinMaxBBoxPointDist2($node->{left}, $x, $minmax);
			}
		}
	}
	return ($minmax);
}
#######################################################################
# ClosestBBoxesToPoint - find a minimum set of boxes, of which one of
# them will include the nearest object to a point
#######################################################################
sub ClosestBBoxesToPoint
{
	my ($node, $x, $func, $bypass, $minmax)=@_;
	my ($min, $max);

#	print"ClosestBBoxesToPoint:\n";
	if($node->{left} eq "NULL"){
		# AH: TODO
		 #if($node->{right} eq "NULL"){die}
#		print"ClosestBBoxesToPoint: x[0]=$x->[0] x[1]=$x->[1]\n";
		$func->($node->{bb}->{obj}, $bypass);
		return;
	}
	BBoxPointDistance2($node->{left}->{bb}, $x, \$min, \$max);
	if($min < $$minmax){
		ClosestBBoxesToPoint($node->{left}, $x , $func, $bypass, $minmax);
	}
	BBoxPointDistance2($node->{right}->{bb}, $x, \$min, \$max);
	if($min < $$minmax){
		ClosestBBoxesToPoint($node->{right}, $x , $func, $bypass, $minmax);
	}

}
#######################################################################
# BBoxForBBoxes - creates a new bbox containing all n bboxes
#######################################################################
sub BBoxForBBoxes
{
	my ($bboxes, $n)=@_;
	my ($bbox,$i,$j);

	if($n<1){return}
	$bbox = BBT_NewBBox($DIM,$bboxes->[0]->{ll}, $bboxes->[0]->{ur},"NULL");
	for($i=1; $i<$n; $i++){
		for($j=0; $j<$DIM; $j++){
			if($bboxes->[$i]->{ll}->[$j] < $bbox->{ll}->[$j]){
				$bbox->{ll}->[$j] = $bboxes->[$i]->{ll}->[$j] }
			if($bboxes->[$i]->{ur}->[$j] > $bbox->{ur}->[$j]){
				$bbox->{ur}->[$j] = $bboxes->[$i]->{ur}->[$j] }
		}
	}
	return($bbox);	
}

#######################################################################
# NewNode - creates a new bounding box tree node
#######################################################################
sub NewNode
{
	my ($bbox)=@_;
	my %newNode;

	$newNode{left} = $newNode{right} = "NULL";
	$newNode{bb} = $bbox;
	return(\%newNode);

}
#######################################################################
# BuildTree - builds a new bounding box tree
#
# PARAMETERS: 
# bboxes - array of pointers to bounding boxes to build tree from 
# n		 - #bounding boxes
#######################################################################
sub BuildTree
{
	my ($bboxes,$n)=@_;
	my ($bbox,$node,@lbboxes,@rbboxes,$ext,$maxext,$cut,$i,$dir,$nl,$nr);
	my ($hrbboxes,$hlbboxes);
	my (@testr,@testl);

	#print"BUILDTREE:\n";
	#@testr=@$bboxes;
	#print"n=$n \$#testr=$#testr\n";
	if($n<1){return}		
	if(!defined($bboxes)){die "bounding boxes not defined \n $!"}
	if(!defined($bboxes->[0])){die "bounding box bbox[0] not defined \n $!"}
	if($n==1){
		return (NewNode($bboxes->[0]));}

	$bbox = BBoxForBBoxes($bboxes, $n);
	$node = NewNode($bbox);

	## split set of bounding boxes by location of the center points
	$maxext = 0.0;
	for($i=0; $i<$DIM; $i++){
		$ext = $bbox->{ur}->[$i] - $bbox->{ll}->[$i];
		if($ext > $maxext){
			$maxext = $ext;
			$dir = $i;
		}
	}
	if($maxext <= 0.0){die "maxext=$maxext must be > 0 $!\n";}
	$cut = ($bbox->{ur}->[$dir] + $bbox->{ll}->[$dir]) / 2.0;
	if($DEBUG[1]){
		printf("BuildTree: %f\n",$cut);}

	$nr = $nl = 0;
	for($i=0; $i<$n; $i++){
		if((($bboxes->[$i]->{ll}->[$dir] + $bboxes->[$i]->{ur}->[$dir])/2.0) > $cut){
			$rbboxes[$nr++] = $bboxes->[$i];
		}else{
			$lbboxes[$nl++] = $bboxes->[$i];
		}
	}
	$hrbboxes = \@rbboxes;
	$hlbboxes = \@lbboxes;
	#print"nl=$nl nr=$nr\n";
	if($nr == 0){
		#print"NR==0\n";
		$nl = sprintf("%.0f",$nl/2);
		$nr = $n - $nl;
		@rbboxes = popn(\@lbboxes,$nr);
		#@testr=@$hlbboxes;
		#print"\$#lbboxes=@lbboxes \$#rbboxes=$#rbboxes \$#testr=$#testr nl=$nl\n";
	}elsif($nl == 0){
		#print"NL==0\n";
		$nr = sprintf("%.0f",$nr/2);
		$nl = $n - $nr;
		@rbboxes = popn(\@lbboxes,$nl);
	}
	#@testr=@$hrbboxes;
	#@testl=@$hlbboxes;
	
	# extend the tree 
	if($DEBUG[1]){
		#@testl=@$hlbboxes;
		#printf("left nl=%d nr=%d \$#testl=%d\n",$nl,$nr,$#testl);
		}	
	$node->{left} = BuildTree($hlbboxes, $nl);
	if($DEBUG[1]){
		#@testr=@$hrbboxes;
		#printf("right nr=%d \$#testr=%d\n",$nr,$#testr);
		}	
	$node->{right} = BuildTree($hrbboxes, $nr);
	
	return ($node);
}

sub TreePointDistanceCallback
{
	my ($obj,$bypass) = @_;
	my ($d,$dist);	

	## function outside tree.pm
	$dist = $bypass->{dist};
	$d = $dist->($bypass->{x},$obj);
	if($d <= $bypass->{min}){
		$bypass->{min} = $d;
		$bypass->{obj} = $obj;
	}
}
#######################################################################
# BBT_NewTree - builds a new bounding box tree from an array of 
# 				pointers to bounding boxes
#
# PARAMETERS: 
# bboxes - array of pointers to bounding boxes to build tree from 
# n		 - #bounding boxes
# dim 	 - dimension	
#######################################################################
sub BBT_NewTree
{
	my ($bboxes, $n, $dim)=@_;	
	my %newTree;

	$newTree{dim} = $dim;
	$newTree{n} = 0;
	$newTree{root} = BuildTree($bboxes, $n);
	return (\%newTree);
}
#######################################################################
# BBT_NewBBox - Creates a new bounding box of dimension dim
#
# PARAMETERS: 
# dim - dimension	
# ll	- lower left corner
# ur	- upper right corner
# obj	- contained object			
#######################################################################
sub BBT_NewBBox
{
	my ($dim,$ll,$ur,$obj)=@_;
	my (@myll,@myur,%newBB,$i);

	$DIM = $dim;
	
	for($i=0; $i<$dim; $i++){
		$myll[$i] = $ll->[$i];
		$myur[$i] = $ur->[$i];
		if($DEBUG[0]){
			print"BBT_NewBBox: $ll->[$i] $ur->[$i]\n";} 
	}
	
	$newBB{dim} = $dim;
	$newBB{obj} = $obj;
	$newBB{ll} = \@myll;
	$newBB{ur} = \@myur;
	
	return (\%newBB);
}
#######################################################################
# BBT_TreePointDistance - call function for a small number of 	
# boxes, of which at least one contains the nearest object to the
# given point x
#
# PARAMETERS:
# tree - bounding box tree
# x	   - point
#######################################################################
sub BBT_TreePointDistance
{
	my ($tree,$x,$obj,$dist,$HUGE)=@_;
	my ($minmax,%bypass);

	if(!defined($tree)){die "Structure tree is not defined $!\n"}
	$DIM = $tree->{dim};
	$minmax = MinMaxBBoxPointDist2($tree->{root},$x, $HUGE); 
	#printf("minmax=%f\n",sqrt($minmax));
	$bypass{min} = $HUGE;
	## Callback function outside tree.pm
	$bypass{dist} = $dist;
	$bypass{x} = $x;
	$bypass{obj} = "NULL";
	ClosestBBoxesToPoint($tree->{root}, $x, \&TreePointDistanceCallback, \%bypass, \$minmax);
#	printf "BBT_TreePointDistance: obj=$bypass{obj}\n";
	$$obj = $bypass{obj};
	return $bypass{min}
}

sub max
{
	my ($a,$b)=@_;
	if($a>$b){return $a}
	else{return $b}
}
sub min
{
	my ($a,$b)=@_;
	if($a<$b){return $a}
	else{return $b}
}
sub popn
{
	my ($a,$n)=@_;
	return(splice(@$a,-$n));
}


1;



