BEGIN { 

IdentHdrBegin = "IdentHdr:"
IdentHdrEnd = "ProcList:"
ORS = " "
n = 0;
}


/IdentObjectHdr:/ {

	# remove pvm output prefix
	#sub (^[ ]*\[.*\][ ]*,"");

	# remove processor number
	#sub (^[ ]*[0-9]*:[ ]*,"");

	# split one line using whitespace as separator
	split ($0,LINE)

	# search start of IdentHdr field
	for (POS=1; POS<=NF; POS++)
		if (LINE[POS] == IdentHdrBegin)
			break;

	# search end of IdentHdr field
	POS++;
	for (POSEND=POS; POSEND<=NF; POSEND++)
		if (LINE[POSEND] == IdentHdrEnd)
			break;

	# sort IdentHdr entries
	min_sort(POS,POSEND);

	# sort ProcList entries
	min_sort(POSEND+1,NF);

	# remember last line
	split (LINE[0],PREVLINE);

	# print LINE 
	print_LINE();

	# increment line counter
	n++;
}

END { }

# functions

function min_sort (from, to)
{
	if (from > to) {
		return 1;
	}

	if (from == to) return 0;

	# search the minimum
	for (sortpos=from; sortpos<to; sortpos++) {
		min = sortpos;
		for (pos=sortpos; pos<=to; pos++) {
			if (LINE[min] > LINE[pos])	
				min = pos;
		}
		if (min != sortpos) {
			dummy = LINE[sortpos];
			LINE[sortpos] = LINE[min];
			LINE[min] = dummy;
		}
	}
}

function print_LINE ()
{
	for (L=1; L<=NF; L++)
		print LINE[L];
	printf("\n");
}
