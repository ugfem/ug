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

	sort_numberlist_LINE();

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

function sort_numberlist_LINE ()
{

	for (L=1; L<=NF; L++) {

		# alpha or num
		alpha = match(LINE[L] , "[a-z|A-Z|:|=]");

		# start of a number list
		if (alpha == 0) {

			# get end of number list
			for (STOP=L+1; (STOP<=NF) && (alpha==0); STOP++) {

				alpha = match(LINE[STOP] , "[a-z|A-Z|:|=]");

			}

			min_sort(L,STOP-1);
		}
	}
}

function print_LINE ()
{
	for (L=1; L<=NF; L++)
		print LINE[L];
	printf("\n");
}
