########################################################################
#
#   bcheck_fix1.awk
#
#   this awk-script replaces spaces in the VALUES part of State-lines
#   by underscores in order to allow correct version consistence checking
#   by bcheck.awk.
#
#   960626 kb  created
#
#
########################################################################

BEGIN {
	# split lines at '=' characters
	FS = "=";
}


########################################################################

/Header/ { print }

/State/ {

	# handle all fields
	for(i=1; i<=NF; i++)
	{
		# split field at spaces
		k = split($i, kk, " ");

		# normal case is: k==2 (i.e., two subfields)
		if (k>2)
		{
			# exception: more spaces than 1 in current field
			for(j=1; j<k-1; j++)
			{
				# replace by underscores
				printf("%s_",kk[j]);
			}

			# keep last space
			printf("%s %s", kk[j], kk[j+1]);
		}
		else
			printf("%s",$i,k);

		# re-insert separator
		if (i<NF) printf("=");
	}
	# next line
	printf("\n");
}


########################################################################

END { }


