########################################################################
#
#   bcheck.awk
#
#   check versions of binary files for consistency of
#   compile time switches. 'ident' should be used for extracting
#   version information from executables, objects or libraries,
#   the output is piped into this awk-script in order to check for
#   valid configurations.
#
#   the executable must contain the keywords $Header$ and $State$
#   for each object file. the $State$ line must give information with
#   syntax PROPERTY=VALUE, where PROPERTY is the name of a compile time
#   switch (or #define) and VALUE is its value at compile time.
#   sample excerpt from 'ident libddd.a' output:
#
#   $Header$
#   $State$
#   $Header$
#   $State$
#   [...]
#
#   this awk-script builds up a database for this P-V-pairs, and 
#   derives consistency information from it.
#
#   960625 kb  created
#
#
########################################################################

BEGIN {
	last_filename = DUMMY;
	expecting = "header";
	no_headers = 0;
	no_states = 0;
	no_aliens = 0;
}


########################################################################

/Header/ && ($0 !~ header_condition) {
	# count number of Header-lines which don't contain header_condition
	# as a reg-expr.
	no_aliens++;
}


/Header/ && ($0 ~ header_condition) {

	# reset first and last entry
	$1 = ""
	$NF = ""


	if (expecting!="header")
	{
		missing[headers[last_filename]] = no_headers;
	}
	expecting = "state";

	no_headers++;

	# remember current filename
	last_filename = $2;

	# remember RCS header for current filename
	headers[last_filename] = $0;
}



/State/ {
	expecting = "header";
	no_states++;

	# loop through items of type "PROPERTY=VALUE"
	for(i=2; i<NF; i++)
	{
		# ii[1] is property, ii[2] is value
		split($i,ii,"=");
		p = ii[1]; v = ii[2];

		# ignore item, if value equals property
		# (i.e., no expansion has occured)
		if (p!=v)
		{
			# count number of entries per property
			properties[p] += 1;

			# pv_count[property,value] counts the number of
			# occurences of every property/value pair
			pv_count[p,v] += 1;

			# build up database, for given (filename,property) the 
			# corresponding value is stored.
			values[last_filename,p] = v;
		}
	}
}


########################################################################

END {
	if (no_headers==0)
	{
		printf("%s: %s doesn't contain id keywords. aborted.\n",
			toolname, filename);
	}
	else
	{
		if (no_states==0)
		{
			printf("%s: %s doesn't contain 'State' keywords.\nsorry, no information available yet.\n",
				toolname, filename);
		}
		else
		{
			# some postprocessing
			compute_max();
			check_properties(prop_ok);

			# automatic report generation
			if (verbose==1)
				print_verbose_report();
			else
				print_error_report();
		}
	}
}



########################################################################
#
#    subroutines
#
########################################################################


# this is the report generator for 'verbose' case

function print_verbose_report ()
{
	inconsistent_values = 0;

	print "--------------------------------------------------";
	for (prop in properties)
	{
		# evaluate consistency
		if (prop_ok[prop])
		{
			result = "consistent.";
		}
		else
		{
			result = "detected inconsistencies!";
			inconsistent_values++;
		}
		print_prop_header(prop, result);

		if (prop_ok[prop])
		{
			printf("VALUE    : %39s\n", prop_vmax_val[prop]);
		}
		else
		{
			printf("DETAILS  :\n");
			print_details(prop);
		}

		print "--------------------------------------------------";
	}

	if (inconsistent_values==0)
		expertise = sprintf("%s is ok.", filename);
	else
		expertise = sprintf("%s is inconsistent.", filename);

	print "--------------------------------------------------";
	printf("SUMMARY  : %39s\n", expertise);

	if (inconsistent_values>0)
		printf("inconsistent property values: %20d\n",
			inconsistent_values);

	if (no_states<no_headers)
		print_missing_states();

	if (no_aliens>0)
		printf("lines containing Header-keyword, but without \"%s\": %d\n",
			header_condition, no_aliens);


	print "--------------------------------------------------";
}


########################################################################


# this is the report generator for 'not verbose' case

function print_error_report ()
{
	inconsistent_values = 0;
	for (prop in properties)
	{
		# evaluate consistency
		if (! prop_ok[prop])
		{
			inconsistent_values++;
			printf("%s: %s has inconsistent property %s (",
				toolname, filename, prop);
			print_wrong_pvs(prop);
			printf(").\n");

			print_details(prop);
		}
	}

	if (no_states<no_headers)
		print_missing_states();
	else
	{
		if (inconsistent_values==0)
			printf("%s: %s is consistent.\n",
				toolname, filename);
	}
}


########################################################################

function print_missing_states ()
{
	if (no_states<no_headers)
	{
		printf("%s: %s has %d Headers, but only %s States! Missing:\n",
			toolname, filename, no_headers, no_states);

		for(m in missing)
		{
			print "      ", m;
		}
	}
}


########################################################################

function print_prop_header (prop, result)
{
	printf("PROPERTY : %39s\n", prop);
	printf("ENTRIES  : %39d\n", properties[prop]);
	printf("RESULT   : %39s\n", result);
}


########################################################################

function print_overview ()
{
	for (prop in properties)
	{
		print "PROPERTY", prop, ",", properties[prop], "entries."

		for (idx in pv_count)
		{
			split(idx,ii,SUBSEP);
			if (ii[1]==prop)
			{
				print "    ",pv_count[idx],"x", ii[2]
			}
		}
	}
}


########################################################################

function print_details (prop)
{
	for (idx in pv_count)
	{
		# ii[1] is property, ii[2] is value
		split(idx,ii,SUBSEP);

		# look for given property
		if (ii[1]==prop)
		{
			if (ii[2] != prop_vmax_val[prop])
			{
				printf("%4d x '%s', in files:\n", pv_count[idx], ii[2]);
				print_pv_headers(prop, ii[2]);
			}
			else
			{
				printf("%4d x '%s' (correct value?)\n", pv_count[idx], ii[2]);
			}
		}
	}
}


########################################################################


function compute_max ()
{
	for (idx in pv_count)
	{
		# ii[1] is property, ii[2] is value
		split(idx,ii,SUBSEP);

		if (prop_vmax[ii[1]] < pv_count[idx])
		{
			prop_vmax[ii[1]] = pv_count[idx];
			prop_vmax_val[ii[1]] = ii[2];
		}
	}
}


########################################################################

function print_max ()
{
	for (prop in prop_vmax)
	{
		print "prop",prop,"has max", prop_vmax[prop];
	}
}



########################################################################

function print_pv_headers (property, value)
{
	for (idx in values)
	{
		# ii[1] is filename, ii[2] is property
		split(idx,ii,SUBSEP);

		if (ii[2]==property && values[idx]==value)
		{
			print "      ", headers[ii[1]];
		}
	}
}



########################################################################

function check_properties (results)
{
	for (prop in properties)
	{
		# a property is ok, iff all values are equal
		results[prop] = (prop_vmax[prop] == properties[prop]);
	}
}


########################################################################

function print_wrong_pvs (prop)
{
	for (idx in pv_count)
	{
		# ii[1] is property, ii[2] is value
		split(idx,ii,SUBSEP);

		if (ii[1]==prop && ii[2]!=prop_vmax_val[prop])
			wrong_pvs+=pv_count[idx];
	}

	printf("%d of %d", wrong_pvs, properties[prop]);
}



########################################################################

