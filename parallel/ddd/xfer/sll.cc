/****************************************************************************/
/*                                                                          */
/* File:      sll.ct                                                        */
/*                                                                          */
/* Purpose:   template routines for linked lists with freelist              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960826 kb  created                                            */
/*            970303 kb  added memory management in segms                   */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* segment of items */
typedef struct aSegm(T)
{
	struct aSegm(T) *next;
	int              nItems;

	T                item[SEGM_SIZE];
} Segm(T);


/* linked list of items */
T *list(T);


/* number of overall items */
int n(T);



/****************************************************************************/
/*                                                                          */
/* definition of local variables                                            */
/*                                                                          */
/****************************************************************************/

/* linked list of segms of items */
static Segm(T) *segms(T);


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static Segm(T) *NewSegm(T) (void)
{
	Segm(T) *segm;

	segm = (Segm(T) *) OO_Allocate (sizeof(Segm(T)));
	if (segm==NULL)
	{
		DDD_PrintError('F', 6060, STR_NOMEM " during XferEnd()");
		return(NULL);
	}

	segm->next   = segms(T);
	segms(T)     = segm;
	segm->nItems = 0;
	
	return(segm);
}


static void FreeSegms(T) (void)
{
	Segm(T) *segm = segms(T);
	Segm(T) *next = NULL;

	while (segm!=NULL)
	{
		next = segm->next;
		OO_Free (segm /*,sizeof(Segm(T))*/ );

		segm = next;
	}

	segms(T) = NULL;
}



T *New(T) (SLLNewArgProtos)
{
	Segm(T) *segm = segms(T);
	T        *item;

	if (segm==NULL || segm->nItems==SEGM_SIZE)
	{
		segm = NewSegm(T) ();
		if (segm==NULL)
			return(NULL);
	}

	item = &(segm->item[segm->nItems++]);


	/* insert item into linked list and count it */
	item->sll_next = list(T);
	list(T) = item;
	n(T)++;

	#ifdef SLL_WithOrigOrder
		/* insert unique counter */
		item->sll_n = n(T);
	#endif

#	ifdef SLL_DebugNew
		strncpy(item->sll_file, file, SLL_NAMELEN);
		item->sll_line = line;
#	endif

	return(item);
}




/*
	create pointer array from linked list and sort it
	according to given comparison function compar().
*/

T **SortedArray(T) (int (*compar) (const void *, const void *))
{
	T **array, *item;
	int  i;

	if (n(T)>0)
	{
		/* alloc array */
		array = (T **) OO_Allocate(sizeof(T *) * n(T));
		if (array==NULL)
		{
			DDD_PrintError('F', 6061, STR_NOMEM " during XferEnd()");
			return(NULL);
		}

		/* fill array with pointer */
		for(item=list(T), i=0; i<n(T); item=item->sll_next, i++)
		{
			array[i] = item;
		}

		/* sort by using compar function */
		if (n(T)>1)
			qsort(array, n(T), sizeof(T *), compar);
	}
	else
	{
		array = NULL;
	}

	return(array);
}



/****************************************************************************/

#ifdef SLL_WithOrigOrder

/*
	sort array of items into order of their New(T) command
	execution. the counter-component T.n is used for doing
	this.
*/

static int sort_OrigOrder(T) (const void *e1, const void *e2)
{
	T *item1 = *((T **)e1);
	T *item2 = *((T **)e2);

	if (item1->sll_n < item2->sll_n) return(-1);
	if (item1->sll_n > item2->sll_n) return(1);
	return(0);
}


void OrigOrder(T) (T **array, int n)
{
	qsort(array, n, sizeof(T *), sort_OrigOrder(T));
}

#endif



/****************************************************************************/


/*
	unify array of items. the array is compressed, the
	resulting number of valid items is returned.
	compar() is a comparison function which returns
	whether two items are equal or not.

	compar(a,b) should return 
		FALSE  if a should be skipped and eventually b could be chosen
		TRUE   if a must be taken.

	the array of items must be sorted in order to
	allow compar() to decide correctly.
*/

int Unify(T) (T **array, int (*compar) (T **, T **))
{
	int  i, cntValid;

	for(i=0, cntValid=0; i<n(T)-1; i++)
	{
		/* test if unique */
		if (compar(&array[i],&array[i+1]))
		{
			/* choose item */
			array[cntValid] = array[i];
			cntValid++;
		}
		/* else: skip item */
	}

	/* always choose last item */
	if (n(T)>0)
	{
		array[cntValid] = array[n(T)-1];
		cntValid++;
	}

	return(cntValid);
}



/****************************************************************************/


/*
	init linked list
*/
void Init(T) (void)
{
	list(T) = NULL;
	n(T) = 0;

	segms(T) = NULL;

	/*
	free(T) = NULL;
	*/
}



/*
	free all items
*/
void FreeAll(T) (void)
{
	list(T) = NULL;
	n(T) = 0;

	FreeSegms(T) ();
}


/*
	get quantitative resource usage
*/
void GetSizes(T) (int *nSegms, int *nItems, size_t *alloc_mem, size_t *used_mem)
{
	size_t   allocated=0, used=0;
	int      ns=0, ni=0;
	Segm(T)  *segm;

	for (segm=segms(T); segm!=NULL; segm=segm->next)
	{
		/* count number of segments and number of items */
		ns++;
		ni+=segm->nItems;

		/* compute memory usage */
		allocated += sizeof(Segm(T));
		used += (sizeof(Segm(T)) - (sizeof(T)*(SEGM_SIZE-segm->nItems)));
	}

	*nSegms    = ns;
	*nItems    = ni;
	*alloc_mem = allocated;
	*used_mem  = used;
}


/****************************************************************************/

#ifdef SLL_WithOrigOrder
#undef SLL_WithOrigOrder
#endif

/****************************************************************************/

