/****************************************************************************/
/*                                                                          */
/* File:      dlmgr.t                                                       */
/*                                                                          */
/* Purpose:   defines for dynamic linked list management                    */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: stefan@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7003                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   960915 sl  start of dynamic list management					*/
/*                                                                          */
/* Remarks:                                                                 */
/*            Management of dynamic linked lists, which consist of          */
/*            several parts. An object can be mapped to a part of the list  */
/*            by its priority using PRIO2LISTPART().                        */
/*            The formulation on object basis allows for management of      */
/*            elements, nodes, vectors and vertices.                        */
/*            A list has the form:                                          */
/*                p0first-p0last->p1first-p1last->...->pnfirst..pnlast      */
/*            where each part p0,p1,..,pn is limited by two pointers,       */
/*            pxfirst and pxlast, x=0..n. The part numbers are ordered      */
/*            increasingly and connected in a manner that one can run       */
/*            through the whole list in increasing order (SUCC) but only    */
/*            through one listpart in decreasing order (PRED).              */
/*            Linking/Unlinking of objects in a list part is done in a      */
/*            way that preserves these conventions.                         */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP

/*
void CAT(GRID_UNLINK_, OTYPE ) (GRID *Grid, OTYPE *Object)
*/
UNLINK(OTYPE)
	{
		INT Prio = DDD_InfoPriority( CAT(HDR(OTYPE),(Object)) );
		INT listpart = PRIO2LISTPART( CAT(OTYPE,_LIST) ,Prio);
		INT listpart1 = listpart;
		OTYPE *Object1 = NULL;

		IFDEBUG(gm,1) 
			printf("%d: GRID_UNLINK_" STR(OTYPE) "():" STR(OTYPE) 
				" has listpart=%d for prio=%dn",me,listpart,Prio);
			fflush(stdout);
		ENDDEBUG 

		if (listpart<0 || listpart>LASTPART_OF_LIST(OTYPE)) {
			printf("%d: GRID_UNLINK_" STR(OTYPE) "(): ERROR " STR(OTYPE)
				" has no valid listpart=%d for prio=%d\n",me,listpart,Prio);
			fflush(stdout);
			ASSERT(0);
		}

		switch (listpart) {

			case FIRSTPART_OF_LIST:

				if (PRED(Object)!=NULL)
					SUCC(PRED(Object)) = SUCC(Object);

				if (CAT(LISTPART_LAST,OTYPE(Grid,listpart)) != Object) {

					if (CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) == Object)
						CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) = SUCC(Object);

					if (SUCC(Object)!=NULL) 
						PRED(SUCC(Object)) = PRED(Object);
				}
				else {

					if (CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) == Object)
						CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) = NULL;

					CAT(LISTPART_LAST,OTYPE(Grid,listpart)) = PRED(Object);
				}

				break;
	
			case LASTPART_OF_LIST(OTYPE):

				if (PRED(Object)!=NULL) 
					SUCC(PRED(Object)) = SUCC(Object);
				else {
					CAT(LISTPART_FIRST,OTYPE(Grid,LASTPART_OF_LIST(OTYPE))) = SUCC(Object);

					do {
						listpart1--;
						Object1 = CAT(LISTPART_LAST,OTYPE(Grid,listpart1));
					}
					while (listpart1>FIRSTPART_OF_LIST && Object1==NULL);

					if (Object1!=NULL)
						SUCC(Object1) = SUCC(Object);
				}
				if (SUCC(Object)!=NULL) 
					PRED(SUCC(Object)) = PRED(Object);
				else {
					CAT(LISTPART_LAST,OTYPE(Grid,LASTPART_OF_LIST(OTYPE))) = PRED(Object);
					if (PRED(Object) != NULL) 
						SUCC(PRED(Object)) = NULL;
				}
				break;
	
			default:

				/* unlink in middle of list */
				if (PRED(Object)!=NULL) 
					SUCC(PRED(Object)) = SUCC(Object);
				else {

					if (SUCC(Object)!=NULL) 
						PRED(SUCC(Object)) = NULL;

					do {
						listpart1--;
						Object1 = CAT(LISTPART_LAST,OTYPE(Grid,listpart1));
					}
					while (listpart1>FIRSTPART_OF_LIST && Object1==NULL);

					if (Object1!=NULL)
						SUCC(Object1) = SUCC(Object);
				}
				if (CAT(LISTPART_LAST,OTYPE(Grid,listpart)) != Object) {

					if (CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) == Object)
						CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) = SUCC(Object);

					if (SUCC(Object)!=NULL) 
						PRED(SUCC(Object)) = PRED(Object);
				}
				else {

					if (CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) == Object)
						CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) = NULL;

					CAT(LISTPART_LAST,OTYPE(Grid,listpart)) = PRED(Object);
				}

				break;
		}
		SUCC(Object) = PRED(Object) = NULL;

	/* debug loop in list */
	{
		INT n = 0;
		for (Object1 = CAT(LISTPART_FIRST,OTYPE(Grid,FIRSTPART_OF_LIST));
			 Object1 != NULL;
			 Object1 = SUCC(Object1))
		{
			n++;	
			if (n>10000) {
				printf("%d: GRID_UNLINK_" STR(OTYPE) "():" STR(OTYPE) 
					" has loop listpart=%d for prio=%d\n",me,listpart,Prio);
				fflush(stdout);
				assert(0);
			}
		
		}

	}
}

#else

/*
void CAT(GRID_UNLINK_,OTYPE) (GRID *Grid, OTYPE *Object)
*/
UNLINK(OTYPE)
{                                                           
	if (PRED(Object)!=NULL)                                 
		SUCC(PRED(Object)) = SUCC(Object);                  
	else {                                                  
		CAT(FIRST,OTYPE(Grid)) = SUCC(Object);                
		if (SUCC(Object)!=NULL) PRED(SUCC(Object)) = NULL;  
	}                                                       
	if (SUCC(Object)!=NULL)                                 
		PRED(SUCC(Object)) = PRED(Object);                  
	else {                                                  
		CAT(LAST,OTYPE(Grid)) = PRED(Object);                 
	    if (PRED(Object)!=NULL) SUCC(PRED(Object)) = NULL;  
	}                                                       
}
#endif

#ifdef ModelP

/*
void CAT(GRID_LINK_,OTYPE) (GRID *Grid, OTYPE *Object, INT Prio)
*/
LINK(OTYPE)
	{
		INT listpart = PRIO2LISTPART(CAT(OTYPE,_LIST),Prio);
		INT listpartprev = listpart;
		INT listpartnext = listpart;
		OTYPE *Object1 = NULL;

		ASSERT(Grid != NULL);
		ASSERT(Object != NULL);
		ASSERT(Prio >= 0);

		IFDEBUG(gm,1) 
			printf("%d: GRID_LINK_" STR(OTYPE) "():" STR(OTYPE) 
				" has listpart=%d for prio=%d %xn",me,listpart,Prio,Object);
			fflush(stdout);
		ENDDEBUG 

		if (listpart<0 || listpart>LASTPART_OF_LIST(OTYPE)) {
			printf("%d: GRID_LINK_" STR(OTYPE) "(): ERROR " STR(OTYPE) 
				" has no valid listpart=%d for prio=%d\n",me,listpart,Prio);
			fflush(stdout);
			ASSERT(0);
		}

		PRED(Object) = SUCC(Object) = NULL;

		switch  (listpart) {

			case FIRSTPART_OF_LIST:

				Object1 = CAT(LISTPART_FIRST,OTYPE(Grid,FIRSTPART_OF_LIST));
				PRED(Object) = NULL;
				CAT(LISTPART_FIRST,OTYPE(Grid,FIRSTPART_OF_LIST)) = Object;
				if (Object1==NULL) {
					CAT(LISTPART_LAST,OTYPE(Grid,FIRSTPART_OF_LIST)) = Object;
					do {
						listpartnext++;
						CAT(Object1=LISTPART_FIRST,OTYPE(Grid,listpartnext));
					}
					while (listpartnext<LASTPART_OF_LIST(OTYPE) && Object1==NULL);
					SUCC(Object) = Object1;
				}
				else {
					SUCC(Object) = Object1;
					PRED(Object1) = Object;
				}
				break;

			case LASTPART_OF_LIST(OTYPE): 

				Object1 = CAT(LISTPART_LAST,OTYPE(Grid,LASTPART_OF_LIST(OTYPE)));
				SUCC(Object) = NULL;
				CAT(LISTPART_LAST,OTYPE(Grid,LASTPART_OF_LIST(OTYPE))) = Object;
				if (Object1==NULL)
				{
					PRED(Object) = NULL;
					CAT(LISTPART_FIRST,OTYPE(Grid,LASTPART_OF_LIST(OTYPE))) = Object;

					do {
						listpartprev--;
						Object1=CAT(LISTPART_LAST,OTYPE(Grid,listpartprev));
					}
					while (listpartprev>FIRSTPART_OF_LIST && Object1==NULL);

					if (Object1!=NULL)
						SUCC(Object1) = Object;
				}
				else
				{
					PRED(Object) = Object1;
					SUCC(Object1) = Object;
				}
				break;

			default: 

				/* link in middle of list */
				Object1 = CAT(LISTPART_FIRST,OTYPE(Grid,listpart));

				CAT(LISTPART_FIRST,OTYPE(Grid,listpart)) = Object;
				SUCC(Object) = Object1;
				PRED(Object) = NULL;

				/* empty list? */
				if (Object1 == NULL) {
					CAT(LISTPART_LAST,OTYPE(Grid,listpart)) = Object;
					do {
						listpartnext++;
						Object1=CAT(LISTPART_FIRST,OTYPE(Grid,listpartnext));
					}
					while (listpartnext<LASTPART_OF_LIST(OTYPE) && Object1==NULL);
					SUCC(Object) = Object1;
				}
				else 
					PRED(Object1) = Object;

				do {
					listpartprev--;
					Object1=CAT(LISTPART_LAST,OTYPE(Grid,listpartprev));
				}
				while (listpartprev>FIRSTPART_OF_LIST && Object1==NULL);

				if (Object1 != NULL)
					SUCC(Object1) = Object;
				break;
		}
	/* debug loop in list */
	{
		INT n = 0;
		for (Object1 = CAT(LISTPART_FIRST,OTYPE(Grid,FIRSTPART_OF_LIST));
			 Object1 != NULL;
			 Object1 = SUCC(Object1))
		{
			n++;	
			if (n>10000) {
				printf("%d: GRID_LINK_" STR(OTYPE) "():" STR(OTYPE) 
					" has loop listpart=%d for prio=%d\n",me,listpart,Prio);
				fflush(stdout);
				assert(0);
			}
		
		}

	}
}

#else

/*
void CAT(GRID_LINK_,OTYPE) (GRID *Grid, OTYPE *Object, INT Prio)
*/
LINK(OTYPE)
	{
		OTYPE *after;

		after = CAT(LAST,OTYPE(Grid));
		SUCC(Object) = NULL;
		if (after==NULL) {
			PRED(Object) = NULL;
			CAT(LAST,OTYPE(Grid)) = Object;
			CAT(FIRST,OTYPE(Grid)) = Object;
		}
		else {
			PRED(Object) = after;
			CAT(LAST,OTYPE(Grid)) = Object;
			SUCC(after) = Object;
		}
	}
#endif

#ifdef ModelP
/*
void CAT3(GRID_INIT_,OTYPE,_LIST(GRID *Grid)) 
*/
INIT(OTYPE)
	{
		INT i;
		for (i=0; i<=LASTPART_OF_LIST(OTYPE); i++){
			CAT(LISTPART_FIRST,OTYPE(Grid,i)) = NULL;
			CAT(LISTPART_LAST,OTYPE(Grid,i)) = NULL;
		}
	}
#else
/*
void CAT3(GRID_INIT_,OTYPE,_LIST(GRID *Grid))
*/
INIT(OTYPE)
    {
		CAT(FIRST,OTYPE(Grid)) = CAT(LAST,OTYPE(Grid)) = NULL;
	}
#endif
			
#ifdef ModelP
/* TODO: define this */
#define GRID_CHECK_OBJECT_LIST(Grid,OTYPE)
#endif

