/****************************************************************************/
/*                                                                          */
/* File:      ifcmd.ct                                                      */
/*                                                                          */
/* Purpose:   template routines for interface functions                     */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   970110 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* derived macros                                                           */
/*                                                                          */
/****************************************************************************/

#define _CAT(a,b)   a ## b
#define CAT(a,b)    _CAT(a,b)


#ifdef IF_WITH_ATTR
#define PART ifAttr
#else
#define PART ifHead
#endif


#ifdef IF_WITH_XARGS
	#define EXEC_LOOP(a,b,c)     NS_DIM_PREFIX IFExecLoopCplX(a,b,c)
	#define COMM_LOOP(a,b,c,d,e) NS_DIM_PREFIX IFCommLoopCplX(a,b,c,d,e)
	#define D_AB(p)              (p->cplAB)
	#define D_BA(p)              (p->cplBA)
	#define D_ABA(p)             (p->cplABA)
#else
	#define EXEC_LOOP(a,b,c)     NS_DIM_PREFIX IFExecLoopObj(a,b,c)
	#ifdef CPP_FRONTEND
	  #ifdef IF_CBR
	    #define COMM_LOOP(a,b,c,d,e) \
			CAT(NS_DIM_PREFIX IFCommLoopObj,a) (GatherScatter,b,c,d,e)
	  #else
	    #define COMM_LOOP(a,b,c,d,e) \
			CAT(NS_DIM_PREFIX IFCommLoopObj,a) (*GatherScatter,b,c,d,e)
	  #endif
	#else
	#define COMM_LOOP(a,b,c,d,e) NS_DIM_PREFIX IFCommLoopObj(a,b,c,d,e)
	#endif
	#define D_AB(p)              (p->objAB)
	#define D_BA(p)              (p->objBA)
	#define D_ABA(p)             (p->objABA)
#endif


#ifdef CPP_FRONTEND
#define IF_FUNCNAME   CAT(DDD_Interface::,IF_NAME)
#else
#define IF_FUNCNAME   CAT(DDD_IF,IF_NAME)
#endif

#ifdef CPP_FRONTEND
	#ifdef IF_CBR
		#define REFDECL &
	#else
		#define REFDECL *
	#endif
#endif


/****************************************************************************/
/*                                                                          */
/* auxiliary macros                                                         */
/*                                                                          */
/****************************************************************************/

#ifndef IF_AUX_MACROS
#define IF_AUX_MACROS

#define FIND_ATTR \
	ifAttr = ifHead->ifAttr;                          \
	while ((ifAttr!=NULL) && (ifAttr->attr!=aAttr))   \
	{                                                 \
		ifAttr = ifAttr->next;                        \
	}                                                 \
	if (ifAttr==NULL) continue;


#define CONCAT(txt,fct)    txt FUNCNAME_STR(fct)
#define FUNCNAME_STR(f)   #f


#endif


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/


#ifdef IF_EXECLOCAL
	/**
	Local execution on \ddd{Interface}.
	This function executes a given handler for each object inside
	a \ddd{Interface}.

	}*/
#else
	#if defined(IF_ONEWAY)
		/**
		Oneway-communication across a \ddd{Interface}.
		This function initiates a communication across a previously
		defined DDD interface. It has to be issued on all processors.

		The communication is carried out one-directional, dependent
		on parameter {\em aDir}. In case of #IF_FORWARD# information
		will flow from object set A to object set B, in case of #IF_BACKWARD#
		it will be vice versa (see \funk{IFDefine}).
		}*/
	#else
		/**
		Bidirectional communication across a \ddd{Interface}.
		This function initiates a communication across a previously
		defined DDD interface. It has to be issued on all processors.

		The communication is carried out in both directions, i.e.~from
		object set A to object set B and vice versa (see \funk{IFDefine}).
		}*/
	#endif

	/*{
	The user-defined function {\em gather\_proc} is called for each object
	inside the interface. It has to pack all necessary information into
	a special memory location. Afterwards the actual communication based
	on the current processor topology happens. After all data had been
	scheduled onto the destination processors, it has to be
	incorporated into the destination part of the DDD interface via
	the user-defined function {\em scatter\_proc}.

	Both communication procedures should have the
	following declaration:

	#int (*ComProcPtr) (DDD_OBJ obj, void *data)#
	}*/
#endif

	/*{
	\todo{further description}
	}*/
	#ifdef IF_WITH_XARGS
		/*{
		\todo{Doku extended args}
		}*/
	#endif

/*{
@param aIF   the \ddd{Interface} to be used.
}*/
#ifdef IF_WITH_ATTR
	/*{
	@param aAttr  the \ddd{Attribute} which specifies the sub-interface.
	}*/
#endif
#ifdef IF_EXECLOCAL
	/*{
	@param ExecProc  handler which is executed for each object in the interface.
	}*/
#else
	#ifdef IF_ONEWAY
		/*{
		@param aDir  the communication direction (#IF_FORWARD# or #IF_BACKWARD#)
		}*/
	#endif
	/*{
	@param aSize  amount of data for each object in the interface (number of bytes)
	}*/
#endif
/*{
*/


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
void IF_FUNCNAME (
#ifdef C_FRONTEND
	NS_DIM_PREFIX DDD_IF aIF,
#endif
	#ifdef IF_WITH_ATTR
		NS_DIM_PREFIX DDD_ATTR aAttr,
	#endif
	#ifdef IF_EXECLOCAL
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX ExecProcXPtr ExecProc)
		#else
			NS_DIM_PREFIX ExecProcPtr ExecProc)
		#endif
	#else
		#ifdef IF_ONEWAY
			NS_DIM_PREFIX DDD_IF_DIR aDir,
		#endif
		size_t aSize,
#ifdef C_FRONTEND
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX ComProcXPtr Gather, NS_DIM_PREFIX ComProcXPtr Scatter)
		#else
			NS_DIM_PREFIX ComProcPtr Gather, NS_DIM_PREFIX ComProcPtr Scatter)
		#endif
#endif
#ifdef CPP_FRONTEND
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX DDD_GatherScatterX REFDECL GatherScatter)
		#else
			NS_DIM_PREFIX DDD_GatherScatter  REFDECL GatherScatter)
		#endif
#endif
	#endif
{
	#ifdef CPP_FRONTEND
		DDD_IF aIF     = _id;
		//ComProcPtr Scatter = GatherScatter->Scatter;
		//ComProcPtr Gather  = GatherScatter->Gather;
	#endif
#endif

#ifdef F_FRONTEND
void IF_FUNCNAME (NS_DIM_PREFIX DDD_IF *_aIF,
	#ifdef IF_WITH_ATTR
		NS_DIM_PREFIX DDD_ATTR *_aAttr,
	#endif
	#ifdef IF_EXECLOCAL
		#ifdef IF_WITH_XARGS
			ExecProcXPtr ExecProc)
		#else
			ExecProcPtr ExecProc)
		#endif
	#else
		#ifdef IF_ONEWAY
			int *_aDir,
		#endif
		size_t *_aSize,
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX ComProcXPtr Gather, NS_DIM_PREFIX ComProcXPtr Scatter)
		#else
			NS_DIM_PREFIX ComProcPtr Gather, NS_DIM_PREFIX ComProcPtr Scatter)
		#endif
	#endif
{
	/* in F_FRONTEND, copy call-by-ref params into variables */
	DDD_IF aIF     = *_aIF;
	#ifdef IF_WITH_ATTR
	DDD_ATTR aAttr   = *_aAttr;
	#endif
	#ifndef IF_EXECLOCAL
		#ifdef IF_ONEWAY
		int    aDir      = *_aDir;
		#endif
		size_t aSize = *_aSize;
	#endif
#endif  /* F_FRONTEND */

	NS_DIM_PREFIX IF_PROC		  *ifHead;
	unsigned long tries;
	int           recv_mesgs;


	/* prohibit using standard interface (IF0) */
	if (aIF==NS_DIM_PREFIX STD_INTERFACE)
	{
		DDD_PrintError('E', 4300,
			CONCAT("cannot use standard interface in ",IF_FUNCNAME));
		HARD_EXIT;
	}


	/* shortcuts can only be used without extended handler arguments */
	#ifndef IF_WITH_XARGS
		/* if shortcut tables are invalid -> recompute */
		NS_DIM_PREFIX IFCheckShortcuts(aIF);
	#endif


#ifdef IF_EXECLOCAL

	ForIF(aIF,ifHead)
	{
		#ifdef IF_WITH_ATTR
			NS_DIM_PREFIX IF_ATTR	    *ifAttr;
			FIND_ATTR; /* find ifAttr */
		#endif

		EXEC_LOOP(ExecProc, D_BA(PART),  PART->nBA);
		EXEC_LOOP(ExecProc, D_AB(PART),  PART->nAB);
		EXEC_LOOP(ExecProc, D_ABA(PART), PART->nABA);
	}

#else /* ! IF_EXECLOCAL */

	/*
	STAT_ZEROTIMER;
	STAT_RESET1;
	STAT_SETCOUNT(20,0);
	STAT_SETCOUNT(21,0);
	STAT_SETCOUNT(22,0);
	*/

	/* allocate storage for in and out buffers */
	ForIF(aIF,ifHead)
	{
		#ifdef IF_WITH_ATTR
			NS_DIM_PREFIX IF_ATTR	    *ifAttr;
		#endif
		#ifdef IF_ONEWAY
			int      nOut, nIn;
		#endif

		#ifdef IF_WITH_ATTR
			/* reset buffers */
			BufferReset(ifHead->bufIn);
			BufferReset(ifHead->bufOut);
			FIND_ATTR; /* find ifAttr */
			/* if no ATTR can be found, buffers will be empty */
		#endif

		#ifdef IF_EXCHANGE
			NS_DIM_PREFIX IFGetMem(ifHead, aSize, PART->nItems, PART->nItems);
		#endif

		#ifdef IF_ONEWAY
			if (aDir==NS_DIM_PREFIX IF_FORWARD) {
				nOut = PART->nAB; nIn = PART->nBA;
			}
			else {
				nOut = PART->nBA; nIn = PART->nAB;
			}
			IFGetMem(ifHead, aSize, nIn+PART->nABA, nOut+PART->nABA);
		#endif

		/* STAT_SETCOUNT(20,STAT_GETCOUNT(20)+ifHead->nItems);
		STAT_SETCOUNT(21,STAT_GETCOUNT(21)+2);
		STAT_SETCOUNT(22,STAT_GETCOUNT(22)+
			BufferLen(ifHead->bufIn) + BufferLen(ifHead->bufOut));
		*/
	}
	/*
	STAT_TIMER1(63);
	*/


	/* init communication, initiate receives */
	recv_mesgs = NS_DIM_PREFIX IFInitComm(aIF);


	/* build messages using gather-handler and send them away */
	ForIF(aIF,ifHead)
	{
		char   *buffer;
		#ifdef IF_WITH_ATTR
			NS_DIM_PREFIX IF_ATTR *ifAttr;
		#endif
		#ifdef IF_ONEWAY
			int      nOut;
			#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX COUPLING **datOut;
			#else
			NS_DIM_PREFIX IFObjPtr  *datOut;
			#endif
		#endif

		#ifdef IF_WITH_ATTR
		FIND_ATTR; /* find ifAttr */
		#endif

		buffer = BufferMem(ifHead->bufOut);

		#ifdef IF_EXCHANGE
			/* exchange BA and AB during send */
			buffer= COMM_LOOP(Gather, D_BA(PART), buffer, aSize, PART->nBA);
			buffer= COMM_LOOP(Gather, D_AB(PART), buffer, aSize, PART->nAB);
		#endif

		#ifdef IF_ONEWAY
			if (aDir==NS_DIM_PREFIX IF_FORWARD) {
				nOut = PART->nAB;  datOut = D_AB(PART);
			}
			else {
				nOut = PART->nBA;  datOut = D_BA(PART);
			}

			buffer= COMM_LOOP(Gather, datOut, buffer, aSize, nOut);
		#endif

		buffer= COMM_LOOP(Gather, D_ABA(PART), buffer, aSize, PART->nABA);

		NS_DIM_PREFIX IFInitSend(ifHead);
	}



	/* poll receives and scatter data */
	for(tries=0; tries<MAX_TRIES && recv_mesgs>0; tries++)
	{
		/* poll receive calls */
		ForIF(aIF,ifHead)
		{
			if ((! BufferIsEmpty(ifHead->bufIn)) && ifHead->msgIn!=NO_MSGID)
			{
				int error = InfoARecv(ifHead->vc, ifHead->msgIn);
				if (error==-1)
				{
					sprintf(NS_DIM_PREFIX cBuffer,
							"PPIF's InfoARecv() failed for recv to proc=%d in IF-Comm",
						ifHead->proc);
					DDD_PrintError('E', 4221, NS_DIM_PREFIX cBuffer);
					HARD_EXIT;
				}

				if (error==1)
				{
					char     *buffer;
					#ifdef IF_ONEWAY
						int      nIn;
						#ifdef IF_WITH_XARGS
						NS_DIM_PREFIX COUPLING **datIn;
						#else
						NS_DIM_PREFIX IFObjPtr  *datIn;
						#endif
					#endif
					#ifdef IF_WITH_ATTR
						NS_DIM_PREFIX IF_ATTR  *ifAttr;
					#endif

					#ifdef CtrlTimeoutsDetailed
					printf("%4d: IFCTRL %02d received msg    from "
						"%4d after %10ld, size %ld\n",
						me, aIF, ifHead->proc,
						(unsigned long)tries,
						(unsigned long)BufferLen(ifHead->bufIn));
					#endif

					recv_mesgs--;
					ifHead->msgIn=NO_MSGID;

					#ifdef IF_WITH_ATTR
					FIND_ATTR; /* find ifAttr */
					#endif

					buffer = BufferMem(ifHead->bufIn);

					/* get data using scatter-handler */
					#ifdef IF_EXCHANGE
						buffer = COMM_LOOP(Scatter,
									D_AB(PART), buffer, aSize, PART->nAB);
						buffer = COMM_LOOP(Scatter,
									D_BA(PART), buffer, aSize, PART->nBA);
					#endif

					#ifdef IF_ONEWAY
						if (aDir==NS_DIM_PREFIX IF_FORWARD) {
							nIn  = PART->nBA;  datIn = D_BA(PART);
						}
						else {
							nIn  = PART->nAB;  datIn = D_AB(PART);
						}

						buffer = COMM_LOOP(Scatter, datIn, buffer, aSize, nIn);
					#endif

					buffer = COMM_LOOP(Scatter,
						D_ABA(PART), buffer, aSize, PART->nABA);
				}
			}
		}
	}

	/* finally poll send completion */
	if (recv_mesgs>0)
	{
		sprintf(NS_DIM_PREFIX cBuffer, CONCAT("receive-timeout for IF %02d in ",IF_FUNCNAME), aIF);
		DDD_PrintError('E', 4200, NS_DIM_PREFIX cBuffer);

		ForIF(aIF,ifHead)
		{
			if ((! BufferIsEmpty(ifHead->bufIn)) && ifHead->msgIn!=NO_MSGID)
			{
				sprintf(NS_DIM_PREFIX cBuffer,
					"  waiting for message (from proc %d, size %ld)",
					ifHead->proc, (unsigned long)BufferLen(ifHead->bufIn));
				DDD_PrintError('E', 4201, NS_DIM_PREFIX cBuffer);
			}
		}
	}
	else
	{
		#ifdef CtrlTimeouts
		printf("%4d: IFCTRL %02d received msg  all after %10ld tries\n",
			me, aIF, (unsigned long)tries);
		#endif

		if (! NS_DIM_PREFIX IFPollSend(aIF))
		{
			sprintf(NS_DIM_PREFIX cBuffer, CONCAT("send-timeout for IF %02d in ",IF_FUNCNAME), aIF);
			DDD_PrintError('E', 4210, NS_DIM_PREFIX cBuffer);

			ForIF(aIF,ifHead)
			{
				if ((! BufferIsEmpty(ifHead->bufOut)) && ifHead->msgOut!=NO_MSGID)
				{
					sprintf(NS_DIM_PREFIX cBuffer,
						"  waiting for send completion (to proc %d, size %ld)",
						ifHead->proc, (unsigned long)BufferLen(ifHead->bufOut));
					DDD_PrintError('E', 4211, NS_DIM_PREFIX cBuffer);
				}
			}
		}
	}


	/* free memory */
	NS_DIM_PREFIX IFExitComm(aIF);

	/*STAT_TIMER1(60);*/

#endif /* IF_EXECLOCAL */
}


/****************************************************************************/

/* undefs for all input macros and defines */

#undef IF_NAME
#undef IF_FUNCNAME

#ifdef IF_ONEWAY
#undef IF_ONEWAY
#endif

#ifdef IF_EXCHANGE
#undef IF_EXCHANGE
#endif

#ifdef IF_EXECLOCAL
#undef IF_EXECLOCAL
#endif


#ifdef IF_WITH_ATTR
#undef IF_WITH_ATTR
#endif

#ifdef IF_WITH_XARGS
#undef IF_WITH_XARGS
#endif

#ifdef IF_CBR
#undef IF_CBR
#endif



/* undefs for derived macros */

#undef PART
#undef EXEC_LOOP
#undef COMM_LOOP
#undef D_AB
#undef D_BA
#undef D_ABA


#ifdef CPP_FRONTEND
#undef REFDECL
#endif


#undef _CAT
#undef CAT



/****************************************************************************/

