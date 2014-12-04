/****************************************************************************/
/*                                                                          */
/* File:      ifstd.ct                                                      */
/*                                                                          */
/* Purpose:   template routines for standard interface functions            */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   980825 kb  begin                                              */
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


#define PART ifHead


#ifdef IF_WITH_XARGS
	#define EXEC_LOOP(a,b,c)     IFExecHdrLoopCplX(a,b,c)
	#define COMM_LOOP(a,b,c,d,e) IFCommHdrLoopCplX(a,b,c,d,e)
#else
	#define EXEC_LOOP(a,b,c)     IFExecHdrLoopCpl(a,b,c)
	#define COMM_LOOP(a,b,c,d,e) IFCommHdrLoopCpl(a,b,c,d,e)
#endif

#define D_ABA(p)             (p->cpl)



#define IF_FUNCNAME   CAT(ddd_StdIF,IF_NAME)


/****************************************************************************/
/*                                                                          */
/* auxiliary macros                                                         */
/*                                                                          */
/****************************************************************************/

#ifndef IF_AUX_MACROS
#define IF_AUX_MACROS

#define CONCAT(txt,fct)    txt FUNCNAME_STR(fct)
#define FUNCNAME_STR(f)   #f


#endif


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



void IF_FUNCNAME (
	#ifdef IF_EXECLOCAL
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX ExecProcHdrXPtr ExecProc)
		#else
			NS_DIM_PREFIX ExecProcHdrPtr ExecProc)
		#endif
	#else
		size_t aSize,
		#ifdef IF_WITH_XARGS
			NS_DIM_PREFIX ComProcHdrXPtr Gather, NS_DIM_PREFIX ComProcHdrXPtr Scatter)
		#else
			NS_DIM_PREFIX ComProcHdrPtr Gather, NS_DIM_PREFIX ComProcHdrPtr Scatter)
		#endif
	#endif
{
	NS_DIM_PREFIX DDD_IF        aIF = NS_DIM_PREFIX STD_INTERFACE;
	NS_DIM_PREFIX IF_PROC		  *ifHead;
	unsigned long tries;
	int           recv_mesgs;


#ifdef IF_EXECLOCAL

	ForIF(aIF,ifHead)
	{
		EXEC_LOOP(ExecProc, D_ABA(PART), PART->nItems);
	}

#else /* ! IF_EXECLOCAL */

	/* allocate storage for in and out buffers */
	ForIF(aIF,ifHead)
	{

		#ifdef IF_EXCHANGE
			NS_DIM_PREFIX IFGetMem(ifHead, aSize, PART->nItems, PART->nItems);
		#endif
	}

	/* init communication, initiate receives */
	recv_mesgs = NS_DIM_PREFIX IFInitComm(aIF);


	/* build messages using gather-handler and send them away */
	ForIF(aIF,ifHead)
	{
		char   *buffer;
		buffer = BufferMem(ifHead->bufOut);

		buffer= COMM_LOOP(Gather, D_ABA(PART), buffer, aSize, PART->nItems);

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

					#ifdef CtrlTimeoutsDetailed
					printf("%4d: IFCTRL %02d received msg    from "
						"%4d after %10ld, size %ld\n",
						me, aIF, ifHead->proc,
						(unsigned long)tries,
						(unsigned long)BufferLen(ifHead->bufIn));
					#endif

					recv_mesgs--;
					ifHead->msgIn=NO_MSGID;

					buffer = BufferMem(ifHead->bufIn);

					/* get data using scatter-handler */
					buffer = COMM_LOOP(Scatter,
						D_ABA(PART), buffer, aSize, PART->nItems);
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

#ifdef IF_EXCHANGE
#undef IF_EXCHANGE
#endif

#ifdef IF_EXECLOCAL
#undef IF_EXECLOCAL
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


#undef _CAT
#undef CAT



/****************************************************************************/

