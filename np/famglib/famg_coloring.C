/****************************************************************************/
/*																			*/
/* File:      famg_coloring.C												*/
/*																			*/
/* Purpose:   parallel graph coloring functions for FAMG					*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Wissenschaftliches Rechnen						*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  internet: Christian.Wrobel@iwr.uni-heidelberg.de				*/
/*																			*/
/*																			*/
/* History:   February 99 begin, Stuttgart									*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#include <stdlib.h> 		// for (s)rand
#include <string.h> 		// for memset
#include <limits.h> 		// for INT_MAX
#include <iostream.h>

#include "gm.h"
#include "ddd.h"
#include "famg_coloring.h"

/* RCS_ID
$Header$
*/

FAMGColor FAMGMyColor;		// the color of this PE

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static DDD_PROC Nb[FAMGColorMaxNb];	// list of all neighbor PE's
static int NrNb = 0;				// number of valid entries in Nb
static int *helpNbPtr= NULL;

static int DetermineNbs( DDD_OBJ obj)
{
	VECTOR *vec = (VECTOR *)obj;
	int *proclist, i;
	
	proclist = DDD_InfoProcList(PARHDR(vec));
	for( i = 2; proclist[i] != -1	; i += 2 )
		helpNbPtr[proclist[i]] = 1;		// perhaps to much neighbors because of counting ghost-copies
	return 0;
}


int ConstructColoringGraph( DDD_ATTR grid_attr)
{
	int i;
	int helpNb[FAMGColorMaxProcs];
	
	assert(FAMGColorMaxProcs>=procs);	// otherwise increase the constant FAMGColorMaxProcs 
	
	// helpNb[] = 0
	memset( helpNb, 0, sizeof(int)*FAMGColorMaxProcs );
	helpNbPtr = helpNb;
	
	// determine the neighboring PE's
	DDD_IFAExecLocal( BorderVectorSymmIF, grid_attr, DetermineNbs );
	
	NrNb = 0;
	for( i = 0; i < FAMGColorMaxProcs; i++ )
		if( helpNb[i] != 0 )
			Nb[NrNb++] = i;
	
	assert(FAMGColorMaxNb>=NrNb);	// otherwise increase the constant FAMGColorMaxNb 

	return 0;
}

static int ColorCompare( const void *c1, const void *c2 )
{
    if ((*(FAMGColor *)c1) < (*(FAMGColor *)c2)) return(-1);
    if ((*(FAMGColor *)c1) > (*(FAMGColor *)c2)) return(1);

    return(0);
}

int ConstructColoring( int OrderingFunctionType )
// according to 
//		Robert K. Gjertsen jr., Mark T. Jones and Paul E. Plassmann
//		Parallel Heuristics for Improved, Balanced Graph Colorings
//		to appear in journal of Parallel and Distributed Computing
{
	int i, nbpe, res, j;	
	VChannelPtr NbCh[FAMGColorMaxNb];		// communication channels to the neighbor Pe's
	VChannelPtr ch;
	double MyWeight;						// weight of this Pe (see PLF algorithm); must be >= 0
	int SendQueue[FAMGColorMaxNb];			// list of the index in Nb of the neighbors that get a color-message from me 
	int Recv[FAMGColorMaxNb];				// 1 for neighbors from which to receive a weight-message, 
											// 2 for neighbors from which to receive a color-message
	int NrSend = 0;							// number of entries in SendQueue (i.e. neighbors with a smaller weight than me)
	int NrWait = 0;							// number of neighbors with a larger weight than me (from which I receive their color-messages)
	FAMGColor NbColor[FAMGColorMaxNb];		// list of colors of the neighbors with a larger weight than me
	double NbWeight[FAMGColorMaxNb];		// list of weights of the neighbors
	msgid MsgOutId[FAMGColorMaxNb];			// id of async send's
	msgid MsgInId[FAMGColorMaxNb];			// id of async recv's
	
	switch( OrderingFunctionType )
	{
		case 1:
			FAMGMyColor = me;
			return 0;
		case 2:
			srand(me);
			MyWeight = rand() / (double)INT_MAX;
			break;
		case 3:
			srand(me);
			MyWeight = NrNb + rand() / (double)INT_MAX;
			break;
		default:
			cout << "ConstructColoring(): unknown ordering function type" << endl << fflush;
			abort();
	}
	
	if( NrNb > FAMGColorMaxNb )
	{
		cout << "ConstructColoring(): error Number of neighbors ("<<NrNb<<") larger than maximum <<FAMGColorMaxNb<<. Increase FAMGColorMaxNb"<<endl<<fflush;
		abort();
	}
	
	//
	// construct the communication channels
	//
	for( i = 0; i < NrNb; i++ )
		NbCh[i] = ConnASync( Nb[i], 7643 );		// just a silly number
	
	//
	// communicate weights with all neighbors
	//
	for( i = 0; i < NrNb; i++ )
	{
		ch = NbCh[i];
		
		// check success of constructing the channel
		while( (res=InfoAConn(ch)) == 0)
			;	// wait until completion
		if( res != 1 )
		{
			cout << "ConstructColoring(): error "<<res<<" during channel construction"<<endl<<fflush;
			abort();
		}
		
		// now ready to send
		MsgOutId[i] = SendASync( ch, &MyWeight, sizeof MyWeight, &res );
		if( res != 0 )
		{
			cout << "ConstructColoring(): error "<<res<<" in SendASync myweight for PE "<<Nb[i]<<endl<<fflush;
			abort();
		}
		assert(MsgOutId[i]!=-1);
		
		// receive 
		MsgInId[i] = RecvASync( ch, NbWeight+i, sizeof(double), &res );
		if( res != 0 )
		{
			cout << "ConstructColoring(): error "<<res<<" in RecvASync weights for PE "<<Nb[i]<<endl<<fflush;
			abort();
		}	
		
		Recv[i] = 1;	// receive a weight message from this neighbor
	}
	
	//
	// received weights from all neighbors? put receive calls for color.
	//
	j = NrNb;
	while( j>0 )
		for( i = 0; i < NrNb; i++ )
		{
			if( Recv[i] == 1 ) // pending receive weight
			{
				ch = NbCh[i];
				res = InfoARecv( ch, MsgInId[i] );
				
				if( res == 1 )
				{	// message arrived
					j--;
					Recv[i] = 0;	// reset
					
					if( MyWeight < NbWeight[i] )
					{
						// receive color
						MsgInId[i] = RecvASync( ch, NbColor+i, sizeof(FAMGColor), &res );
						if( res != 0 )
						{
							cout << "ConstructColoring(): error "<<res<<" in RecvASync color for PE "<<Nb[i]<<endl<<fflush;
							abort();
						}	
						assert(MsgInId[i]!=-1);
						NrWait++;
						Recv[i] = 2;	// receive a color message from this neighbor
					}
					else
					#ifdef Debug
					if( MyWeight==NbWeight[i] )
					{
						cout << "ConstructColoring(): error: received weight "<<NbWeight[i]<<" from PE "<<Nb[i]<<". Same as my weight!"<<endl<<fflush;
						abort();
					}
					else
					#endif						
						SendQueue[NrSend++] = i;
					
					if( j==0 )
						break;					
				}
				else if( res != 0 )
				{
					cout << "ConstructColoring(): error: return status "<<res<<" in InfoARecv for weight from PE "<<Nb[i]<<". Same as my weight!"<<endl<<fflush;
					abort();
				}
			}

		}
	
	//
	// received color from all neighbors?
	//
	j = NrWait;
	while( j>0 )
		for( i = 0; i < NrNb; i++ )
		{
			if( Recv[i] == 2 ) // pending receive color
			{
				res = InfoARecv( NbCh[i], MsgInId[i] );
				
				if( res == 1 )
				{	// message arrived
					j--;
					Recv[i] = 3;	// reset this flag and set flag for next step
					if( j==0 )
						break;
				}
				else if( res != 0 )
				{
					cout << "ConstructColoring(): error: return status "<<res<<" in InfoARecv for color from PE "<<Nb[i]<<endl<<fflush;
					abort();
				}
			}

		}
	
	//
	// rearrange color array
	//
	
	// for every i with Recv[i]==2 NbColor has a valid color entry
	// compress the color array such that the colors are stored 
	// in the first consecutive array positions
	i = j = 0;
	for( i = 0; i < NrNb; i++ )
		if( Recv[i]==3 )
			NbColor[j++] = NbColor[i];
	assert(j == NrWait);
	
	// sort the colors
	qsort( (void*)NbColor, NrWait, sizeof(FAMGColor), ColorCompare );
	
	//
	// determine the smallest unused color for me
	//
	for( i = 0; i < NrWait; i++ )
		if( NbColor[i] != i )
		{
			FAMGMyColor = i;	// use the first gap in the color sequence
			break;
		}
	if( i == NrWait )
		FAMGMyColor = NrWait;	// introduce a new color
	
	//
	// check if first sends are finished
	//
	j = NrNb;
	while( j>0 )
		for( i = 0; i < NrNb; i++ )
		{
			if( MsgOutId[i] != -1 ) // pending send weight
			{
				res = InfoASend( NbCh[i], MsgOutId[i] );
				
				if( res == 1 )
				{	// message sended successfully
					j--;
					MsgOutId[i] = -1;	// reset
					if( j==0 )
						break;
				}
				else if( res != 0 )
				{
					cout << "ConstructColoring(): error: return status "<<res<<" in InfoASend for weight from PE "<<Nb[i]<<endl<<fflush;
					abort();
				}
			}
		}
	
	//
	// send FAMGMyColor to all neighbors with less weight
	//
	for( i = 0; i < NrSend; i++ )
	{
		MsgOutId[i] = SendASync( NbCh[SendQueue[i]], &FAMGMyColor, sizeof FAMGMyColor, &res );
		if( res != 0 )
		{
			cout << "ConstructColoring(): error "<<res<<" in SendASync FAMGMyColor for PE "<<Nb[SendQueue[i]]<<endl<<fflush;
			abort();
		}
		assert(MsgOutId[i]!=-1);
	}
	
	// wait for completion of sending
	j = NrSend;
	while( j>0 )
		for( i = 0; i < NrSend; i++ )
		{
			if( MsgOutId[i] != -1 ) // pending send weight
			{
				res = InfoASend( NbCh[SendQueue[i]], MsgOutId[i] );
				
				if( res == 1 )
				{	// message sended successfully
					j--;
					MsgOutId[i] = -1;	// reset
					if( j==0 )
						break;
				}
				else if( res != 0 )
				{
					cout << "ConstructColoring(): error: return status "<<res<<" in InfoASend for weight from PE "<<Nb[SendQueue[i]]<<endl<<fflush;
					abort();
				}
			}
		}
	
	//
	// close channels
	//
	for( i = 0; i < NrNb; i++ )
	{
		res = DiscASync( NbCh[i]);
		if( res != 0 )
		{
			cout << "ConstructColoring(): error "<<res<<" in DiscASync to PE "<<Nb[i]<<endl<<fflush;
			abort();
		}
		MsgOutId[i] = 1;	// set flag
	}
	
	// wait for completion of disconnecting
	j = NrNb;
	while( j>0 )
		for( i = 0; i < NrNb; i++ )
		{
			if( MsgOutId[i] != 0 ) // pending disconnect
			{
				res = InfoADisc( NbCh[i] );
				
				if( res == 1 )
				{	// channeld destructed successfully
					j--;
					MsgOutId[i] = 0;	// reset
					if( j==0 )
						break;
				}
				else if( res != 0 )
				{
					cout << "ConstructColoring(): error "<<res<<" in InfoADisc to PE "<<Nb[i]<<endl<<fflush;
					abort();
				}
			}
		}
	
	return 0;
}
