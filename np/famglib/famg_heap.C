/****************************************************************************/
/*																			*/
/* File:      famg_heap.C													*/
/*																			*/
/* Purpose:   famg heap class functions										*/
/*                                                                          */
/* Author:    Christian Wagner 					                            */
/*            Institut fuer Computeranwendungen  III			            */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27	                                        */
/*            70569 Stuttgart                                               */
/*            internet: chris@ica3.uni-stuttgart.de                         */
/*            											                    */
/*                                                                          */
/* History:   November 97 begin, Stuttgart                                  */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include <iostream.h>
#include <strstream.h>
#include <math.h>
#include <stdlib.h>
#include "famg_misc.h"
#include "famg_system.h"
#include "famg_heap.h"

/* RCS_ID
$Header$
*/

/****************************************************************************/
/*                                                                          */
/* Global Variable                                                          */
/*                                                                          */
/****************************************************************************/

FAMGHeap *famgheapptr;

FAMGHeap::~FAMGHeap()
{
    free(buffer);
}

FAMGHeap::FAMGHeap(unsigned long size)
{
	ntop = nbottom = 0;
    size = FAMGCEIL(size);
    buffer = malloc(size);
	bottom = (unsigned long) buffer;
    if(buffer == NULL)
    {
		ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not allocate " << size << " byte." << endl;
		FAMGError(ostr);
		size = 0;	// induce error in the first GetMem
    }
	top = bottom + size;
	info_max_bottom = 0;
	info_max_top = top;
	info_min_free = size;
	info_size = size;
}


void *FAMGHeap::GetMem(unsigned long size, int mode)
{
	void *ptr;

	if (size<=0)
	{
        return((void *)-1);
    }
    
    size = FAMGCEIL(size);
    if (mode==FAMG_FROM_TOP)
    {
		if (top<size)
			return NULL;
        top -= size;
        if(top < bottom)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory for " << size << " byte." << endl;
            FAMGError(ostr);
            return(NULL);
        }
		info_max_top = MIN(info_max_top,top);
		info_min_free = MIN(info_min_free,top-bottom);
        return((void *)top);
    }
	else 
    {
        ptr = (void *) bottom;
		if ((bottom+size) < bottom)
		{
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": exeeds max. address for " << size << " byte." << endl;
            FAMGError(ostr);
            return(NULL);
		}
        bottom += size;
        if(top < bottom) 
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory for " << size << " byte." << endl;
            FAMGError(ostr);
            return(NULL);
        }
		info_max_bottom = MAX(info_max_bottom,bottom);
		info_min_free = MIN(info_min_free,top-bottom);
        return(ptr);
    }
}

int FAMGHeap::Mark(int mode)
{
	
    if (mode==FAMG_FROM_TOP)
    {
        if (ntop >= FAMGMAXSTACK)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack too small." << endl;
            FAMGError(ostr);
            return(1);
        }
        topstack[ntop++] = top;
    }
	else 
    {
        if (nbottom >= FAMGMAXSTACK)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack too small." << endl;
            FAMGError(ostr);
            return(1);
        }
        bottomstack[nbottom++] = bottom;
    }

	return(0);
}


int FAMGHeap::Release(int mode)
{
	if (mode==FAMG_FROM_TOP)
	{
        if (ntop <= 0)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack error." << endl;
            FAMGError(ostr);
            return(1);
        }
        top = topstack[--ntop];
	}
	else
	{
        if (nbottom <= 0)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack error." << endl;
            FAMGError(ostr);
            return(1);
        }
        bottom = bottomstack[--nbottom];
	}
	return(0);
}

void FAMGHeap::PrintInfo()
// print statistics about the heap usage
{
	unsigned long comm_buffer[6], tmp_buffer[6];

	int l, size = 6*sizeof(unsigned long);	

	comm_buffer[0] = info_max_bottom-(unsigned long)buffer;
	comm_buffer[1] = comm_buffer[0];
	comm_buffer[2] = (unsigned long)buffer+info_size-info_max_top;
	comm_buffer[3] = comm_buffer[2];
	comm_buffer[4] = info_size-info_min_free;
	comm_buffer[5] = comm_buffer[4];
	
#ifdef ModelP
	// global MIN/MAX similar to UG_GlobalMinNINT
	for (l=degree-1; l>=0; l--)
	{
		GetConcentrate(l,tmp_buffer,size);
		comm_buffer[0] = MIN(comm_buffer[0],tmp_buffer[0]);
		comm_buffer[1] = MAX(comm_buffer[1],tmp_buffer[1]);
		comm_buffer[2] = MIN(comm_buffer[2],tmp_buffer[2]);
		comm_buffer[3] = MAX(comm_buffer[3],tmp_buffer[3]);
		comm_buffer[4] = MIN(comm_buffer[4],tmp_buffer[4]);
		comm_buffer[5] = MAX(comm_buffer[5],tmp_buffer[5]);
		
	}
	Concentrate(comm_buffer,size);
	Broadcast(comm_buffer,size);	

	if( me==master )
#endif
		// output of each quantity its min and max value across the pe's
		cout << "FAMGHeap: allocated " << info_size << " byte"
			 << " used " << comm_buffer[4] << ".." << comm_buffer[5] 
			 << " from bottom " << comm_buffer[0] << ".." << comm_buffer[1]
			 << " from top " << comm_buffer[2] << ".." << comm_buffer[3] << endl;
}

void *FAMGGetMem(unsigned long size, int mode)
{
    return famgheapptr->GetMem(size,mode);
}

int FAMGMarkHeap(int mode)
{
    return famgheapptr->Mark(mode);
}

int FAMGReleaseHeap(int mode)
{
    return famgheapptr->Release(mode);
}

void FAMGSetHeap(FAMGHeap *ptr)
{
    famgheapptr = ptr;
}

void FAMGFreeHeap()
{
	famgheapptr->PrintInfo();
	delete famgheapptr;
	FAMGSetHeap( NULL );
}

FAMGHeap *FAMGGetHeap()
{
    return famgheapptr;
}

