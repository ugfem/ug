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
    size = FAMGCEIL(size);
    buffer = malloc(size);
    if(buffer == NULL)
    {
         ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not allocate " << size << " byte." << endl;
         FAMGError(ostr);
    }
    else
    {
        ntop = nbottom = 0;
        bottom = (unsigned long) buffer;
        top = bottom + size;
    }
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
        top -= size;
        if(top < bottom)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory for " << size << " byte." << endl;
            FAMGError(ostr);
            return(NULL);
        }
        return((void *)top);
    }
	else 
    {
        ptr = (void *) bottom;
        bottom += size;
        if(top < bottom) 
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory for " << size << " byte." << endl;
            FAMGError(ostr);
            return(NULL);
        }
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

FAMGHeap *FAMGGetHeap()
{
    return famgheapptr;
}

