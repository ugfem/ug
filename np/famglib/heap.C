/****************************************************************************/
/*																			*/
/* File:      heap.C														*/
/*																			*/
/* Purpose:   cmg heap class functions										*/
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
#include "misc.h"
#include "system.h"
#include "heap.h"

/* RCS_ID
$Header$
*/

/****************************************************************************/
/*                                                                          */
/* Global Variable                                                          */
/*                                                                          */
/****************************************************************************/

CMGHeap *cmgheapptr;

CMGHeap::~CMGHeap()
{
    free(buffer);
}

CMGHeap::CMGHeap(unsigned long size)
{
    size = CMGCEIL(size);
    buffer = malloc(size);
    if(buffer == NULL)
    {
         ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not allocate " << size << " byte." << endl;
         CMGError(ostr);
    }
    else
    {
        ntop = nbottom = 0;
        bottom = (unsigned long) buffer;
        top = bottom + size;
    }
}


void *CMGHeap::GetMem(unsigned long size, int mode)
{
	void *ptr;

	if (size<=0)
	{
        return((void *)-1);
    }
    
    size = CMGCEIL(size);
    if (mode==CMG_FROM_TOP)
    {
        top -= size;
        if(top < bottom)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory for << size << byte." << endl;
            CMGError(ostr);
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
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": not enough memory." << endl;
            CMGError(ostr);
            return(NULL);
        }
        return(ptr);
    }
}

int CMGHeap::Mark(int mode)
{
	
    if (mode==CMG_FROM_TOP)
    {
        if (ntop >= CMGMAXSTACK)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack too small." << endl;
            CMGError(ostr);
            return(1);
        }
        topstack[ntop++] = top;
    }
	else 
    {
        if (nbottom >= CMGMAXSTACK)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack too small." << endl;
            CMGError(ostr);
            return(1);
        }
        bottomstack[nbottom++] = bottom;
    }

	return(0);
}


int CMGHeap::Release(int mode)
{
	if (mode==CMG_FROM_TOP)
	{
        if (ntop <= 0)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack error." << endl;
            CMGError(ostr);
            return(1);
        }
        top = topstack[--ntop];
	}
	else
	{
        if (nbottom <= 0)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": stack error." << endl;
            CMGError(ostr);
            return(1);
        }
        bottom = bottomstack[--nbottom];
	}
	return(0);
}

void *CMGGetMem(unsigned long size, int mode)
{
    return cmgheapptr->GetMem(size,mode);
}

int CMGMarkHeap(int mode)
{
    return cmgheapptr->Mark(mode);
}

int CMGReleaseHeap(int mode)
{
    return cmgheapptr->Release(mode);
}

void CMGSetHeap(CMGHeap *ptr)
{
    cmgheapptr = ptr;
}

CMGHeap *CMGGetHeap()
{
    return cmgheapptr;
}

