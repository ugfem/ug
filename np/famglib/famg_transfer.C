/****************************************************************************/
/*																			*/
/* File:      famg_transfer.C												*/
/*																			*/
/* Purpose:   famg transfer classes functions								*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   November 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#include "famg_transfer.h"
#include "famg_heap.h"
#include "famg_misc.h"
#include "famg_grid.h"

/* RCS_ID
$Header$
*/

// Class FAMGTransferEntry

void FAMGTransferEntry::Init(int i)
{
    id.f0 = i;
    next = NULL;
    data = 0.0;

    return;
}        

FAMGTransferEntry *FAMGTransferEntry::GetEntry(int j)
{
    FAMGTransferEntry *transij;

    for(transij = this; transij != NULL; transij = transij->GetNext())
    {
        if(transij->GetId() == j) return &(*transij);
    }

    return NULL;
}



FAMGTransferEntry *FAMGTransferEntry::NewEntry(FAMGTransferEntry* rowj)
{
    FAMGTransferEntry *transij, *transji, *ptr;
    int i,j;

    // Allocate transij and transji at once. That's not nice but it saves
    // an pointer.

    ptr = (FAMGTransferEntry *) FAMGGetMem(sizeof(FAMGTransferEntry[2]),FAMG_FROM_TOP);
    if(ptr == NULL) return(NULL);

    i = id.f0;
    j = rowj->GetId();
    ptr->SetRev(0);
    (ptr+1)->SetRev(1);
    transij = ptr;
    transji = ptr+1;

    transij->SetData(0.0);
    transji->SetData(0.0);
    transij->SetId(j);
    transji->SetId(i);

    transij->SetNext(next);
    transji->SetNext(rowj->GetNext());
    next = transij;
    rowj->SetNext(transji);

    return transij;
}
    
int FAMGTransferEntry::SaveEntry(FAMGTransferEntry* rowj, double val)
{
    FAMGTransferEntry *transij;
    int j;

    // test 
    if(Abs(val) < 1e-20) return 0;

    j = rowj->GetId();
    if(id.f0 == j) // diagonal entry
    {
        data += val;
        return 0;
    }
    transij = GetEntry(j);
    if(transij == NULL)
    {
        transij = NewEntry(rowj);
        if(transij == NULL) return 1;
    }

    transij->AddData(val);
    return 0;
}

int FAMGTransferEntry::SaveEntry(FAMGTransferEntry* rowj, double val,FAMGTransferEntry** p)
{
    FAMGTransferEntry *transij;
    int j;

    *p = NULL;    
    // test 
    if(Abs(val) < 1e-20) return 0;

    j = rowj->GetId();
    if(id.f0 == j) // diagonal entry
    {
        data += val;
        return 0;
    }
    transij = GetEntry(j);
    if(transij == NULL)
    {
        transij = NewEntry(rowj);
        if(transij == NULL) return 1;
    }

    transij->AddData(val);
    *p = transij;
    return 0;
}

      

// Class FAMGTransfer


int FAMGTransfer::Init(FAMGGrid *grid) 
{
    FAMGTransferEntry *rowi;
    int i;

    n = grid->GetN();

    row = (FAMGTransferEntry*) FAMGGetMem(n*sizeof(FAMGTransferEntry), FAMG_FROM_TOP);
    if (row == NULL) return(1);
    for(i = 0; i < n; i++)
    {
        rowi = row+i;
        rowi->Init(i);
        rowi->SetData(1.0);
    }
    
    return(0);
} 

int FAMGTransfer::Order(int *mapping)
{
    FAMGTransferEntry *helptrans, *transij;
    int i;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helptrans = (FAMGTransferEntry *) FAMGGetMem(n*sizeof(FAMGTransferEntry), FAMG_FROM_TOP);
    if (helptrans == NULL) return 1;
    for(i = 0; i < n; i++) helptrans[i] = row[i];
    for(i = 0; i < n; i++) row[mapping[i]] = helptrans[i];
    for(i = 0; i < n; i++)
    {
        for(transij = row+i; transij != NULL; transij = transij->GetNext())
        {
            transij->SetId(mapping[transij->GetId()]);
        }
    }
    FAMGReleaseHeap(FAMG_FROM_TOP);
    
    return 0;
}

int FAMGTransfer::Reorder(int *mapping)
{
    FAMGTransferEntry *helptrans, *transij, *transji;
    int i;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helptrans = (FAMGTransferEntry *) FAMGGetMem(n*sizeof(FAMGTransferEntry), FAMG_FROM_TOP);
    if (helptrans == NULL) return 1;
    for(i = 0; i < n; i++) helptrans[i] = row[i];
    for(i = 0; i < n; i++) row[i] = helptrans[mapping[i]];
    for(i = 0; i < n; i++)
    {
        (row+i)->SetId(i);
        for(transij = (row+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            transji = transij->GetReverse();
            transji->SetId(i);
        }
    }
    FAMGReleaseHeap(FAMG_FROM_TOP);
    
    return 0;
}
    
