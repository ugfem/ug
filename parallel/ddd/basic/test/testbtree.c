// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* some macros for customizing oopp */
#define _NEWPARAMS
#define _NEWPARAMS_OR_VOID    void

#define __INDENT(n)   { int i; for(i=0; i<n; i++) fputs("   ",fp);}
#define _PRINTPARAMS  , int indent, FILE *fp
#define _PRINTPARAMS_DEFAULT  ,0,stdout
#define _INDENT       __INDENT(indent)
#define _INDENT1      __INDENT(indent+1)
#define _INDENT2      __INDENT(indent+2)
#define _PRINTNEXT    , indent+1, fp
#define _PRINTSAME    , indent, fp

#define _CHECKALLOC(ptr)   assert(ptr!=NULL)
#define OO_Allocate        malloc
#define OO_Free            free

#include "../oopp.h"    /* for object-orientated style via preprocessor */

/****************************************************************************/
/* Tree element:                                                            */
/****************************************************************************/

#define ClassName TestTreeElement
Class_Data_Begin
int value;
Class_Data_End
void Method(Print)   (DefThis _PRINTPARAMS);
int  Method(Compare) (ClassPtr, ClassPtr);

#undef ClassName


/* define container class */
#ifndef SetOf  /* necessary for inline documentation only */
#define SetOf          TestTreeElement
#define Set_SegmSize   256
#define Set_BTreeOrder 32
#ifdef XferMemFromHeap
#define ArrayAllocate  xfer_AllocHeap
#define NoArrayFree
#endif
#endif

// Create method definitions, not just declarations
#define ContainerImplementation

#include "../ooppcc.h"

void TestTreeElement_Print(_TestTreeElement* element, int indent, FILE* fp)
{
  printf("TestTreeElement: %d\n", element->value);
}

int TestTreeElement_Compare(_TestTreeElement* a, _TestTreeElement* b)
{
  if (a->value < b->value)
    return -1;
  else if (a->value > b->value)
    return 1;

  return 0;
}

int main(int argc, char** argv)
{
  TestTreeElementBTree* foo = New_TestTreeElementBTree();

  for (int i=0; i<40; i+=2) {

    TestTreeElement* newItem = new TestTreeElement;
    newItem->value = i;

    TestTreeElementBTree_Insert(foo, newItem);
  }

  for (int i=1; i<40; i+=2) {

    TestTreeElement* newItem = new TestTreeElement;
    newItem->value = i;

    TestTreeElementBTree_Insert(foo, newItem);
  }

  TestTreeElementBTree_Print(foo, 0, stdout);

  TestTreeElementBTree_Reset(foo);

  return 0;
}
