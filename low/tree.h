// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      tree.h                                                        */
/*                                                                          */
/* Purpose:   header file for quad- and oct-tree                            */
/*                                                                          */
/* Author:      Carsten Schwarz                                             */
/*              Institut fuer Hydromechanik und Wasserwirtschaft            */
/*              ETH Hoenggerberg                                            */
/*              8093 Zuerich                                                */
/*                                                                          */
/* History:   07.04.97 begin, ug version 3.4                                */
/*                                                                          */
/* Revision:  07.04.97                                                      */
/*                                                                          */
/****************************************************************************/
/* RCS_ID
   $Header$
 */
/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __TREE__
#define __TREE__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __FIFO__
#include "fifo.h"
#endif

#ifndef __HEAPS__
#include "heaps.h"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXTREEDIM DIM

#define TREE_CHANGED 1
#define TREE_SEARCH  2
#define TREE_INVALID 255

#define TREELEAF 1
#define TREENODE   2

#define TREEROOT(p)             ((p)->root)
#define TREEPOS(p,i,j)          ((p)->posrange[(i) + (p)->dim * (j)])
#define TREESEARCHLL(p,i)       TREEPOS((p),(i),2)
#define TREESEARCHUR(p,i)       TREEPOS((p),(i),3)
#define TNODETYPE(p)            ((p)->etype)
#define TNODEFATHER(p)          ((p)->tnode.father)
#define TNODESON(p)             ((p)->tnode.son)
#define TNODENEXT(p)            ((p)->tnode.next)
#define TNODEPOS(p,i,j,d)       ((p)->tnode.boxcorners[(i) + (d)*(j)])
#define TNODEOBJ(p)             ((p)->tleaf.obj)
#define TLEAFPOS(p,i)           ((p)->tleaf.pos[(i)])


/****************************************************************************/
/*                                                                            */
/* data structures exported by the corresponding source file                */
/*                                                                            */
/****************************************************************************/

struct tree_node {
  INT etype;
  union tree_entry *father;
  union tree_entry *son;
  union tree_entry *next;
  DOUBLE boxcorners[1];
};

typedef struct tree_node TREE_NODE;

struct tree_leaf {
  INT etype;
  union tree_entry *father;
  void *obj;
  DOUBLE pos[1];
};

typedef struct tree_leaf TREE_LEAF;

union tree_entry {
  INT etype;
  TREE_NODE tnode;
  TREE_LEAF tleaf;
};

typedef union tree_entry TREE_ENTRY;

typedef struct {
  INT status;
  size_t fifo_max_mem;
  HEAP *heap;
  FIFO *fifo;
  INT dim;
  TREE_ENTRY *root;
  DOUBLE posrange[1];
} TREE;


/****************************************************************************/
/*                                                                            */
/* function declarations                                                    */
/*                                                                            */
/****************************************************************************/

TREE *CreateTree (HEAP *theHeap, INT dim, DOUBLE *posrange);
INT DeleteTree(TREE *theTree);
INT InsertinTree (TREE *theTree, DOUBLE *Position, void *obj);
void *DeleteObjinTree (TREE *theTree, DOUBLE *Position);
TREE_ENTRY *GetFirstLeafinQuader(TREE *theTree, DOUBLE *ll, DOUBLE *ur);
TREE_ENTRY *GetNextLeafinQuader(TREE *theTree);
#endif
