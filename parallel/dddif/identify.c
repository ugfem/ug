// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  identify.c													*/
/*																			*/
/* Purpose:   identification of distributed ug objects             			*/
/*																			*/
/* Author:	  Stefan Lang                    		 						*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   26.11.96 begin, first version extracted from refine.c			*/
/*																			*/
/* Remarks: 																*/
/*																			*/
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

#include "general.h"
#include "compiler.h"
#include "debug.h"
#include "gm.h"
#include "rm.h"
#include "refine.h"
#include "ddd.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* flags for identification */
#define CLEAR 0
#define IDENT 1

/* maximum count of objects for identification */
#define MAX_OBJECT	3

/* maximum count of tokens for identification */
#define MAX_TOKEN	10

/* determine the ddd header for identification of a node */
#define GET_IDENT_HDR(node) ( (NTYPE(node) == CORNER_NODE) ? \
								PARHDR(NFATHER(node)) : PARHDR(node) )

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/* temp node flag for Identification */
static INT ce_NIDENT;
#define NIDENT_LEN                    1
#define NIDENT(p)                     CW_READ(p,ce_NIDENT)
#define SETNIDENT(p,n)                CW_WRITE(p,ce_NIDENT,n)

#ifdef __THREEDIM__
/* temp edge flag for Identification */
static INT ce_EDIDENT;
#define EDIDENT_LEN                   1
#define EDIDENT(p)                    CW_READ(p,ce_EDIDENT)
#define SETEDIDENT(p,n)               CW_WRITE(p,ce_EDIDENT,n)
#endif

/* this function is called for low level identification */
static INT (*Ident_FctPtr) (DDD_HDR *IdentObjectHdr, INT nobject, 
			int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident) = NULL;

#ifdef Debug
static INT debug = 0;
static INT identlevel = 0;
#endif

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   name - short_description
   
   SYNOPSIS:
   
   PARAMETERS:
.  par - meaning
   
   DESCRIPTION:
   
   RETURN VALUE:
   
   SEE ALSO:
D*/
/****************************************************************************/

INT compare_gid (const void *e0, const void *e1)
{
	INT num0, num1;

	num0 = DDD_InfoGlobalId(*(DDD_HDR *)e0);
	num1 = DDD_InfoGlobalId(*(DDD_HDR *)e1);

	if (num0 < num1) return(1);
	if (num0 > num1) return(-1);
	return(0);
}

static void ResetIdentFlags (GRID *UpGrid)
{
	NODE *theNode;

	for (theNode=FIRSTNODE(UpGrid); theNode!=NULL; theNode=SUCCN(theNode)) {

		#ifdef Debug
		if (NIDENT(theNode) == IDENT)
			printf("%d: NGID=%d has NIDENT flag set\n",
				me,DDD_InfoGlobalId(PARHDR(theNode)));
		#endif
		SETNIDENT(theNode,CLEAR);
	}

}

#ifdef Debug
static INT Print_Identify_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
	INT i;

	ASSERT(nobject>0);
	ASSERT(nident>0);
	ASSERT(*proclist!=-1);

	/* print the interesting call parameters */
	PrintDebug("%d:    Print_Identify_ObjectList(): nobject=%d nident=%d"
				" skiptag=%d\n",me,nobject,nident,skiptag);

	/* print the objects to identify */
	PrintDebug("%d: l=%d   IdentObjectHdr:",me,identlevel);
	for (i=0; i<nobject; i++) {
		PrintDebug(" %d",DDD_InfoGlobalId(IdentObjectHdr[i]));
	}

	/* print the objects used for identify */
	PrintDebug("    IdentHdr:");
	for (i=0; i<nident; i++) {
		PrintDebug(" %d",DDD_InfoGlobalId(IdentHdr[i]));
	}

	/* print the proclist to identify to */
	PrintDebug("    ProcList: %d",me);
	while (*proclist != -1) {
		PrintDebug(" %d",*proclist);
		proclist += 2;
	}
	PrintDebug("\n");

	return;

}
#endif

static INT Identify_by_ObjectList (DDD_HDR *IdentObjectHdr, INT nobject, int *proclist, int skiptag, DDD_HDR *IdentHdr, INT nident)
{
	INT i,j,n;
	
	ASSERT(nobject>0);
	ASSERT(nident>0);
	ASSERT(*proclist!=-1);

	n = 0;
	while (*proclist != -1) {
		ASSERT(n<procs);

		if (*(proclist+1) == skiptag) {
			proclist += 2;
			continue;
		}

		/* identify the object */
		for (j=0; j<nobject; j++) {
			for (i=0; i<nident; i++) {

				PRINTDEBUG(dddif,5,("%d: Identify_by_ObjectList(): identcall" \
						" IdentObjectHdrGID=%d proclist=%d IdentHdr=%d\n",me,
						DDD_InfoGlobalId(IdentObjectHdr[j]),
						*proclist,
						DDD_InfoGlobalId(IdentHdr[i])));
				
				/* hand identification hdr to ddd */
				DDD_IdentifyObject(IdentObjectHdr[j], *proclist, IdentHdr[i]);
			}
		}

		n++;
		assert(n<procs);
		proclist += 2;
	}

	ASSERT(n>0);
}

#ifdef __THREEDIM__
static INT IdentifySideVector (ELEMENT* theElement, ELEMENT *theNeighbor, ELEMENT *Son, INT SonSide) 
{
	INT k,nident;
	DDD_HDR IdentHdr[MAX_TOKEN];
	int *proclist;
	NODE *theNode;

	nident = 0;

	proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

	/* identify using corner nodes */
	for (k=0; k<CORNERS_OF_SIDE(Son,SonSide); k++) { 
		theNode = CORNER(Son,CORNER_OF_SIDE(Son,SonSide,k));
		IdentHdr[nident++] = PARHDR(theNode);
	}

	proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

	Ident_FctPtr(&(PARHDR(SVECTOR(Son,SonSide))), 1, 
		proclist+2, PrioGhost, IdentHdr, nident); 

}
#endif

static void IdentifyNode (ELEMENT *theNeighbor, NODE *theNode, NODE *Nodes[MAX_SIDE_NODES], INT node, INT ncorners, INT Vec)
{
	INT nobject,nident;
	DDD_HDR IdentObjectHdr[MAX_OBJECT];
	DDD_HDR IdentHdr[MAX_TOKEN];

	nobject = nident = 0;

	/* is this node identified? */
	#ifdef Debug 
	if (debug == 1) {
		if (NIDENT(theNode) == CLEAR) return;
	}
	else
	#endif
	/* return if already identified */
	if (NIDENT(theNode) == IDENT) return;

	switch (NTYPE(theNode)) {
		int *proclist;

		case (CORNER_NODE): 

			PRINTDEBUG(dddif,1,("%d: Identify CORNERNODE gid=%08x node=%d vec=%d\n",
					me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

			IdentObjectHdr[nobject++] = PARHDR(theNode);
			if (Vec)
				IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

			/* identify to proclist of node */
			proclist = DDD_InfoProcList(PARHDR(NFATHER(theNode)));

			/* identify using father node */
			IdentHdr[nident++] = PARHDR(NFATHER(theNode));

			Ident_FctPtr(IdentObjectHdr, nobject, 
				proclist+2, PrioGhost, IdentHdr, nident); 

			break;

		case (MID_NODE): {

			#ifdef __TWODIM__
			NODE **EdgeNodes;
			EdgeNodes = Nodes;
			#endif

			#ifdef __THREEDIM__
			NODE *EdgeNodes[MAX_SIDE_NODES];
			EDGE *theEdge;

			EdgeNodes[0] = Nodes[node-ncorners];
			EdgeNodes[1] = Nodes[(node-ncorners+1)%ncorners];
			EdgeNodes[2] = theNode;
			#endif

			ASSERT(EdgeNodes[0]!=NULL);
			ASSERT(EdgeNodes[1]!=NULL);
			ASSERT(EdgeNodes[2]!=NULL);

			PRINTDEBUG(dddif,1,("%d: Identify MIDNODE gid=%08x node=%d Vec=%d\n",
					me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

			/* identify midnode, vertex, vector */
			IdentObjectHdr[nobject++] = PARHDR(theNode);
			IdentObjectHdr[nobject++] = PARHDRV(MYVERTEX(theNode));
			if (Vec)
				IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

			#ifdef __TWODIM__
			/* 2D: identify to proclist of neighbor element */
			proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
			#endif

			#ifdef __THREEDIM__
			/* 3D: identify to proclist of edge */
			theEdge = GetEdge(NFATHER(EdgeNodes[0]),NFATHER(EdgeNodes[1]));
			ASSERT(theEdge!=NULL);

			proclist = DDD_InfoProcList(PARHDR(theEdge));
			#endif

			/* identify using edge nodes */
			IdentHdr[nident++] = PARHDR(NFATHER(EdgeNodes[0]));
			IdentHdr[nident++] = PARHDR(NFATHER(EdgeNodes[1]));

			Ident_FctPtr(IdentObjectHdr, nobject, 
				proclist+2, PrioGhost, IdentHdr, nident); 

			break;
		}

		#ifdef __THREEDIM__
		case (SIDE_NODE): {
			
			INT  i;

			PRINTDEBUG(dddif,1,("%d: Identify SIDENODE gid=%08x node=%d Vec=%d\n",
					me, DDD_InfoGlobalId(PARHDR(theNode)), node, Vec));

			/* identify sidenode, vertex and vector */ 
			IdentObjectHdr[nobject++] = PARHDR(theNode);
			IdentObjectHdr[nobject++] = PARHDRV(MYVERTEX(theNode));
			if (Vec)
				IdentObjectHdr[nobject++] = PARHDR(NVECTOR(theNode));

			/* identify to proclist of neighbor element */
			proclist = DDD_InfoProcList(PARHDRE(theNeighbor));

			/* identify using corner nodes of side */
			for (i=0; i<ncorners; i++)
				IdentHdr[nident++] = PARHDR(NFATHER(Nodes[i]));

			/* identify side node */
			Ident_FctPtr(IdentObjectHdr, nobject, 
				proclist+2, PrioGhost, IdentHdr, nident); 

			break;
		}
		#endif
		default:
			ASSERT(0);
			break;
	}

	#ifdef Debug 
	if (debug == 1) {
		SETNIDENT(theNode,CLEAR);
	}
	else
	#endif
	/* lock this node for identification */
	SETNIDENT(theNode,IDENT);

	return;
}

#ifdef __THREEDIM__
EDGE *FatherEdge(SideNodes,Nodes)
{
	INT sonedge,i;
	EDGE *fatherEdge = NULL;

	ASSERT(Nodes[0]!=NULL);
	ASSERT(Nodes[1]!=NULL);

	for (i=0; i<MAX_SIDE_NODES; i++) {
		if (SideNodes[i] == Nodes[0])
			break;
	}
	ASSERT(i<MAX_SIDE_NODES)

	switch (NTYPE[Nodes[0]) {

		case (CORNER_NODE):

			if ( (SideNodes[(i+1)%ncorners] == Nodes[1]) ||
			     (SideNodes[(i+1)%ncorners+ncorners] == Nodes[1]) ) {
				
				fatherEdge = GetEdge(NFATHER(Nodes[0],NFATHER(SideNodes[(i+1)%ncorners]));
				ASSERT(fatherEdge!=NULL);
			}

			if ( (SideNodes[(i-1+ncorners)%ncorners] == Nodes[1]) ||
			     (SideNodes[(i-1+ncorners)%ncorners+ncorners] == Nodes[1]) ) {

				fatherEdge = GetEdge(NFATHER(Nodes[0],
									NFATHER(SideNodes[(i-1+ncorners)%ncorners]));
				ASSERT(fatherEdge!=NULL);
			}
			
			break;

		case (MID_NODE):

			if (SideNodes[(i+1)%ncorners] == Nodes[1]) {

				fatherEdge = GetEdge(NFATHER(Nodes[0],NFATHER(SideNodes[(i+1)%ncorners]));
				ASSERT(fatherEdge!=NULL);
			}

			if (SideNodes[(i-1+ncorners)%ncorners] == Nodes[1]) ) {

				fatherEdge = GetEdge(NFATHER(Nodes[0],
									NFATHER(SideNodes[(i-1+ncorners)%ncorners]));
				ASSERT(fatherEdge!=NULL);
			}

			break;

		case (SIDE_NODE):

			break;

		default:
			assert(0);
			break;
	}
	
	return(fatherEdge);
}
#endif
		
static INT IdentifyEdge (ELEMENT *theElement, ELEMENT *theNeighbor, NODE **SideNodes, INT ncorners, ELEMENT *Son, INT SonSide, INT k, INT Vec)
{
	NODE *Nodes[2];
	EDGE *theEdge;
	VECTOR *theVector;
	INT nobject,nident;
	INT *proclist;
	DDD_HDR IdentObjectHdr[MAX_OBJECT];
	DDD_HDR IdentHdr[MAX_TOKEN];

	nobject = nident = 0;

	#ifdef __TWODIM__
	Nodes[0] = CORNER(Son,CORNER_OF_EDGE(Son,SonSide,0));
	Nodes[1] = CORNER(Son,CORNER_OF_EDGE(Son,SonSide,1));
	#endif

	#ifdef __THREEDIM__
	Nodes[0] = CORNER(Son,CORNER_OF_EDGE(Son,EDGE_OF_SIDE(Son,SonSide,k),0));
	Nodes[1] = CORNER(Son,CORNER_OF_EDGE(Son,EDGE_OF_SIDE(Son,SonSide,k),1));
	#endif

	ASSERT(Nodes[0]!=NULL);
	ASSERT(Nodes[1]!=NULL);

	theEdge = GetEdge(Nodes[0],Nodes[1]);
	ASSERT(theEdge!=NULL);

	PRINTDEBUG(dddif,1,("%d: Identify EDGE nedge=%d pe=%08x/%x eID=%d"
		" ntype0=%d  ntype1=%d\n",me,k,
		DDD_InfoGlobalId(PARHDRE(Son)),Son,ID(Son),
		NTYPE(Nodes[0]), NTYPE(Nodes[1])))

	#ifdef __THREEDIM__
	IdentObjectHdr[nobject++] = PARHDR(theEdge);
	#endif
	if (Vec)
		IdentObjectHdr[nobject++] = PARHDR(EDVECTOR(theEdge));

	#ifdef __TWODIM__
	/* identify to proclist of neighbor */
	proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
	#endif

	#ifdef __THREEDIM__
	/* identify to proclist of edge */
	{
	EDGE *fatherEdge = NULL;

	/* check whether edge inside the side of the element */
	fatherEdge = FatherEdge(SideNodes,Nodes);
	
	if (fatherEdge != NULL)		
		proclist = DDD_InfoProcList(PARHDR(fatherEdge));
	else
		proclist = DDD_InfoProcList(PARHDRE(theNeighbor));
	}
	#endif

	IdentHdr[nident++] = PARHDR(Nodes[0]);
	IdentHdr[nident++] = PARHDR(Nodes[1]);

	Ident_FctPtr(IdentObjectHdr, nobject, 
		proclist+2, PrioGhost, IdentHdr, nident); 
}


static INT IdentifyObjectsOfElementSide(GRID *theGrid, ELEMENT *theElement, INT i, ELEMENT *theNeighbor) 
{
	INT nodes,j;
	NODE *SideNodes[MAX_SIDE_NODES];
	INT  ncorners;
	NODE *theNode;

	GetSonSideNodes(theElement,i,&nodes,SideNodes);
	ncorners = CORNERS_OF_SIDE(theElement,i);

	PRINTDEBUG(dddif,1,("IdentifyObjectsOfElementSide():identify NODES ncorners=%d nodes=%d\n",
		ncorners,nodes));

	/* identify nodes, vertices and node vectors of son elements */
	for (j=0; j<nodes; j++) {
		INT prio;

		theNode = SideNodes[j];
		ASSERT(theNode != NULL);

		/* identify new node including its vector and vertex 	*/
		IdentifyNode(theNeighbor, theNode, SideNodes, j, ncorners,
			TYPE_DEF_IN_GRID(theGrid,NODEVECTOR));
	}

	/* identify edge vectors (2D); edges, edge and side vectors (3D) */
	if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3) {

		ELEMENT *SonList[MAX_SONS];
		INT SonsOfSide,SonSides[MAX_SONS];
		INT j;

		PRINTDEBUG(dddif,1,("%d: IdentifyObjectsOfElementSide(): identify EDGES and VECTORS\n",me));

		if (Get_Sons_of_ElementSide(theElement,i,&SonsOfSide,
				SonList,SonSides,1)!=GM_OK) 
			RETURN(GM_FATAL);
				
		for (j=0; j<SonsOfSide; j++) {

			if (TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR) || DIM==3) {

				INT k;
				INT nedges = EDGES_OF_SIDE(SonList[j],SonSides[j]);
			
				/* identify the edge and vector */
				for (k=0; k<nedges; k++) {
					IdentifyEdge(theElement,theNeighbor,SideNodes,ncorners,
								 SonList[j],SonSides[j],k,
								 TYPE_DEF_IN_GRID(theGrid,EDGEVECTOR));
				}
			}

			#ifdef __THREEDIM__
			if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR)) {

				IdentifySideVector(theElement,theNeighbor,SonList[j],SonSides[j]);
			}
			#endif
		}
	}
}

INT	IdentifyDistributedObjects (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
	INT l,i,j,prio;
	ELEMENT *theElement,*theNeighbor;
	NODE *theNode;
	GRID *theGrid;
	
	PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): FromLevel=%d ToLevel=%d\n",
				me,FromLevel,ToLevel));

	/* identify distributed objects */
	for (l=FromLevel; l<ToLevel; l++) {

		PRINTDEBUG(dddif,1,("%d: IdentifyDistributedObjects(): identification level=%d\n",me,l));

		theGrid = GRID_ON_LEVEL(theMG,l);

		#ifdef Debug
		identlevel = l;
		#endif

		/* reset control word flags for ident on upper level */
		/* TODO: is resetting needed? */
		ResetIdentFlags(GRID_ON_LEVEL(theMG,l+1));

		for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement)) { 

			if (!IS_REFINED(theElement) ||
				(prio = DDD_InfoPriority(PARHDRE(theElement))) == PrioGhost) {
				continue;
			}

			PRINTDEBUG(dddif,0,("%d: element loop\n",me));

			for (i=0; i<SIDES_OF_ELEM(theElement); i++) {

				theNeighbor = NBELEM(theElement,i);
				if (theNeighbor == NULL) continue;

				PRINTDEBUG(dddif,0,("%d: side loop prio=%d nsons=%d\n",me,
					DDD_InfoPriority(PARHDRE(theNeighbor)),
					NSONS(theNeighbor)));

				if ((prio = DDD_InfoPriority(PARHDRE(theNeighbor))) != PrioGhost
					 || NSONS(theNeighbor)!=0) 
					 continue;

				PRINTDEBUG(dddif,1,("%d: Identify element: pe=%08x/%x eID=%d side=%d\n",me,
							DDD_InfoGlobalId(PARHDRE(theElement)),theElement,
							ID(theElement),i));

				IdentifyObjectsOfElementSide(theGrid,theElement,i,theNeighbor);
			}
		}
	}
	
	return(GM_OK);
}


INT	IdentifyGridLevels (MULTIGRID *theMG, INT FromLevel, INT ToLevel)
{
	#ifdef Debug
	debug = 0;
	#endif

	/* allocate a control word entry to lock nodes */
	if (AllocateControlEntry(NODE_CW,1,&ce_NIDENT) != GM_OK)
		assert(0);

	#ifdef __THREEDIM__
	/* allocate a control word entry to lock edges */
	if (AllocateControlEntry(EDGE_CW,1,&ce_EDIDENT) != GM_OK)
		assert(0);
	#endif

	/* set Ident_FctPtr to identification mode */
	Ident_FctPtr = Identify_by_ObjectList;

	/* identify new created objects */
	DDD_IdentifyBegin();

	/* give identification calls to DDD */
	IdentifyDistributedObjects(theMG,FromLevel,ToLevel);

	/* start identification process */
	DDD_IdentifyEnd();

	/* only for debugging */
	IFDEBUG(dddif,1)

	debug = 1;

	/* set Ident_FctPtr to print mode */
	Ident_FctPtr = Print_Identify_ObjectList;

	PrintDebug("AFTER Identify\n");
	IdentifyDistributedObjects(theMG,FromLevel,ToLevel);

	ENDDEBUG

	FreeControlEntry(ce_NIDENT);

	#ifdef __THREEDIM__
	FreeControlEntry(ce_EDIDENT);
	#endif
}

#endif /* end ModelP */
