/****************************************************************************/
/*																			*/
/* File:	  famg_graphics.c												*/
/*																			*/
/* Purpose:   graphic functions             								*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   August 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
 

extern "C"
{
#include "wop.h"
#include "wpm.h"
#include "misc.h"
#include "evm.h"
#include "cw.h"
#include "graph.h"
#include "gm.h"
#include "commands.h"
}
    
#include "famg_uginterface.h"

/* RCS_ID
$Header$
*/

static long BlackColor; /* Black  */
static long RedColor; /* Red  */
static VECTOR *GlobalVec1;
static VECTOR *GlobalVec2;
static INT GlobalIds;

    
struct FAMGPlotObject
{
    struct PlotObjHead theHead;    /* the head */
    int level;
    int ids;
};


static INT SetFAMGGraph (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
	BVP_DESC *theBVPDesc;
    struct FAMGPlotObject *theObj;
    int l,i;
    
	theObj = (struct FAMGPlotObject *) &(thePlotObj->theExternObject);
	theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
	V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
	PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

    if (ReadArgvINT("l",&l,argc,argv)) l = 999;
    theObj->level = l;
    if (ReadArgvINT("i",&i,argc,argv)) i = 1;
    if(i > 0) theObj->ids = 1;
    else theObj->ids = 0;
   
	return (ACTIVE);
}
    
static INT  DisplayFAMGGraph (PLOTOBJ *thePlotObj)
{
    struct FAMGPlotObject *theObj;
    
	theObj = (struct FAMGPlotObject *) &(thePlotObj->theExternObject);

    UserWriteF("level: %d \n",theObj->level); 
	return (0);
}

static INT PreProcessFAMGGraph (PICTURE *thePicture, WORK *theWork)
{
	OUTPUTDEVICE *theOD;
    struct FAMGPlotObject *theObj;
    int level,maxlevel;
    GRID *grid;
    MULTIGRID *mg;
    
	mg = GetCurrentMultigrid();
	if (mg==NULL) return(OKCODE);

	theObj = (struct FAMGPlotObject *) &(PIC_PO(thePicture)->theExternObject);
	theOD  = PIC_OUTPUTDEV(thePicture);

	BlackColor = theOD->black;
	RedColor = theOD->red;
    if(theObj->level == 999) 
    {
        level =  CURRENTLEVEL(mg);
    }
    else
    {
        level = theObj->level;
    }
    grid =  GRID_ON_LEVEL(mg,level);
    if(grid == NULL)
    {
        GlobalVec1 = NULL;
        GlobalVec2 = NULL;
        return (0);
    }

    GlobalVec1 = FIRSTVECTOR(grid);
    GlobalVec2 = FIRSTVECTOR(grid);
    GlobalIds = theObj->ids;

    
    return(0);
}


static INT EvalFAMGGraph1 (DRAWINGOBJ *theDO, VECTOR *vec)
{
    VERTEX *vertex, *nbvertex;
    DOUBLE_VECTOR mypos,nbpos;
    int j;
    VECTOR *nbvec;
    MATRIX *mat, *imat;

    UgSetLineWidth(1);


    if(vec == NULL) return (0);

    vertex = MYVERTEX(VMYNODE(vec));
    V_DIM_COPY(CVECT(vertex),mypos);


    /* plot  marker */
    if(VCCOARSE(vec))
    {
        DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO); 
        DO_2c(theDO) = 1; DO_inc(theDO); 
        DO_2l(theDO) = BlackColor; DO_inc(theDO);
        DO_2s(theDO) = FILLED_CIRCLE_MARKER; DO_inc(theDO); 
        DO_2s(theDO) = 8; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }
    else
    {
        DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO); 
        DO_2c(theDO) = 1; DO_inc(theDO); 
        DO_2l(theDO) = BlackColor; DO_inc(theDO);
        DO_2s(theDO) = EMPTY_CIRCLE_MARKER; DO_inc(theDO); 
        DO_2s(theDO) = 8; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }
    

    if(GlobalIds)
    {
        /* print id */
        DO_2c(theDO) = DO_TEXT; DO_inc(theDO);
        DO_2l(theDO) = BlackColor; DO_inc(theDO);
        DO_2c(theDO) = TEXT_REGULAR; DO_inc(theDO); 
        DO_2c(theDO) = TEXT_NOT_CENTERED; DO_inc(theDO); 
        DO_2s(theDO) = 6; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        sprintf(DO_2cp(theDO),"%d",VINDEX(vec));
        DO_inc_str(theDO);
    }

    /* plot matrix */
    for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
    {
        nbvec = MDEST(mat);
        nbvertex = MYVERTEX(VMYNODE(nbvec));
        V_DIM_COPY(CVECT(nbvertex),nbpos);
        DO_2c(theDO) = DO_LINE; DO_inc(theDO); 
        DO_2l(theDO) =BlackColor; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }

    DO_2c(theDO) = DO_NO_INST;

            
    return(0);
    
}
static INT EvalFAMGGraph2 (DRAWINGOBJ *theDO, VECTOR *vec)
{
    VERTEX *vertex, *nbvertex;
    DOUBLE_VECTOR mypos,nbpos;
    int j;
    VECTOR *nbvec;
    MATRIX *mat, *imat;

    UgSetLineWidth(1);

    if(vec == NULL) return (0);

    vertex = MYVERTEX(VMYNODE(vec));
    V_DIM_COPY(CVECT(vertex),mypos);



    /* plot transfer */
    for (imat=VISTART(vec); imat!=NULL; imat=MNEXT(imat))
    {
        nbvec = MDEST(imat);
        nbvertex = MYVERTEX(VMYNODE(nbvec));
        V_DIM_COPY(CVECT(nbvertex),nbpos);
        DO_2c(theDO) = DO_LINE; DO_inc(theDO); 
        DO_2l(theDO) =RedColor; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2); 
     }

     DO_2c(theDO) = DO_NO_INST;

    return(0);
    
}


static INT EvalFAMGGraph (DRAWINGOBJ *theDO, INT *end)
{

    if(GlobalVec2 == NULL) { *end = 1; return(0);}
    else *end = 0;
    if((GlobalVec2 != NULL) && (GlobalVec1 == NULL))
    {
        EvalFAMGGraph2(theDO,GlobalVec2);
        GlobalVec2 = SUCCVC(GlobalVec2);
    }
    else if(GlobalVec1 != NULL)
    {
        EvalFAMGGraph1(theDO,GlobalVec1);
        GlobalVec1 = SUCCVC(GlobalVec1);
    }

    if(GlobalVec2 == NULL) *end = 1;
    else *end = 0;
            
    return(0);
    
}

    
static INT PostProcessFAMGGraph (PICTURE *thePicture, WORK *theWork)
{
    return(0);
}


INT InitFAMGGraph (void)
{
	PLOTOBJHANDLING *thePOH;
	WORKPROCS *theWP;
    EXTERNWORK *theEXW;
	PLOTOBJTYPE *thePOT;

    /* create WorkHandling for 'FAMGGraph' */
    if ((thePOH=CreatePlotObjHandling ("FAMGGraph")) 	== NULL) 
		return (1);
  
    POH_DYNAMIC_INFO(thePOH) = NULL;
    POH_CLICKACTION(thePOH)  = NULL;
		
    /* draw work */
    POH_NBCYCLES(thePOH,DRAW_WORK) = 1;
		
    theWP = POH_WORKPROGS(thePOH,DRAW_WORK,0);
    WP_WORKMODE(theWP) = EXTERN;
    theEXW = WP_EXTERNWISE(theWP);
    theEXW->EXT_PreProcessProc = PreProcessFAMGGraph;
    theEXW->EXT_EvaluateProc   = EvalFAMGGraph;
    theEXW->EXT_ExecuteProc	   = Draw2D;
    theEXW->EXT_PostProcessProc	= NULL; /* PostProcessFAMGGraph; */

    if ((thePOT=GetPlotObjType("FAMGGraph"))    == NULL) 
		return (1);
    thePOT->Dimension = TYPE_2D;
    thePOT->SetPlotObjProc = SetFAMGGraph;
    thePOT->DispPlotObjProc = DisplayFAMGGraph;

    return(0);
}

