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
#include "npscan.h"

#ifdef ModelP
#include "parallel.h"
#include "pargm.h"
#endif

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
static INT LineWidth;
static INT VecCoordComp;
static long ColorTable[8];

#ifdef ModelP
static int DrawBorderVec;
static long GreenColor; /* Green for non master vectors */
static VECDATA_DESC *ConsVector;

static int Gather_CoordVectorComp (DDD_OBJ obj, void *data)
{
	DOUBLE *vecdata = VVALUEPTR((VECTOR *)obj,VecCoordComp);
								
	V_DIM_COPY(vecdata,(DOUBLE*)data);	// set data[DIM] := vec_coord
		
	return NUM_OK;
}
 
static int Scatter_CoordVectorComp (DDD_OBJ obj, void *data)
{
	DOUBLE *vecdata = VVALUEPTR((VECTOR *)obj,VecCoordComp);
								
	V_DIM_COPY((DOUBLE*)data, vecdata);	// set vec_coord := data[DIM]

	return NUM_OK;
}

INT l_coord_project (GRID *g, const VECDATA_DESC *x)
{
	int n;
	
	if( g==NULL )
		return NUM_OK;
	
	ConsVector = (VECDATA_DESC *)x;

	DDD_IFAOneway(BorderVectorIF, GRID_ATTR(g), IF_BACKWARD, DIM * sizeof(DOUBLE),
				  Gather_CoordVectorComp, Scatter_CoordVectorComp);

	DDD_IFAOneway(VectorVIF, GRID_ATTR(g), IF_BACKWARD, DIM * sizeof(DOUBLE),
				  Gather_CoordVectorComp, Scatter_CoordVectorComp);

	return NUM_OK;
}
#endif


struct FAMGPlotObject
{
	struct PlotObjHead theHead;    /* the head */
	int level;
	int ids;
	int LineWidth;					// LineWidth used for plotting: for X11 uses 1 (default), for postscript uses 2
	VECDATA_DESC *CoordVec;			// vector to hold the coordinate info for the node
#ifdef ModelP
	int DrawBorder;					// if 1, also border vectors are plotted
#endif
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
   
	if (ReadArgvINT("linewidth",&i,argc,argv))
		theObj->LineWidth = 1;
	else
		theObj->LineWidth = i;
	
#ifdef ModelP
    if (ReadArgvINT("b",&i,argc,argv))
    	theObj->DrawBorder = 0;
	else
    	theObj->DrawBorder = i;
#endif
	
    theObj->CoordVec = ReadArgvVecDesc(PO_MG(thePlotObj),"coordvec",argc,argv);
	if( theObj->CoordVec != NULL )
		if (VD_ncmps_in_otype(theObj->CoordVec,NODEVEC) < DIM)
		{
			PrintErrorMessage('E',"plot FAMG Graph","coordinate vector has too few components");
			return(NOT_ACTIVE);
		}

	return (ACTIVE);
}
    
static INT  DisplayFAMGGraph (PLOTOBJ *thePlotObj)
{
    struct FAMGPlotObject *theObj;
    
	theObj = (struct FAMGPlotObject *) &(thePlotObj->theExternObject);

    UserWriteF("level: %d \n",theObj->level);
	
	if (theObj->CoordVec != NULL) 
		UserWriteF(DISPLAY_NP_FORMAT_SS,"CoordVec",ENVITEM_NAME(theObj->CoordVec));
	
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
#ifdef ModelP
	GreenColor = theOD->green;
	DrawBorderVec = theObj->DrawBorder;
#endif
	
	ColorTable[0] = theOD->blue;
	ColorTable[1] = theOD->magenta;
	ColorTable[2] = theOD->green;
	ColorTable[3] = theOD->orange;
	ColorTable[4] = theOD->cyan;
	ColorTable[5] = theOD->gray;
	ColorTable[6] = theOD->yellow;
	ColorTable[7] = theOD->red;
	
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
    LineWidth = theObj->LineWidth;

	VecCoordComp = -1;		// dummy to indicate "no vector for coordinates"
    if( (theObj->CoordVec != NULL) && ( level < 0 ) ) 
    {
		INT i, n;
	    VERTEX *vertex;
		NODE *node;
		DOUBLE *vertex_coord, *vector_coord, *vcoarse_coord;
		MATRIX *im;
		VECTOR *vec;
		
		// copy coordinate info into coordinate vector
	
		// initialize coordinate vector with the vertex-coordinate info
		grid = GRID_ON_LEVEL(mg,0);
		VecCoordComp = VD_ncmp_cmpptr_of_otype(theObj->CoordVec,NODEVEC,&n)[0];
	    assert(n>=DIM);
		for( node=PFIRSTNODE(grid); node!=NULL; node=SUCCN(node) )
		{
		    vertex_coord = CVECT(MYVERTEX(node));
			vector_coord = VVALUEPTR(NVECTOR(node),VecCoordComp);
			V_DIM_COPY(vertex_coord,vector_coord);	// set vector_coord := vertex_coord
		}
		
		// now propagate the coordinate info to the algebraic levels
		for( i=0; i>=level; i-- )
		{
			grid = GRID_ON_LEVEL(mg,i);
			for( vec=PFIRSTVECTOR(grid); vec!=NULL; vec=SUCCVC(vec) )
			{
				im = VISTART(vec);
				if( im!=NULL && MNEXT(im)==NULL )
				{	// this vector has exactly 1 interpolation matrix entry; thus it is a coarse grid vector and its coord value must be restricted to its coarse grid instance
					vector_coord = VVALUEPTR(vec,VecCoordComp);
					vcoarse_coord = VVALUEPTR(MDEST(im),VecCoordComp);
					V_DIM_COPY(vector_coord,vcoarse_coord);	// set vcoarse_coord := vector_coord
				}
			}
			
			#ifdef ModelP
			// communicate the coord info to copies
			l_coord_project (DOWNGRID(grid), theObj->CoordVec);
			#endif
		}
	}
	
    return(0);
}


static INT EvalFAMGGraph1 (DRAWINGOBJ *theDO, VECTOR *vec)
{
    VERTEX *vertex, *nbvertex;
    DOUBLE_VECTOR mypos,nbpos;
    int j, CircleSize;
    VECTOR *nbvec;
    MATRIX *mat, *imat;
	long VectorColor;

    if(vec == NULL)
		goto EvalFAMGGraph1_finish;

    UgSetLineWidth(LineWidth);

#ifdef ModelP
	if( !IS_FAMG_MASTER(vec) )
	{
		if( !DrawBorderVec )
			goto EvalFAMGGraph1_finish;
		VectorColor = GreenColor;
	}
	else
	{
		if( DrawBorderVec )
			VectorColor = BlackColor;
		else
			VectorColor = ColorTable[me%8];
	}
	CircleSize = 12;
#else
	VectorColor = BlackColor;
	CircleSize = 8;
#endif
	
	if( VecCoordComp == -1 )
	{
		vertex = MYVERTEX(VMYNODE(vec));				// take coord from vertex
		V_DIM_COPY(CVECT(vertex),mypos);
	}
	else
	{
		V_DIM_COPY(VVALUEPTR(vec,VecCoordComp),mypos);	// take coord from special vector
	}

    /* plot  marker */
    if(VCCOARSE(vec))
    {
        DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO); 
        DO_2c(theDO) = 1; DO_inc(theDO); 
        DO_2l(theDO) = VectorColor; DO_inc(theDO);
        DO_2s(theDO) = FILLED_CIRCLE_MARKER; DO_inc(theDO); 
        DO_2s(theDO) = CircleSize; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }
    else
    {
        DO_2c(theDO) = DO_POLYMARK; DO_inc(theDO); 
        DO_2c(theDO) = 1; DO_inc(theDO); 
        DO_2l(theDO) = VectorColor; DO_inc(theDO);
        DO_2s(theDO) = EMPTY_CIRCLE_MARKER; DO_inc(theDO); 
        DO_2s(theDO) = CircleSize; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }
    

    if(GlobalIds)
    {
        /* print id */
        DO_2c(theDO) = DO_TEXT; DO_inc(theDO);
        DO_2l(theDO) = VectorColor; DO_inc(theDO);
        DO_2c(theDO) = TEXT_REGULAR; DO_inc(theDO); 
        DO_2c(theDO) = TEXT_NOT_CENTERED; DO_inc(theDO); 
        DO_2s(theDO) = CircleSize + 2; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        sprintf(DO_2cp(theDO),"%d",VINDEX(vec));
        DO_inc_str(theDO);
    }

    /* plot matrix */
    for (mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
    {
        nbvec = MDEST(mat);
		if( VecCoordComp == -1 )
		{
			nbvertex = MYVERTEX(VMYNODE(nbvec));
			V_DIM_COPY(CVECT(nbvertex),nbpos);					// take coord from vertex
		}
		else
		{
			V_DIM_COPY(VVALUEPTR(nbvec,VecCoordComp),nbpos);	// take coord from special vector
		}
        DO_2c(theDO) = DO_LINE; DO_inc(theDO); 
        DO_2l(theDO) = BlackColor; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
    }

EvalFAMGGraph1_finish:	
	#ifdef ModelP
	WOP_DObjPnt = theDO;
	#endif
            
    return(0);
}

static INT EvalFAMGGraph2 (DRAWINGOBJ *theDO, VECTOR *vec)
{
    VERTEX *vertex, *nbvertex;
    DOUBLE_VECTOR mypos,nbpos;
    int j;
    VECTOR *nbvec;
    MATRIX *mat, *imat;

    if(vec == NULL)
		goto EvalFAMGGraph2_finish;

    UgSetLineWidth(LineWidth);

	if( VecCoordComp == -1 )
	{
		vertex = MYVERTEX(VMYNODE(vec));				// take coord from vertex
		V_DIM_COPY(CVECT(vertex),mypos);
	}
	else
	{
		V_DIM_COPY(VVALUEPTR(vec,VecCoordComp),mypos);	// take coord from special vector
	}

    /* plot transfer */
    for (imat=VISTART(vec); imat!=NULL; imat=MNEXT(imat))
    {
        nbvec = MDEST(imat);
		if( VecCoordComp == -1 )
		{
			nbvertex = MYVERTEX(VMYNODE(nbvec));
			V_DIM_COPY(CVECT(nbvertex),nbpos);					// take coord from vertex
		}
		else
		{
			V_DIM_COPY(VVALUEPTR(nbvec,VecCoordComp),nbpos);	// take coord from special vector
		}
        DO_2c(theDO) = DO_LINE; DO_inc(theDO); 
        DO_2l(theDO) = RedColor; DO_inc(theDO);
        V2_COPY(mypos,DO_2Cp(theDO)); DO_inc_n(theDO,2);
        V2_COPY(nbpos,DO_2Cp(theDO)); DO_inc_n(theDO,2); 
     }

EvalFAMGGraph2_finish:	
	#ifdef ModelP
	WOP_DObjPnt = theDO;
	#endif
            
    return(0);
    
}


static INT EvalFAMGGraph (DRAWINGOBJ *theDO, INT *end)
{

    if(GlobalVec2 == NULL)
	{
		*end = 1;
		
		#ifdef ModelP
		WOP_DObjPnt = theDO;
		#endif
	
		return(0);
	}
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

