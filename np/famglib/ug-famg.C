/****************************************************************************/
/*																			*/
/* File:      ug-famg.C														*/
/*																			*/
/* Purpose:   ug - famg interface											*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   Aug 1997 begin, ug version 3.7								*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <strstream.h>

extern "C"
{
#include <math.h>    

/* ug library */
#include "gm.h"        /* for data structure               */
#include "ugm.h"       /* for CreateNewLevelAMG            */
#include "evm.h"       /* for data structure               */
#include "ugdevices.h" /* for UserWrite, PrintErrorMessage */ 
#include "np.h"        /* for CreateNumProc,VECDATA_DESC   */
#include "debug.h"
#include "ugstruct.h"
#include "iter.h"
#include "disctools.h"   // for AssembleDirichletBoundary

/* test */
#include "wpm.h"
#include "wop.h"
#include "connectuggrape.h"
#include "uginterface.h"

#ifdef USE_UG_DS
#include "amgtransfer.h"
#include "npscan.h"
#endif

// give C-linkage to these functions
INT FAMGRestrictDefect (NP_TRANSFER *theNP, INT level,
						   VECDATA_DESC *to, VECDATA_DESC *from, 
						   MATDATA_DESC *A, VEC_SCALAR damp,
						   INT *result);

INT FAMGInterpolateCorrection (NP_TRANSFER *theNP, INT level,
								  VECDATA_DESC *to, VECDATA_DESC *from, 
								  MATDATA_DESC *A, VEC_SCALAR damp,
								  INT *result);
INT InitFAMG (void);

} // extern "C"

#include "ug-famg.h"

#ifdef USE_UG_DS
#include "famg_ugalgebra.h"
#else
#include "famg_arrayalgebra.h"
#endif

#include "famg_uginterface.h"
#include "famg_misc.h"
#include "famg_heap.h"

#include "famg_grid.h" // nur fuer printm fuers debuggen. WEG!
#include "famg_sparse.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/
#define DISPLAY_NP_FORMAT_SE			"%-16.13s = %-7.4E\n"


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
REP_ERR_FILE;

static struct FAMG_Interface famg_interface;
static struct FAMGParameter_ug famg_parameter;

/* RCS_ID
$Header$
*/

/****************************************************************************/
/*                                                                          */
/* functions                                								*/
/*                                                                          */
/****************************************************************************/

extern "C" INT InitFAMG ();
INT InitFAMGGraph (void);

static void FAMGReadArgvParameter(INT argc, char **argv)
{
	if (ReadArgvINT("h",&(famg_parameter.heap),argc,argv))
        famg_parameter.heap = (int)1e+7;
	if (ReadArgvINT("n1",&(famg_parameter.n1),argc,argv))
		famg_parameter.n1 = 1;
	if (ReadArgvINT("n2",&(famg_parameter.n2),argc,argv))
		famg_parameter.n2 = 1;
	if (ReadArgvINT("g",&(famg_parameter.gamma),argc,argv))
		famg_parameter.gamma = 1;
	if (ReadArgvINT("cgn",&(famg_parameter.cgnodes),argc,argv))
		famg_parameter.cgnodes = 1;
	if (ReadArgvINT("cgl",&(famg_parameter.cglevels),argc,argv))
		famg_parameter.cglevels = 100;
}
	
static void FAMGReadStringParameter(void)
{
    char *str;
	
    famg_parameter.ilut = 1e+10;
    GetStringValueDouble(":famg:ilut",&(famg_parameter.ilut));

    famg_parameter.cgilut = 0.0;
    GetStringValueDouble(":famg:cgilut",&(famg_parameter.cgilut));

    famg_parameter.conloops = 0;
    GetStringValueInt(":famg:conloops",&(famg_parameter.conloops));

    famg_parameter.mincoarse = 0.8;
    GetStringValueDouble(":famg:mincoarse",&(famg_parameter.mincoarse));

    famg_parameter.type = 0;
    GetStringValueInt(":famg:type",&(famg_parameter.type));

    famg_parameter.stv = 0;
    GetStringValueInt(":famg:stv",&(famg_parameter.stv));

    famg_parameter.tol = 0.95;
    GetStringValueDouble(":famg:tol",&(famg_parameter.tol));

    famg_parameter.sigma = 0.45;
    GetStringValueDouble(":famg:sigma",&(famg_parameter.sigma));

    famg_parameter.omegar = 1.0;
    GetStringValueDouble(":famg:omegar",&(famg_parameter.omegar));

    famg_parameter.omegal = 1.0;
    GetStringValueDouble(":famg:omegal",&(famg_parameter.omegal));

    famg_parameter.error1 = 1e-6;
    GetStringValueDouble(":famg:error1",&(famg_parameter.error1));

    famg_parameter.error2 = 1.0;
    GetStringValueDouble(":famg:error2",&(famg_parameter.error2));

    famg_parameter.maxit = 100;
    GetStringValueInt(":famg:maxit",&(famg_parameter.maxit));

    famg_parameter.alimit = 1e-14;
    GetStringValueDouble(":famg:alimit",&(famg_parameter.alimit));

    famg_parameter.rlimit = 1e-10;
    GetStringValueDouble(":famg:rlimit",&(famg_parameter.rlimit));

    famg_parameter.divlimit= 10.0;
    GetStringValueDouble(":famg:divlimit",&(famg_parameter.divlimit));
    
    famg_parameter.reduction = 1.0;
    GetStringValueDouble(":famg:reduction",&(famg_parameter.reduction));

    strcpy(famg_parameter.solver,"linit");
    str = GetStringVar(":famg:solver");
    if(str != NULL) strcpy(famg_parameter.solver,str);

    strcpy(famg_parameter.presmoother,"fgs");
    str = GetStringVar(":famg:presmoother");
    if(str != NULL) strcpy(famg_parameter.presmoother,str);

    strcpy(famg_parameter.postsmoother,"bgs");
    str = GetStringVar(":famg:postsmoother");
    if(str != NULL) strcpy(famg_parameter.postsmoother,str);

    strcpy(famg_parameter.cgsmoother,"ilut");
    str = GetStringVar(":famg:cgsmoother");
    if(str != NULL) strcpy(famg_parameter.cgsmoother,str);

    famg_parameter.coloringmethod = 3;
    GetStringValueInt(":famg:coloringmethod",&(famg_parameter.coloringmethod));

}


static INT FAMGPreProcess  (MULTIGRID *mg, INT *mark_key, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
							INT *result)
{
    GRID *grid, *amggrid;
    VECTOR *vec, *newvec, *w, *fw;
	MATRIX *m, *newmat;
	NODE *node;
	INT nrVec = 0, nrLinks = 0, lev, ll, found, i;
	SHORT mc,xmask,bmask;	
	// WEG SHORT xc,bc;	
    
	#ifdef ModelP
	assert(0); // not for parallel
	#endif
	
    MarkTmpMem(MGHEAP(mg),mark_key); /* release in PostProcess */
	
    if (MD_IS_SCALAR(A) && VD_IS_SCALAR(x) && VD_IS_SCALAR(b))
	{
		// WEG xc    = VD_SCALCMP(x);
		mc    = MD_SCALCMP(A);
		// WEG bc    = VD_SCALCMP(b);
		xmask  = VD_SCALTYPEMASK(x);
		bmask  = VD_SCALTYPEMASK(b);
    }
    else
    {
        UserWrite("Not a scalar equation. \n");
        REP_ERR_RETURN(1); 
    }

	amggrid = CreateNewLevelAMG(mg);
	if( amggrid == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create coarser grid " << endl;
		FAMGError(ostr);
		assert(0);
	}
	assert(GLEVEL(amggrid)==-1);	


    // create vectors on first AMG level
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
            {
				if (CreateVector(amggrid,VOTYPE(vec),VOBJECT(vec),&newvec))
			    {
					ostrstream ostr;
					ostr << __FILE__ << __LINE__  << "can not create vector on algebraic level" << endl;
					FAMGError(ostr);
					assert(0);
				}
				VINDEX(newvec) = nrVec++;
				SETVCLASS(newvec,3);
				SETVNCLASS(newvec,VCLASS(vec));
				SETNEW_DEFECT(newvec,1);
				SETFINE_GRID_DOF(newvec,0);
				SETPRIO(newvec,PRIO(vec));
				VECSKIP(newvec)=VECSKIP(vec);
				SETVCCOARSE(newvec,0);
#ifdef DYNAMIC_MEMORY_ALLOCMODEL
				VSTART(newvec) = NULL;
				VISTART(newvec) = NULL;
#endif

				// store a temporary link from the geometric-level vector to its amg-level copy
				// use therefore the geom_object field
				vec->object = (union geom_object *)newvec;
            }
        }
    }

    grid =  GRID_ON_LEVEL(mg,level);
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
        {
			if (CreateVector(amggrid,VOTYPE(vec),VOBJECT(vec),&newvec))
		    {
				ostrstream ostr;
				ostr << __FILE__ << __LINE__  << "can not create vector on algebraic level" << endl;
				FAMGError(ostr);
				assert(0);
			}
			VINDEX(newvec) = nrVec++;
			SETVCLASS(newvec,3);
			SETVNCLASS(newvec,VCLASS(vec));
			SETNEW_DEFECT(newvec,1);
			SETFINE_GRID_DOF(newvec,0);
			SETPRIO(newvec,PRIO(vec));
			VECSKIP(newvec)=VECSKIP(vec);
			SETVCCOARSE(newvec,0);
#ifdef DYNAMIC_MEMORY_ALLOCMODEL
			VSTART(newvec) = NULL;
			VISTART(newvec) = NULL;
#endif

			// store a temporary link from the geometric-level vector to its amg-level copy
			// use therefore the geom_object field
			vec->object = (union geom_object *)newvec;
        }
    }

    // copy matrix entries to first AMG level
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
            {
				assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
				newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(vec->object));
				if( newmat == NULL )
				{
			        ostrstream ostr;
			   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
			       	FAMGError(ostr);
					assert(0);
			   	}
				m = VSTART(vec);
				MVALUE(newmat,mc) = MVALUE(m,mc);
				MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
				nrLinks++;
                for (m=MNEXT(m); m!=NULL; m=MNEXT(m))
                {
                    w = MDEST(m);
                    found = 0;
                    if(FINE_GRID_DOF(w))
                    {
                        if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                        {
							assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
							assert(OBJT((VECTOR*)(w->object))==VEOBJ);
							newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(w->object));
							if( newmat == NULL )
							{
						        ostrstream ostr;
						   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
						       	FAMGError(ostr);
								assert(0);
						   	}
							MVALUE(newmat,mc) = MVALUE(m,mc);
							MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
							nrLinks++;
                        }
                        found = 1;
                    }
                    if(!found)
                    {
                        node = VMYNODE(w);
                        while(CORNERTYPE(node))
                        {
                            node = (NODE *)NFATHER(node);
                            fw = NVECTOR(node);
                            if(FINE_GRID_DOF(fw))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
									assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
									assert(OBJT((VECTOR*)(w->object))==VEOBJ);
									newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(w->object));
									if( newmat == NULL )
									{
								        ostrstream ostr;
						   	    		ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
								       	FAMGError(ostr);
										assert(0);
								   	}
									MVALUE(newmat,mc) = MVALUE(m,mc);
									MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
									nrLinks++;
                                }
                                found = 1;
                                break;
                            }
                        }
                    }   
                    if(!found)
                    {
                        node = (NODE *)SONNODE(VMYNODE(w));
                        ll = lev; ll++;
                        while(node != NULL)
                        {
                            fw = NVECTOR(node);
                            if((FINE_GRID_DOF(fw) && (ll < level))
                               || (NEW_DEFECT(fw) && (ll == level)))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
									assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
									assert(OBJT((VECTOR*)(w->object))==VEOBJ);
									newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(w->object));
									if( newmat == NULL )
									{
								        ostrstream ostr;
						   	    		ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
								       	FAMGError(ostr);
										assert(0);
								   	}
									MVALUE(newmat,mc) = MVALUE(m,mc);
									MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
									nrLinks++;
                                }
                                found = 1;
                                break;
                            }
                            node = (NODE *)SONNODE(node);
                        }
                    }
                    if(!found)
                    {
                        UserWrite("error in FAMGSolve \n");
                    }
                }
            }
        }
    }

    grid =  GRID_ON_LEVEL(mg,level);
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
        {
			newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(vec->object));
			if( newmat == NULL )
			{
		        ostrstream ostr;
		   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
		       	FAMGError(ostr);
				assert(0);
		   	}
			m = VSTART(vec);
			MVALUE(newmat,mc) = MVALUE(m,mc);
			MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
			nrLinks++;
            for (m=MNEXT(m); m!=NULL; m=MNEXT(m))
            {
                w = MDEST(m);
                found = 0;
                if(NEW_DEFECT(w))
                {
                    if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                    {
						assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
						assert(OBJT((VECTOR*)(w->object))==VEOBJ);
						newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(w->object));
						if( newmat == NULL )
						{
					        ostrstream ostr;
					   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
					       	FAMGError(ostr);
							assert(0);
					   	}
						MVALUE(newmat,mc) = MVALUE(m,mc);
						MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
						nrLinks++;
                    }
                    found = 1;
                }
                if(!found)
                {
                    node = VMYNODE(w);
                    while(CORNERTYPE(node))
                    {
                        node = (NODE *)NFATHER(node);
                        fw = NVECTOR(node);
                        if(FINE_GRID_DOF(fw))
                        {
                            if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                            {
								assert(OBJT((VECTOR*)(vec->object))==VEOBJ);
								assert(OBJT((VECTOR*)(w->object))==VEOBJ);
								newmat = CreateConnection(amggrid,(VECTOR*)(vec->object), (VECTOR*)(w->object));
								if( newmat == NULL )
								{
							        ostrstream ostr;
							   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
						    	   	FAMGError(ostr);
									assert(0);
							   	}
								MVALUE(newmat,mc) = MVALUE(m,mc);
								MVALUE(MADJ(newmat),mc) = MVALUE(MADJ(m),mc);
								nrLinks++;
                            }
                            found = 1;
                            break;
                        }
                    }
                }   
                if(!found)
                {
                    UserWrite("error in FAMGSolve \n");
                }                   
            }
        }
    }

    // restore the geom_object info in the vectors on geometric levels
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
				VOBJECT(vec) = VOBJECT((VECTOR *)(vec->object));

    grid =  GRID_ON_LEVEL(mg,level);
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
				VOBJECT(vec) = VOBJECT((VECTOR *)(vec->object));

	
	famg_interface.gridvector = new FAMGugGridVector(amggrid);
	if( famg_interface.gridvector == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create gridvector" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_DEFECT] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, b );
	if( famg_interface.vector[FAMG_DEFECT] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_DEFECT" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_UNKNOWN] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, x );
	if( famg_interface.vector[FAMG_UNKNOWN] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 0;
	}

	for(i = 0; i < FAMG_NVECTORS; i++)
	{
		switch(i)
		{
			case FAMG_DEFECT: 
				continue;	// alreday allocated
				
			case FAMG_UNKNOWN:
				continue;	// alreday allocated
				
			default:
				famg_interface.vector[i] = famg_interface.vector[FAMG_UNKNOWN]->create_new();
		}
		
		if( famg_interface.vector[i] == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector nr. " << i << endl;
			FAMGError(ostr);
			return 0;
		}
	}
	
	if( (famg_interface.matrix = new FAMGugMatrix( amggrid, A, nrVec, nrLinks )) == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
		FAMGError(ostr);
		return 0;
	}

	// init testvectors 
	*famg_interface.vector[FAMG_TVA] = 1.0;
	*famg_interface.vector[FAMG_TVB] = 1.0;

	FAMGConstructParameter(&famg_parameter);

    FAMGConstruct(famg_interface.gridvector, famg_interface.matrix, famg_interface.matrix, famg_interface.vector);
	
	result[0]=0;

	return NUM_OK;
}



#ifdef FAMG_SPARSE_BLOCK
static INT FAMGPreProcessForCoarseGridSolver  (MULTIGRID *mg, INT *mark_key, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, MATDATA_DESC *ACons, MATDATA_DESC *D, VECDATA_DESC *tv, VECDATA_DESC *tvT,
							INT *result)
// in this case, use grid on level 0 as the fine grid for FAMG and
// do not copy the start grid into level -1
{
	GRID *grid;
    VECTOR *vec;
	MATRIX *m;
	INT nrVec = 0, nrLinks = 0, i;

    grid = GRID_ON_LEVEL(mg,0);
	assert(grid!=NULL);


    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
		VINDEX(vec) = nrVec;
		nrVec++;
        // if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
        if(!VSKIPME(vec,0))
		{
			for( m=VSTART(vec); m!=NULL; m = MNEXT(m) )
				nrLinks++;
		}
		else
		{
			// check whether a dirichlet vector has exactly 1 matrix entry (the main diagonnal)
            
            // todo: check does not work for systems

			m = VSTART(vec);
			nrLinks++;
			// assert(fabs(MVALUE(m,mcd))>=1e-8);
			for( m=MNEXT(m); m!=NULL; m = MNEXT(m) )
			{
				nrLinks++;
				// assert(fabs(MVALUE(m,mc))<1e-8);
			}
		}
    }

	famg_interface.gridvector = new FAMGugGridVector(grid);
	if( famg_interface.gridvector == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create gridvector" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_DEFECT] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, b );
	if( famg_interface.vector[FAMG_DEFECT] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_DEFECT" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_UNKNOWN] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, x );
	if( famg_interface.vector[FAMG_UNKNOWN] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 0;
	}

    
    // sparse vector structure of the test vectors
    // todo: create this automatically from the script
    // test vectros gets the same structure as the off-diagonal of A
    SPARSE_MATRIX *spma = A->sm[MTP(0,0)];

    short ncomp = spma->nrows;
    short *compmap, cmpm, j, jj;
    compmap = new short[ncomp];
    short *diagoff = new short[spma->N];
    short ndiagoff = 0;

    for(i = 0; i < spma->nrows; i++)
    {
        for(jj = spma->row_start[i]; jj < spma->row_start[i+1]; jj++)
        {
            if (i == spma->col_ind[jj])
            {
                diagoff[i] = spma->offset[jj];
                ndiagoff++;
            }
        }
    }

    if(ncomp != ndiagoff)
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 1;
	}

    for(i = 0; i < ncomp; i++) compmap[i] = -1;
    cmpm = 0;
    for(i = 0; i < ncomp; i++)
    {
        if(compmap[i] >= 0) continue;
        compmap[i] = cmpm;
        for(j = i+1; j < ncomp; j++)
        {
            if(diagoff[i] == diagoff[j]) compmap[j] = cmpm;
        }
        cmpm++;
    }
        

	 famg_interface.vector[FAMG_TVA] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector,tv,compmap,ncomp);
	if( famg_interface.vector[FAMG_TVA] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 0;
	}

	 famg_interface.vector[FAMG_TVB] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector,tvT,compmap,ncomp);
	if( famg_interface.vector[FAMG_TVB] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 0;
	}

    delete diagoff;
    delete compmap;

    for(i = 0; i < FAMG_NVECTORS; i++)
	{
		switch(i)
		{
			case FAMG_DEFECT: 
				continue;	// alreday allocated
				
			case FAMG_UNKNOWN:
				continue;	// alreday allocated
				
			case FAMG_TVA:
				continue;	// alreday allocated
				
			case FAMG_TVB:
				continue;	// alreday allocated
				
			default:
				famg_interface.vector[i] = famg_interface.vector[FAMG_UNKNOWN]->create_new();
		}
		
		if( famg_interface.vector[i] == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector nr. " << i << endl;
			FAMGError(ostr);
			return 0;
		}
	}
	
	if( (famg_interface.matrix = new FAMGugMatrix( grid, A, nrVec, nrLinks )) == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
		FAMGError(ostr);
		return 0;
	}

	if( A!=ACons )
	{
		if( (famg_interface.Consmatrix = new FAMGugMatrix( grid, ACons, nrVec, nrLinks )) == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
			FAMGError(ostr);
			return 0;
		}
	}
	else
		famg_interface.Consmatrix = famg_interface.matrix;

    // todo: I don't know how to alloc temporarily a sparse matrix.
    //       So D has to be created explicitly in the script. 
    //       I do not know how to check this.

    if(D != NULL)
    {
		if( (famg_interface.diagmatrix = new FAMGugMatrix( grid, D, nrVec, nrLinks)) == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
			FAMGError(ostr);
			return 0;
		}
    }

        
	
	// init testvectors
    // test
    double *val = new double[3];
    double *valT = new double[3];
    val[0] = 1.0; val[1] = -0.5; val[2] = 1.0;
    valT[0] = -1.5; valT[1] = 1.0; valT[2] = 1.0;
 	SetValueSkip(*famg_interface.vector[FAMG_TVA],val);
 	SetValueSkip(*famg_interface.vector[FAMG_TVB],valT);
    delete val;
    // 	SetValueSkip(*famg_interface.vector[FAMG_TVA],1.0);
 	// SetValueSkip(*famg_interface.vector[FAMG_TVB],1.0);

	FAMGConstructParameter(&famg_parameter);

    FAMGConstruct(famg_interface.gridvector, famg_interface.matrix, famg_interface.Consmatrix, famg_interface.diagmatrix, famg_interface.vector);

	result[0]=0;

	return NUM_OK;
}
#else
static INT FAMGPreProcessForCoarseGridSolver  (MULTIGRID *mg, INT *mark_key, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, MATDATA_DESC *ACons,
							INT *result)
// in this case, use grid on level 0 as the fine grid for FAMG and
// do not copy the start grid into level -1
{
	GRID *grid;
    VECTOR *vec;
	MATRIX *m;
	INT nrVec = 0, nrLinks = 0, i;
	SHORT mc,mcd,bmask;	
	// WEG SHORT xc,bc,xmask;	
    
    MarkTmpMem(MGHEAP(mg),mark_key); /* release in PostProcess */
	
	    
    if (MD_IS_SCALAR(A) && MD_IS_SCALAR(ACons) && VD_IS_SCALAR(x) && VD_IS_SCALAR(b))
	{
		// WEG xc    = VD_SCALCMP(x);
		mc    = MD_SCALCMP(A);
		// WEG bc    = VD_SCALCMP(b);
		// WEG xmask  = VD_SCALTYPEMASK(x);
		bmask  = VD_SCALTYPEMASK(b);
    }
    else
    {
        UserWrite("Not a scalar equation. \n");
        REP_ERR_RETURN(1);
    }
    mcd = mc;

    grid = GRID_ON_LEVEL(mg,0);
	assert(grid!=NULL);

	if (AssembleDirichletBoundary (grid,A,x,b))
       NP_RETURN(1,result[0]); 

    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
		VINDEX(vec) = nrVec;
		nrVec++;
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
		{
			for( m=VSTART(vec); m!=NULL; m = MNEXT(m) )
				nrLinks++;
		}
		else
		{
			// check whether a dirichlet vector has exactly 1 matrix entry (the main diagonnal)
			m = VSTART(vec);
			nrLinks++;
			assert(fabs(MVALUE(m,mcd))>=1e-8);
			for( m=MNEXT(m); m!=NULL; m = MNEXT(m) )
			{
				nrLinks++;
				assert(fabs(MVALUE(m,mc))<1e-8);
			}
		}
    }

	famg_interface.gridvector = new FAMGugGridVector(grid);
	if( famg_interface.gridvector == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create gridvector" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_DEFECT] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, b );
	if( famg_interface.vector[FAMG_DEFECT] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_DEFECT" << endl;
		FAMGError(ostr);
		return 0;
	}
	
	famg_interface.vector[FAMG_UNKNOWN] = new FAMGugVector( *(FAMGugGridVector*)famg_interface.gridvector, x );
	if( famg_interface.vector[FAMG_UNKNOWN] == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector FAMG_UNKNOWN" << endl;
		FAMGError(ostr);
		return 0;
	}

	for(i = 0; i < FAMG_NVECTORS; i++)
	{
		switch(i)
		{
			case FAMG_DEFECT: 
				continue;	// alreday allocated
				
			case FAMG_UNKNOWN:
				continue;	// alreday allocated
				
			default:
				famg_interface.vector[i] = famg_interface.vector[FAMG_UNKNOWN]->create_new();
		}
		
		if( famg_interface.vector[i] == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector nr. " << i << endl;
			FAMGError(ostr);
			return 0;
		}
	}
	
	if( (famg_interface.matrix = new FAMGugMatrix( grid, A, nrVec, nrLinks )) == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
		FAMGError(ostr);
		return 0;
	}

	if( A!=ACons )
	{
		if( (famg_interface.Consmatrix = new FAMGugMatrix( grid, ACons, nrVec, nrLinks )) == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create matrix" << endl;
			FAMGError(ostr);
			return 0;
		}
	}
	else
		famg_interface.Consmatrix = famg_interface.matrix;
	

	// init testvectors
	*famg_interface.vector[FAMG_TVA] = 1.0;
	*famg_interface.vector[FAMG_TVB] = 1.0;

 	// SetValueSkip(*famg_interface.vector[FAMG_TVA],1.0);  
 	// SetValueSkip(*famg_interface.vector[FAMG_TVB],1.0);

	FAMGConstructParameter(&famg_parameter);

    FAMGConstruct(famg_interface.gridvector, famg_interface.matrix, famg_interface.Consmatrix, famg_interface.vector);

	result[0]=0;

	return NUM_OK;
}
#endif

/////////////////////////////////////////////////////////////////////////////////
//
// famg as a complete solver
//
/////////////////////////////////////////////////////////////////////////////////

static INT FAMGIterInit (NP_BASE *theNP, INT argc, char **argv)
{
	NP_FAMG_ITER *np;
	
	np = (NP_FAMG_ITER *) theNP;

	FAMGReadArgvParameter(argc, argv);
	FAMGReadStringParameter();

    np->heap = famg_parameter.heap;
    np->n1 = famg_parameter.n1;
    np->n2 = famg_parameter.n2;
    np->gamma = famg_parameter.gamma;
    np->cgnodes = famg_parameter.cgnodes;
    np->cglevels = famg_parameter.cglevels;
    np->maxit = famg_parameter.maxit;
    np->alimit = famg_parameter.alimit;
    np->rlimit = famg_parameter.rlimit;
    np->divlimit = famg_parameter.divlimit;
    np->reduction = famg_parameter.reduction;

	return (NPIterInit(&np->iter,argc,argv));
}

static INT FAMGIterDisplay (NP_BASE *theNP)
{
	NP_FAMG_ITER *np;	

	np = (NP_FAMG_ITER *) theNP;

	NPIterDisplay(&np->iter);

	UserWrite("configuration parameters:\n");	
	UserWriteF(DISPLAY_NP_FORMAT_SI,"h",(int)np->heap);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->n1);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->n2);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"cgn",(int)np->cgnodes);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"cgl",(int)np->cglevels);
	UserWriteF(DISPLAY_NP_FORMAT_SI,"maxit",(int)np->maxit);
	UserWriteF(DISPLAY_NP_FORMAT_SE,"alimit",(double)np->alimit);
	UserWriteF(DISPLAY_NP_FORMAT_SE,"rlimit",(double)np->rlimit);
	UserWriteF(DISPLAY_NP_FORMAT_SE,"divlim",(double)np->divlimit);
	UserWriteF(DISPLAY_NP_FORMAT_SE,"red",(double)np->reduction);

	return (0);
}

#ifdef USE_UG_DS
static INT FAMGIterPreProcess  (NP_ITER *theNP, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
							INT *baselevel, INT *result)
{
	MULTIGRID *mg;
    NP_FAMG_ITER *np;

	np = (NP_FAMG_ITER *) theNP;
	mg = NP_MG(theNP);
	
	#ifdef ModelP
	assert(0);	// not implemented for parallel
	#endif
	
	return FAMGPreProcess( mg, &np->famg_mark_key, level, x, b, A, result);
}
#else
static INT FAMGIterPreProcess  (NP_ITER *theNP, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
							INT *baselevel, INT *result)
{
	MULTIGRID *mg;
    GRID *grid;
    NODE *node;
    VECTOR *vec,*w,*fw;
	MATRIX *m;
	SHORT xc,bc,mc,xmask,bmask;
	INT i,j,lev,found,ll, nnb, offset;
    DOUBLE d, sum;
	
    int n, nl, nv, *index, *start;
    double *entry, *vector[FAMG_NVECTORS];
    void **extra;
    NP_FAMG_ITER *np;	

	#ifdef ModelP
	assert(0);	// not implemented for parallel
	#endif
	
	np = (NP_FAMG_ITER *) theNP;
    
	mg = NP_MG(theNP);
    MarkTmpMem(MGHEAP(mg),&np->famg_mark_key); /* release in PostProcess */
	
	    
    if (MD_IS_SCALAR(A) && VD_IS_SCALAR(x) && VD_IS_SCALAR(b))
	{
		xc    = VD_SCALCMP(x);
		mc    = MD_SCALCMP(A);
		bc    = VD_SCALCMP(b);
		xmask  = VD_SCALTYPEMASK(x);
		bmask  = VD_SCALTYPEMASK(b);
    }
    else
    {
        UserWrite("Not a scalar equation. \n");
        REP_ERR_RETURN(1);
    }


    /* count unknowns */
    n = 0;
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (FINE_GRID_DOF(vec)))
            {
                VINDEX(vec) = n;
                n++;
            }
        }
    }
		
    grid =  GRID_ON_LEVEL(mg,level);  
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (NEW_DEFECT(vec)))
        {
            VINDEX(vec) = n;
            n++;
        }
    }
        
    
    /* ug node information */
    extra = (void **) GetTmpMem(MGHEAP(mg),n*sizeof(void *),np->famg_mark_key);
    if (extra == NULL)
    {
        UserWrite("FAMGCreateSystem: not enough memory. \n");
        REP_ERR_RETURN(1);
    }
    
    i = 0;
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (FINE_GRID_DOF(vec)))
            {
                extra[i] = (void *) MYVERTEX(VMYNODE(vec));
                i++;
            }
        }
    }
		
    grid =  GRID_ON_LEVEL(mg,level);  
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (NEW_DEFECT(vec)))
        {
            extra[i] = (void *) MYVERTEX(VMYNODE(vec));
            i++;
        }
    }

    /* allocate row/column start array */
    start = (int *) GetTmpMem(MGHEAP(mg),(n+1)*sizeof(int),np->famg_mark_key);
    if (start == NULL)
    {
        UserWrite("ug - famg: not enough memory. \n");
        REP_ERR_RETURN(1);
    }

    /* allocate vectors */
    for(j = 0; j < FAMG_NVECTORS; j++)
    {
        vector[j] = (DOUBLE *) GetTmpMem(MGHEAP(mg),n*sizeof(DOUBLE),np->famg_mark_key);
        if (vector[j] == NULL)
        {
            UserWrite("ug - famg: not enough memory. \n");
            REP_ERR_RETURN(1);
        }
    }

    /* copy UG system matrix  into interface data structure */
    /* first step: count links */
    i = 0;
    start[0] = 0;
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
            {
                nnb = 1;
                for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
                {
                    w = MDEST(m);
                    found = 0;
                    if(FINE_GRID_DOF(w))
                    {
                        if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                        {
                            nnb++;
                        }
                        found = 1;
                    }
                    if(!found)
                    {
                        node = VMYNODE(w);
                        while(CORNERTYPE(node))
                        {
                            node = (NODE *)NFATHER(node);
                            fw = NVECTOR(node);
                            if(FINE_GRID_DOF(fw))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
                                    nnb++;
                                }
                                found = 1;
                                break;
                            }
                        }
                    }   
                    if(!found)
                    {
                        node = (NODE *)SONNODE(VMYNODE(w));
                        ll = lev; ll++;
                        while(node != NULL)
                        {
                            fw = NVECTOR(node);
                            if((FINE_GRID_DOF(fw) && (ll < level))
                               || (NEW_DEFECT(fw) && (ll == level)))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
                                    nnb++;
                                }
                                found = 1;
                                break;
                            }
                            node = (NODE *)SONNODE(node);
                            ll++;
                        }
                    }
                    if(!found)
                    {
                        UserWrite("error in FAMGSolve \n");
                    }
                }
                i++;
                start[i] = nnb+start[i-1];
            }
        }
    }
    grid =  GRID_ON_LEVEL(mg,level);
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
        {
            nnb = 1;
            for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
            {
                w = MDEST(m);
                found = 0;
                if(NEW_DEFECT(w))
                {
                    if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                    {
                        nnb++;
                    }
                    found = 1;
                }
                if(!found)
                {
                    node = VMYNODE(w);
                    while(CORNERTYPE(node))
                    {
                        node = (NODE *)NFATHER(node);
                        fw = NVECTOR(node);
                        if(FINE_GRID_DOF(fw))
                        {
                            if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                            {
                                nnb++;
                            }
                            found = 1;
                            break;
                        }
                    }
                }   
                if(!found)
                {
                    UserWrite("error in FAMGSolve \n");
                }                   
            }
            i++;
            start[i] = nnb+start[i-1];
        }
    }

    if(i != n)
    {
        UserWrite("error in FAMGIterPreProcess. \n");
        REP_ERR_RETURN(1);
    }

    /* allocate index and matrix array */
    nl = start[n];
    index = (int *) GetTmpMem(MGHEAP(mg),nl*sizeof(int),np->famg_mark_key);
    if (index == NULL)
    {
        UserWrite("ug - famg: not enough memory. \n");
        REP_ERR_RETURN(1);
    }
    entry = (double *) GetTmpMem(MGHEAP(mg),nl*sizeof(double),np->famg_mark_key);
    if (entry == NULL)
    {
        UserWrite("ug - famg: not enough memory. \n");
        REP_ERR_RETURN(1);
    }

    /* second step: save matrix entries */
    i = 0; offset = 0;
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
            {
                if(offset != start[i])
                {
                    UserWrite("error in FAMGIterPreProcess. \n");
                    REP_ERR_RETURN(1);
                }
                entry[offset] = MVALUE(VSTART(vec),mc);
                index[offset] = i;
                offset++;
                for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
                {
                    w = MDEST(m);
                    found = 0;
                    if(FINE_GRID_DOF(w))
                    {
                        if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                        {
                            entry[offset] = MVALUE(m,mc);
                            index[offset] = VINDEX(w);
                            offset++;
                        }
                        found = 1;
                    }
                    if(!found)
                    {
                        node = VMYNODE(w);
                        while(CORNERTYPE(node))
                        {
                            node = (NODE *)NFATHER(node);
                            fw = NVECTOR(node);
                            if(FINE_GRID_DOF(fw))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
                                    entry[offset] = MVALUE(m,mc);
                                    index[offset] = VINDEX(fw);
                                    offset++;
                                }
                                found = 1;
                                break;
                            }
                        }
                    }   
                    if(!found)
                    {
                        node = (NODE *)SONNODE(VMYNODE(w));
                        ll = lev; ll++;
                        while(node != NULL)
                        {
                            fw = NVECTOR(node);
                            if((FINE_GRID_DOF(fw) && (ll < level))
                               || (NEW_DEFECT(fw) && (ll == level)))
                            {
                                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                                {
                                    entry[offset] = MVALUE(m,mc);
                                    index[offset] = VINDEX(fw);
                                    offset++;
                                }
                                found = 1;
                                break;
                            }
                            node = (NODE *)SONNODE(node);
                            ll++;
                        }
                    }
                    if(!found)
                    {
                        UserWrite("error in FAMGSolve \n");
                    }
                }
                i++;
            }
        }
    }

    grid =  GRID_ON_LEVEL(mg,level);
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
        {
            if(offset != start[i])
            {
                UserWrite("error in FAMGIterPreProcess. \n");
                REP_ERR_RETURN(1);
            }
            entry[offset] = MVALUE(VSTART(vec),mc);
            index[offset] = i;
            offset++;
            for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
            {
                w = MDEST(m);
                found = 0;
                if(NEW_DEFECT(w))
                {
                    if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
                    {
                        entry[offset] = MVALUE(m,mc);
                        index[offset] = VINDEX(w);
                        offset++;
                    }
                    found = 1;
                }
                if(!found)
                {
                    node = VMYNODE(w);
                    while(CORNERTYPE(node))
                    {
                        node = (NODE *)NFATHER(node);
                        fw = NVECTOR(node);
                        if(FINE_GRID_DOF(fw))
                        {
                            if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                            {
                                entry[offset] = MVALUE(m,mc);
                                index[offset] = VINDEX(fw);
                                offset++;
                            }
                            found = 1;
                            break;
                        }
                    }
                }   
                if(!found)
                {
                    UserWrite("error in FAMGSolve \n");
                }                   
            }
            i++;
        }
    }

    famg_interface.n = n;
    famg_interface.nl = nl;
    famg_interface.nv = FAMG_NVECTORS;
    famg_interface.start = start;
    famg_interface.index = index;
    famg_interface.entry = entry;
    famg_interface.extra = extra;
    for(j = 0; j < FAMG_NVECTORS; j++)
    {
        famg_interface.vector[j] = vector[j];
    }
        
	FAMGConstructParameter(&famg_parameter);

    FAMGConstruct((FAMG_Interface*)&famg_interface);

    return(0);
}
#endif

static INT FAMGIterSolve (NP_ITER *theNP, INT level,
				 VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
				 INT *result)
{
	MULTIGRID *mg;
    GRID *grid;
    NODE *snode;
    VECTOR *vec,*svec;
	SHORT cc,bc,cmask;
	INT i, n, lev;
    FAMGVector *unknown, *rhs, *tv, *tvT, *defect;

	mg = NP_MG(theNP);
    grid = GRID_ON_LEVEL(mg,level);

    if (MD_IS_SCALAR(A) && VD_IS_SCALAR(c) && VD_IS_SCALAR(b))
	{
		cc    = VD_SCALCMP(c);
		// WEG mc    = MD_SCALCMP(A);
		bc    = VD_SCALCMP(b);
		cmask  = VD_SCALTYPEMASK(c);
    }
    else
    {
        UserWrite("Not a scalar equation. \n");
        REP_ERR_RETURN(1);
    }
    
    unknown = famg_interface.vector[FAMG_UNKNOWN];
    rhs = famg_interface.vector[FAMG_RHS];
    defect = famg_interface.vector[FAMG_DEFECT];
    tv = famg_interface.vector[FAMG_TVA];
    tvT = famg_interface.vector[FAMG_TVB];
    
	FAMGVectorIter viter(*rhs);
	FAMGVectorEntry ve;
	n=0;
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&cmask) && (FINE_GRID_DOF(vec)))
            {
				viter(ve);
				(*rhs)[ve] = VVALUE(vec,bc);
				(*tv)[ve] = 1.0;
				(*tvT)[ve] = 1.0;
				n++;
            }
        }
    }
		
    grid =  GRID_ON_LEVEL(mg,level);  
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&cmask) && (NEW_DEFECT(vec)))
        {
			viter(ve);
			(*rhs)[ve] = VVALUE(vec,bc);
			(*tv)[ve] = 1.0;
			(*tvT)[ve] = 1.0;
			n++;
        }
    }
	
	viter(ve);
    if((*rhs).is_valid(ve))
    {
        UserWrite("error in FAMGIterSolve: now more entries than initialized\n");
        REP_ERR_RETURN(1);
    }
    UserWriteF("unknowns: %d \n",n);

    /* solve */
    FAMGSolveSystem(&famg_interface);
        
    i = 0;
	viter.reset();
    for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
    {    
        for (vec=PFIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
        {
            if((VDATATYPE(vec)&cmask) && (FINE_GRID_DOF(vec)))
            {
                if(!VSKIPME(vec,0))
                {
					viter(ve);
                    VVALUE(vec,cc) = (*unknown)[ve];
                    VVALUE(vec,bc) = (*defect)[ve];
                    i++;
                }
                else
                {
                    VVALUE(vec,cc) = 0.0;
                    VVALUE(vec,bc) = 0.0;
                }
                /* make solution consistent on higher levels */
                /* we are not responsible for consistency on lower levels */
                snode = (NODE *)SONNODE(VMYNODE(vec));
                while (snode != NULL)
                {
                    svec = NVECTOR(snode);
                    if (VDATATYPE(svec)&cmask)
                    {
                       VVALUE(svec,cc) = VVALUE(vec,cc); 
                       VVALUE(svec,bc) = VVALUE(vec,bc);
                    }
                    snode = (NODE *)SONNODE(snode);
                }
            }
        }
    }
        
		
    grid =  GRID_ON_LEVEL(mg,level);  
    for (vec=PFIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
    {
        if((VDATATYPE(vec)&cmask) && (NEW_DEFECT(vec)))
        {
            if(!VSKIPME(vec,0))
            {
				viter(ve);
                VVALUE(vec,cc) = (*unknown)[ve];
                VVALUE(vec,bc) = (*defect)[ve];
                i++;
            }
            else
            {
                VVALUE(vec,cc) = 0.0;
                VVALUE(vec,bc) = 0.0;
            }
        }
    }
    if(i != n)
    {
        UserWrite("error in FAMGIterSolve: vector length differ\n");
        REP_ERR_RETURN(1);
    }

               
       
    return(0);
}


static INT FAMGIterPostProcess (NP_ITER *theNP, INT level,
							VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
							INT *result)
{
	MULTIGRID *mg;
    NP_FAMG_ITER *np;	
	int i;
	
	np = (NP_FAMG_ITER *) theNP;

	if( famg_interface.Consmatrix != famg_interface.matrix )
		delete famg_interface.Consmatrix;
		
	delete famg_interface.matrix;
	
	for(i = 0; i < FAMG_NVECTORS; i++)
		delete famg_interface.vector[i];
	
	delete famg_interface.gridvector;

	FAMGDeconstruct();
	
	FAMGDeconstructParameter();
	
	mg = NP_MG(theNP);
	FAMGFreeHeap();
    ReleaseTmpMem(MGHEAP(mg),np->famg_mark_key); /* mark in PreProcess */

    return (0);
}

// famg as a complete solver 
static INT FAMGConstructIterNP (NP_BASE *theNP)
{
    NP_ITER *np;
	
	theNP->Init = FAMGIterInit;
	theNP->Display = FAMGIterDisplay;
	theNP->Execute = NPIterExecute;

	np = (NP_ITER *) theNP;
	np->PreProcess = FAMGIterPreProcess;
	np->Iter = FAMGIterSolve;
	np->PostProcess = FAMGIterPostProcess;

    return(0);
}

/////////////////////////////////////////////////////////////////////////////////
//
// famg as a pure transfer num proc
//
/////////////////////////////////////////////////////////////////////////////////

INT FAMGTransferInit (NP_BASE *theNP, INT argc, char **argv)
{
	NP_FAMG_TRANSFER *famgtrans = (NP_FAMG_TRANSFER *)theNP;
	
	FAMGReadArgvParameter(argc, argv);
	FAMGReadStringParameter();

	if (ReadArgvINT("coarsegridsolver",&(famgtrans->coarsegridsolver),argc,argv))
        famgtrans->coarsegridsolver = 1;
	
	if (ReadArgvINT("coarsegridagglo",&(famgtrans->coarsegridagglo),argc,argv))
        famgtrans->coarsegridagglo = 0;
	
#ifdef FAMG_SPARSE_BLOCK
	famgtrans->D = ReadArgvMatDesc(theNP->mg, "D", argc, argv);
	famgtrans->tv = ReadArgvVecDesc(theNP->mg, "tv", argc, argv);
	famgtrans->tvT = ReadArgvVecDesc(theNP->mg, "tvT", argc, argv);
#endif

	// famgtrans->ConsMat = ReadArgvMatDesc(famgtrans->amg_trans.transfer.base.mg,"ConsMat",argc,argv);

	famgtrans->smooth_sol = NULL;	// default to detect errors
	famgtrans->smooth_def = NULL;	// default to detect errors
	
	//return AMGTransferInit (theNP, argc, argv); is not good because we can't provide the necessary parameters
	return NP_EXECUTABLE;
}


INT FAMGTransferPreProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
						   VECDATA_DESC *x, VECDATA_DESC *b, 
						   MATDATA_DESC *A, INT *result)
{
	MULTIGRID *mg;
    NP_FAMG_TRANSFER *np;
	INT res;
	MATDATA_DESC *ACons;
	double StartTimeTotal = CURRENT_TIME_LONG;
	
	np = (NP_FAMG_TRANSFER *) theNP;
	mg = NP_MG(theNP);
	
	if( GRID_ON_LEVEL(mg,-1) != NULL )		// remove AMG grids if not done
		if( DisposeAMGLevels(mg) )
			NP_RETURN(1,result[0]);
	
	FAMGReadStringParameter();	// reread to be able to configure each solver call individually
	
#ifdef ModelP
	// check assumptions for IS_FAMG_MASTER and IS_FAMG_GHOST
	assert(PrioMaster>PrioBorder);
	assert(PrioHGhost<PrioBorder);
	assert(PrioVGhost<PrioBorder);
	assert(PrioVHGhost<PrioBorder);
	
	// a consistent copy of the stiffmat is needed
	// in order to be able to construct a consistent interpolationmatrix
	
	ACons = np->ConsMat;
	if (ACons==NULL)
	{
		if( AllocMDFromMD(NP_MG(theNP),tl,tl,A,&ACons) || (ACons==NULL) )	
		{
			ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not read ConsMat symbol" << endl;
			FAMGError(ostr);
			assert(0);
		}
		np->ConsMat = ACons;
		np->ConsMatTempAllocated = 1;
	}
	else 
		np->ConsMatTempAllocated = 0;
#else
	ACons = A;
#endif
	
	if( np->coarsegridsolver )
    {
#ifdef FAMG_SPARSE_BLOCK
		res = FAMGPreProcessForCoarseGridSolver( mg, &np->famg_mark_key, *fl, x, b, A, ACons, np->D, np->tv, np->tvT, result);
#else
		res = FAMGPreProcessForCoarseGridSolver( mg, &np->famg_mark_key, *fl, x, b, A, ACons, result);
#endif
    }
	else
	{
		#ifdef ModelP
		assert(0);	// not implemented for parallel
		#endif
		res = FAMGPreProcess( mg, &np->famg_mark_key, *fl, x, b, A, result);
	}

	/* we set the baselevel for the following cycle!! */
	*fl = mg->bottomLevel;

#ifdef ModelP
	// coarse grid agglomeration
//prm(mg->bottomLevel,0);
	if( np->coarsegridagglo )
	{
		AMGAgglomerate(mg);
		l_amgmatrix_collect(GRID_ON_LEVEL(mg,mg->bottomLevel),A);
		UserWrite("coarse grid agglomerated\n");
		printf("%d: coarse grid agglomerated\n", me);
	}
//prm(mg->bottomLevel,0);
#endif
	
	cout << me << ": total time for constructing FAMG Transfer = " << CURRENT_TIME_LONG - StartTimeTotal << endl;
	
	return res;
}

static INT FAMGTransferPostProcess (NP_TRANSFER *theNP, INT *fl, INT tl,
								   VECDATA_DESC *x, VECDATA_DESC *b, 
								   MATDATA_DESC *A, INT *result)
{
	MULTIGRID *theMG;
	NP_FAMG_TRANSFER *np;
	INT i;
	
	result[0]=0;
	np = (NP_FAMG_TRANSFER *) theNP;
	theMG = NP_MG(theNP);

	FAMGDeconstruct();
	
	FAMGDeconstructParameter();
	
	if( famg_interface.Consmatrix != famg_interface.matrix )
		delete famg_interface.Consmatrix;
		
	delete famg_interface.matrix;
	
	for(i = 0; i < FAMG_NVECTORS; i++)
		delete famg_interface.vector[i];	// free temp allocated symbols
	
	delete famg_interface.gridvector;
	
	ReleaseTmpMem(MGHEAP(theMG),np->famg_mark_key); /* mark in PreProcess */
	// all memory released? TODO
	
	if (np->ConsMatTempAllocated)
	{
		if( FreeMD(theMG,*fl,tl,np->ConsMat))
			NP_RETURN(1,result[0]);
		np->ConsMat = NULL;
		np->ConsMatTempAllocated = 0;
	}
	
	return 0;
}


/* actions for FAMG defect restriction:
UG	t := 0				new [ug/np/procs/iter.c/Lmgc()]
	t += M^{-1} * b		fine grid Jacobi smoothing [FAMGGrid::Restriction]
						t carries the solution update back to ug
	d -= K * t			update defect [FAMGGrid::Restriction]
						anticipates the following c += t in ug
	d_{l+1} = R * d_l	restrict defect [FAMGGrid::Restriction]
UG	c += t				new [ug/np/procs/iter.c/Lmgc()]
*/
INT FAMGRestrictDefect (NP_TRANSFER *theNP, INT level,
						   VECDATA_DESC *to, VECDATA_DESC *from, 
						   MATDATA_DESC *A, VEC_SCALAR damp,
						   INT *result)
{
	INT famglevel;
    NP_FAMG_TRANSFER *np;

	np = (NP_FAMG_TRANSFER *) theNP;
	if( np->coarsegridsolver )
		famglevel = -level;
	else
		famglevel = -1-level;
	
#ifdef ModelP
	// TODO: sollte eigentlich ueberfluessig sein: defect sollte auf border vectoren eh schon 0 sein!
	if (l_vector_collect(GRID_ON_LEVEL(np->amg_trans.transfer.base.mg,level),from)!=NUM_OK) 
		NP_RETURN(1,result[0]);
#endif
	
	result[0] = FAMG_RestrictDefect(famglevel, to, from, np->smooth_sol, np->smooth_def);
	return result[0];
}


/* actions for FAMG correction interpolation
UG	t := 0				new [ug/np/procs/iter.c/Lmgc()]
	t_l += P * c_{l+1}	prolong correction [FAMGGrid::Restriction]
	d -= K * t			update defect [FAMGGrid::Restriction]
	c += t				new, only for transfer [FAMGGrid::Restriction]
	t := 0				prepare for the followign smoothing [FAMGGrid::Restriction]
	t += M^{-1} * b		fine grid Jacobi smoothing [FAMGGrid::Restriction]
						t carries the solution update back to ug
UG	c += t				[ug/np/procs/iter.c/Lmgc()]
UG	d -= K * t			update defect [ug/np/procs/iter.c/Lmgc()]; not in [FAMGGrid::Restriction]
*/
INT FAMGInterpolateCorrection (NP_TRANSFER *theNP, INT level,
								  VECDATA_DESC *to, VECDATA_DESC *from, 
								  MATDATA_DESC *A, VEC_SCALAR damp,
								  INT *result)
{
	INT famglevel;
    NP_FAMG_TRANSFER *np;

	np = (NP_FAMG_TRANSFER *) theNP;	
	if( np->coarsegridsolver )
	{
		famglevel = -level;
	}
	else
		famglevel = -1-level;
	
	result[0] = FAMG_ProlongCorrection(famglevel, to, from, np->smooth_sol, np->smooth_def);
	return result[0];
}


static INT FAMGConstructTransferNP (NP_BASE *theNP)
{
    NP_TRANSFER *np;
    NP_AMG_TRANSFER *amg_np;

	AMGTransferConstruct(theNP);

	amg_np = (NP_AMG_TRANSFER *)theNP;
	amg_np->AMGtype = FAMG;
	amg_np->Coarsen = NULL;
	amg_np->SetupIR = NULL;
	
	// change default settings
	theNP->Init = FAMGTransferInit;
	np = (NP_TRANSFER *)theNP;
	np->PreProcess = FAMGTransferPreProcess;
	np->PostProcess = FAMGTransferPostProcess;
	np->RestrictDefect = FAMGRestrictDefect; 
	np->InterpolateCorrection = FAMGInterpolateCorrection;

    return(0);
}


INT InitFAMG (void)
{	
	if (CreateClass(ITER_CLASS_NAME ".famg",sizeof(NP_FAMG_ITER),FAMGConstructIterNP))
		REP_ERR_RETURN (__LINE__);

	if (CreateClass(TRANSFER_CLASS_NAME ".famgTransfer",sizeof(NP_FAMG_TRANSFER),FAMGConstructTransferNP))
		REP_ERR_RETURN (__LINE__);

	if (InitFAMGGraph())
		REP_ERR_RETURN (__LINE__);

	return(0);
}
    
