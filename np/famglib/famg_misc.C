/****************************************************************************/
/*																			*/
/* File:      famg_misc.C													*/
/*																			*/
/* Purpose:   misc functions												*/
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

#include <math.h>
#include "famg_misc.h"

#ifdef UG_DRAW
extern "C"
{
#include <assert.h>
#include "ugdevices.h"
}
#define FAMG_IOBUFFFER_LEN 10000
static char testf[FAMG_IOBUFFFER_LEN];
#endif

#ifdef ModelP
#include "parallel.h"
#include "pargm.h"
#endif

/* RCS_ID
$Header$
*/

void FAMGError(ostrstream &OutputString)
{
#ifdef UG_DRAW
	ostrstream ostr(testf,FAMG_IOBUFFFER_LEN);
	
	assert(OutputString.pcount()<FAMG_IOBUFFFER_LEN-20);
	ostr << "Error !" << OutputString.rdbuf() << '\0';
	UserWrite( ostr.str() );
#else	
    cerr << "Error ! " << OutputString.rdbuf() << flush;
#endif
}

void FAMGWarning(ostrstream &OutputString)
{
#ifdef UG_DRAW
	ostrstream ostr(testf,FAMG_IOBUFFFER_LEN);
	
	assert(OutputString.pcount()<FAMG_IOBUFFFER_LEN-20);
	ostr << "Warning !" << OutputString.rdbuf() << '\0';
	UserWrite( ostr.str() );
#else	
    cerr  << "Warning !" << OutputString.rdbuf() << flush;
#endif
}

void FAMGWrite(ostrstream &OutputString)
{
#ifdef UG_DRAW
	ostrstream ostr(testf,FAMG_IOBUFFFER_LEN);
	
	assert(OutputString.pcount()<FAMG_IOBUFFFER_LEN);
	ostr << OutputString.rdbuf() << '\0';
	UserWrite( ostr.str() );
#else	
	cout << OutputString.rdbuf() << flush;
#endif
}

double FAMGNorm(const int n, const double *v)
{
    double norm = 0.0;
    int i;

    for(i = 0; i < n; i++) norm += v[i]*v[i];
    norm = sqrt(norm);

    return norm;
}



int LR_Decomp (const short n, double *decomp, short *pivotmap)
{
	register double dinv, piv, val, factor;
	register short i, j, k, offset_i, offset_j;
	
	for (i = 0; i < n; i++) pivotmap[i] = i;

	/* LR decomposition */
	for (i = 0; i < n; i++)
	{
		/* pivot search */
		k = i;
		piv = fabs(decomp[pivotmap[i]*n+i]);
		for (j = i+1; j < n; j++)
		{
			val = fabs(decomp[pivotmap[j]*n+i]);
			if (val > piv)
			{
				k = j;
				piv = val;
			}
		}
	   
		/* reorder  */
		if (k != i)
		{
			j = pivotmap[k];
			pivotmap[k] = pivotmap[i];
			pivotmap[i] = j;
		}
		
		offset_i = pivotmap[i]*n;
		dinv = decomp[offset_i+i];
		if (fabs(dinv) < 1e-15) 
        {
            cout << "Warning! LR_Decomp: pivot too small," << endl << flush;
            return(1);
        }
            

		dinv = decomp[offset_i+i] = 1.0/dinv; // store inverse of the diagonal
		for (j = i+1; j < n; j++)
		{
			offset_j = pivotmap[j]*n;
			factor = decomp[offset_j+i] * dinv; 
            for (k = i+1; k < n; k++) decomp[offset_j+k] -= decomp[offset_i+k] * factor;
            decomp[offset_j+i] = factor;

		}
	}
	
	return(0);
}

int LR_Decomp (const short n, double *decomp)
{
    // without pivot search
	register double dinv, factor;
	register short i, j, k, offset_i, offset_j;
	

	/* LR decomposition */
	for (i = 0; i < n; i++)
	{
		offset_i = i*n;
		dinv = decomp[offset_i+i];
		if (fabs(dinv) < 1e-15) 
        {
            cout << "Warning! LR_Decomp: pivot too small." << endl << flush;
            dinv = 1e-15;
        }

		dinv = decomp[offset_i+i] = 1.0/dinv; // store inverse of the diagonal
		for (j = i+1; j < n; j++)
		{
			offset_j = j*n;
			factor = decomp[offset_j+i] * dinv; 
            for (k = i+1; k < n; k++) decomp[offset_j+k] -= decomp[offset_i+k] * factor;
            decomp[offset_j+i] = factor;

		}

	}
        	
	return(0);
}


int LR_Solve (const short n, const double *decomp, const short *pivotmap, double *x, const double *b)
{
	register short i, j, offset_i;
	register double sum;

    // foreward subsitution 
	for (i = 0; i < n; i++)
	{
		sum = b[pivotmap[i]];
        offset_i = pivotmap[i]*n;
		for (j = 0; j < i; j++)
		{
			sum -= decomp[offset_i+j] * x[j];
		}
		x[i] = sum;
	}

    // backward subsitution 
	for (i=n-1; i>=0; i--)
	{
		offset_i = pivotmap[i]*n;
		sum = x[i];
		for (j = i+1; j < n; j++) sum -= decomp[offset_i+j] * x[j];
		x[i] = sum * decomp[offset_i+i];	// the inverse of the diagonal is stored
	}

	return (0);
}

int LR_Solve (const short n, const double *decomp, double *x, const double *b)
{
    // without pivot search
	register short i, j, offset_i ;
	register double sum;

    // foreward subsitution 
	for (i = 0; i < n; i++)
	{
		sum = b[i];
        offset_i = i*n;
		for (j = 0; j < i; j++)
		{
			sum -= decomp[offset_i+j] * x[j];
		}
		x[i] = sum;
	}

    // backward subsitution 
	for (i = n-1; i >= 0; i--)
	{
		offset_i = i*n;
		sum = x[i];
		for (j = i+1; j < n; j++) sum -= decomp[offset_i+j] * x[j];
		x[i] = sum * decomp[offset_i+i];	// the inverse of the diagonal is stored
	}

	return (0);
}

int LR_SolveT (const short n, const double *decomp, double *x, const double *b)
{
    // without pivot search
	register short i, j;
	register double sum;

    // foreward subsitution 
	for (i = 0; i < n; i++)
	{
		sum = b[i];
		for (j = 0; j < i; j++)
		{
			sum -= decomp[j*n+i] * x[j];
		}
		x[i] = sum*decomp[i*n+i]; // the inverse of the diagonal is stored
	}

    // backward subsitution 
	for (i = n-1; i >= 0; i--)
	{
		sum = x[i];
		for (j = i+1; j < n; j++)
        {
            sum -= decomp[j*n+i] * x[j];
        }
		x[i] = sum ;	
	}

	return (0);
}

void FAMGSetVector(const int n, double *v, const double val)
{
    int i;

    for(i = 0; i < n; i++) v[i] = val;
}

void FAMGCopyVector(const int n, double *v1, const double *v2)
{
    int i;

    for(i = 0; i < n; i++) v1[i] = v2[i];
}

void FAMGCopyScaledVector(const int n, double *v1, const double *v2, const double factor)
{
    int i;

    for(i = 0; i < n; i++) v1[i] = factor*v2[i];
}

void FAMGSubVector(const int n, double *v1, const double *v2)
{
    int i;

    for(i = 0; i < n; i++) v1[i] -= v2[i];
}

void FAMGAddVector(const int n, double *v1, const double *v2)
{
    int i;

    for(i = 0; i < n; i++) v1[i] += v2[i];
}

void FAMGAddVector(const int n, double *v1, const double *v2, const double factor)
{
    int i;

    for(i = 0; i < n; i++) v1[i] += v2[i]*factor;
}

void FAMGAddVector(const int n, double *v1, const double factor, const double *v2)
{
    int i;

    for(i = 0; i < n; i++) v1[i] = v1[i]*factor+v2[i];
}


void FAMGMultVector(const int n, double *v1, const double factor)
{
    int i;

    for(i = 0; i < n; i++) v1[i] = v1[i]*factor;
}

void FAMGSetSubVector(const int n, double *v1, const double *v2, const double *v3)
{
    int i;

    for(i = 0; i < n; i++) v1[i] = v2[i] - v3[i];
}

double FAMGSum(const int n, const double *v1)
{
    double sum=0.0;
    int i;

    for(i = 0; i < n; i++) sum += v1[i];

    return sum;
}

double FAMGScalProd(const int n, const double *v1, const double *v2)
{
    double sum=0.0;
    int i;

    for(i = 0; i < n; i++) sum += v1[i]*v2[i];

    return sum;
}


void FAMGEigenVector(int n, double *a, double *b, double *e)
{
    int i, z;
    double norm;

    e[n-1] = 1.0;
    for(i = 0; i < n-2; i++) e[i] = 0.0;
    
    for(z = 0; z < 100; z++)
    {
        for(i = 1; i < n; i++) e[i] += e[i-1];
        for(i = 0; i < n; i++) e[i] = a[i]*e[i];
        for(i = n-2; i >= 0; i--) e[i] += b[i]*e[i+1];

        norm = 0.0;
        for(i = 0; i < n; i++) norm += e[i]*e[i];
        norm = sqrt(norm);
        for(i = 0; i < n; i++) e[i] = e[i]/norm;
    }
    

    return;
}

double FAMGTimeVar;

void PrintTIME( double time, char *text )
{
#ifdef ModelP
	int maxpe = -1;
	double maxt = time;
	maxt = UG_GlobalMaxDOUBLE( (DOUBLE)time );

	if( maxt == time )
		maxpe = me;

	maxpe = UG_GlobalMaxINT( (INT)maxpe );

	if( me == 0 )
		printf(PFMT" PrintTIME %s %g maxt %g mpe %d\n",me, text, time, maxt, maxpe );
#else
	printf("  0: PrintTIME %s mt %g mpe 0\n",me, text, time );
#endif
	return;
}
