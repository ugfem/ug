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
        
