/****************************************************************************/
/*																			*/
/* File:      famg.C														*/
/*																			*/
/* Purpose:   famg file interface											*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   April 98 begin, Stuttgart										*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#include <strstream.h>
#include <iostream.h>
#include <fstream.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "famg_interface.h"
#include "famg_system.h"

/* RCS_ID
$Header$
*/

int firsti;

static int ReadRHS(double *&rhs, int n, int argc, char **argv)
{
    double val;
    int i, res;
    char rhsname[32];

    rhs = (double*) new double[n];
    if(rhs == NULL) return 1;


    if(argc > 2)
    {
	    res= sscanf(argv[2],"%s",rhsname);
        if (res< 1) 
        {
            cout << "matrix not found." <<endl;
            return 1;
        }
        ifstream rhsfile(rhsname,ios::in);
        if (!rhsfile) 
        {
            cout << "right hand side not found." << endl;
            return 1;
        }
    
        while(rhsfile >> i >> val)
        {
			i -= firsti;
            rhs[i] = val;
            if((i+1)==n) break;
        }
    }
    else
    {
        for(i = 0; i < n; i++) rhs[i] = 1.0;
    }


    return 0;
            
}

static int ReadMatrix(double *&entry, int *&index, int *&start, int &n, int
&nl, int argc, char **argv)
{
    double m;
    int i,j,sym,res,offset,oldi;
    char mname[32];

    if(argc < 2)
    {
        cout << "matrix not found." <<endl;
        return 1;
    }

	res=sscanf(argv[1],"%s",mname);
	if (res< 1)
    {
        cout << "matrix not found." <<endl;
        return 1;
    }

    ifstream mfile(mname,ios::in);
	if (!mfile)
    {
        cout << "matrix not found." <<endl;
        return 1;
    }

	n = 0;
    mfile >> n >> sym;

    start = new int[n+1];
    if (start == NULL) return 1;

    // count matrix entries
    mfile >> i >> j >> m;
    firsti = i;
    offset = 1;
    while(mfile >> i >> j >> m)
    {
        offset++;
    }

    nl = offset;
    index = new int[nl];
    if (index == NULL) return 1;
    entry = new double[nl];
    if (entry == NULL) return 1;


    // read matrix

    mfile.clear();
    mfile.seekg(0, ios::beg);
    mfile >> n >> sym;

    offset = 0; oldi = firsti-1;
    while(mfile >> i >> j >> m)
    {
        i -= firsti;
        index[offset] = j;
        entry[offset] = m;
        if(i != oldi)
        {
            if (i != (oldi+1))
            {
                cout << __FILE__ << "  " << "ReadMatrix: wrong matrix format."
<< endl;
                return 1;
            }
            start[i] = offset;
            oldi = i;
        }
        offset++;
    }
    start[n] = offset;
    if(offset != nl)
    {
        cout << "wrong matrix format." << endl;
        return 1;
    }


    return 0;
}

static int WriteMatrix(double *&entry, int *&index, int *&start, int &n, int &nl, double *rhs)
{
    int ii, i;

    ofstream mfile("mat.out",ios::out);
	if (!mfile) 
    {
        cout << "matrix not found." <<endl;
        return 1;
    }
	
    mfile << n << nl << endl;

    for(i = 0; i < n; i++)
    {
        mfile << i+1 << "\t" << rhs[i] << endl;
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            mfile << i+1 << "\t" << index[ii] << "\t" << entry[ii] << endl;
        }
    }
        
    return 0;
}

static int isEmpty(char * str, char comment, int maxlen)
{
    for(int i=0;i<maxlen;i++)
	if(str[i]==comment || str[i]=='\0' || str[i]=='\n')return 1;
	else if(str[i]=='\t' || str[i]==' ')continue;
	else return 0;
    return 1;
}

int FAMGParameter::Read()
{
    heap = 100e+6;
    nv = 16;
    gamma = 1;
    n1 = 0;
    n2 = 1;
    ilut = 1e+10;
    cgilut = 0.0;
    cgnodes = 100;
    mincoarse = 0.8;
    conloops = 1;
    type = 0;
    stv = 0;
    tol = 0.95;
    sigma = 0.1;
    omegar = 1.0;
    omegal = 1.0;
    error1 = 1e-5;
    error2 = 10.0;
    maxit = 100;
    alimit = 1e-14;
    rlimit = 1e-10;
    divlimit = 1e+10;
    reduction = 1.0;
    strcpy(solver,"linit");
    strcpy(presmoother,"fgs");
    strcpy(postsmoother,"bgs");
    strcpy(cgsmoother,"ilut");
	coloringmethod = 3;
	
	const int bufLen=80;
	char buf[bufLen], pastr[30];

    ifstream infile("famg.in",ios::in);
	if (!infile) 
    {
        cout << "parameter input file not found." <<endl;
        return 1;
    }

    while(!infile.eof())
    {
        for(buf[0]='#';isEmpty(buf, '#', bufLen);)
        {
            if(infile.eof()) break;
            infile.getline(buf,bufLen);
        }
        
        istrstream ist(buf,bufLen);
        ist.precision(15);
        if (ist>>pastr)
        {
            if(strcmp(pastr,"heap") == 0)
            {
                double paheap;
                if(!(ist >> paheap)) heap = 120e+6;
                else heap = ceil(paheap); 
            }
            else if(strcmp(pastr,"nv") == 0)
            {
                if(!(ist >> nv)) nv = 16;
            }
            else if(strcmp(pastr,"gamma") == 0)
            {
                if(!(ist >> gamma)) gamma = 1;
            }
            else if(strcmp(pastr,"ilut") == 0)
            {
                if(!(ist >> ilut)) ilut = 1e+10;;
            }
            else if(strcmp(pastr,"cgilut") == 0)
            {
                if(!(ist >> cgilut)) cgilut = 0.0;
            }
            else if(strcmp(pastr,"n1") == 0)
            {
                if(!(ist >> n1)) n1 = 0;
            }
            else if(strcmp(pastr,"n2") == 0)
            {
                if(!(ist >> n2)) n2 = 1;
            }
            else if(strcmp(pastr,"cgnodes") == 0)
            {
                if(!(ist >> cgnodes)) cgnodes = 1;
            }
            else if(strcmp(pastr,"mincoarse") == 0)
            {
                if(!(ist >> mincoarse)) mincoarse = 0.8;
            }
            else if(strcmp(pastr,"conloops") == 0)
            {
                if(!(ist >> conloops)) conloops = 1;
            }
            else if(strcmp(pastr,"type") == 0)
            {
                if(!(ist >> type)) type = 0;
            }
            else if(strcmp(pastr,"stv") == 0)
            {
                if(!(ist >> stv)) stv = 0;
            }
            else if(strcmp(pastr,"tol") == 0)
            {
                if(!(ist >> tol)) tol = 0;
            }
            else if(strcmp(pastr,"sigma") == 0)
            {
                if(!(ist >> sigma)) sigma = 0.2;
            }
            else if(strcmp(pastr,"omegar") == 0)
            {
                if(!(ist >> omegar)) omegar = 1.0;
            }
            else if(strcmp(pastr,"omegal") == 0)
            {
                if(!(ist >> omegal)) omegal = 1.0;
            }
            else if(strcmp(pastr,"error1") == 0)
            {
                if(!(ist >> error1)) error1 = 1e-5;
            }
            else if(strcmp(pastr,"error2") == 0)
            {
                if(!(ist >> error2)) error2 = 100.0;
            }
            else if(strcmp(pastr,"maxit") == 0)
            {
                if(!(ist >> maxit)) maxit = 100;
            }
            else if(strcmp(pastr,"alimit") == 0)
            {
                if(!(ist >> alimit)) alimit = 1e-14;
            }
            else if(strcmp(pastr,"rlimit") == 0)
            {
                if(!(ist >> rlimit)) rlimit = 1e-10;
            }
            else if(strcmp(pastr,"divlimit") == 0)
            {
                if(!(ist >> divlimit)) divlimit = 1e+10;
            }
            else if(strcmp(pastr,"reduction") == 0)
            {
                if(!(ist >> reduction)) reduction = 1.08;
            }
            else if(strcmp(pastr,"solver") == 0)
            {
                if(!(ist >> solver)) strcpy(solver,"linit");
            }
            else if(strcmp(pastr,"presmoother") == 0)
            {
                if(!(ist >> presmoother)) strcpy(presmoother,"fgs");
            }
            else if(strcmp(pastr,"postsmoother") == 0)
            {
                if(!(ist >> postsmoother)) strcpy(postsmoother,"bgs");
            }
            else if(strcmp(pastr,"cgsmoother") == 0)
            {
                if(!(ist >> cgsmoother)) strcpy(cgsmoother,"ilut");
            }
            else if(strcmp(pastr,"coloringmethod") == 0)
            {
                if(!(ist >> coloringmethod)) coloringmethod = 3;
            }
       }
    }

    return 0;
}    



main(int argc, char **argv)
{
    double *entry, *rhs;
    int *index, *start, n, nl, status, i;
	struct FAMG_Interface my_famg_parameter;
	
    if (ReadMatrix(entry,index,start,n,nl,argc,argv)) return 1;
    if (ReadRHS(rhs,n,argc,argv)) return 1;

    // read parameter 
    FAMGParameter parameter;
    if(parameter.Read()) return 1; 

    // and here we go ....
    
    double *unknown = new double[n];
   double *defect = new double[n];


   if(FAMGConstructParameter(&parameter)) return 1;

	my_famg_parameter.n = n;
	my_famg_parameter.nl = nl;
	my_famg_parameter.extra = extra;
	my_famg_parameter.entry = entry;
	my_famg_parameter.index = index;
	my_famg_parameter.start = start;
	my_famg_parameter.nv = FAMG_NVECTORS;
	for( i=0; i<FAMG_NVECTORS ; i++ )
		my_famg_parameter.vector[i] = NULL;
	if(FAMGConstruct(&my_famg_parameter)) return 1;

   status = FAMGSolve(rhs,defect,unknown);

   FAMGDeconstruct();

   FAMGDeconstructParameter();



    delete [] entry; 
    delete [] index; 
    delete [] start; 
    delete [] rhs; 
    delete [] unknown; 
    delete [] defect; 


    return status;
}
