// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Orthogonalize to all one's */
void orthog1(x,beg,end)
double *x; 
int beg, end;
{
        int i;
        double *pntr;
        double sum;
        int len;
 
        len = end - beg + 1;
        sum = 0.0;
        pntr = x + beg;
        for (i = len; i; i--)  {
                sum += *pntr++;
        }
        sum /= len;
        pntr = x + beg;
        for (i = len; i; i--)  {
                *pntr++ -= sum;
                
        }
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Orthogonalize to all one's */
void orthog1(x,beg,end)
double *x; 
int beg, end;
{
        int i;
        double *pntr;
        double sum;
        int len;
 
        len = end - beg + 1;
        sum = 0.0;
        pntr = x + beg;
        for (i = len; i; i--)  {
                sum += *pntr++;
        }
        sum /= len;
        pntr = x + beg;
        for (i = len; i; i--)  {
                *pntr++ -= sum;
                
        }
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Copy a range of a double vector*/
void cpvec(copy,beg,end,vec)
double *copy;
int beg, end;
double *vec;
{
	int i;

	copy = copy + beg;
	vec = vec + beg;
	for (i = end - beg + 1; i; i--)  {
		(*copy++) = *vec++;
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "../main/defs.h"

/* Finds needed eigenvalues of tridiagonal T using either the QL 
   algorithm or Sturm sequence bisection, whichever is predicted 
   to be faster based on a simple complexity model. */
void get_ritzvals(alpha,beta,j,Anorm,workj,ritz,d,left_goodlim,right_goodlim,eigtol) 
double *alpha;                  /* vector of Lanczos scalars */
double *beta;                   /* vector of Lanczos scalars */
int j;                          /* number of Lanczos iterations taken */
double Anorm;			/* Gershgorin estimate */
double *workj;                  /* work vector for Sturm sequence */
double *ritz;			/* array holding evals */
int d;				/* problem dimension = num. eigenpairs needed */
int left_goodlim;		/* number of ritz pairs checked on left end */
int right_goodlim;		/* number of ritz pairs checked on right end */
double eigtol;			/* tolerance on eigenpair */
{
	extern double BISECTION_SAFETY;	/* bisection tolerance function divided by this */
	extern int DEBUG_EVECS;		/* debug flag for eigen computation */
	int nvals_left;			/* numb. evals to find on left end of spectrum */
	int nvals_right;		/* numb. evals to find on right end of spectrum */
	double bisection_tol;		/* width of interval bisection should converge to */
	double predicted_steps;		/* predicts number of required bisection steps */
	void ql(), shell_sort(), bisect();

	/* Determine number of ritzvals to find on left and right ends */
	nvals_left = max(d,left_goodlim);
	nvals_right = min(j-nvals_left,right_goodlim);

	/* Predict work for bisection and ql assuming bisection takes roughly 5j flops per
           step, ql takes roughly 30j^j flops per call. (Ignored sort following ql.) */
	bisection_tol = eigtol * eigtol / BISECTION_SAFETY;
	predicted_steps = (nvals_left + nvals_right)*
			  (log10(Anorm/bisection_tol)/log10(2.0));

	if (5*predicted_steps < 30*j*j) { 
		bisect(alpha,beta-1,j,Anorm,workj,ritz,nvals_left,nvals_right,bisection_tol);
		if (DEBUG_EVECS > 2) {char buf[150]; sprintf(buf,"  bisection\n");UserWrite(buf);}}
	}
	else {
		if (DEBUG_EVECS > 2) {char buf[150]; sprintf(buf,"  ql\n");UserWrite(buf);}}
		ql(ritz,workj-1,j);
		shell_sort(j,ritz); 
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"

/* Finds some eigenvalues of T using Sturm sequence bisection. 
   Based on Wilkinson's algorithm, AEP, p.302. Beta must be offset 
   as it is passed in so that indices are consistent with AEP. */
void bisect(alpha,beta,j,Anorm,workj,ritz,nevals_left, nevals_right,tol) 
double *alpha;                  /* vector of Lanczos scalars */
double *beta;                   /* vector of Lanczos scalars */
int j;                          /* number of Lanczos iterations taken */
double Anorm;			/* Gershgorin estimate */
double *workj;                  /* work vector for Sturm sequence */
int nevals_left;		/* number of evals on right to find */
int nevals_right;		/* number of evals on left to find */
double *ritz;			/* array holding evals */
double tol;			/* tolerance on bracket width */
{
	extern int DEBUG_EVECS;	/* debug flag for eigen computation */
	extern double DOUBLE_MAX;	/* largest double value */
	int index;		/* index of sturm polynomial */
 	int i;			/* loop index */
	double *pntr;		/* pntr to double array */
	double x1,x2;		/* the bracketing interval */
	int x1cnt, x2cnt;	/* Sturm counts at x1 and x2 */
	double x;		/* the inserted point */
	int xcnt;		/* the Sturm count at x */
	int steps;		/* number of bisection steps for a Ritzval */
	int tot_steps;		/* number of bisection steps for all Ritzvals */
	int numbracketed;	/* number of evals between x1 and x2 */
 	int sturmcnt();

	/* Initialize portion of ritz we will use (use max double so scanmin 
           will work properly when called later on) */ 
	pntr = &ritz[1];
	for (i=j; i; i--) {
		*pntr++ = DOUBLE_MAX;
	}

	tot_steps = 0;

	/* find evals on left in decreasing index order */
	x2 = Anorm;
	x2cnt = j;
	numbracketed = j;
	for (index = nevals_left; index >= 1; index--) {
		x1 = 0;
		x1cnt = 0;
		steps = 0;
		while ((x2 - x1) > tol || numbracketed > 1) {
			x = 0.5*(x1 + x2);
			xcnt = sturmcnt(alpha,beta,j,x,workj);
			if (xcnt >= index) { 
				x2 = x;
				x2cnt = xcnt;
			 }
			else {
				x1 = x;
				x1cnt = xcnt;
			}
			numbracketed = x2cnt - x1cnt;
			steps++;
		}
		ritz[index] = 0.5*(x1 + x2);
		if (DEBUG_EVECS > 2) {
		{char buf[150]; sprintf(buf,"index %d, bisection steps %d, root %20.16f\n", index,steps,ritz[index]);UserWrite(buf);}
		}
		tot_steps += steps;
	}

	/* find evals on right in increasing index order */
	x1 = 0;
	x1cnt = 0;
	for (index = j - nevals_right + 1; index <= j; index++) {
		x2 = Anorm;
		x2cnt = j;
		steps = 0;
		while ((x2 - x1) > tol || numbracketed > 1) {
			x = 0.5*(x1 + x2);
			xcnt = sturmcnt(alpha,beta,j,x,workj);
			if (xcnt >= index) { 
				x2 = x;
				x2cnt = xcnt;
			 }
			else {
				x1 = x;
				x1cnt = xcnt;
			}
			numbracketed = x2cnt - x1cnt;
			steps++;
		}
		ritz[index] = 0.5*(x1 + x2);
		if (DEBUG_EVECS > 2) {
		{char buf[150]; sprintf(buf,"index %d, bisection steps %d, root %20.16f\n", index,steps,ritz[index]);UserWrite(buf);}
		}
		tot_steps += steps;
	}
	if (DEBUG_EVECS > 2) {
	{char buf[150]; sprintf(buf,"Total number of bisection steps %d\n",tot_steps);UserWrite(buf);}	
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Eigensolution of real symmetric tridiagonal matrix. */
/* After Numerical Recipies p. 380, removed eigenvector calc. */
/* CAREFUL: e is used as workspace, d returns eigenvals. */

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

void ql(d,e,n)
double d[],e[];
int n;

{
        int m,l,iter,i;
        double s,r,p,g,f,dd,c,b;
        void nrerror();
	double fabs();
	double sqrt();

        /* off diagonal data starts at second vector location, so renumber */
        for (i = 2; i <= n; i++) e[i-1] = e[i]; 
        e[n] = 0.0;

        for (l = 1; l <= n; l++)  {
                iter = 0;
                do {
                        for (m = l; m <= n-1; m++)  {
                                dd = fabs(d[m]) + fabs(d[m+1]);
                                if (fabs(e[m]) + dd == dd) break;
                        }
                        if (m != l)  {
                                if (iter++ == 50) nrerror("ql: Too many iterations (max 50).");
                                g = (d[l+1] - d[l])/(2.0 * e[l]);
                                r = sqrt((g*g) + 1.0);
                                g = d[m] - d[l] + e[l]/(g + SIGN(r,g));
                                s = c = 1.0;
                                p = 0.0;
                                for (i = m-1; i >= l; i--)  {
                                        f = s * e[i];
                                        b = c * e[i];
                                        if (fabs(f) >= fabs(g))  {
                                                c = g/f;
                                                r = sqrt((c*c) + 1.0);
                                                e[i+1] = f * r;
                                                c *= (s = 1.0/r);
                                        }
                                        else  {
                                                s = f/g;
                                                r = sqrt((s*s) + 1.0);
                                                e[i+1] = g * r;
                                                s *= (c = 1.0/r);
                                        }
                                        g = d[i+1] - p;
                                        r = (d[i] - g) * s + 2.0*c*b;
                                        p = s * r;
                                        d[i+1] = g + p;
                                        g = c*r - b;
                                }
                                d[l] = d[l] - p;
                                e[l] = g;
                                e[m] = 0.0;
                        }
                }
                while (m != l);

        }
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"

/* Return minimum of vector over range */
void scanmin(vec,beg,end,scanlist)
double *vec;			/* vector to scan */
int beg,end;                    /* index range */
struct scanlink **scanlist;	/* pntr to list holding results of scan */
{
	extern double DOUBLE_MAX;
        struct scanlink *top;
        struct scanlink *curlnk;
        struct scanlink *prevlnk;
	double val;
	int i;
 
	curlnk = *scanlist;
	while (curlnk != NULL) {
		curlnk->indx = 0;
		curlnk->val = DOUBLE_MAX;
		curlnk = curlnk->pntr;
	}

/* Note: Uses current top link (which would need to be deleted anyway) each time
	 an insertion to the list is required. */

	for (i=beg; i<=end; i++) { 
		/* consider each element for insertion */
		top = *scanlist;
		val = vec[i];
		if (val < top->val) { 
			if (top->pntr == NULL) {
				/* the list is only one long, so just replace */
				top->val = val;
				top->indx = i;
			}
			else {
				/* beats top element; scan for insertion point */
				if (val < (top->pntr)->val) {  
					/* 2nd link becomes list pntr; otherwise stays same */
					*scanlist = top->pntr;  
				}
       				prevlnk = curlnk = top;
        			while ((val < curlnk->val)&&(curlnk->pntr != NULL)) {
					prevlnk = curlnk;
					curlnk = curlnk->pntr;
				}
				if (val < curlnk->val) {
					/* got to end of list; add top to bottom */
					curlnk->pntr = top;
					top->val = val;
					top->indx = i;
					top->pntr = NULL;
				}
				else { 
					/* stopped within list; insert top here */
					prevlnk->pntr = top;
					top->val = val;
					top->indx = i;
					top->pntr = curlnk;
				}
			}
		}
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Copy a range of a double vector*/
void cpvec(copy,beg,end,vec)
double *copy;
int beg, end;
double *vec;
{
	int i;

	copy = copy + beg;
	vec = vec + beg;
	for (i = end - beg + 1; i; i--)  {
		(*copy++) = *vec++;
	}
}
/* An array of these stores all the data for the graph/matrix. */
/* PERHAPS DON'T NEED GNEDGES AND VNUM */
struct vtx_data {
	/*int vnum;*/		/* global number of this vertex */
	int vwgt;		/* weight of vertex */
	int nedges;		/* number of neighbors of vertex in subgraph */
	/*int gnedges;*/	/* number of neighbors in original graph */
				/* Note: above fields always include self-edge */
	int *edges;		/* neighbor list in subgraph numbering scheme */
	float *ewgts;		/* weights of all the edges */
				/* Note: above 2 fields have self-edge first */
};


/* Doubly linked list to hold coarsened graphs */
struct graphlist {
	struct vtx_data **graph;	/* current graph */
	int nvtxs;			/* number of vertices in graph */
	int nedges;			/* number of edges in graph */
	double *vwsqrt;			/* square roots of vertex weights */
	double **evecs;			/* lowest eigenvectors of graph */
	double *evals;			/* corresponding eigenvalues of graph */
	int *mflag;			/* matching data for contracted edges */
	int *v2cv;			/* mapping from fine to coarse vtxs */
	struct ipairs *merged;		/* pairs of merged vertices */
	struct fpairs *reduction;	/* values used to reduce edge weights */
	struct graphlist *finer;	/* pointer to finer graph */
	struct graphlist *coarser;	/* pointer to coarser graph */
};


/* Used in root finding in the Lanczos code */
struct arglist {
	double zleft;		/* left zero of rational func */
	double pole;		/* pole of rational func */
	double zright;		/* rigth zero of rational func */
	double *alpha;		/* Lanczos constants */
	double *beta;		/* Lanczos constants */
	int j;			/* step in the Lanczos process */ 
};


/* A linked list of these stores the selective orthogonalization set */
struct orthlink {
	int depth;		/* bottom of list is 0, previous is 1 etc */
	int index;		/* position in list of ritz vals (i index) */
	double ritzval;		/* good ritz value */		
	double betaji;		/* residual bound on good ritz pair */ 
	double tau;		/* from orthogonality recursion */
	double prevtau;		/* from orthogonality recursion */
	double *vec;		/* vector to orthogonalize against */
	struct orthlink *pntr;	/* pointer to next link */
};


/* A linked list of these stores the minimum elements of a vector */
struct scanlink {
	double val;		/* value of vector entry */
	int indx;		/* index of corresponding entry */
	struct scanlink *pntr;	/* pointer to next link */
};

/* These store the phantom edges needed to keep a subgraph connected */
struct edgeslist {
	 int vtx1;               /* first vertex in edge */
	 int vtx2;               /* second vertex in edge */
	 struct edgeslist *next; /* pointer to next element in list */
 };

/* These store all the data needed to modify edges for connectivity. */
struct connect_data {
	struct ilists *old_edges;	/* overwritten old edges */
	struct flists *old_ewgts;	/* overwritten old weights */
	struct edgeslist *new_edges;	/* list of new edges */
	int old_nedges;		/* original number of edges in graph */
};

/* Linked list stuff for various uses */
struct list {				/* linked list of integers */
	int num;			/* element number */
	struct list *next;		/* ptr to next element in list */
};

struct lists {				/* linked list of lists */
	struct list *begin;		/* pointer to list */
	struct lists *nextlist;		/* next list header */
};

struct bilist {				/* bidirectional list */
	struct bilist *prev;		/* pointer to previous element */
	struct bilist *next;		/* ptr to next element in list */
};

struct ipairs {				/* stores pairs of integers */
	int val1;
	int val2;
};

struct dpairs {				/* stores pairs of doubles */
	double val1;
	double val2;
};

struct fpairs {				/* stores pairs of floats */
	double val1;
	double val2;
};

struct ilists {				/* linked list of integer lists */
	int *list;			
	struct ilists *next;
};

struct flists {				/* linked list of floating lists */
	float *list;			
	struct flists *next;
};
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Set a double precision vector to constant over range. */
void setvec(vec,beg,end,setval)
double vec[];
int beg, end;
double setval;
{
        int i;

	vec = vec + beg;
        for(i = end - beg + 1; i; i--)   {
		(*vec++) = setval;
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Finds eigenvector s of T and returns residual norm. */
double Tevec(alpha,beta,j,ritz,s) 
double *alpha;                          /* vector of Lanczos scalars */
double *beta;                           /* vector of Lanczos scalars */
int j;                                  /* number of Lanczos iterations taken */
double ritz;                            /* approximate eigenvalue  of T */
double *s;                              /* approximate eigenvector of T */
{
        int i;                          /* index */
        double residual;                /* how well recurrence gives eigenvector */
 
        double normalize();             /* normalizes vector, returns norm */
	double fabs();			/* intrinsic absolute value */
        
        s[j] = 1.0;

        if (j == 1)  {
                residual = (alpha[1]-ritz);
        }
 
        if (j >= 2)  {
                s[j-1] = - (alpha[j]-ritz)/beta[j];
                for (i=j; i>=3; i--)  {
                        s[i-2] = -((alpha[i-1]-ritz)*s[i-1] + beta[i]*s[i])/beta[i-1];
                }
                residual = (alpha[1]-ritz)*s[1] + beta[2]*s[2];
        }
 
        residual = fabs(residual)/normalize(s,1,j);
        return(residual);
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include  <stdio.h>

/* Print a double precesion number with filtering format. */
void doubleout(number,mode)
double number;
int mode;
{
	double fabs();		/* intrinsic absolute value function */

	if (mode == 1) {
		if (fabs(number) < 100) {
		{char buf[150]; sprintf(buf,"  %19.16f",number);UserWrite(buf);}
		}
		else {
		{char buf[150]; sprintf(buf,"  %19g",number);UserWrite(buf);}
		}
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Scaled add - fills vec1 with vec1 + alpha*vec2 over range*/
void scadd(vec1,beg,end,fac,vec2)
double *vec1;
int beg, end;
double fac;
double *vec2;
{
	int i;

	vec1 = vec1 + beg;
	vec2 = vec2 + beg;
	for (i = end - beg + 1; i; i--)  {
		(*vec1++) += fac * (*vec2++);
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Check an eigenpair of A by direct multiplication.  */
double checkeig(err,A,y,n,lambda,vwsqrt,work)
double *err;
struct vtx_data **A;
double *y;
int n;
double lambda;
double *vwsqrt;
double *work;
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	double resid;
	double normy;
	double norm();
	void splarax(), scadd();

	splarax(err, A, n, y, vwsqrt, work);
    if (!MEM_OK) return;
	scadd(err,1,n,-lambda,y);
	normy = norm(y,1,n);
	resid = norm(err,1,n)/normy;
	return(resid);
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>

/* Returns 2-norm of an n-vector over range. */
double norm(vec,beg,end)
double *vec;
int beg, end;
{
	double dot();

	return(sqrt(dot(vec,beg,end,vec)));
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

void solistout(solist,n,ngood,j)
struct orthlink **solist;      	/* vector of pntrs to orthlnks */
int n;                        	/* length of vecs to orth. against */
int ngood;			/* number of good vecs on list*/
int j;				/* current number of Lanczos steps */
{
	int i;				/* index */
	extern int DEBUG_EVECS;		/* debugging output level for eigen computations */

	/* to placate alint */
 	n = n;
 
	for (i=1; i<=ngood; i++) {
		if ( (solist[i])->index <= (int)(j/2) ) {
		{char buf[150]; sprintf(buf,".");UserWrite(buf);}
		} 
		else {
		{char buf[150]; sprintf(buf,"+");UserWrite(buf);}
		}
		/*  Really detailed output:
	 {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
	{char buf[150]; sprintf(buf,"depth %d\n",(solist[i])->depth);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"index %d\n",(solist[i])->index);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"ritzval %g\n",(solist[i])->ritzval);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"betaji %g\n",(solist[i])->betaji);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"tau %g\n",(solist[i])->tau);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"prevtau %g\n",(solist[i])->prevtau);UserWrite(buf);}
		vecout((solist[i])->vec,1,n,"vec");
		*/
	}
{char buf[150]; sprintf(buf,"%d\n",ngood);UserWrite(buf);}

	if (DEBUG_EVECS > 2) {
	{char buf[150]; sprintf(buf,"  actual indicies: ");UserWrite(buf);}
		for (i=1; i<=ngood; i++) {
			{char buf[150]; sprintf(buf," %2d",solist[i]->index);UserWrite(buf);}
		}
	{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"

struct scanlink *mkscanlist(depth)
int depth;
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	struct scanlink *prevlnk;
	struct scanlink *newlnk;
	int i;


	prevlnk = (struct scanlink *) (MEM_OK = smalloc(sizeof(struct scanlink));
    if (!MEM_OK) return;
	prevlnk->pntr = NULL;
	newlnk = prevlnk; /* in case the list is one long */
	for (i=1; i<=(depth-1); i++) {
		newlnk = (struct scanlink *) (MEM_OK = smalloc(sizeof(struct scanlink));
        if (!MEM_OK) return;
		newlnk->pntr = prevlnk;
		prevlnk = newlnk;
	}
	return(newlnk);
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/defs.h"
#include "../main/structs.h"


/* Sparse linked A(matrix) times x(vector). */
void splarax(result, mat, n, vec, vwsqrt, work)
double *result;			/* result of matrix vector multiplication */
struct vtx_data **mat;		/* graph data structure */
int n;				/* number of rows/columns in matrix */
double *vec;			/* vector being multiplied by matrix */
double *vwsqrt;			/* square roots of vertex weights */
double *work;			/* work vector from 1-n */
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	extern int PERTURB;		/* perturb matrix? */
	extern int NPERTURB;		/* if so, number of edges to perturb */
	extern double PERTURB_MAX;	/* maximum value of perturbation */
	struct vtx_data *mat_i;		/* an entry in "mat" */
	double sum;			/* sums inner product of matrix-row & vector */
	int *colpntr;			/* loops through indices of nonzeros in a row */
	float *wgtpntr;			/* loops through values of nonzeros */
	int i, j;			/* loop counters */
	double *wrkpntr;		/* loops through indices of work vector */
	double *vwsqpntr;		/* loops through indices of vwsqrt */
	double *vecpntr;		/* loops through indices of vec */
	double *respntr;		/* loops through indices of result */
	void perturb();

	if (vwsqrt == NULL) {	/* No vertex weights */
	    if (mat[1]->ewgts == NULL) {    /* No edge weights */
		respntr = result;
	        for (i=1; i<=n; i++) {
		    mat_i = mat[i]; 
		    colpntr = mat_i->edges;
		    sum = (mat_i->nedges-1)*vec[*colpntr++];
		    for (j=mat_i->nedges-1; j; j--) {
		        sum -= vec[*colpntr++];
		    }
		    *(++respntr) = sum;
	        }
	    }
	    else {    /* Edge weights */
		respntr = result;
	        for (i=1; i<=n; i++) {
		    mat_i = mat[i]; 
	    	    colpntr = mat_i->edges;
		    wgtpntr = mat_i->ewgts;
		    sum = 0.0;
		    for (j=mat_i->nedges; j; j--) {
		        sum -= *wgtpntr++ * vec[*colpntr++];
		    }
		    *(++respntr) = sum;	/* -sum if want -Ax */
	        }
	    }
 	    if (PERTURB && NPERTURB > 0  && PERTURB_MAX > 0.0) {perturb(result, vec);
                                                            if (!MEM_OK) return;
                                                           }
	}
	else {    /* Vertex weights */
	    if (vwsqrt != NULL) {
		wrkpntr = work;
		vecpntr = vec;
		vwsqpntr = vwsqrt;
	        for (i=n; i; i--) {
	    	    *(++wrkpntr) = *(++vecpntr) / *(++vwsqpntr);
	        }
	    }

	    if (mat[1]->ewgts == NULL) {    /* No edge weights. */
		respntr = result;
	        for (i=1; i<=n; i++) {
		    mat_i = mat[i]; 
		    colpntr = mat_i->edges;
		    sum = (mat_i->nedges-1)*work[*colpntr++];
		    for (j=mat_i->nedges-1; j; j--) {
		        sum -= work[*colpntr++];
		    }
		    *(++respntr) = sum;
	        }
	    }
	    else{    /* Edge weights. */
		respntr = result;
	        for (i=1; i<=n; i++) {
		    mat_i = mat[i]; 
		    colpntr = mat_i->edges;
		    wgtpntr = mat_i->ewgts;
		    sum = 0.0;
		    for (j=mat_i->nedges; j; j--) {
		        sum -= *wgtpntr++ * work[*colpntr++];
		    }
		    *(++respntr) = sum;	/* -sum if want -Ax */
	        }
	    }
 	    if (PERTURB && NPERTURB > 0  && PERTURB_MAX > 0.0) {perturb(result, work);
                                                            if (!MEM_OK) return;
                                                           }

	    if (vwsqrt != NULL) {
		respntr = result;
		vwsqpntr = vwsqrt;
	        for (i=n; i; i--) {
		    *(++respntr) /= *(++vwsqpntr);
	        }
	    }
	}
}
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <math.h>
#include "../main/structs.h"


void makevwsqrt(vwsqrt, graph, nvtxs)
/* Make vector of square roots of vertex weights. */
double *vwsqrt;			/* vector returned */
struct vtx_data **graph;	/* graph data structure */
int nvtxs;			/* number of vertices in graph */
{
   extern int NSQRTS;		/* number of sqrts already computed */
   extern double *SQRTS;	/* values computed */
   short vwgt;			/* vertex weight */
   int i;			/* loop counter */

   for (i=1; i<= nvtxs; i++) {
      vwgt = graph[i]->vwgt;
      if (vwgt <= NSQRTS) vwsqrt[i] = SQRTS[(int) vwgt];
      else vwsqrt[i] = sqrt((double) vwgt);
   }
}


/* Extract the subgraph vwsqrt values */
void make_subvector(vec, subvec, subnvtxs, loc2glob)
double *vec;			/* vector for all vertices */
double *subvec;			/* vector for vertices in subgraph */
int subnvtxs;			/* number of vtxs in subgraph */
int *loc2glob;			/* subgraph -> graph numbering map */
{
   int i;

   for (i=1; i<=subnvtxs; i++) {
      ++subvec;
      (*subvec) = vec[loc2glob[i]];
   }
}
