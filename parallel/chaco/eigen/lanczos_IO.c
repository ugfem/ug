// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "../main/structs.h"
#include "../main/defs.h"

/* Comments from parent version, Lanczos_FO. */
/* Lanczos iteration with FULL orthogonalization.*/
/* Finds the lowest (zero) eigenvalue, but does not print it 
   or compute corresponding eigenvector. */
/* Does NOT find distinct eigenvectors corresponding to multiple
   eigenvectors and hence will not lead to good partitioners
   for symmetric graphs. Will be useful for graphs known
   (believed) not to have multiple eigenvectors, e.g. graphs
   in which a random set of edges have been added or perturbed. */
/* May fail on small graphs with high multiplicities (e.g. k5) if Ritz 
   pairs converge before there have been as many iterations as the number
   of eigenvalues sought. This is rare and a different random number
   seed will generally alleviate the problem. */
/* Convergence check uses Paige bji estimate over the whole
   spectrum of T. This is a lot of work, but we are trying to be 
   extra safe. Since we are orthogonalizing fully, we assume the 
   bji esitmates are very good and don't provide a contingency for 
   when they don't match the residuals. */
/* A lot of the time in this routine (say half) is spent in ql finding
   the evals of T on each iteration. This could be reduced by only using
   ql say every 10 steps. This might require that the orthogonalization
   be done with Householder (rather than Gram-Scmidt as currently) 
   to avoid a problem with zero beta values causing subsequent breakdown.
   But this routine is really intended to be used for smallish problems
   where the Lanczos runs will be short, e.g. when we are using the inverse
   operator method. In the inverse operator case, it's really important to
   stop quickly to avoid additional back solves. If the ql work is really 
   a problem then we should be using a selective othogonalization algorithm. 
   This routine provides a convnenient reference point for how well those 
   routines are functioning since it has the same basic structure but just 
   does more orthogonalizing. */
/* The algorithm orthogonalizes the starting vector and the residual vectors
   against the vector of all ones since we know that is the null space of
   the Laplacian. This generally saves net time because Lanczos tends to
   converge faster. */ 

/* Comments for this version, Lanczos with inverted operator. */
/* Used Symmlq for the back solve since already maintaining that code for
   the RQI/Symmlq multilevel method. Straight CG would only be marginally
   faster. */
/* The opthogonalization against the vector of all ones in Symmlq is not
   as efficient as possible - could in principle use orthog1 method, but
   don't want to disrupt the RQI/Symmlq code. Also, probably would only
   need to orthogonlize out a given mode periodically. But we want to be
   extra robust numerically. */

void lanczos_IO(A, n, d, y, lambda, bound, eigtol, vwsqrt)
struct vtx_data **A;		/* graph data structure */
int n;				/* number of rows/colums in matrix */
int d;				/* problem dimension = # evecs to find */
double **y;			/* columns of y are eigenvectors of A  */
double *lambda;			/* ritz approximation to eigenvals of A */
double *bound;			/* on ritz pair approximations to eig pairs of A */
double eigtol;			/* tolerance on eigenvectors */
double *vwsqrt;			/* square root of vertex weights */

{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	extern int DEBUG_EVECS;		/* print debugging output? */
	extern int WARNING_EVECS;	/* print warning messages? */
	extern double WARNING_ORTHTOL;	/* Warning on modest level of orthogonality loss */
	extern double WARNING_MISTOL;	/* Warning on serious level of orthogonality loss */
	extern double WARNING_SRESTOL;	/* Warning on inaccurate computation of evec of T */
	extern double DOUBLE_MAX;	/* Warning on inaccurate computation of evec of T */
	extern double splarax_time;	/* time matvecs */
	extern double orthog_time;	/* time orthogonalization work */
	extern double tevec_time;	/* time tridiagonal eigvec work */
	extern double evec_time;	/* time to generate eigenvectors */
	extern double ql_time;		/* time tridiagonal eigval work */
	extern double blas_time;	/* time for blas (not assembly coded) */
	extern double init_time;	/* time for allocating memory, etc. */
	extern double scan_time;	/* time for scanning bounds list */
	extern double debug_time;	/* time for debug computations and output */
	int i,j,k;			/* indicies */
	int maxj;			/* maximum number of Lanczos iterations */
	double *u,*r;			/* Lanczos vectors */
	double *Aq;			/* sparse matrix-vector product vector */
	double *alpha, *beta;  		/* the Lanczos scalars from each step */
	double *ritz;			/* copy of alpha for tqli */
	double *workj;  		/* work vector (eg. for tqli) */
	double *workn;  		/* work vector (eg. for checkeig) */
	double *s;			/* eigenvector of T */
	double **q;			/* columns of q = Lanczos basis vectors */
	double *bj;			/* beta(j)*(last element of evecs of T) */
	double Sres;			/* how well Tevec calculated eigvecs */
	double Sres_max;		/* Maximum value of Sres */
	double *Ares;			/* how well Lanczos calculated each eigpair */
	double *inv_lambda;		/* eigenvalues of inverse operator */
	int *index;			/* the Ritz index of an eigenpair */
	struct orthlink *orthlist;	/* vectors to orthogonalize against in Lanczos */
	struct orthlink *orthlist2;	/* vectors to orthogonalize against in Symmlq */
	struct orthlink *temp;		/* for expanding orthogonalization list */
	double *ritzvec;		/* ritz vector for current iteration */
	double *zeros;			/* vector of all zeros */
	double *ones;			/* vector of all ones */
	struct scanlink *scanlist;	/* list of fields for min ritz vals */
	struct scanlink *curlnk;	/* for traversing the scanlist */
	double bji_tol;			/* tol on bji estimate of A e-residual */
	int converged;			/* has the iteration converged? */
	double time;			/* current clock time */
	int warning1;			/* is warning1 cond. (eigtol not achieved) true? */
	int warning2;			/* is warning2 cond. (premature orth. loss) true? */
	int warning3;			/* is warning3 cond. (suspected misconvergence) true? */
	double shift, rtol;		/* symmlq input */
	long precon, goodb, nout;	/* symmlq input */
	long checka, intlim;		/* symmlq input */
	double anorm, acond;		/* symmlq output */
	double rnorm, ynorm;		/* symmlq output */
	long istop, itn;		/* symmlq output */
	double macheps;			/* machine precision calculated by symmlq */
	double normxlim;		/* a stopping criteria for symmlq */
	long itnmin;			/* enforce minimum number of iterations */
	int symmlqitns;			/* # symmlq itns */
	double *wv1,*wv2, *wv3;		/* Symmlq work space */
	double *wv4, *wv5, *wv6;	/* Symmlq work space */
	long long_n;			/* long int copy of n for symmlq */

	double *mkvec();		/* allocates space for a vector */
	double dot();			/* standard dot product routine */
	struct orthlink *makeorthlnk();	/* make space for entry in orthog. set */ 
	double norm();			/* vector norm */
	double Tevec();			/* calc evec of T by linear recurrence */
	double checkeig();		/* calculate residual of eigenvector of A */
	struct scanlink *mkscanlist();	/* make scan list for min ritz vecs */
	double seconds();		/* current clock timer */
	int symmlq_();
	void setvec(), scale(), update(), ql(), scadd(), vecran();
	void cpvec(), scanmax(), frvec(), orthogonalize(), orthog1(), shell_sort();
	void doubleout(), orthogvec(), exit();

	if (DEBUG_EVECS > 0) {
	{char buf[150]; sprintf(buf,"Full orthogonalization Lanczos, inverted operator, matrix size = %d\n",n);UserWrite(buf);}
	}

	/* Initialize time. */
	time = seconds();

	if (n < d+1) {
	{char buf[150]; sprintf(buf,"ERROR: System too small for number of eigenvalues requested.\n");UserWrite(buf);}
		exit(0);
		/* d+1 since don't use zero eigenvalue pair */
	}

	/* Allocate Lanczos space. */
	u = mkvec(1,n);
    if (!MEM_OK) return;
	r = mkvec(1,n);
    if (!MEM_OK) return;
	Aq = mkvec(1,n);
    if (!MEM_OK) return;
	ritzvec = mkvec(1,n);
    if (!MEM_OK) return;
	zeros = mkvec(1,n); 
    if (!MEM_OK) return;
	setvec(zeros,1,n,0.0);
	workn = mkvec(1,n);
    if (!MEM_OK) return;
	Ares = mkvec(1,d);
    if (!MEM_OK) return;
	inv_lambda = mkvec(1,d);
    if (!MEM_OK) return;
	index = (int *) (MEM_OK = smalloc((unsigned) (d+1) * sizeof(int));
    if (!MEM_OK) return;
	maxj = n;
	alpha = mkvec(1,maxj);
    if (!MEM_OK) return;
	beta = mkvec(1,maxj+1);
    if (!MEM_OK) return;
	ritz = mkvec(1,maxj);
    if (!MEM_OK) return;
	s = mkvec(1,maxj);
    if (!MEM_OK) return;
	bj = mkvec(1,maxj);
    if (!MEM_OK) return;
	workj = mkvec(1,maxj+1);
    if (!MEM_OK) return;
	q = (double **) (MEM_OK = smalloc((unsigned) (maxj+1) * sizeof(double *)); 
    if (!MEM_OK) return;
	scanlist = mkscanlist(d);
    if (!MEM_OK) return;

	/* Allocate Symmlq space all in one chunk. */
	wv1 = (double *) (MEM_OK = smalloc((unsigned) 6*(n+1)*sizeof(double));
    if (!MEM_OK) return;
	wv2 = &wv1[(n+1)];
	wv3 = &wv1[2*(n+1)];
	wv4 = &wv1[3*(n+1)];
	wv5 = &wv1[4*(n+1)];
	wv6 = &wv1[5*(n+1)];

	/* Set invariant symmlq parameters */
        precon = FALSE;		/* FALSE until we figure out a good way */
        goodb = FALSE; 		/* should be FALSE for this application */
        checka = FALSE;		/* if don't know by now, too bad */
        intlim = n;		/* set to enforce a maximum number of Symmlq itns */
	itnmin = 0;		/* set to enforce a minimum number of Symmlq itns */
	shift = 0.0;		/* since just solving rather than doing RQI */
	symmlqitns = 0;	 	/* total number of Symmlq iterations */
        nout = 0;		/* Effectively disabled - see notes in symmlq.f */
	rtol = 1.0e-5; 		/* requested residual tolerance */
	normxlim = DOUBLE_MAX;	/* Effectively disables ||x|| termination criterion */ 	
	long_n = n;		/* copy to long for linting */

	/* Initialize. */ 
	vecran(r,1,n); 
	if (vwsqrt == NULL) {
		/* whack one's direction from initial vector */
		orthog1(r,1,n);

		/* list the ones direction for later use in Symmlq */
		orthlist = NULL;
		orthlist2 = makeorthlnk();
        if (!MEM_OK) return;
		ones = mkvec(1,n);
        if (!MEM_OK) return;
		setvec(ones,1,n,1.0);
		orthlist2->vec = ones;
		orthlist2->pntr = NULL;
	}
	else {
		/* whack vwsqrt direction from initial vector */
		orthogvec(r,1,n,vwsqrt);

		/* list the vwsqrt direction for later use in Symmlq */
		orthlist = NULL;
		orthlist2 = makeorthlnk();
        if (!MEM_OK) return;
		orthlist2->vec = vwsqrt;
		orthlist2->pntr = NULL;
	}
	beta[1] = norm(r,1,n);
	q[0] = zeros;
	bji_tol = eigtol;
	Sres_max = 0.0;
	init_time += seconds() - time;

	/* Main Lanczos loop. */
	j = 1; 
	converged = FALSE;
	while((j <= maxj) && (converged == FALSE)) {
		time = seconds();
		q[j] = mkvec(1,n);
        if (!MEM_OK) return;
		scale(q[j],1,n,1.0/beta[j],r);
		blas_time += seconds() -time;
		time = seconds();
        	symmlq_(&long_n,&(q[j][1]),&wv1[1],&wv2[1],&wv3[1],&wv4[1],&Aq[1],&wv5[1],&wv6[1],
			&checka,&goodb,&precon,&shift,&nout,
        		&intlim,&rtol,&istop,&itn,&anorm,&acond,
	  		&rnorm,&ynorm,(double *)A,vwsqrt,(double *)orthlist2,
			&macheps,&normxlim,&itnmin);    
            if (!MEM_OK) return;
		symmlqitns += itn;
		if (DEBUG_EVECS > 1) {
		{char buf[150]; sprintf(buf,"Symmlq report:      rtol %g\n",rtol);UserWrite(buf);}
		{char buf[150]; sprintf(buf,"  system norm %g, solution norm %g\n",anorm,ynorm);UserWrite(buf);}
		{char buf[150]; sprintf(buf,"  system condition %g, residual %g\n",acond,rnorm);UserWrite(buf);}
		{char buf[150]; sprintf(buf,"  termination condition %2d, iterations %3d\n",istop,itn);UserWrite(buf);}
		}
		splarax_time += seconds() - time;
		time = seconds();
		update(u,1,n,Aq,-beta[j],q[j-1]);
		alpha[j] = dot(u,1,n,q[j]);
		update(r,1,n,u,-alpha[j],q[j]);
		blas_time += seconds() - time;
		time = seconds();
		if (vwsqrt == NULL) {
			orthog1(r,1,n);
		}
		else {
			orthogvec(r,1,n,vwsqrt);
		}
		orthogonalize(r,n,orthlist);
		temp = orthlist;
		orthlist = makeorthlnk();
        if (!MEM_OK) return;
		orthlist->vec = q[j];
		orthlist->pntr = temp;
		beta[j+1] = norm(r,1,n);
		orthog_time += seconds() - time;

		time = seconds();
		cpvec(ritz,1,j,alpha);
		cpvec(workj,1,j+1,beta);
		ql(ritz,workj,j);
		shell_sort(j,ritz);
		ql_time += seconds() - time;

		/* Convergence check using Paige bji estimates. */
		time = seconds();
		for (i=1; i<=j; i++)  {
			Sres = Tevec(alpha,beta,j,ritz[i],s);
			if (Sres > Sres_max) {Sres_max = Sres;}
			bj[i] = s[j]*beta[j+1];
/*{char buf[150]; sprintf(buf,"bj[%d] %20.16f    ritz[%d] %20.16f     1/ritz  %20.16f\n",i,bj[i],i,ritz[i],1.0/ritz[i]);UserWrite(buf);} */
		}
		tevec_time += seconds() - time;


		time = seconds();
		scanmax(ritz,1,j,&scanlist);
		converged = TRUE;
		if (j <= d) converged = FALSE; 
		else {
			curlnk = scanlist;
			while (curlnk != NULL) {
				if (bj[curlnk->indx] > bji_tol) {
					converged = FALSE;
				}
				curlnk = curlnk->pntr;
			}
		}
		scan_time += seconds() - time; 
		j++;
	}
	j--;

	/* Collect eigenvalue and bound information. */
	time = seconds();
	i = d; 
	curlnk = scanlist;
	while (curlnk != NULL) {
		inv_lambda[i] = curlnk->val;
		lambda[i] = 1.0/curlnk->val;
		bound[i] = bj[curlnk->indx];
		index[i] = curlnk->indx;
		curlnk = curlnk->pntr;
		i--;
	}
	scan_time += seconds() - time; 

	/* Compute eigenvectors. */
	time = seconds();
	for (i=1; i<=d; i++) {
		Sres = Tevec(alpha,beta,j,inv_lambda[i],s);
		if (Sres > Sres_max) {Sres_max = Sres;}
		setvec(y[i],1,n,0.0);
		for (k = 1; k <= j; k++) {
			scadd(y[i],1,n,s[k],q[k]);
		}
	}
	evec_time += seconds() - time;

	time = seconds();
	if (DEBUG_EVECS > 0 || WARNING_EVECS > 0) {
		for (i=1; i<=d; i++) {
			Ares[i] = checkeig(workn,A,y[i],n,lambda[i],vwsqrt,u);
            if (!MEM_OK) return;
		}
		/* NB lambda[i] and inv_lambda[i] have the same eigenvector */
	}

	if (DEBUG_EVECS > 0) {
	   {char buf[150]; sprintf(buf,"Lanczos itns. = %d\n",j);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"order     lambda               Ares est.              Ares           index\n");UserWrite(buf);}
		for (i=1; i<=d; i++) {
		{char buf[150]; sprintf(buf,"%2d.",i);UserWrite(buf);}
			doubleout(lambda[i],1);
			doubleout(bound[i],1);
			doubleout(Ares[i],1);
		{char buf[150]; sprintf(buf,"   %3d\n",index[i]);UserWrite(buf);}
		}
	{char buf[150]; sprintf(buf,"\nTotal Symmlq iterations %3d\n",symmlqitns);UserWrite(buf);}
	}

	if (WARNING_EVECS > 0) {
		warning1 = FALSE;
		warning2 = FALSE;
		warning3 = FALSE;
		for (i=1; i<=d; i++) {
			if (Ares[i] > eigtol) {warning1 = TRUE;}
			if (Ares[i]/bound[i] > WARNING_ORTHTOL && Ares[i] > 1.0e-14) {
				warning2 = TRUE;
			}
			if (Ares[i]/bound[i] > WARNING_MISTOL && Ares[i] > 1.0e-14) {
				warning3 = TRUE;
			}
			/* if Ares ~ 1.0e-14 its probably just round-off, so don't alarm user */
		}
		if (warning1) {
		{char buf[150]; sprintf(buf,"WARNING: Eigen pair tolerance (%g) not acheived.\n",eigtol);UserWrite(buf);}
		}
		if (warning2 && ! warning3) {
		{char buf[150]; sprintf(buf,"WARNING: Minor loss of orthogonality (Ares/est. > %g).\n",WARNING_ORTHTOL);UserWrite(buf);}
		}
		if (warning3) {
		{char buf[150]; sprintf(buf,"WARNING: Substantial loss of orthogonality (Ares/est. > %g).\n",WARNING_MISTOL);UserWrite(buf);}
		}
		if (j == maxj) {
		{char buf[150]; sprintf(buf,"WARNING: Maximum number of Lanczos iterations reached.\n");UserWrite(buf);}
		}
	}

	if (WARNING_EVECS > 1) {
		if (warning1 || warning2 || warning3) {
			if (DEBUG_EVECS <= 0) {
			{char buf[150]; sprintf(buf,"order     lambda               Ares est.              Ares           index\n");UserWrite(buf);}
				for (i=1; i<=d; i++) {
				{char buf[150]; sprintf(buf,"%2d.",i);UserWrite(buf);}
					doubleout(lambda[i],1);
					doubleout(bound[i],1);
					doubleout(Ares[i],1);
				{char buf[150]; sprintf(buf,"   %3d\n",index[i]);UserWrite(buf);}
				}
			}
			/* otherwise gets printed above */
		}
		if (WARNING_EVECS > 2 && Sres_max > WARNING_SRESTOL) {
		{char buf[150]; sprintf(buf,"WARNING: High maximum residual in computing eigenvector of T %g.\n",Sres_max);UserWrite(buf);}
		}
	}
	debug_time += seconds() - time;

	/* Free any memory allocated in this routine. */
	time = seconds();
	frvec(u,1);
	frvec(r,1);
	frvec(Aq,1);
	frvec(ritzvec,1);
	frvec(zeros,1);
	if (vwsqrt == NULL) {
		frvec(ones,1);
	}
	frvec(workn,1);
	frvec(Ares,1);
	frvec(inv_lambda,1);
	sfree((char *) index);
	frvec(alpha,1);
	frvec(beta,1);
	frvec(ritz,1);
	frvec(s,1);
	frvec(bj,1);
	frvec(workj,1);
	frvec(wv1,0);
	while (scanlist != NULL) {
		curlnk = scanlist->pntr;
		sfree((char *) scanlist);
		scanlist = curlnk;
	}
	for (i=1; i<=j; i++)  {
		frvec(q[i],1);
	}
	while (orthlist != NULL) {
		temp = orthlist->pntr;
		sfree((char *) orthlist);
		orthlist = temp;
	}
	while (orthlist2 != NULL) {
		temp = orthlist2->pntr;
		sfree((char *) orthlist2);
		orthlist2 = temp;
	}
	sfree((char *)q);	
	init_time += seconds() - time;
}
