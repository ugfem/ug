// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <math.h>
#include "../main/structs.h"
#include "../main/defs.h"

/* These comments are for version 1: */
/* This version of selective orthogonalization largely follows that described in 
   Pareltt and Scott, "The Lanczos Algorithm with Selective Orthogonalization",
   Math Comp v33 #145, 1979. Different heuristics are used to control loss of
   orthogonality. Specifically, the pauses to check for convergence of Ritz pairs
   at first come at a small but increasing interval as the computation builds up.
   Later they come at a regular, pre-set interval, e.g. every 10 steps. Hence this 
   is similar to Grear's periodic reorthogonalization scheme, but with an adaptive 
   period. A small number of ritz pairs at both ends of the spectrum are monitored
   for convergence at each pause. The number monitored gradually increases as more 
   Ritz pairs converge so that the number monitored at each pause is always at least
   2 greater than the number that have previously converged. This is a fairly 
   conservative strategy, at least in the context of the typical Laplacian graph 
   matrices studied in connection with load balancing. But it may occasionally 
   result in premature loss of orthogonality. If the loss of orthogonality is
   severe, Lanczos may fail to meet the set eigen tolerance or mis-converge to
   the wrong eigenpair. These conditions will usually be detected by the algorithm
   and a warning issued. If the graph is small, full orthogonalization is
   then a fall back option. If the graph is large, one of the multi-level schemes
   should be tried. (Qualitatively, the multi-level schemes appear *less* susceptible 
   to misconvergence when a small eigen pair is sought.). The algorithm orthogonalizes
   the starting vector and each residual vector against the vector of all ones
   since we know that is an eigenvector of the Laplacian (which we don't want to see).
   This is accomplished through calls to orthog1. These can be removed, but then
   we must compute an extra eigenvalue and the vector of all ones shows up as
   a Ritz pair we need to orthogonalize against in order to maintain the basis.
   This is more expensive then simply taking it out of the residual on each step,
   and doing the latter also reduces the number of iterations. The Ritz values
   being monitored at each step are computed using the QL algorithm or bisection
   on the Sturm sequence, whichever is faster based on a simple complexity model.
   In practice the time spent computing Ritz values is only a small portion of
   the total Lanczos time, so more complex schemes based on interpolating the
   bottom pivot function represent only a very marginal savings in execution
   time at the expense of either some robustness or substantially increased
   code complexity. */

/* These comments are for version 2: */
/* These comments are cumulative from version 1. Modified the heuristics to
   make it faster. Now only checks left end of spectrum for converging Ritz pairs.
   Only checking the left end of the spectrum for converging Ritz pairs is in principle
   not as good as checking both ends, but in practice seems to be no worse and sometimes
   better than checking both ends. Since the distribution of eigenvalues for Laplacian
   graphs of interest seems to generally result in much more rapid convergence on the
   right end of the spectrum, ignoring the converging Ritz pairs there saves a lot of
   work. */ 

/* This option hasn't performed that well, so I've rescinded it from the menu. */  
/* These comments are for version 3: */
/* Uses Paige suggestion of monitoring dot product of current Lanczos vector 
   against first to sense loss of orthogonality. This can save time by avoiding
   calls to the QL or bisection routines for finding Ritz values. */

void lanczos_SO_no_time(A, n, d, y, lambda, bound, eigtol, vwsqrt, maxdeg, version)
struct vtx_data **A;		/* sparse matrix in row linked list format */
int n;				/* problem size */
int d;				/* problem dimension = number of eigvecs to find */
double **y;			/* columns of y are eigenvectors of A  */
double *lambda;			/* ritz approximation to eigenvals of A */
double *bound;			/* on ritz pair approximations to eig pairs of A */
double eigtol;			/* tolerance on eigenvectors */
double *vwsqrt;			/* square roots of vertex weights */
double maxdeg;			/* maximum degree of graph */
int version;			/* flags which version of sel. orth. to use */
{
    extern Heap   *heap;     /* pointer to heap of multigrid */
    extern double *MEM_OK;   /* variable for memory overflow exeception */
	extern int LANCZOS_SO_INTERVAL;	/* interval between orthogonalizations */
	extern int DEBUG_EVECS;		/* print debugging output? */
	extern int WARNING_EVECS;	/* print warning messages? */
	extern double WARNING_ORTHTOL;	/* Warning on modest level of orthogonality loss */
	extern double WARNING_MISTOL;	/* Warning on serious level of orthogonality loss */
	extern double WARNING_SRESTOL;	/* Warning on inaccurate recurrence for evec of T */
	extern double DOUBLE_EPSILON;	/* machine precision */
	extern double DOUBLE_MAX;	/* largest double value */
	int i,j,k;			/* indicies */
	int maxj;			/* maximum number of Lanczos iterations */
	double *u,*r;			/* Lanczos vectors */
	double *alpha, *beta;  		/* the Lanczos scalars from each step */
	double *ritz;			/* copy of alpha for ql */
	double *workj;  		/* work vector, e.g. copy of beta for ql */
	double *workn;  		/* work vector, e.g. product Av for checkeig */
	double *s;			/* eigenvector of T */
	double **q;			/* columns of q are Lanczos basis vectors */
	double *bj;			/* beta(j)*(last el. of corr. eigvec s of T) */
	double Sres;			/* how well Tevec calculated eigvec s */
	double Sres_max;		/* Max value of Sres */
	double *Ares;			/* how well Lanczos calc. eigpair lambda,y */
	int *index;			/* the Ritz index of an eigenpair */
	struct orthlink **solist;	/* vec. of structs with vecs. to orthog. against */
	struct scanlink *scanlist;	/* linked list of fields to do with min ritz vals */
	struct scanlink *curlnk;	/* for traversing the scanlist */
	double bji_tol;			/* tol on bji est. of eigen residual of A */
	int converged;			/* has the iteration converged? */
	int warning1;			/* is warning1 cond. (eigtol not achieved) true? */
	int warning2;			/* is warning2 cond. (premature orth. loss) true? */
	int warning3;			/* is warning3 cond. (suspected misconvergence) true? */
	double goodtol;			/* error tolerance for a good Ritz vector */
	int ngood;			/* total number of good Ritz pairs at current step */
	int maxngood;			/* biggest val of ngood through current step */
	int left_ngood;			/* number of good Ritz pairs on left end */
	int right_ngood;		/* number of good Ritz pairs on right end */
	int lastpause;			/* Most recent step with good ritz vecs */
	int interval;			/* number of steps between pauses */
	int left_goodlim;		/* number of ritz pairs checked on left end */
	int right_goodlim;		/* number of ritz pairs checked on right end */
	double Anorm;			/* Norm estimate of the Laplacian matrix */
	int pausemode;			/* which Lanczos pausing criterion to use */
	int pause;			/* whether to pause */
	int temp;			/* used to prevent redundant index computations */

	double *mkvec();		/* allocates space for a vector */
	double dot();			/* standard dot product routine */
	struct orthlink *makeorthlnk();	/* makes space for new entry in orthog. set */ 
	double norm();			/* vector norm */
	double Tevec();			/* calc eigenvector of T by linear recurrence*/
	double checkeig();		/* calculate residual of eigenvector of A */
	struct scanlink *mkscanlist();	/* makes initial scan list for min ritz vecs */
	int lanpause();
	void setvec(), scale(), splarax(), update(), sorthog(), exit();
	void scanmin(), frvec(), scadd(), cpvec(), orthog1(), vecran(); 
	void solistout(), doubleout(), orthogvec(), get_ritzvals();

	if (DEBUG_EVECS > 0) {
	{char buf[150]; sprintf(buf,"Selective orthogonalization Lanczos (v. %d), matrix size = %d.\n",version,n);UserWrite(buf);}
	}

	if (n < d+1) {
	{char buf[150]; sprintf(buf,"ERROR: System too small for number of eigenvalues requested.\n");UserWrite(buf);}
		exit(0);
		/* d+1 since don't use zero eigenvalue pair */
	}

	/* Allocate space. */
	u = mkvec(1,n);
    if (!MEM_OK) return;
	r = mkvec(1,n);
    if (!MEM_OK) return;
	workn = mkvec(1,n);
    if (!MEM_OK) return;
	Ares = mkvec(0,d);
    if (!MEM_OK) return;
	index = (int *) (MEM_OK = smalloc((unsigned) (d+1) * sizeof(int));
    if (!MEM_OK) return;
	maxj = n;
	alpha = mkvec(1,maxj);
    if (!MEM_OK) return;
	beta = mkvec(0,maxj);
    if (!MEM_OK) return;
	ritz = mkvec(1,maxj);
    if (!MEM_OK) return;
	s = mkvec(1,maxj);
    if (!MEM_OK) return;
	bj = mkvec(1,maxj);
    if (!MEM_OK) return;
	workj = mkvec(0,maxj);
    if (!MEM_OK) return;
	q = (double **) (MEM_OK = smalloc((unsigned) (maxj+1) * sizeof(double *)); 
    if (!MEM_OK) return;
	solist = (struct orthlink **) (MEM_OK = smalloc((unsigned) (maxj+1) * sizeof(struct orthlink *));
    if (!MEM_OK) return;
	scanlist = mkscanlist(d);
    if (!MEM_OK) return;
	ngood = 0;
	maxngood = 0;
	bji_tol = eigtol;

	/* Set some constants govering the orthogonalization heuristic. */
	Anorm = 2 * maxdeg; 			/* Gershgorin estimate for ||A|| */
	goodtol = Anorm * sqrt(DOUBLE_EPSILON);	/* Parlett & Scott's bound, p.224 of
						   "Lanczos Algs. w/ Sel. Orth.", Math. Comp. */ 
	/* The second term is empirical and conservative. */
        interval = (int) min(LANCZOS_SO_INTERVAL-1,n/(2*LANCZOS_SO_INTERVAL)) + 1; 

	if (DEBUG_EVECS > 0) {
	{char buf[150]; sprintf(buf,"  maxdeg %g\n",maxdeg);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"  goodtol %g\n",goodtol);UserWrite(buf);}
        {char buf[150]; sprintf(buf,"  interval %d\n",interval);UserWrite(buf);}
	}

	/* Initialize space. */
	vecran(r,1,n); 
	if (vwsqrt == NULL) {
		orthog1(r,1,n);
	}
	else {
		orthogvec(r,1,n,vwsqrt);
	}
	beta[0] = norm(r,1,n);
	q[0] = mkvec(1,n);
    if (!MEM_OK) return;
	setvec(q[0],1,n,0.0);
	setvec(bj,1,maxj,DOUBLE_MAX);

	/* Main Lanczos loop. */
	j = 1; 
	lastpause = 0;
	pausemode = 1;
	left_ngood = 0;
	right_ngood = 0;
	left_goodlim = 0;
	right_goodlim = 0;
	converged = FALSE;
	Sres_max = 0.0;
	while((j <= maxj) && (!converged)) {
		q[j] = mkvec(1,n);	 		
        if (!MEM_OK) return;
		scale(q[j],1,n,1.0/beta[j-1],r);	
		splarax(u, A, n, q[j], vwsqrt, workn);
        if (!MEM_OK) return;
		update(r,1,n,u,-beta[j-1],q[j-1]);
		alpha[j] = dot(r,1,n,q[j]);
		update(r,1,n,r,-alpha[j],q[j]);		
		if (vwsqrt == NULL) {
			orthog1(r,1,n);
		}
		else {
			orthogvec(r,1,n,vwsqrt);
		}
		if ( (j==(lastpause+1)) || (j== (lastpause+2)) ) {
			sorthog(r,n,solist,ngood);
		}
		beta[j] = norm(r,1,n);
		pause = lanpause(j,lastpause,interval,q,n,&pausemode,version);
		if (pause) {
			lastpause = j;			 

			/* Compute limits for checking Ritz pair convergence. */
			if (version == 1) {
				if (left_ngood + 2> left_goodlim) {
					left_goodlim = left_ngood + 2;
				}
				if (right_ngood + 3> right_goodlim) {
					right_goodlim = right_ngood + 3;
				}
			}
			if (version == 2) {
				if (left_ngood + 2> left_goodlim) {
					left_goodlim = left_ngood + 2;
				}
				right_goodlim = 0; 
			}

			/* Special case: need at least d Ritz vals on left. */
			left_goodlim = max(left_goodlim,d);

			/* Special case: can't find more than j total Ritz vals. */
			if (left_goodlim + right_goodlim > j) {
				left_goodlim = min(left_goodlim,j);
				right_goodlim = j - left_goodlim;
			}

			/* Find Ritz vals using faster of Sturm bisection or ql. */
			cpvec(ritz,1,j,alpha);
			cpvec(workj,0,j,beta);		 
			get_ritzvals(alpha,beta,j,Anorm,workj,ritz,d,
				left_goodlim,right_goodlim,eigtol);

			/* Scan for minimum evals of tridiagonal. */
			scanmin(ritz,1,j,&scanlist);

			/* Compute Ritz pair bounds at left end. */
			setvec(bj,1,j,0.0);
			for (i=1; i<=left_goodlim; i++)  { 
				Sres = Tevec(alpha,beta-1,j,ritz[i],s);
				if (Sres > Sres_max) {Sres_max = Sres;}
				bj[i] = s[j]*beta[j];
			}

			/* Compute Ritz pair bounds at right end. */
			for (i=j; i>j-right_goodlim; i--)  {
				Sres = Tevec(alpha,beta-1,j,ritz[i],s);
				if (Sres > Sres_max) {Sres_max = Sres;}
				bj[i] = s[j]*beta[j];
			}

			/* Show the portion of the spectrum checked for convergence. */
			if (DEBUG_EVECS > 2) {
			{char buf[150]; sprintf(buf,"\nindex         Ritz vals            bji bounds\n");UserWrite(buf);}
				for (i=1; i<=left_goodlim; i++) {
				{char buf[150]; sprintf(buf,"  %3d",i);UserWrite(buf);}  
					doubleout(ritz[i],1);
					doubleout(bj[i],1);
				{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
				}
			{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
				curlnk = scanlist;
				while (curlnk != NULL) {
					temp = curlnk->indx;
					if ((temp > left_goodlim)&&(temp < j-right_goodlim)) {
					{char buf[150]; sprintf(buf,"  %3d",temp);UserWrite(buf);}  
						doubleout(ritz[temp],1);
						doubleout(bj[temp],1);
					{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}  
					}
					curlnk = curlnk->pntr;
				}
			{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
				for (i=j-right_goodlim+1; i<=j; i++) {
				{char buf[150]; sprintf(buf,"  %3d",i);UserWrite(buf);}  
					doubleout(ritz[i],1);
					doubleout(bj[i],1);
				{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}  
				}
			{char buf[150]; sprintf(buf,"                            -------------------\n",goodtol);UserWrite(buf);}
			{char buf[150]; sprintf(buf,"                goodtol:    %19.16f\n\n",goodtol);UserWrite(buf);}
			}

			/* Check convergence of the predicted residual bounds. */
			converged = TRUE;
			if (j < d) converged = FALSE;
			else {
				curlnk = scanlist;
				while (curlnk != NULL) {
					if (bj[curlnk->indx] > bji_tol) {
						converged = FALSE;
					}
					curlnk = curlnk->pntr;
				}
			}

			/* Show current estimates of evals and bounds (for help in tuning) */
			if (DEBUG_EVECS > 2 && !converged) {
			        /* Collect eigenvalue and bound information for display, return. */
			        i = d;
			        curlnk = scanlist;
			        while (curlnk != NULL) {
			                lambda[i] = curlnk->val;
			                bound[i] = bj[curlnk->indx];
					index[i] = curlnk->indx;
			                curlnk = curlnk->pntr;
       		        		i--;
       		 		}
 
			        /* Compute eigenvectors and display associated info. */
			{char buf[150]; sprintf(buf,"          lambda                Ares est.               Ares          index\n");UserWrite(buf);}
			        for (i=1; i<=d; i++) {
			                Sres = Tevec(alpha,beta-1,j,lambda[i],s);
					if (Sres > Sres_max) {Sres_max = Sres;}
			                setvec(y[i],1,n,0.0);
			                for (k = 1; k <= j; k++) {
			                        scadd(y[i],1,n,s[k],q[k]);
			                }
					Ares[i] = checkeig(workn,A,y[i],n,lambda[i],vwsqrt,u);
                    if (!MEM_OK) return;
				{char buf[150]; sprintf(buf,"%2d.",i);UserWrite(buf);}
					doubleout(lambda[i],1);
					doubleout(bound[i],1);
					doubleout(Ares[i],1);
				{char buf[150]; sprintf(buf,"   %3d\n",index[i]);UserWrite(buf);}
			        }
			{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
			}

			/* Could check convergence of actual residual bounds against 
                           predicted bounds here and possibly fix-up if orthogonality 
                           is eroding. But how to fix things up without starting over? */

			if (!converged) {	
				ngood = 0;
				left_ngood = 0;		/* for setting left_goodlim on next loop */
				right_ngood = 0;	/* for setting right_goodlim on next loop */

				/* Compute converged Ritz pairs on left end */
				for (i=1; i<=left_goodlim; i++)  {
					if (bj[i] <= goodtol) {
						ngood += 1;
						left_ngood += 1;
						if (ngood > maxngood) {
							maxngood = ngood;
							solist[ngood] = makeorthlnk();
                            if (!MEM_OK) return;
							(solist[ngood])->vec = mkvec(1,n);
                            if (!MEM_OK) return;
						}
						(solist[ngood])->index = i;
						Sres = Tevec(alpha,beta-1,j,ritz[i],s);
						if (Sres > Sres_max) {Sres_max = Sres;}
						setvec((solist[ngood])->vec,1,n,0.0);
						for (k = 1; k <= j; k++) {
							scadd((solist[ngood])->vec,1,n,s[k],q[k]);
						}
					}
				}

				/* Compute converged Ritz pairs on right end */
				for (i=j; i>j-right_goodlim; i--)  {
					if (bj[i] <= goodtol) {
						ngood += 1;
						right_ngood += 1;
						if (ngood > maxngood) {
							maxngood = ngood;
							solist[ngood] = makeorthlnk();
                            if (!MEM_OK) return;
							(solist[ngood])->vec = mkvec(1,n);
                            if (!MEM_OK) return;
						}
						(solist[ngood])->index = i;
						Sres = Tevec(alpha,beta-1,j,ritz[i],s);
						if (Sres > Sres_max) {Sres_max = Sres;}
						setvec((solist[ngood])->vec,1,n,0.0);
						for (k = 1; k <= j; k++) {
							scadd((solist[ngood])->vec,1,n,s[k],q[k]);
						}
					}
				}

				if (DEBUG_EVECS > 2) {
				{char buf[150]; sprintf(buf,"  j %3d; goodlim lft %2d, rgt %2d; list ", j,left_goodlim,right_goodlim);UserWrite(buf);}
					solistout(solist,n,ngood,j);
				}
			}
		}
		j++;
	}
	j--;

	/* Collect eigenvalue and bound information. Only compute and display info for 
	   the eigpairs actually used in the partitioning since don't want to spend the 
	   time or space to compute the null-space of the Laplacian. */
	i = d; 
	curlnk = scanlist;
	while (curlnk != NULL) {
		lambda[i] = curlnk->val;
		bound[i] = bj[curlnk->indx];
		index[i] = curlnk->indx;
		curlnk = curlnk->pntr;
		i--;
	}

	/* Compute eigenvectors. */
	for (i=1; i<=d; i++) {
		Sres = Tevec(alpha,beta-1,j,lambda[i],s);
		if (Sres > Sres_max) {Sres_max = Sres;}
		setvec(y[i],1,n,0.0);
		for (k = 1; k <= j; k++) {
			scadd(y[i],1,n,s[k],q[k]);
		}
	}

	if (DEBUG_EVECS > 0 || WARNING_EVECS > 0) {
		for (i=1; i<=d; i++) {
			Ares[i] = checkeig(workn,A,y[i],n,lambda[i],vwsqrt,u);
            if (!MEM_OK) return;
		}
	}

	if (DEBUG_EVECS > 0) {
	   {char buf[150]; sprintf(buf,"Lanczos itns. = %d\n",j);UserWrite(buf);}
	{char buf[150]; sprintf(buf,"          lambda                Ares est.               Ares          index\n");UserWrite(buf);}
		for (i=1; i<=d; i++) {
		{char buf[150]; sprintf(buf,"%2d.",i);UserWrite(buf);}
			doubleout(lambda[i],1);
			doubleout(bound[i],1);
			doubleout(Ares[i],1);
		{char buf[150]; sprintf(buf,"   %3d\n",index[i]);UserWrite(buf);}
		}
	{char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
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
		}
		if (warning1) {
		{char buf[150]; sprintf(buf,"WARNING: Eigen pair tolerance (%g) not acheived.\n",eigtol);UserWrite(buf);}
		}
		if (warning2 && ! warning3) {
		{char buf[150]; sprintf(buf,"WARNING: Minor loss of orthogonality (Ares/est. > %g).\n", WARNING_ORTHTOL);UserWrite(buf);}
		}
		if (warning3) {
		{char buf[150]; sprintf(buf,"WARNING: Substantial loss of orthogonality (Ares/est. > %g).\n", WARNING_MISTOL);UserWrite(buf);}
		}
		if (j == maxj) {
		{char buf[150]; sprintf(buf,"WARNING: Maximum number of Lanczos iterations reached.\n");UserWrite(buf);}
		}
	}

	if (WARNING_EVECS > 1) {
		if (warning1 || warning2 || warning3) {
			if (DEBUG_EVECS <= 0) {
			{char buf[150]; sprintf(buf,"          lambda                Ares est.               Ares          index\n");UserWrite(buf);}
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
	}
	if (WARNING_EVECS > 2) {
		if (Sres_max > WARNING_SRESTOL) {
		{char buf[150]; sprintf(buf,"WARNING: High maximum residual in computing eigenvector of T %g.\n", Sres_max);UserWrite(buf);}
		}
	}

	/* free up memory */
	frvec(u,1);
	frvec(r,1);
	frvec(workn,1);
	frvec(Ares,0);
	sfree((char *) index);
	frvec(alpha,1);
	frvec(beta,0);
	frvec(ritz,1);
	frvec(s,1);
	frvec(bj,1);
	frvec(workj,0);
	for (i=0; i<=j; i++)  {
		frvec(q[i], 1);
	}
	sfree((char *)q);	
	while(scanlist != NULL) {
		curlnk = scanlist->pntr;
		sfree((char *) scanlist);
		scanlist = curlnk;
	}
	for (i=1; i<=maxngood; i++) {
	   frvec((solist[i])->vec, 1);
	   sfree((char *) solist[i]);
	}
	sfree((char *) solist);
}
