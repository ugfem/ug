// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	<string.h>
#include	"../main/defs.h"
#include	"../main/params.h"


void input_geom(fingeom, geomname, nvtxs, igeom, x, y, z)
FILE *fingeom;			/* geometry input file */
char *geomname;			/* name of geometry file */
int nvtxs;			/* number of coordinates to read */
int *igeom;			/* dimensionality of geometry */
float **x, **y, **z;		/* coordiates of vertices */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_INPUT;	/* echo that read was successful? */
   float xc, yc, zc;		/* first x, y, z coordinate */
   float temp;			/* temporary value */
   char line[LINE_LENGTH+1];	/* line of date from file */
   int nread;			/* number of coordinates read */
   int flag;			/* Have I skipped past blank lines? */
   int nvals;			/* number of values in an input line */
   int i;			/* loop counter */
   void exit();

   flag = FALSE;
   while (!flag && fgets(line, LINE_LENGTH, fingeom) != NULL) {
      nvals = sscanf(line, "%f%f%f", &xc, &yc, &zc, &temp);
      if (nvals != 0) flag = TRUE;
   }

   if (nvals > 3) {
     {char buf[150]; sprintf(buf," Too many values on geometry file input line\n");UserWrite(buf);}
      exit(1);
   }

   *igeom = nvals;

   *y = *z = NULL;

   *x = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
   if (!MEM_OK) return;
   (*x)[1] = xc;
   if (nvals > 1) {
      *y = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
      if (!MEM_OK) return;
      (*y)[1] = yc;
   }
   if (nvals > 2) {
      *z = (float *) (MEM_OK = smalloc((unsigned) (nvtxs+1)*sizeof(float));
      if (!MEM_OK) return;
      (*z)[1] = zc;
   }

   for (nread=2; nread<=nvtxs; nread++) {
      if (nvals == 1) {
	 i = fscanf(fingeom, "%f", &((*x)[nread]));
      }
      else if (nvals == 2) {
	 i = fscanf(fingeom, "%f%f", &((*x)[nread]), &((*y)[nread]));
      }
      else if (nvals == 3) {
	 i = fscanf(fingeom, "%f%f%f", &((*x)[nread]), &((*y)[nread]),
		    &((*z)[nread]));
      }

      if (i == EOF) {
	{char buf[150]; sprintf(buf,"Too few values in geometry file; nvtxs=%d, but nread=%d\n", nvtxs, nread);UserWrite(buf);}
	 exit(1);
      }
      else if (i != nvals) {
	{char buf[150]; sprintf(buf,"Wrong number of values in line of geometry file %s.\n", geomname);UserWrite(buf);}
	 exit(1);
      }
   }

   fclose(fingeom);

   if (DEBUG_INPUT > 0) {
     {char buf[150]; sprintf(buf,"Finished reading geometry file %s; dimension = %d.\n", geomname, nvals);UserWrite(buf);}
   }
}
