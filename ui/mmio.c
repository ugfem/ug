// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *   Matrix Market I/O library for ANSI C
 *
 *   See http://math.nist.gov/MatrixMarket for details.
 *
 *
 */

#include <stdio.h>
#include <string.h>
#ifndef __MACOSXSERVER__
#ifndef __MWCW__
#include <malloc.h>
#endif
#endif

#include <ctype.h>

#include "mmio.h"

#include "algebra.h"
#include "cmdline.h"  /* for CreateCommand */
#include "commands.h"
#include "debug.h"
#include "gm.h"
#include "npscan.h"
#include "scan.h"
#include "udm.h"

REP_ERR_FILE;

/* MacOS doesn't support the non-standard (!) strdup function */
#ifdef __MWCW__
char *strdup(char *text);  /* forward declaration to make ANSI compilers happy */

char *strdup(const char *text)
{
  char *c;
  c = (char *)malloc( strlen(text)+1 );
  return c;
}
#endif


int mm_is_valid(MM_typecode matcode)
{
  if (!mm_is_matrix(matcode)) return 0;
  if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
  if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
  if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
                                 mm_is_skew(matcode))) return 0;
  return 1;
}

int mm_read_banner(FILE *f, MM_typecode *matcode)
{
  char line[MM_MAX_LINE_LENGTH];
  char banner[MM_MAX_TOKEN_LENGTH];
  char mtx[MM_MAX_TOKEN_LENGTH];
  char crd[MM_MAX_TOKEN_LENGTH];
  char data_type[MM_MAX_TOKEN_LENGTH];
  char storage_scheme[MM_MAX_TOKEN_LENGTH];
  char *p;

  mm_clear_typecode(matcode);

  if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
    return MM_PREMATURE_EOF;

  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
             storage_scheme) != 5)
    return MM_PREMATURE_EOF;

  for (p=mtx; *p!='\0'; *p=tolower(*p),p++) ;   /* convert to lower case */
  for (p=crd; *p!='\0'; *p=tolower(*p),p++) ;
  for (p=data_type; *p!='\0'; *p=tolower(*p),p++) ;
  for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++) ;

  /* check for banner */
  if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
    return MM_NO_HEADER;

  /* first field should be "mtx" */
  if (strcmp(mtx, MM_MTX_STR) != 0)
    return MM_UNSUPPORTED_TYPE;
  mm_set_matrix(matcode);


  /* second field describes whether this is a sparse matrix (in coordinate
          storgae) or a dense array */


  if (strcmp(crd, MM_SPARSE_STR) == 0)
    mm_set_sparse(matcode);
  else
  if (strcmp(crd, MM_DENSE_STR) == 0)
    mm_set_dense(matcode);
  else
    return MM_UNSUPPORTED_TYPE;


  /* third field */

  if (strcmp(data_type, MM_REAL_STR) == 0)
    mm_set_real(matcode);
  else
  if (strcmp(data_type, MM_COMPLEX_STR) == 0)
    mm_set_complex(matcode);
  else
  if (strcmp(data_type, MM_PATTERN_STR) == 0)
    mm_set_pattern(matcode);
  else
  if (strcmp(data_type, MM_INT_STR) == 0)
    mm_set_integer(matcode);
  else
    return MM_UNSUPPORTED_TYPE;


  /* fourth field */

  if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
    mm_set_general(matcode);
  else
  if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
    mm_set_symmetric(matcode);
  else
  if (strcmp(storage_scheme, MM_HERM_STR) == 0)
    mm_set_hermitian(matcode);
  else
  if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
    mm_set_skew(matcode);
  else
    return MM_UNSUPPORTED_TYPE;


  return 0;
}

int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
  if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
    return MM_COULD_NOT_WRITE_FILE;
  else
    return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

  /* set return null parameter values, in case we exit with errors */
  *M = *N = *nz = 0;

  /* now continue scanning until you reach the end-of-comments */
  do
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
      return MM_PREMATURE_EOF;
  } while (line[0] == '%');

  /* line[] is either blank or has M,N, nz */
  if (sscanf(line, "%d %d %d", M, N, nz) == 3)
    return 0;

  else
    do
    {
      num_items_read = fscanf(f, "%d %d %d", M, N, nz);
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);

  return 0;
}


int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

  /* set return null parameter values, in case we exit with errors */
  *M = *N = 0;

  /* now continue scanning until you reach the end-of-comments */
  do
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
      return MM_PREMATURE_EOF;
  } while (line[0] == '%');

  /* line[] is either blank or has M,N, nz */
  if (sscanf(line, "%d %d %d", M, N) == 2)
    return 0;

  else   /* we have a blank line */
    do
    {
      num_items_read = fscanf(f, "%d %d %d", M, N);
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);

  return 0;
}

int mm_write_mtx_array_size(FILE *f, int M, int N)
{
  if (fprintf(f, "%d %d\n", M, N) != 2)
    return MM_COULD_NOT_WRITE_FILE;
  else
    return 0;
}



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
                         double val[], MM_typecode matcode)
{
  int i;
  if (mm_is_complex(matcode))
  {
    for (i=0; i<nz; i++)
      if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
          != 4) return MM_PREMATURE_EOF;
  }
  else if (mm_is_real(matcode))
  {
    for (i=0; i<nz; i++)
    {
      if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
          != 3) return MM_PREMATURE_EOF;

    }
  }

  else if (mm_is_pattern(matcode))
  {
    for (i=0; i<nz; i++)
      if (fscanf(f, "%d %d", &I[i], &J[i])
          != 2) return MM_PREMATURE_EOF;
  }
  else
    return MM_UNSUPPORTED_TYPE;

  return 0;

}

int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
                          double *real, double *imag, MM_typecode matcode)
{
  if (mm_is_complex(matcode))
  {
    if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
        != 4) return MM_PREMATURE_EOF;
  }
  else if (mm_is_real(matcode))
  {
    if (fscanf(f, "%d %d %lg\n", I, J, real)
        != 3) return MM_PREMATURE_EOF;

  }

  else if (mm_is_pattern(matcode))
  {
    if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
  }
  else
    return MM_UNSUPPORTED_TYPE;

  return 0;

}


/************************************************************************
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J,
                    double **val, MM_typecode *matcode)
{
  int ret_code;
  FILE *f;

  if (strcmp(fname, "stdin") == 0) f=stdin;
  else
  if ((f = fopen(fname, "r")) == NULL)
    return MM_COULD_NOT_READ_FILE;


  if ((ret_code = mm_read_banner(f, matcode)) != 0)
    return ret_code;

  if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
        mm_is_matrix(*matcode)))
    return MM_UNSUPPORTED_TYPE;

  if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
    return ret_code;


  *I = (int *)  malloc(*nz * sizeof(int));
  *J = (int *)  malloc(*nz * sizeof(int));
  *val = NULL;

  if (mm_is_complex(*matcode))
  {
    *val = (double *) malloc(*nz * 2 * sizeof(double));
    ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                    *matcode);
    if (ret_code != 0) return ret_code;
  }
  else if (mm_is_real(*matcode))
  {
    *val = (double *) malloc(*nz * sizeof(double));
    ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                    *matcode);
    if (ret_code != 0) return ret_code;
  }

  else if (mm_is_pattern(*matcode))
  {
    ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                    *matcode);
    if (ret_code != 0) return ret_code;
  }

  if (f != stdin) fclose(f);
  return 0;
}

int mm_write_banner(FILE *f, MM_typecode matcode)
{
  char *str = mm_typecode_to_str(matcode);
  int ret_code;

  ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
  free(str);
  if (ret_code !=2 )
    return MM_COULD_NOT_WRITE_FILE;
  else
    return 0;
}

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
                     double val[], MM_typecode matcode)
{
  FILE *f;
  int i;

  if (strcmp(fname, "stdout") == 0)
    f = stdout;
  else
  if ((f = fopen(fname, "w")) == NULL)
    return MM_COULD_NOT_WRITE_FILE;

  /* print banner followed by typecode */
  fprintf(f, "%s ", MatrixMarketBanner);
  fprintf(f, "%s\n", mm_typecode_to_str(matcode));

  /* print matrix sizes and nonzeros */
  fprintf(f, "%d %d %d\n", M, N, nz);

  /* print values */
  if (mm_is_pattern(matcode))
    for (i=0; i<nz; i++)
      fprintf(f, "%d %d\n", I[i], J[i]);
  else
  if (mm_is_real(matcode))
    for (i=0; i<nz; i++)
      fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
  else
  if (mm_is_complex(matcode))
    for (i=0; i<nz; i++)
      fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i],
              val[2*i+1]);
  else
  {
    if (f != stdout) fclose(f);
    return MM_UNSUPPORTED_TYPE;
  }

  if (f !=stdout) fclose(f);

  return 0;
}


char  *mm_typecode_to_str(MM_typecode matcode)
{
  char buffer[MM_MAX_LINE_LENGTH];
  char *types[4];
  int error =0;

  /* check for MTX type */
  if (mm_is_matrix(matcode))
    types[0] = MM_MTX_STR;
  else
    error=1;

  /* check for CRD or ARR matrix */
  if (mm_is_sparse(matcode))
    types[1] = MM_SPARSE_STR;
  else
  if (mm_is_dense(matcode))
    types[1] = MM_DENSE_STR;
  else
    return NULL;

  /* check for element data type */
  if (mm_is_real(matcode))
    types[2] = MM_REAL_STR;
  else
  if (mm_is_complex(matcode))
    types[2] = MM_COMPLEX_STR;
  else
  if (mm_is_pattern(matcode))
    types[2] = MM_PATTERN_STR;
  else
  if (mm_is_integer(matcode))
    types[2] = MM_INT_STR;
  else
    return NULL;


  /* check for symmetry type */
  if (mm_is_general(matcode))
    types[3] = MM_GENERAL_STR;
  else
  if (mm_is_symmetric(matcode))
    types[3] = MM_SYMM_STR;
  else
  if (mm_is_hermitian(matcode))
    types[3] = MM_HERM_STR;
  else
  if (mm_is_skew(matcode))
    types[3] = MM_SKEW_STR;
  else
    return NULL;

  sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
  return strdup(buffer);

}


/****************************************************************************/
/*D
   readMM - reads a matrix in MatrixMarket format

   DESCRIPTION:
   Reads a matrix in MatrixMarket format.
   (ug interface to mmio.c from http://math.nist.gov/MatrixMarket).

   'readMM <name> $A MATDATA_DESC [$blocked]'

   .  <name> - file name

   KEYWORDS:
   data, io

   SEE ALSO:
   D*/
/****************************************************************************/

static ReadMMCommand(INT argc, char **argv)
{
  MULTIGRID *theMG;
  HEAP *theHeap;
  GRID *theGrid;
  NODE *AdamNode;
  VECTOR **vec_array;
  MATRIX *mat, *mat1;
  MATDATA_DESC *theMD;

  MM_typecode matcode;
  FILE *stream;
  char filename[NAMESIZE];

  INT comp, i, k, M, N, nz, nvec, dof, dofsq, err_code;
  int from, to, from_vec, to_vec;
  INT blocked;
  INT MarkKey;
  double value;
  double *mptr;

  if ( (theMG=GetCurrentMultigrid())==NULL )
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: no current multi grid */

  if (TOPLEVEL(theMG)!=0)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: MG not empty */

  theGrid = GRID_ON_LEVEL(theMG, 0);
  AdamNode = FIRSTNODE(theGrid);
  if (AdamNode==NULL)
    REP_ERR_RETURN(CMDERRORCODE);              /* Error: MG has no nodes */

  /* read filename: */
  if (sscanf(argv[0],expandfmt(CONCAT3("readMM %",NAMELENSTR,"[ -~]")),filename)!=1)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: could not read name of the file */

  if ((stream = fopen(filename, "r")) == NULL)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: Could not open file */

  theMD = ReadArgvMatDesc(theMG,"A",argc,argv);
  if (theMD == NULL)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: could not read mat-data descriptor */

  dof = MD_ROWS_IN_MTYPE(theMD, 0);
  if (dof != MD_COLS_IN_MTYPE(theMD, 0))
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: md should be symmetric */

  dofsq = dof*dof;

  if ((MD_SUCC_COMP(theMD))==0)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: cannot handle sparse formats */

  comp  = MD_MCMP_OF_MTYPE(theMD, 0, 0);

  blocked = 0;
  if (ReadArgvOption("blocked",argc,argv)==1)
    blocked = 1;

  if (mm_read_banner(stream, &matcode) != 0)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: Could not process Matrix Market banner. */

  /* get size of matrix */
  if (mm_read_mtx_crd_size(stream, &M, &N, &nz) !=0)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: Could not read size of sparse matrix! */

  if (M!=N)
    REP_ERR_RETURN(PARAMERRORCODE);              /* Error: Matrix must be symmetric. */

  if (blocked)
    nvec = M;
  else
  {
    if (M%dof!=0)
      REP_ERR_RETURN(PARAMERRORCODE);                    /* Error: nr of vects wrong */
    nvec = M/dof;
  }

  theHeap=MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);

  /* reserve memory for array */
  vec_array = (VECTOR **) GetTmpMem(theHeap,sizeof(VECTOR *) * nvec, MarkKey);
  if (vec_array==NULL)
    REP_ERR_RETURN(CMDERRORCODE);              /* Error: Could not allocate array. */

  /* allocate nvec vectors */
  for (i=0; i<nvec; i++)
  {
    err_code = CreateVector (theGrid, NODEVEC, (GEOM_OBJECT *) AdamNode, vec_array+i);
    if (err_code)
      REP_ERR_GOTO(, error);                    /* Error: could not create vector on new level. */

    VINDEX(vec_array[i]) = i;
    SETNEW_DEFECT(vec_array[i],1);
    SETFINE_GRID_DOF(vec_array[i],0);
  }

  /* read in matrix */
  if (blocked)
  {
    for (i=0; i<nz; i++)
    {
      fscanf(stream, "%d %d", &from, &to);
      from--; to--;                    /* to C-style */
      if (from<0 || from>=M || to<0 || to>=M)
        REP_ERR_GOTO(, error);                          /* Error: index out of range. */

      mat = GetMatrix(vec_array[from], vec_array[to]);
      if (mat==NULL)
      {
        mat = (MATRIX *) CreateConnection(theGrid, vec_array[from], vec_array[to]);
        if (mat == NULL)
          REP_ERR_GOTO(, error);                                /* Error: Could not create connection. */

        /* here we have only to clear the transpose,
           since it might not be set */
        if (!MDIAG(mat))
        {
          mat1 = CMATRIX1(mat);
          for (k=0; k<dofsq; k++)
            MVALUE(mat1, comp+k) = 0.0;
        }
      }

      mptr = MVALUEPTR(mat, comp);
      for (k=0; k<dofsq; k++)
        fscanf(stream, "%lg", mptr+k);

      fscanf(stream, "\n");
    }
  }
  else
  {
    for (i=0; i<nz; i++)
    {
      fscanf(stream, "%d %d %lg \n", &from, &to, &value);
      from--; to--;                    /* to C-style */
      if (from<0 || from>=M || to<0 || to>=M)
        REP_ERR_GOTO(, error);                          /* Error: index out of range. */

      from_vec = from/dof; to_vec = to/dof;
      mat = GetMatrix(vec_array[from_vec], vec_array[to_vec]);
      if (mat==NULL)
      {
        mat = (MATRIX *) CreateConnection(theGrid, vec_array[from_vec], vec_array[to_vec]);
        if (mat == NULL)
          REP_ERR_GOTO(, error);                                /* Error: Could not create connection. */

        /* we have to clear the newly generated matrices */
        for (k=0; k<dofsq; k++)
          MVALUE(mat, comp+k) = 0.0;
        if (!MDIAG(mat))
        {
          mat1 = CMATRIX1(mat);
          for (k=0; k<dofsq; k++)
            MVALUE(mat1, comp+k) = 0.0;
        }
      }

      MVALUE(mat, comp + (from%dof)*dof + (to%dof) ) = value;
    }
  }
  fclose(stream);
  ReleaseTmpMem(theHeap,MarkKey);
  return 0;

error:
  fclose(stream);
  ReleaseTmpMem(theHeap,MarkKey);
  REP_ERR_RETURN(CMDERRORCODE);
}


int InitMMIO (void)
{
  if (CreateCommand("readMM",     ReadMMCommand)==NULL)
    return (__LINE__);

  return(0);
}
