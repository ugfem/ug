// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      blasm.c                                                       */
/*                                                                          */
/* Purpose:   basic linear algebra routines                                 */
/*            working on sparse matrix-vector structure                     */
/*                                                                          */
/* Author:    Nicolas Neuss                                                 */
/*            email: Nicolas.Neuss@IWR.Uni-Heidelberg.De                    */
/*                                                                          */
/* History:   02.01.98 begin sparse matrix routines                         */
/*            20.01.98 end of implementation phase                          */
/*            28.01.98 scalar matrix operations work                        */
/*            xx.02.98 matrix operations work                               */
/*                                                                          */
/* Note:      The files sm.[ch], blasm.[ch] may be obtained also in a       */
/*            standalone form under the GNU General Public License.         */
/*            The use inside of UG under the actual UG license              */
/*            is allowed.                                                   */
/*                                  HD, 13.7.99,  Nicolas Neuss.            */
/*                                                                          */
/****************************************************************************/

#include <stddef.h>

#ifdef _2
#define __UG__
#endif

#ifdef _3
#define __UG__
#endif

#ifndef __UG__
   #include "blasm.h"
   #include "sm.h"
   #include "udm.h"

   #define ERR_FILE_ID 7
   #define MGFORMAT(mg) mg->format
   #define VCFLAG(vect) FLAG(vect)
   #define SETVCFLAG(vect,n) SETFLAG(vect,n)
#else /* __UG__ */

#include "gm.h"
#include "sm.h"
#include "blasm.h"
#include "udm.h"
#include "debug.h"

#define ERR_RETURN(x) REP_ERR_RETURN(x)
REP_ERR_FILE;

#define BOTTOM_GRID(mg) GRID_ON_LEVEL(mg,BOTTOMLEVEL(mg))
#define FINER(grid) UPGRID(grid)

#define GORDERED(grid) GSTATUS(grid,GSTATUS_ORDERED)
#define SETGORDERED(grid,n) grid->status=(grid->status&~GSTATUS_ORDERED)+n*GSTATUS_ORDERED

#define FINE_GRID_DOF_MASK (1<<FINE_GRID_DOF_SHIFT)
#define NEW_DEFECT_MASK (1<<NEW_DEFECT_SHIFT)
#define VACTIVE_MASK (1<<VACTIVE_SHIFT)
#define VTYPE_MASK (((1<<VTYPE_LEN)-1)<<VTYPE_SHIFT)

#define MDIAG_MASK (1<<MDIAG_SHIFT)
#define MLOWER_MASK (1<<MLOWER_SHIFT)
#define MUPPER_MASK (1<<MUPPER_SHIFT)
#define MACTIVE_MASK (1<<MACTIVE_SHIFT)

/****************************************************************************/
/*  this routine is an import from sm/algebra.c                             */
/****************************************************************************/
/*D
   Mark_and_Sort_Matrix - mark and sort the matrix graph

   SYNOPSIS:
   INT Mark_and_Sort_Matrix (GRID *grid, int operation)

   DESCRIPTION:
   These routines marks the links in the matrix graph as lower or upper.
   If operation==1, it additionally sorts the links (not yet used and tested).
   It should be called by all routines that are changing the order of vectors
   and/or matrices (if the correct setting can't be done locally)!
   D*/
/****************************************************************************/

INT Mark_and_Sort_Matrix (GRID *grid, int operation)
{
  int i;
  VECTOR *vect, *vect2;
  MATRIX *mat, *matD;

  /* first reset FLAG bits (and set index field for compatibility) */
  i = 1;
  for (vect=FIRSTVECTOR(grid); vect!=NULL; vect=SUCCVC(vect))
  {
    SETVCFLAG(vect,0);
    vect->index = i++;
  }

  /* then loop through the vector list once more sorting the matrices */
  for (vect=FIRSTVECTOR(grid); vect!=NULL; vect=SUCCVC(vect))
  {
    /* Mark the vector as visited, i.e. it is in the lower triangular part */
    /* of its neighbors. This may be done at the beginning, */
    /* since the diagonal is handled separately. */
    SETVCFLAG(vect,1);

    matD = VSTART(vect);
    if (matD==NULL)
      continue;
    /* ERR_RETURN(-1);  // Error (internal): cannot handle matrix without diagonal */
    if (MDEST(matD)!=vect)
      ERR_RETURN(-1);                    /* Error (internal): diagonal matrix not at beginning of list */

    /* further we mark vect as inactive, if class(vect)<ACTIVE_CLASS */
    /* matrices are marked inactive, if they point to inactive vectors */
    if (VCLASS(vect)<ACTIVE_CLASS)
    {
      SETVACTIVE(vect,0);
      SETMACTIVE(matD,0);
    }
    else
    {
      SETVACTIVE(vect,1);
      SETMACTIVE(matD,1);
    }

    /* set diagonal to be neither lower or upper */
    SETMLOWER(matD,0);
    SETMUPPER(matD,0);

    switch (operation)
    {
    case 0 :            /* without matrix reordering */
      for (mat=MNEXT(matD); mat!=NULL; mat=MNEXT(mat))
      {
        vect2 = MDEST(mat);
        if (vect2==NULL)
          ERR_RETURN(-1);                                /* Error (internal): matrix has no destination vector */
        if (vect2==vect)
          ERR_RETURN(-1);                                /* Error (internal): two diagonal matrices */

        /* matrix is inactive if it points to an inactive vector, see above */
        if (VCLASS(vect2)<ACTIVE_CLASS)
          SETMACTIVE(mat,0);
        else
          SETMACTIVE(mat,1);

        if (VCFLAG(vect2))
        {
          SETMLOWER(mat,1);
          SETMUPPER(mat,0);
        }
        else
        {
          SETMLOWER(mat,0);
          SETMUPPER(mat,1);
        }
      }
      break;

    case 1 :             /* with matrix reordering */
    {
      MATRIX *next_mat, *last_mat, *L_first, *L_last, *U_first, *U_last;
      L_first = L_last = NULL;
      U_first = U_last = NULL;
      mat = MNEXT(matD);
      while (mat!=NULL)
      {
        next_mat = MNEXT(mat);
        vect2 = MDEST(mat);
        if (vect2==NULL)
          ERR_RETURN(-1);                                /* Error (internal): matrix has no destination vector */
        if (vect2==vect)
          ERR_RETURN(-1);                                /* Error (internal): two diagonal matrices */

        /* matrix is inactive if it points to an inactive vector, see above */
        if (VCLASS(vect2)<ACTIVE_CLASS)
          SETMACTIVE(mat,0);
        else
          SETMACTIVE(mat,1);

        if (VCFLAG(vect2))
        {
          /* dest vector comes before vect */
          MNEXT(mat) = L_first;
          L_first = mat;
          if (L_last==NULL)
            L_last=mat;
          SETMLOWER(mat,1);
          SETMUPPER(mat,0);
        }
        else
        {
          /* vector comes after vect */
          MNEXT(mat) = U_first;
          U_first = mat;
          if (U_last==NULL)
            U_last=mat;
          SETMLOWER(mat,0);
          SETMUPPER(mat,1);
        }

        mat = next_mat;
      }

      /* then couple the lists in the prescribed order */
      last_mat = matD;
      switch (operation)
      {
      case 0 :                   /* first lower, then upper part */
        if (L_first!=NULL)
        {
          MNEXT(last_mat) = L_first;
          last_mat = L_last;
        }
        if (U_first!=NULL)
        {
          MNEXT(last_mat) = U_first;
          last_mat = U_last;
        }
        break;

      case 1 :                   /* first upper, then lower part */
        if (U_first!=NULL)
        {
          MNEXT(last_mat) = U_first;
          last_mat = U_last;
        }
        if (L_first!=NULL)
        {
          MNEXT(last_mat) = L_first;
          last_mat = L_last;
        }
        break;
      default :
        ERR_RETURN(-1);                          /* Error (internal): operation not implemented */
      }
      break;
    }
    }
  }

  SETGORDERED(grid,1);

  return(0);
}

/****************************************************************************/
/* the following routines are imports from sm/udm.c                         */
/****************************************************************************/
/* VV_Compatible, MM_Compatible, MV_Compatible                              */
/****************************************************************************/

static INT MV_Compatible (const FORMAT *format, const MATDATA_DESC *A, const VECDATA_DESC *y)
{
  int type, vtype;

  if (format==NULL)
    ERR_RETURN(-1);              /* Error: check to make compiler silent (format now used) */

  for (type=0; type<NMATTYPES; type++)
  {
    SPARSE_MATRIX *sm = MD_SM(A,type);
    if (sm != NULL)
    {
      vtype = MTYPE_CT(type);
      if (sm->ncols!=VD_NCMPS_IN_TYPE(y,vtype))
        return (1);                          /* Result: A can't be applied to y */
    }
  }

  return (0);
}

static INT VM_Compatible (const FORMAT *format, const VECDATA_DESC *x, const MATDATA_DESC *A)
{
  int type, vtype;

  if (format==NULL)
    ERR_RETURN(-1);              /* Error: check to make compiler silent (format now used) */

  for (type=0; type<NMATTYPES; type++)
  {
    SPARSE_MATRIX *sm = MD_SM(A,type);
    if (sm != NULL)
    {
      vtype = MTYPE_RT(type);
      if (sm->nrows!=VD_NCMPS_IN_TYPE(x,vtype))
        return (1);                          /* Result: x is not row compatible with A */
    }
  }

  return (0);
}

static INT MM_Compatible (const FORMAT *format, const MATDATA_DESC *A, const MATDATA_DESC *B)
{
  int type;

  if (format==NULL)
    ERR_RETURN(-1);              /* Error: check to make compiler silent (format now used) */

  for (type=0; type<NMATTYPES; type++)
    if (MD_SM(A,type) != NULL)
    {
      if (MD_SM(B,type)==NULL)
        return (1);
      else
        return(SM_Compare(MD_SM(A,type), MD_SM(B,type)));
    }
    else
    if (MD_SM(B,type)!=NULL)
      return (1);

  return (0);
}

#endif /* __UG__ */

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   Matrix_Loop_? - loop through the grid performing Block_Loop_?

   SYNOPSIS:
   static INT Matrix_Loop_? ( parameters )

   DESCRIPTION:
   These routines loop through the matrix graph. For each link they perform
   a sparse block operation Block_Loop_?, which is included as a macro.
   D*/
/****************************************************************************/

#define Block_Loop_M(operation,additional_increment)       \
  operation;                                                 \
  for (j=nm1; j!=0; j--)                                     \
  {                                                          \
    M_ptr = (DOUBLE *) ((char *) M_ptr + *M_diff_ptr++);   \
    additional_increment;                                  \
    operation;                                             \
  }

static INT Matrix_Loop_M(GRID *grid, unsigned int _vroot_mask_, unsigned int _vroot_value_,
                         unsigned int _mat_mask_, unsigned int _mat_value_,
                         int n, int n_D,
                         ptrdiff_t M_start_off, const ptrdiff_t *M_diff,
                         ptrdiff_t M_start_off_D, const ptrdiff_t *M_diff_D,
                         int mode, const DOUBLE *value, const DOUBLE *value_D )
{
  /* most important in registers */
  int j;
  DOUBLE *M_ptr;
  const DOUBLE *v_ptr;
  const ptrdiff_t *M_diff_ptr;

  /* next important in registers */
  int nm1;
  int operation;
  MATRIX *mat;
  unsigned int mat_mask  = _mat_mask_;
  unsigned int mat_value = _mat_value_;
  VECTOR *vect;
  VECTOR *vect2;
  const DOUBLE *v_start_ptr = value;
  ptrdiff_t M_start_offset = M_start_off;
  const ptrdiff_t *M_diff_start = M_diff;

  /* not important in registers */
  unsigned int vroot_mask = _vroot_mask_;
  unsigned int vroot_value = _vroot_value_;

  operation = (mode&BLAS_OP_MASK)>>BLAS_OP_SHIFT;

  for (vect=FIRSTVECTOR(grid); vect!=NULL; vect=SUCCVC(vect))
    if ((CTRL(vect)&vroot_mask) == vroot_value)
    {
      mat = VSTART(vect);
      if (n_D != 0)
      {
        if ((CTRL(mat)&mat_mask) == mat_value)
        {
          /* handle the diagonal matrix (always at beginning!) */
          nm1 = n_D-1;
          vect2 = MDEST(mat);
          if (vect2!=vect)
            ERR_RETURN (-1);                              /* Error: matrix has no diagonal connection */

          M_ptr = (DOUBLE *) ((char *) mat + M_start_off_D);
          M_diff_ptr = M_diff_D;
          switch (operation)
          {
          case BLAS_M_CLEAR :
            Block_Loop_M(*M_ptr=0.0,;);
            break;
          case BLAS_M_SET :
            v_ptr = value_D;
            Block_Loop_M(*M_ptr=*v_ptr,v_ptr++);
            break;
          default :
            ERR_RETURN (-1);                              /* Error: no such Block_Loop_M-operation */
          }
        }
      }

      if (n != 0)
      {
        /* now the rest of the matrices ... */
        nm1 = n-1;
        for (mat=MNEXT(mat); mat!=NULL; mat = MNEXT(mat))
          if ((CTRL(mat)&mat_mask) == mat_value)
          {
            M_ptr = (DOUBLE *) ((char *) mat + M_start_offset);
            M_diff_ptr = M_diff_start;
            switch (operation)
            {
            case BLAS_M_CLEAR :
              Block_Loop_M(*M_ptr=0.0,;);
              break;
            case BLAS_M_SET :
              v_ptr = v_start_ptr;
              Block_Loop_M(*M_ptr=*v_ptr,v_ptr++);
              break;
            default :
              ERR_RETURN (-1);                                    /* Error: no such Block_Loop_M-operation */
            }
          }
      }

    }      /* end loop vect */

  return (0);

}

#define Block_Loop_MN(operation,additional_increment)      \
  operation;                                                 \
  for (j=nm1; j!=0; j--)                                     \
  {                                                          \
    M_ptr = (DOUBLE *) ((char *) M_ptr + *M_diff_ptr++);   \
    N_ptr = (DOUBLE *) ((char *) N_ptr + *N_diff_ptr++);   \
    additional_increment;                                  \
    operation;                                             \
  }

static INT Matrix_Loop_MN(GRID *grid, unsigned int _vroot_mask_, unsigned int _vroot_value_,
                          unsigned int _mat_mask_, unsigned int _mat_value_,
                          int n, int n_D,
                          const ptrdiff_t M_start_off,   const ptrdiff_t *M_diff,
                          const ptrdiff_t M_start_off_D, const ptrdiff_t *M_diff_D,
                          const ptrdiff_t N_start_off,   const ptrdiff_t *N_diff,
                          const ptrdiff_t N_start_off_D, const ptrdiff_t *N_diff_D,
                          int mode, const DOUBLE *value, const DOUBLE *value_D)
{
  /* most important in registers */
  int j;
  DOUBLE *M_ptr;
  const DOUBLE *N_ptr;
  const ptrdiff_t *M_diff_ptr;
  const ptrdiff_t *N_diff_ptr;
  const DOUBLE *v_ptr;

  /* next important in registers */
  int nm1;
  int operation;
  MATRIX *mat;
  unsigned int mat_mask  = _mat_mask_;
  unsigned int mat_value = _mat_value_;
  VECTOR *vect;
  VECTOR *vect2;
  const DOUBLE *v_start_ptr = value;
  ptrdiff_t M_start_offset = M_start_off;
  const ptrdiff_t *M_diff_start = M_diff;
  ptrdiff_t N_start_offset = N_start_off;
  const ptrdiff_t *N_diff_start = N_diff;

  /* not important in registers */
  unsigned int vroot_mask = _vroot_mask_;
  unsigned int vroot_value = _vroot_value_;

  operation = (mode&BLAS_OP_MASK)>>BLAS_OP_SHIFT;

  for (vect=FIRSTVECTOR(grid); vect!=NULL; vect=SUCCVC(vect))
    if ((CTRL(vect)&vroot_mask) == vroot_value)
    {
      mat = VSTART(vect);
      if (n_D != 0)
      {
        if ((CTRL(mat)&mat_mask) == mat_value)
        {
          /* handle the diagonal matrix (always at beginning!) */
          nm1 = n_D-1;
          vect2 = MDEST(mat);
          if (vect2!=vect)
            ERR_RETURN (-1);                              /* Error: matrix has no diagonal connection */

          M_ptr = (DOUBLE *) ((char *) mat + M_start_off_D);
          M_diff_ptr = M_diff_D;
          N_ptr = (const DOUBLE *) ((char *) mat + N_start_off_D);
          N_diff_ptr = N_diff_D;
          switch (operation)
          {
          case BLAS_M_COPY :
            Block_Loop_MN(*M_ptr=*N_ptr,;);
            break;
          case BLAS_M_ADD1 :
            Block_Loop_MN(*M_ptr+=*N_ptr,;);
            break;
          case BLAS_M_MINUS1 :
            Block_Loop_MN(*M_ptr-=*N_ptr,;);
            break;
          case BLAS_M_SCALMUL :
            v_ptr = value_D;
            Block_Loop_MN(*M_ptr=(*v_ptr)*(*N_ptr),v_ptr++);
            break;
          default :
            ERR_RETURN (-1);                              /* Error: no such Block_Loop_M-operation */
          }
        }
      }

      if (n != 0)
      {
        /* now the rest of the matrices ... */
        nm1 = n-1;
        for (mat=MNEXT(mat); mat!=NULL; mat = MNEXT(mat))
          if ((CTRL(mat)&mat_mask) == mat_value)
          {
            M_ptr = (DOUBLE *) ((char *) mat + M_start_offset);
            M_diff_ptr = M_diff_start;
            N_ptr = (const DOUBLE *) ((char *) mat + N_start_offset);
            N_diff_ptr = N_diff_start;
            switch (operation)
            {
            case BLAS_M_COPY :
              Block_Loop_MN(*M_ptr = *N_ptr; , );
              break;
            case BLAS_M_ADD1 :
              Block_Loop_MN(*M_ptr += *N_ptr; , );
              break;
            case BLAS_M_MINUS1 :
              Block_Loop_MN(*M_ptr -= *N_ptr; , );
              break;
            case BLAS_M_SCALMUL :
              v_ptr = v_start_ptr;
              Block_Loop_MN(*M_ptr = (*v_ptr) * (*N_ptr); , v_ptr++);
              break;
            default :
              ERR_RETURN (-1);                                    /* Error: no such Block_Loop_MN-operation */
            }
          }
      }

    }      /* end loop vect */

  return (0);

}


/* Important: at the end of the j-loop the diff-fields should have put */
/* M_ptr, y_ptr back to their initial position! */
/* x_ptr can only advance simply for the moment (perhaps for always, */
/* if no necessity arises). */
/* Note also, that the M difference field is different to the one */
/* needed in the routines above. */
#define Block_Loop_Mxy(pre_inner,operation,post_inner)                     \
  {                                                                                                                          \
    for (i=nr; i>0; i--)                                                                       \
    {                                                                                                                  \
      pre_inner;                                                                                         \
      for (j=*ncol_ptr; j>0; j--)                                                \
      {                                                                                                          \
        operation;                                                                                 \
        M_ptr = (DOUBLE *) ((char *) M_ptr + *M_diff_ptr++);   \
        y_ptr = (DOUBLE *) ((char *) y_ptr + *y_diff_ptr++);   \
      }                                                                                                          \
      post_inner;                                                                                        \
                                                                                                                           \
      ncol_ptr++;                                                                                        \
      x_ptr++;                                                                                           \
    }                                                                                                                  \
  }

static INT Matrix_Loop_Mxy(GRID *grid, unsigned int _vroot_mask_, unsigned int _vroot_value_,
                           unsigned int _mat_mask_, unsigned int _mat_value_,
                           int n, int n_D, int _nr_, const int *ncols, const int *ncols_D,
                           const ptrdiff_t M_start_off, const ptrdiff_t *M_diff,
                           const ptrdiff_t M_start_off_D, const ptrdiff_t *M_diff_D,
                           const ptrdiff_t y_start_off, const ptrdiff_t *y_diff,
                           const ptrdiff_t y_start_off_D, const ptrdiff_t *y_diff_D,
                           const ptrdiff_t y_0_off_D, const ptrdiff_t x_start_off,
                           int mode, /*const DOUBLE *value_start,*/ DOUBLE *result)
{
  /* most important in registers */
  int j;
  const DOUBLE *M_ptr;
  const DOUBLE *y_ptr;
  const ptrdiff_t *M_diff_ptr;
  const ptrdiff_t *y_diff_ptr;
  DOUBLE sum;
  /* const DOUBLE *v_ptr; */

  /* next important in registers */
  const int nr = _nr_;
  const int *ncol_ptr;
  DOUBLE *x_ptr;
  int i;

  /* not very important in registers */
  int operation;
  MATRIX *mat;
  unsigned int mat_mask  = _mat_mask_;
  unsigned int mat_value = _mat_value_;
  VECTOR *vect;
  VECTOR *vect2;
  /* const DOUBLE *v_start_ptr = value_start; */
  ptrdiff_t M_start_offset = M_start_off;
  const ptrdiff_t *M_diff_start = M_diff;
  ptrdiff_t y_start_offset = y_start_off;
  const ptrdiff_t *y_diff_start = y_diff;
  ptrdiff_t x_start_offset = x_start_off;
  DOUBLE global_sum;

  /* not important in registers */
  MATRIX *matD=NULL;

#ifndef __UG__
  DOUBLE *LR=NULL, *sum_vec=NULL;
  int *pivot=NULL;
#else /* __UG__ */
  DOUBLE LR[MAX_MAT_COMP], sum_vec[MAX_VEC_COMP];
  int pivot[MAX_VEC_COMP];
#endif /* __UG__ */

  unsigned int vroot_mask = _vroot_mask_;
  unsigned int vroot_value = _vroot_value_;

  operation = (mode&BLAS_OP_MASK)>>BLAS_OP_SHIFT;

  if (operation == BLAS_MV_LGS)
  {
    if (n_D==0)
      ERR_RETURN(-1);                    /* Error: no lgs possible with empty diagonal */

    if (GORDERED(grid)==0)
      ERR_RETURN(-1);                    /* Error: lgs needs ordered matrix graph (command set_index) */

#ifndef __UG__
    /* we have to allocate some memory for solving an nr*nr-system */
    LR = (DOUBLE *) Get_Temporary_Memory_Global_Heap (nr*nr*sizeof(DOUBLE));
    if (LR == NULL)
      ERR_RETURN (-1);                    /* Error: could not allocate LR field */

    sum_vec = (DOUBLE *) Get_Temporary_Memory_Global_Heap (nr*sizeof(DOUBLE));
    if (sum_vec == NULL)
      ERR_RETURN (-1);                    /* Error: could not allocate sum_vec field */

    pivot = (int *) Get_Temporary_Memory_Global_Heap (nr*sizeof(int));
    if (pivot == NULL)
      ERR_RETURN (-1);                    /* Error: could not allocate LR field */
#else /* __UG__ */
    if (nr>MAX_VEC_COMP)
      ERR_RETURN (-1);                    /* Error: nr too large, increase MAX_VEC_COMP (udm.h) */
    if (nr*nr>MAX_MAT_COMP)
      ERR_RETURN (-1);                    /* Error: nr too large, increase MAX_MAT_COMP (udm.h) */
#endif /* __UG__ */
  }

  if (operation == BLAS_MV_BILFORM)
    global_sum = 0.0;

  for (vect=FIRSTVECTOR(grid); vect!=NULL; vect=SUCCVC(vect))
    if ((CTRL(vect)&vroot_mask) == vroot_value)
    {
      /* probably not needed */
      /*if (operation&BLAS_CLEARX) */
      /*{ */
      /*	x_ptr = (DOUBLE *) ((char *) vect + x_start_offset); */
      /*	for (i=0; i<nr; i++) */
      /*		x_ptr++ = 0.0; */
      /*} */

      matD = VSTART(vect);
      if (n_D != 0)
      {
        if ((CTRL(matD)&mat_mask) == mat_value)
        {
          /* handle the diagonal matrix (always at start of list!) */
          vect2 = MDEST(matD);
          if (vect2!=vect)
            ERR_RETURN (-1);                              /* Error: matrix has no diagonal connection */

          M_ptr = (DOUBLE *) ((char *) matD + M_start_off_D);
          M_diff_ptr = M_diff_D;
          y_ptr = (const DOUBLE *) ((char *) vect + y_start_off_D);
          y_diff_ptr = y_diff_D;
          x_ptr = (DOUBLE *) ((char *) vect + x_start_offset);
          ncol_ptr = ncols_D;
          switch (operation)
          {
          case BLAS_MV_MUL :
            /* note that we do "*x_ptr = 0.0, *x_ptr += sum" here, */
            /* thus this will not work, if n_D==0 (-->Fatal error) (NOTE_M_MUL) */
            Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr = sum );
            break;
          case BLAS_MV_MULADD :
            Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr += sum );
            break;
          case BLAS_MV_MULMINUS :
            Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr -= sum );
            break;
          case BLAS_MV_LGS :
            /* here we initialize sum_vec to hold the defect */
            for (i=0; i<nr; i++) sum_vec[i] = x_ptr[i];
            break;
          case BLAS_MV_BILFORM :
            Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), global_sum += *x_ptr * sum);
            break;
          default :
            ERR_RETURN (-1);                              /* Error: no such Block_Loop_M-operation */
          }
        }
        else
        {
          if (operation==BLAS_MV_LGS)
          {
            if (mat_value-(CTRL(matD)&mat_mask)==MACTIVE_MASK)
            {
              /* to be UG compatible, we have to set the */
              /* result to zero on inactive vectors. such */
              /* vectors are recognized as having their */
              /* diagonal matrix inactive */
              x_ptr = (DOUBLE *) ((char *) vect + y_0_off_D);
              for (i=0; i<nr; i++)
                x_ptr[i] = 0.0;
              continue;                                    /* next vector */
            }
          }
        }
      }

      if (n != 0)
      {
        /* now the rest of the matrices ... */
        for (mat=MNEXT(matD); mat!=NULL; mat = MNEXT(mat))
          if ((CTRL(mat)&mat_mask) == mat_value)
          {
            vect2 = MDEST(mat);
            M_ptr = (DOUBLE *) ((char *) mat + M_start_offset);
            M_diff_ptr = M_diff_start;
            y_ptr = (const DOUBLE *) ((char *) vect2 + y_start_offset);
            y_diff_ptr = y_diff_start;
            ncol_ptr = ncols;
            switch (operation)
            {
            case BLAS_MV_MUL :
              /* note that we have to add here, in contrast to the diagonal (NOTE_M_MUL) */
              x_ptr = (DOUBLE *) ((char *) vect + x_start_offset);
              Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr += sum);
              break;
            case BLAS_MV_MULADD :
              x_ptr = (DOUBLE *) ((char *) vect + x_start_offset);
              Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr += sum);
              break;
            case BLAS_MV_MULMINUS :
              x_ptr = (DOUBLE *) ((char *) vect + x_start_offset);
              Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr -= sum);
              break;
            case BLAS_MV_LGS :
              x_ptr = sum_vec;
              Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), *x_ptr -= sum);
              break;
            case BLAS_MV_BILFORM :
              x_ptr = (DOUBLE *) ((char *) vect + x_start_offset);
              Block_Loop_Mxy(sum=0.0, sum += (*M_ptr)*(*y_ptr), global_sum += *x_ptr * sum);
              break;
            default :
              ERR_RETURN (-1);                                    /* Error: no such Block_Loop_M-operation */
            }
          }
      }

      if (operation==BLAS_MV_LGS)
      {
        /* we have to multiply sum_vec with the inverse of the diagonal block */
        /* and store it in the diagonal of the y-field. */
        /* since y_ptr is a const-pointer, we use x_ptr here for the y-field! */
        if (nr==1)
        {
          /* we handle this case separately, since the function calls might */
          /* represent too much overhead. n_D>0 (i.e.==1) was checked above, */
          /* therefore we may access M_ptr */
          M_ptr = (DOUBLE *) ((char *) matD + M_start_off_D);
          x_ptr = (DOUBLE *) ((char *) vect + y_start_offset);
          {
            register DOUBLE diag = *M_ptr;
            if (diag==0.0)
              ERR_RETURN (-1);                                    /* Error: diagonal is not invertible */
            *x_ptr = *sum_vec / diag;
            continue;
          }
        }

        /* perhaps also the case nr==2 should be handled separately ... */

        /* initialize LR field */
        for (i=0; i<nr*nr; i++)
          LR[i] = 0.0;

        /* inject sparse matrix into rectangular field, */
        /* x_ptr is temporarily used as pointer to LR */
        M_ptr = (DOUBLE *) ((char *) matD + M_start_off_D);
        x_ptr = (DOUBLE *) ((char *) LR + (y_start_off_D-y_0_off_D));
        M_diff_ptr = M_diff_D;
        y_diff_ptr = y_diff_D;
        ncol_ptr = ncols_D;
        for (i=nr; i>0; i--)
        {
          for (j=*ncol_ptr; j>0; j--)
          {
            *x_ptr = *M_ptr;
            M_ptr = (DOUBLE *) ((char *) M_ptr + *M_diff_ptr++);
            x_ptr = (DOUBLE *) ((char *) x_ptr + *y_diff_ptr++);
          }
          ncol_ptr++;
          x_ptr += nr;
        }

        /* then decompose the small block */
        if (Decompose_LR_pivot (nr, LR, pivot)!=0)
          ERR_RETURN(-1);                        /* Error: ... could not decompose small block */

        /* and solve the system LR x = sum_vec */
        x_ptr = (DOUBLE *) ((char *) vect + y_0_off_D);
        if (Solve_LR (nr, LR, pivot, x_ptr, sum_vec) <0)
          ERR_RETURN(-1);                        /* Error: ... could not solve small block system */
      }

    }      /* end loop vect */

  if (operation==BLAS_MV_BILFORM)
    *result = global_sum;

  return (0);

}

/* This routine fixes a certain rtype/ctype pair and branches according to
   the loop_type. It would be relatively easy (and probably will have to
   be done for several problems) to shift the type choice inside the grid loop. */
static INT Matrix_Loop (FORMAT *format, GRID *grid,
                        unsigned int _vroot_mask_, unsigned int _vroot_value_,
                        unsigned int _mat_mask_, unsigned int _mat_value_,
                        const int *_n_, const int *_nr_, int **_NC_,
                        const ptrdiff_t *_M_start_off_, ptrdiff_t **_M_diff_,
                        const ptrdiff_t *_N_start_off_, ptrdiff_t **_N_diff_,
                        const ptrdiff_t *_y_start_off_, ptrdiff_t **_y_diff_,
                        const ptrdiff_t *_y_0_off_, const ptrdiff_t *_x_start_off_,
                        int mode, const DOUBLE **_value_,
                        DOUBLE *result)

{
  INT error;
  int rtype, ctype, type, dtype, loop_type, operation;
  const DOUBLE *value, *value_D;
  unsigned int vroot_mask, vroot_value;
  unsigned int mat_mask, mat_value;

#ifdef __UG__
  if (format==NULL)
    ERR_RETURN(-1);              /* Error: check to make compiler silent (format now used) */
#endif

  loop_type = (mode&BLAS_LOOP_MASK)>>BLAS_LOOP_SHIFT;
  operation = (mode&BLAS_OP_MASK)>>BLAS_OP_SHIFT;

  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=0; ctype<NVECTYPES; ctype++)
    {
      /* initialize the actual parameters */
      int n=0, n_D=0, nr=0;
      const int *NC = NULL, *NC_D = NULL;
      ptrdiff_t M_start_off   = 0; const ptrdiff_t *M_diff   = NULL;
      ptrdiff_t M_start_off_D = 0; const ptrdiff_t *M_diff_D = NULL;
      ptrdiff_t N_start_off   = 0; const ptrdiff_t *N_diff   = NULL;
      ptrdiff_t N_start_off_D = 0; const ptrdiff_t *N_diff_D = NULL;
      ptrdiff_t y_start_off   = 0; const ptrdiff_t *y_diff   = NULL;
      ptrdiff_t y_start_off_D = 0; const ptrdiff_t *y_diff_D = NULL;
      ptrdiff_t y_0_off_D = 0, x_start_off = 0;

      type = MTP(rtype,ctype); n = _n_[type];

      n_D = 0;
      if (rtype==ctype) {
        dtype = DMTP(rtype); n_D = _n_[dtype];
      }

      /* note, that n is either Nred/nr depending on loop type. */
      /* But in a check for zero, this does not matter. */
      if ((n==0)&&(n_D==0)) continue;            /* no matrix for this type combination */

      /* inits for M */
      if (n!=0)
      {
        M_start_off = _M_start_off_[type];
        M_diff = _M_diff_[type];
      }

      if (n_D!=0)
      {
        M_start_off_D = _M_start_off_[dtype];
        M_diff_D = _M_diff_[dtype];
      }

      /* inits for value field */
      if (_value_!=NULL)
      {
        if (n!=0) value = _value_[type];
        if (n_D!=0) value_D = _value_[dtype];
      }
      else
        value = value_D = NULL;

      /* now initialize the masks that are needed for our loops */
      vroot_mask  = VTYPE_MASK | _vroot_mask_;
      vroot_value = (rtype << VTYPE_SHIFT) | _vroot_value_;
      mat_mask    = _mat_mask_;
      mat_value   = _mat_value_;

      /* additionally the type criterion must be set */
      mat_value |= (ctype << MDESTTYPE_SHIFT);

      switch (loop_type)
      {
      case BLAS_LOOP_M :
        error = Matrix_Loop_M (grid, vroot_mask, vroot_value,
                               mat_mask, mat_value, n, n_D,
                               M_start_off, M_diff, M_start_off_D, M_diff_D,
                               mode, value, value_D);
        if (error<0)
          ERR_RETURN(-1);                        /* Error: ... error in matrix loop (type M) */

        break;

      case BLAS_LOOP_MN :
        /* inits for N */
        if (n!=0)
        {
          N_start_off = _N_start_off_[type];
          N_diff = _N_diff_[type];
        }
        if (n_D!=0)
        {
          N_start_off_D = _N_start_off_[dtype];
          N_diff_D = _N_diff_[dtype];
        }
        error = Matrix_Loop_MN (grid, vroot_mask, vroot_value,
                                mat_mask, mat_value, n, n_D,
                                M_start_off, M_diff, M_start_off_D, M_diff_D,
                                N_start_off, N_diff, N_start_off_D, N_diff_D,
                                mode, value, value_D);
        if (error<0)
          ERR_RETURN(-1);                        /* Error: ... error in matrix loop (type MN) */

        break;

      case BLAS_LOOP_Mxy :
        /* inits for x, y */
        if (n!=0)
        {
          y_start_off = _y_start_off_[type];
          y_diff      = _y_diff_[type];
          x_start_off = _x_start_off_[type];
          NC          = _NC_[type];
        }
        if (n_D!=0)
        {
          y_start_off_D = _y_start_off_[dtype];
          y_0_off_D     = _y_0_off_[dtype];
          y_diff_D      = _y_diff_[dtype];
          NC_D          = _NC_[dtype];
        }
        else
        {
          if (operation==BLAS_MV_MUL)
            ERR_RETURN(-1);                              /* Error: matmul needs diagonal parts (see NOTE_M_MUL in this file) */
        }
        nr = _nr_[type];
        error = Matrix_Loop_Mxy (grid, vroot_mask, vroot_value,
                                 mat_mask, mat_value, n, n_D, nr, NC, NC_D,
                                 M_start_off, M_diff, M_start_off_D, M_diff_D,
                                 y_start_off, y_diff, y_start_off_D, y_diff_D,
                                 y_0_off_D, x_start_off,
                                 mode, /* value, value_D,*/ result);
        if (error<0)
          ERR_RETURN(-1);                        /* Error: ... error in matrix loop (type Mxy) */

        break;

      default :
        ERR_RETURN(-1);                  /* Error: no such loop_type */
      }
    }

  return(0);
}

/****************************************************************************/
/*D
   MG_Matrix_Loop - loop through types, levels, vectors and matrices

   SYNOPSIS:
   static INT MG_Matrix_Loop(MULTIGRID *mg, INT fl, INT tl, INT mode, INT op,
                             const MATDATA_DESC *M, const MATDATA_DESC *N,
                                                     const VECDATA_DESC *x, const VECDATA_DESC *y,
                             const DOUBLE *value)
   DESCRIPTION:
   Basic routine for the BLAS level 2 operations.

   REMARKS:
   D*/
/****************************************************************************/

INT MG_Matrix_Loop(MULTIGRID *mg, INT fl, INT tl, INT mode,
                   const MATDATA_DESC *M, const MATDATA_DESC *N,
                   const VECDATA_DESC *x, const VECDATA_DESC *y,
                   int N_vals, const DOUBLE *values, DOUBLE *result)
{
  int j, type, vtype;
  SPARSE_MATRIX *sm;
  GRID *grid;
  FORMAT *format;

  /* The M, N, y fields are of size N_mat_types, */
  /* the x-field is of size N_vec_types. */
  /* The diff-fields depend on the sparse matrix and if */
  /* reduced structure or unreduced structure is used. */
  ptrdiff_t *M_start_off = NULL, **M_diff = NULL;
  ptrdiff_t *N_start_off = NULL, **N_diff = NULL;
  ptrdiff_t *y_start_off = NULL, **y_diff = NULL;
  ptrdiff_t *x_start_off = NULL;
  ptrdiff_t *y_0_off = NULL;
  int *n = NULL, *nr = NULL, **NC = NULL;
  const DOUBLE **value_ptr;

#ifdef __UG__ /* in UG dynamic allocation might cause problems */
  ptrdiff_t Ptrdiff[2*NMATTYPES+2*MAX_MAT_COMP];
  void *Ptr[4*NMATTYPES];
  int Int[2*NMATTYPES+NMATTYPES*MAX_VEC_COMP];
  SHORT red_offsets[MAX_MAT_COMP];
#else
  ptrdiff_t *Ptrdiff;
  void **Ptr;
  int *Int;
  SHORT *red_offsets;
#endif /* not __UG__ */

  int pos_pd, pos_int, pos_ptr, pos_v, size_pd, size_ptr, size_int;
  int red_size, red_max;

  unsigned int operation, loop_type, blas_mode;
  unsigned int vroot_mask, vroot_value;
  unsigned int mat_mask, mat_value;

  loop_type = (mode&BLAS_LOOP_MASK)>>BLAS_LOOP_SHIFT;
  operation = (mode&BLAS_OP_MASK)>>BLAS_OP_SHIFT;
  blas_mode = (mode&BLAS_MODE_MASK)>>BLAS_MODE_SHIFT;

  if (mg==NULL) ERR_RETURN(-1);        /* Error: you have to supply an mg */

  format = MGFORMAT(mg);
  if (format==NULL) ERR_RETURN (-1);        /* Error: you have to supply a format */

  /* consistency checks */
  if (M==NULL) ERR_RETURN(-1);        /* Error: you have to supply an MD M */

  /* check further data compatibility */
  /* compute also memory requirements of the (temporary) fields */
  size_pd = 0;
  size_ptr = 0;
  size_int = 0;
  red_max = 0;
  switch (loop_type)
  {
  case BLAS_LOOP_M :
    /* *M_diff */
    for (type=0; type<NMATTYPES; type++)
      if (MD_SM(M,type)!=NULL)
      {
        red_size = SM_Compute_Reduced_Size (MD_SM(M,type));
        if (red_size<0)
          ERR_RETURN (-1);                                /* Error: ... could not compute reduced size of M */
        if (red_size>red_max)
          red_max = red_size;
        size_pd += red_size;
      }
    size_pd  += NMATTYPES;              /* M_start_off */
    size_ptr += NMATTYPES;              /* actual M_diff */
    if (values!=NULL) size_ptr += NMATTYPES;              /*ptrs to value field */
    size_int += NMATTYPES;              /* n */
    break;
  case BLAS_LOOP_MN :
    if (N==NULL) ERR_RETURN(-1);              /* Error: this loop type needs an MD N */
    if (MM_Compatible(format, M, N)!=0)
      ERR_RETURN(-1);                    /* Error: M and N are not compatible */

    /* *M_diff, *N_diff */
    red_max = 0;
    for (type=0; type<NMATTYPES; type++)
      if (MD_SM(M,type)!=NULL)
      {
        red_size = SM_Compute_Reduced_Size (MD_SM(M,type));
        if (red_size<0)
          ERR_RETURN (-1);                                /* Error: ... could not compute reduced size of M */
        if (red_size>red_max)
          red_max = red_size;
        size_pd += 2*red_size;
      }
    size_pd  += 2*NMATTYPES;              /* M_start_off, N_start_off */
    size_ptr += 2*NMATTYPES;              /* M_diff, N_diff */
    if (values!=NULL) size_ptr += NMATTYPES;              /*ptrs to value field */
    size_int +=   NMATTYPES;              /* n */
    break;
  case BLAS_LOOP_Mxy :
    if (x==NULL) ERR_RETURN(-1);              /* Error: this loop type needs a VD x */
    if (y==NULL) ERR_RETURN(-1);              /* Error: this loop type needs a VD y */
    if (VM_Compatible(format, x, M)!=0)
      ERR_RETURN(-1);                    /* Error: x and M are not compatible */
    if (MV_Compatible(format, M, y)!=0)
      ERR_RETURN(-1);                    /* Error: M and y are not compatible */

    for (type=0; type<NMATTYPES; type++)
      if (MD_SM(M,type)!=NULL)
      {
        size_pd  += 2*MD_SM(M,type)->N;                         /* *M_diff, *y_diff */
        size_int += MD_SM(M,type)->nrows;                          /* length of rows in *NC */
      }
    size_pd  += 4*NMATTYPES;              /* M_start_off, y_start_off, y_0_off, x_start_off */
    size_ptr += 3*NMATTYPES;              /* M_diff, y_diff, NC */
    size_int += 2*NMATTYPES;              /* n, nr */
    break;
  default :
    ERR_RETURN(-1);              /* Error: no such loop_type */
  }

        #ifndef __UG__
  /* now allocate the temporary memory */
  Ptrdiff = (ptrdiff_t *) Get_Temporary_Memory_Global_Heap (size_pd*sizeof(ptrdiff_t));
  if (Ptrdiff == NULL)
    ERR_RETURN (-1);              /* Error: could not allocate ptrdiff-fields */

  Ptr = (void **) Get_Temporary_Memory_Global_Heap (size_ptr*sizeof(void *));
  if (Ptr == NULL)
    ERR_RETURN (-1);              /* Error: could not allocate ptr-fields */

  Int = (int *) Get_Temporary_Memory_Global_Heap (size_int*sizeof(int));
  if (Int == NULL)
    ERR_RETURN (-1);              /* Error: could not allocate int-fields */

  red_offsets = (SHORT *) Get_Temporary_Memory_Global_Heap (red_max*sizeof(SHORT));

        #else
  /* check if the arrays are sufficiently large */
  if (size_pd  > 2*NMATTYPES+2*MAX_MAT_COMP)
    ERR_RETURN (-1);              /* Error: ptrdiff-array too small */
  if (size_ptr > 4*NMATTYPES)
    ERR_RETURN (-1);              /* Error: ptr-array too small */
  if (size_int  > 2*NMATTYPES+NMATTYPES*MAX_VEC_COMP)
    ERR_RETURN (-1);              /* Error: int-array too small */
  if (red_max  > MAX_MAT_COMP)
    ERR_RETURN (-1);              /* Error: SHORT-array too small */
        #endif

  /* then set the actual arrays */
  pos_pd  = 0;
  pos_ptr = 0;
  pos_int = 0;
  pos_v   = 0;
  if (values!=NULL)
  {
    value_ptr = (const DOUBLE **) (Ptr + pos_ptr);
    pos_ptr += NMATTYPES;
  }
  else
    value_ptr = NULL;

  switch (loop_type)
  {
  case BLAS_LOOP_M :
    M_diff =    (ptrdiff_t **)    (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    M_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    n           = Int + pos_int;    pos_int += NMATTYPES;
    for (type=0; type<NMATTYPES; type++)
    {
      M_diff[type] = NULL;
      if (values!=NULL)
        value_ptr[type] = NULL;
      M_start_off[type] = 0;
      n[type] = 0;

      sm = MD_SM(M,type);
      if (sm==NULL) continue;

      n[type] = SM_Compute_Reduced_Size(sm);
      M_diff[type] = Ptrdiff + pos_pd; pos_pd += n[type];
      if (values!=NULL)
      {
        value_ptr[type] = values + pos_v;
        pos_v += n[type];
      }
      if (SM_Compute_Reduced_Offsets(sm, red_offsets)<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute reduced offsets for M */

      if (SM_Compute_Diff_From_Offset(n[type], red_offsets, M_diff[type])<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute M_diff-field */

      M_start_off[type] = offsetof(MATRIX, value)
                          + sizeof(DOUBLE) * red_offsets[0];
    }
    break;

  case BLAS_LOOP_MN :
    M_diff =    (ptrdiff_t **)    (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    N_diff =    (ptrdiff_t **)    (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    M_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    N_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    n           = Int + pos_int;    pos_int += NMATTYPES;
    for (type=0; type<NMATTYPES; type++)
    {
      M_diff[type] = N_diff[type] = NULL;
      if (values!=NULL)
        value_ptr[type] = NULL;
      M_start_off[type] = N_start_off[type] = 0;
      n[type] = 0;
      sm = MD_SM(M,type);
      if (sm==NULL) continue;

      n[type] = SM_Compute_Reduced_Size(sm);
      M_diff[type] = Ptrdiff + pos_pd; pos_pd += n[type];
      N_diff[type] = Ptrdiff + pos_pd; pos_pd += n[type];
      if (values!=NULL)
      {
        value_ptr[type] = values + pos_v;
        pos_v += n[type];
      }

      /* set M fields */
      if (SM_Compute_Reduced_Offsets(sm, red_offsets)<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute reduced offsets for M */

      if (SM_Compute_Diff_From_Offset(n[type], red_offsets, M_diff[type])<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute M_diff-field */

      M_start_off[type] = offsetof(MATRIX, value)
                          + sizeof(DOUBLE) * red_offsets[0];

      /* and now the same for N */
      if (SM_Compute_Reduced_Offsets(MD_SM(N,type), red_offsets)<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute reduced offsets for N */

      if (SM_Compute_Diff_From_Offset(n[type], red_offsets, N_diff[type])<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute N_diff-field */

      N_start_off[type] = offsetof(MATRIX, value)
                          + sizeof(DOUBLE) * red_offsets[0];
    }
    break;

  case BLAS_LOOP_Mxy :
    M_diff      = (ptrdiff_t **) (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    y_diff      = (ptrdiff_t **) (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    NC          = (int **)       (Ptr + pos_ptr); pos_ptr += NMATTYPES;
    M_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    y_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    y_0_off     = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    x_start_off = Ptrdiff + pos_pd; pos_pd  += NMATTYPES;
    n           = Int + pos_int;    pos_int += NMATTYPES;
    nr          = Int + pos_int;    pos_int += NMATTYPES;
    for (type=0; type<NMATTYPES; type++)
    {
      M_diff[type] = y_diff[type] = NULL;
      M_start_off[type] = y_start_off[type] = x_start_off[type] = 0;
      NC[type] = NULL;
      n[type] = nr[type] = 0;

      sm = MD_SM(M,type);
      if (sm==NULL) continue;

      n[type] = sm->N;
      nr[type] = sm->nrows;
      M_diff[type] = Ptrdiff + pos_pd; pos_pd  += sm->N;
      y_diff[type] = Ptrdiff + pos_pd; pos_pd  += sm->N;
      NC[type] = Int + pos_int;        pos_int += sm->nrows;

      /* set M fields */
      if (SM_Compute_Diff_From_Offset(sm->N, sm->offset, M_diff[type])<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute M_diff-field */

      M_start_off[type] = offsetof(MATRIX, value)
                          + sizeof(DOUBLE) * sm->offset[0];

      /* set also the number of columns per row */
      for (j=0; j<nr[type]; j++)
        NC[type][j] = sm->row_start[j+1]-sm->row_start[j];

      /* then the y fields */
      vtype = MTYPE_CT(type);
      if (SM_Compute_yDiff_From_Offset(sm->N, sm->col_ind,
                                       &(VD_CMP_OF_TYPE(y,vtype,0)),
                                       y_diff[type])<0)
        ERR_RETURN (-1);                          /* Error: ... could not compute y_diff-field */

      y_start_off[type] = offsetof(VECTOR, value)
                          + sizeof(DOUBLE) * VD_CMP_OF_TYPE(y,vtype,sm->col_ind[0]);

      y_0_off[type] = offsetof(VECTOR, value) +
                      + sizeof(DOUBLE) * VD_CMP_OF_TYPE(y,vtype,0);

      /* and finally x_start_off */
      vtype = MTYPE_RT(type);
      x_start_off[type] = offsetof(VECTOR, value)
                          + sizeof(DOUBLE) * VD_CMP_OF_TYPE(x,vtype,0);
    }
    break;
  }

  /* consistency check */
  if (pos_pd != size_pd)
    ERR_RETURN(-1);              /* Error (internal): ptrdiff field computation inconsistent */
  if (pos_ptr != size_ptr)
    ERR_RETURN(-1);              /* Error (internal): ptr field computation inconsistent */
  if (pos_int != size_int)
    ERR_RETURN(-1);              /* Error (internal): int field computation inconsistent */
  if (pos_v != N_vals)
    ERR_RETURN(-1);              /* Error (internal): size of value field does not fit */

  /******** might be useful, if only off-diagonal matrices occur somewhere
                        see token NOTE_M_MUL
                        if (mode&CLEARX) dset(mg,fl,tl,mode,x,0.0);********/

  /* active criterion for vectors */
  if (blas_mode&BLAS_VACTIVE)
    vroot_mask = vroot_value = VACTIVE_MASK;
  else
    vroot_mask = vroot_value = 0;

  /* active criterion for matrices */
  if (blas_mode&BLAS_MACTIVE)
    mat_mask = mat_value = MACTIVE_MASK;
  else
    mat_mask = mat_value = 0;

  /* further criterion for matrices: diag, lower, upper? */
  switch ((mode&MBLAS_MTYPE_MASK)>>MBLAS_MTYPE_SHIFT)
  {
  case MBLAS_NONE :
    mat_value |= MDIAG_MASK;
    break;
  case MBLAS_DIAG :
    mat_mask  |= MDIAG_MASK;
    mat_value |= MDIAG_MASK;
    break;
  case MBLAS_LOWER :
    mat_mask  |= MLOWER_MASK;
    mat_value |= MLOWER_MASK;
    break;
  case MBLAS_NOT_UPPER :
    mat_mask  |= MUPPER_MASK;
    break;
  case MBLAS_UPPER :
    mat_mask  |= MUPPER_MASK;
    mat_value |= MUPPER_MASK;
    break;
  case MBLAS_NOT_LOWER :
    mat_mask  |= MLOWER_MASK;
    break;
  case MBLAS_NOT_DIAG :
    mat_mask  |= MDIAG_MASK;
    break;
  case MBLAS_ALL :      /* no restriction */
    break;
  default :
    ERR_RETURN(-1);              /* Error (internal): bad MTYPE criterion */
  }

  /* setup phase is over (wow!), we enter the computational phase */
  ASSERT(fl<=tl);
  ASSERT((mode != ON_SURFACE) || (fl <= FULLREFINELEVEL(mg)));
  if (blas_mode&BLAS_SURFACE)
  {
    for (grid=BOTTOM_GRID(mg); grid!=NULL; grid=FINER(grid))
    {
      if (grid->level<tl)
      {
        /* not on the finest level */
        if (grid->level>=FULLREFINELEVEL(mg))
        {
          if (Matrix_Loop (format, grid,
                           vroot_mask|FINE_GRID_DOF_MASK, vroot_mask|FINE_GRID_DOF_MASK,
                           mat_mask, mat_value,
                           n, nr, NC,     M_start_off, M_diff,
                           N_start_off, N_diff, y_start_off, y_diff,
                           y_0_off, x_start_off, mode, value_ptr, result)
              <0)
            ERR_RETURN(-1);                                      /* Error: ... in Matrix_Loop */
        }
      }
      else
      {
        /* on the finest level */
        if (Matrix_Loop (format, grid,
                         vroot_mask|NEW_DEFECT_MASK, vroot_mask|NEW_DEFECT_MASK,
                         mat_mask, mat_value,
                         n, nr, NC,     M_start_off, M_diff,
                         N_start_off, N_diff, y_start_off, y_diff,
                         y_0_off, x_start_off, mode, value_ptr, result)
            <0)
          ERR_RETURN(-1);                                /* Error: ... in Matrix_Loop */
      }
    }
  }
  else
  {
    for (grid=BOTTOM_GRID(mg); grid!=NULL; grid=FINER(grid))
      if ( (grid->level>=fl) && (grid->level<=tl) )
        if (Matrix_Loop (format, grid,
                         vroot_mask, vroot_mask,
                         mat_mask, mat_value,
                         n, nr, NC, M_start_off, M_diff,
                         N_start_off, N_diff, y_start_off, y_diff,
                         y_0_off, x_start_off, mode, value_ptr, result)
            <0)
          ERR_RETURN(-1);                                /* Error: ... in Matrix_Loop */
  }

  return (0);
}
