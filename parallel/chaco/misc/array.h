// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* header file for array.c */

/* exported functions */
void merge_goals(double *merged_goal, double *goal, int ndims, int part_type);
int partition_type(int xlen, int ylen, int ndims, int *array2proc,
                   int *xlens, int *ylens);
int map (int x, int y, int ndims, int part_type,
         int *xlens, int *ylens);
int suba2a(short *assignment, int *loc2glob, int i, short *subassign, int ndims,
           int xlen, int ylen, int set, int *xlens, int *ylens, int part_type);
int compute_sub_ndims (int xlen, int ylen, int ndims);
