// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#define LIST_LEN 200

#define MAXTAG 8
#define TETTAG 4

#define Scalar_float                                    0
#define Scalar_double                                   1
#define Vector_float                                    2
#define Vector_double                                   3

#define N_RULES                                         242
#define MAX_VECTOR_ENTRIES                              10

#define PATHDEPTHMASK 0xF0000000
#define PATHDEPTHSHIFT 28
#define PATHDEPTH(i)                            (((i) & PATHDEPTHMASK)>>PATHDEPTHSHIFT)

/* 2 bits at position n for element side */
#define NEXTSIDEMASK 0x00000003
#define NEXTSIDE(i,n)                           (((i) & (NEXTSIDEMASK<<(2*(n))))>>(2*(n)))
