// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <math.h>


#define SMALL 1e-8
#define MAXDOUBLE 1e200
#define PI 3.141592654

#define Lenght(vec)             sqrt(vec[0]*vec[0]      \
                                     +vec[1]*vec[1]  \
                                     +vec[2]*vec[2])

#define Cross(vec,vec1,vec2)    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]; \
  vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]; \
  vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

#define Minus(sol,vec1,vec2)    sol[0] = vec1[0] - vec2[0];     \
  sol[1] = vec1[1] - vec2[1];     \
  sol[2] = vec1[2] - vec2[2];

#define Plus(sol,vec1,vec2)     sol[0] = vec1[0] + vec2[0];     \
  sol[1] = vec1[1] + vec2[1];     \
  sol[2] = vec1[2] + vec2[2];

#define Mult(vec1,vec2)         ( vec1[0] * vec2[0]     \
                                  + vec1[1] * vec2[1]     \
                                  + vec1[2] * vec2[2])

#define Scale(vec,l)            vec[0] = vec[0] * l;            \
  vec[1] = vec[1] * l;            \
  vec[2] = vec[2] * l;

#define Det(a)                  (       - a[0] * a[4] * a[8]    \
                                        + a[0] * a[5] * a[7]    \
                                        + a[3] * a[1] * a[8]    \
                                        - a[3] * a[2] * a[7]    \
                                        - a[6] * a[1] * a[5]    \
                                        + a[6] * a[2] * a[4])

#define Det3d(a)                (         a[0][0] * a[1][1] * a[2][2]   \
                                          + a[0][1] * a[1][2] * a[2][0]   \
                                          + a[0][2] * a[1][0] * a[2][1]   \
                                          - a[0][0] * a[1][2] * a[2][1]   \
                                          - a[0][1] * a[1][0] * a[2][2]   \
                                          - a[0][2] * a[1][1] * a[2][0])

#define InvMatMult(b,c,a)               b[0] =    ( a[5] * a[7] - a[4] * a[8] ) * c[0]          \
                                               + ( a[1] * a[8] - a[2] * a[7] ) * c[1]          \
                                               + ( a[2] * a[4] - a[1] * a[5] ) * c[2];         \
  b[1] =    ( a[3] * a[8] - a[5] * a[6] ) * c[0]          \
         + ( a[2] * a[6] - a[0] * a[8] ) * c[1]          \
         + ( a[0] * a[5] - a[2] * a[3] ) * c[2];         \
  b[2] =    ( a[4] * a[6] - a[3] * a[7] ) * c[0]          \
         + ( a[0] * a[7] - a[1] * a[6] ) * c[1]          \
         + ( a[1] * a[3] - a[0] * a[4] ) * c[2];         \
  b[0] = b[0] / Det(a);                                   \
  b[1] = b[1] / Det(a);                                   \
  b[2] = b[2] / Det(a);

#define InvMatMult3d(b,c,a)             b[0] =    ( a[1][1] * a[2][2] - a[1][2] * a[2][1] ) * c[0]              \
                                               + ( a[0][2] * a[2][1] - a[0][1] * a[2][2] ) * c[1]              \
                                               + ( a[0][1] * a[1][2] - a[0][2] * a[1][1] ) * c[2];             \
  b[1] =    ( a[1][2] * a[2][0] - a[1][0] * a[2][2] ) * c[0]              \
         + ( a[0][0] * a[2][2] - a[0][2] * a[2][0] ) * c[1]              \
         + ( a[0][2] * a[1][0] - a[0][0] * a[1][2] ) * c[2];             \
  b[2] =    ( a[1][0] * a[2][1] - a[1][1] * a[2][0] ) * c[0]              \
         + ( a[0][1] * a[2][0] - a[0][0] * a[2][1] ) * c[1]              \
         + ( a[0][0] * a[1][1] - a[0][1] * a[1][0] ) * c[2];             \
  b[0] = b[0] / Det3d(a);                                                                                 \
  b[1] = b[1] / Det3d(a);                                                                                 \
  b[2] = b[2] / Det3d(a);



#define Det2d(a)                                ( a[0][0] * a[1][1] - a[1][0] * a[0][1] )

#define InvMatMult2d(b,c,a)             b[0] =    a[1][1] * c[0]                \
                                               - a[0][1] * c[1];               \
  b[1] =  - a[1][0] * c[0]                \
         + a[0][0] * c[1];               \
  b[0] = b[0] / Det2d(a);         \
  b[1] = b[1] / Det2d(a);

#define MatCopy3d(a, b)                 a[0][0] = b[0][0];      a[0][1] = b[0][1];      a[0][2] = b[0][2];      \
  a[1][0] = b[1][0];      a[1][1] = b[1][1];      a[1][2] = b[1][2];      \
  a[2][0] = b[2][0];      a[2][1] = b[2][1];      a[2][2] = b[2][2];

#define VecCopy3d(a, b)                 a[0] = b[0];    \
  a[1] = b[1];    \
  a[2] = b[2];

int GetNormalVector(double *p0, double *p1, double *p2, double *n);
