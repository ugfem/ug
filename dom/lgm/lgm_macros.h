// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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

#define Det2d(a)                                ( a[0] * a[3] - a[1] * a[2] )

#define InvMatMult2d(b,c,a)             b[0] =    a[3] * c[0]           \
                                               - a[1] * c[1];          \
  b[1] =  - a[2] * c[0]           \
         + a[0] * c[1];          \
  b[0] = b[0] / Det2d(a);         \
  b[1] = b[1] / Det2d(a);

#define LGM_VECTOR_PRODUCT(A,B,C)       {       (C)[0] = (A)[1]*(B)[2] - (A)[2]*(B)[1];\
                                                (C)[1] = (A)[2]*(B)[0] - (A)[0]*(B)[2];\
                                                (C)[2] = (A)[0]*(B)[1] - (A)[1]*(B)[0];}

#define LGM_SCALAR_PRODUCT(A,B,c)       (c) = ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]);
