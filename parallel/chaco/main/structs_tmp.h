// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* An array of these stores all the data for the graph/matrix. */
/* PERHAPS DON'T NEED GNEDGES AND VNUM */
struct vtx_data {
  /*int vnum;*/		/* global number of this vertex */
  int vwgt;                     /* weight of vertex */
  int nedges;                   /* number of neighbors of vertex in subgraph */
  /*int gnedges;*/	/* number of neighbors in original graph */
  /* Note: above fields always include self-edge */
  int *edges;                   /* neighbor list in subgraph numbering scheme */
  float *ewgts;                 /* weights of all the edges */
                                /* Note: above 2 fields have self-edge first */
};


/* Doubly linked list to hold coarsened graphs */
struct graphlist {
  struct vtx_data **graph;              /* current graph */
  int nvtxs;                            /* number of vertices in graph */
  int nedges;                           /* number of edges in graph */
  double *vwsqrt;                       /* square roots of vertex weights */
  double **evecs;                       /* lowest eigenvectors of graph */
  double *evals;                        /* corresponding eigenvalues of graph */
  int *mflag;                           /* matching data for contracted edges */
  int *v2cv;                            /* mapping from fine to coarse vtxs */
  struct ipairs *merged;                /* pairs of merged vertices */
  struct fpairs *reduction;             /* values used to reduce edge weights */
  struct graphlist *finer;              /* pointer to finer graph */
  struct graphlist *coarser;            /* pointer to coarser graph */
};


/* Used in root finding in the Lanczos code */
struct arglist {
  double zleft;                 /* left zero of rational func */
  double pole;                  /* pole of rational func */
  double zright;                /* rigth zero of rational func */
  double *alpha;                /* Lanczos constants */
  double *beta;                 /* Lanczos constants */
  int j;                        /* step in the Lanczos process */
};


/* A linked list of these stores the selective orthogonalization set */
struct orthlink {
  int depth;                    /* bottom of list is 0, previous is 1 etc */
  int index;                    /* position in list of ritz vals (i index) */
  double ritzval;               /* good ritz value */
  double betaji;                /* residual bound on good ritz pair */
  double tau;                   /* from orthogonality recursion */
  double prevtau;               /* from orthogonality recursion */
  double *vec;                  /* vector to orthogonalize against */
  struct orthlink *pntr;        /* pointer to next link */
};


/* A linked list of these stores the minimum elements of a vector */
struct scanlink {
  double val;                   /* value of vector entry */
  int indx;                     /* index of corresponding entry */
  struct scanlink *pntr;        /* pointer to next link */
};

/* These store the phantom edges needed to keep a subgraph connected */
struct edgeslist {
  int vtx1;                      /* first vertex in edge */
  int vtx2;                      /* second vertex in edge */
  struct edgeslist *next;        /* pointer to next element in list */
};

/* These store all the data needed to modify edges for connectivity. */
struct connect_data {
  struct ilists *old_edges;             /* overwritten old edges */
  struct flists *old_ewgts;             /* overwritten old weights */
  struct edgeslist *new_edges;          /* list of new edges */
  int old_nedges;               /* original number of edges in graph */
};

/* Linked list stuff for various uses */
struct list {                           /* linked list of integers */
  int num;                              /* element number */
  struct list *next;                    /* ptr to next element in list */
};

struct lists {                          /* linked list of lists */
  struct list *begin;                   /* pointer to list */
  struct lists *nextlist;               /* next list header */
};

struct bilist {                         /* bidirectional list */
  struct bilist *prev;                  /* pointer to previous element */
  struct bilist *next;                  /* ptr to next element in list */
};

struct ipairs {                         /* stores pairs of integers */
  int val1;
  int val2;
};

struct dpairs {                         /* stores pairs of doubles */
  double val1;
  double val2;
};

struct fpairs {                         /* stores pairs of floats */
  double val1;
  double val2;
};

struct ilists {                         /* linked list of integer lists */
  int *list;
  struct ilists *next;
};

struct flists {                         /* linked list of floating lists */
  float *list;
  struct flists *next;
};
