// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ansys2lgm.h                                                   */
/*                                                                          */
/* Purpose:   header file for ansys2UG			                                                */
/*                                                                          */
/* Author:    Dirk Feuchter                                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   27.5.97 begin, ug version 3.7                                     */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

struct struct_DOM_INF_TP {
  INT NumberOfSubdomains;
  INT NumberOfSurfaces;
  INT NumberOfPolyLines;
  INT NumberOfPoints;
} ;
typedef struct struct_DOM_INF_TP DOMAIN_INFO_TYP;
#define NMB_OF_SBDMS(p)                         (p->NumberOfSubdomains)
#define NMB_OF_SFCES(p)                         (p->NumberOfSurfaces)
#define NMB_OF_PLINES(p)                        (p->NumberOfPolyLines)
#define NMB_OF_POINTS(p)                        (p->NumberOfPoints)



struct struct_SFE_KNOTEN_TYP {
  INT nodeid[3];
  struct struct_SFE_KNOTEN_TYP *next;
  struct struct_SFE_KNOTEN_TYP *neighbour[3];
  DOUBLE identifier[2];
  INT thefourthnode[2];
  INT flg;
  INT local_triangle_ID;
} ;
typedef struct struct_SFE_KNOTEN_TYP SFE_KNOTEN_TYP;
#define SFE_NDID_ARRAY(p)               (p->nodeid)
#define SFE_NDID1(p)                    (p->nodeid[0])
#define SFE_NDID2(p)                    (p->nodeid[1])
#define SFE_NDID3(p)                    (p->nodeid[2])
#define SFE_NDID(p,n)                   (p->nodeid[n])
#define SFE_NGHB1(p)                    (p->neighbour[0])
#define SFE_NGHB2(p)                    (p->neighbour[1])
#define SFE_NGHB3(p)                    (p->neighbour[2])
#define SFE_NGHB(p,n)                   (p->neighbour[n])
#define SFE_NEXT(p)                     (p->next)
#define SFE_IDF_1(p)                    (p->identifier[1])
#define SFE_IDF_0(p)                    (p->identifier[0])
#define SFE_IDF_N(p,n)                  (p->identifier[n])
#define SFE_4ND_1(p)                    (p->thefourthnode[1])
#define SFE_4ND_0(p)                    (p->thefourthnode[0])
#define SFE_4ND_N(p,n)                  (p->thefourthnode[n])
#define SFE_ORIENTATION_FLAG(p) (p->flg)
#define SFE_LOCAL_TRI_ID(p)             (p->local_triangle_ID)



struct structIDF_TYP {
  DOUBLE idf_value;
  struct structIDF_TYP *next;
  SFE_KNOTEN_TYP *triangle;
  INT dritteid;
} ;
typedef struct structIDF_TYP IDF_TYP;
#define IDF_VAL(p)                              (p->idf_value)
#define IDF_NXT(p)                              (p->next)
#define IDF_SFE_TRIA(p)                 (p->triangle)
#define IDF_ID3(p)                              (p->dritteid)



struct struct_IDF_SHORT_TYP {
  DOUBLE idf_short_value;
  struct struct_IDF_SHORT_TYP *next;
} ;
typedef struct struct_IDF_SHORT_TYP IDF_SHORT_TYP;
#define IDF_SHORT_VAL(p)                (p->idf_short_value)
#define IDF_SHORT_NXT(p)                (p->next)



struct structLI_KNOTEN_TYP {
  INT nodeid[2];
  struct structLI_KNOTEN_TYP *next;
  IDF_TYP *identifiers;
} ;
typedef struct structLI_KNOTEN_TYP LI_KNOTEN_TYP;
#define LI_NDID(p)                              (p->nodeid)
#define LI_NDID1(p)                             (p->nodeid[0])
#define LI_NDID2(p)                             (p->nodeid[1])
#define LI_NEXT(p)                              (p->next)
#define LI_IDFS(p)                              (p->identifiers)


struct structPL_LINE_TYP {
  struct structPL_LINE_TYP *next;
  LI_KNOTEN_TYP *LineOfPolyline;
} ;
typedef struct structPL_LINE_TYP PL_LINE_TYP;
#define PL_LINES_NXT(p)                 (p->next)
#define PL_LINES_LINE(p)                (p->LineOfPolyline)


/*Typdefinition fuer die Polyline
   Bedeutugngsdefinition einer Polyline in ansys2UG in Worten:
   "...eine Polyline ist die sortierte Vereinigung aller Lines,
   die die selben Surfaceidentifiers besitzen..."*/
struct structPL_TYP {
  IDF_TYP *Characteristic_Identifiers;
  INT NumberOfCharacteristicIDfsOfThePolyline;
  struct structPL_TYP *next;
  PL_LINE_TYP *LinesOfPolyline;
  INT NumberOfPointsOfThePolyline;
} ;
typedef struct structPL_TYP PL_TYP;
#define PL_IDFS(p)                              (p->Characteristic_Identifiers)
#define PL_NMB_OF_CH_IDFS(p)    (p->NumberOfCharacteristicIDfsOfThePolyline)
#define PL_NXT(p)                               (p->next)
#define PL_LINES(p)                             (p->LinesOfPolyline)
#define PL_NMB_OF_POINTS(p)             (p->NumberOfPointsOfThePolyline)        /*nicht unbedingt notwendig*/



struct struct_TRIANGLE_TYP {
  SFE_KNOTEN_TYP *sfepter;
  struct struct_TRIANGLE_TYP *next;
} ;
typedef struct struct_TRIANGLE_TYP TRIANGLE_TYP; /* Typ fuer die Triangles einer Surface */
#define TRIA_SFE_KN(p)                  (p->sfepter)
#define TRIA_NEXT(p)                    (p->next)


struct struct_SFPL_TYP {
  PL_TYP *plline;
  struct struct_SFPL_TYP *next;
} ;
typedef struct struct_SFPL_TYP SFPL_TYP; /* Typ fuer die Polylines einer Surface */
#define SFPL_PL(p)                              (p->plline)
#define SFPL_NEXT(p)                            (p->next)


/* SURFACEDETECTOR . . .  */
struct struct_PLZ_TYP {
  struct struct_PLZ_TYP *next;
  INT plz_nmb_of_pllns;
  SFPL_TYP *plz_polylines;
} ;
typedef struct struct_PLZ_TYP PLZ_TYP; /*Typ fuer die Polylinezyklen */
#define PLZ_NEXT(p)                                     (p->next)
#define PLZ_NMB_OF_POLYLINES(p)         (p->plz_nmb_of_pllns)
#define PLZ_POLYLINES(p)                        (p->plz_polylines)



struct struct_RS_TYP {
  struct struct_RS_TYP *next;
  PLZ_TYP *RealSurface_PolylineZyklen;
  INT rs_nmb_of_plz;
} ;
typedef struct struct_RS_TYP RS_TYP; /*Typ fuer die "realSurfaces" */
#define RS_NEXT(p)                                      (p->next)
#define RS_PL_ZKLN(p)                           (p->RealSurface_PolylineZyklen)
#define RS_NMB_OF_PL_ZKLN(p)            (p->rs_nmb_of_plz)
/* . . . SURFACEDETECTOR */

struct struct_SF_TYP {
  struct struct_SF_TYP *next;
  TRIANGLE_TYP *triangles;
  INT nmboftrias;
  INT nmbofpoints;
  DOUBLE sfcname[2];
  INT rightsubdomain;
  INT leftsubdomain;
  SFPL_TYP *polylines;
  INT nmbofpolylines;
  INT nmbofpolylincycles;
  PLZ_TYP *polylinecycles;
  INT nmbofrealsurfaces;
  RS_TYP *realsurfaces;
} ;
typedef struct struct_SF_TYP SF_TYP; /*Typ fuer die surfaces */
#define SF_NAME1(p)                             (p->sfcname[0])
#define SF_NAME2(p)                             (p->sfcname[1])
#define SF_NEXT(p)                              (p->next)
#define SF_TRIAS(p)                             (p->triangles)
#define SF_NMB_OF_TRIAS(p)              (p->nmboftrias)
#define SF_RIGHT_SBD(p)                 (p->rightsubdomain)
#define SF_LEFT_SBD(p)                  (p->leftsubdomain)
#define SF_NMB_OF_POINTS(p)             (p->nmbofpoints)
#define SF_POLYLINES(p)                 (p->polylines)
#define SF_NMB_OF_POLYLINES(p)          (p->nmbofpolylines)
#define SF_NMB_OF_POLYLI_ZYK(p)         (p->nmbofpolylincycles)
#define SF_NMB_OF_REALSFCS(p)           (p->nmbofrealsurfaces)
#define SF_POLYLI_ZYK(p)            (p->polylinecycles)
#define SF_REALSFCS(p)          (p->realsurfaces)



struct struct_SFC_TYP {
  struct struct_SFC_TYP *next;
  SF_TYP *surface;
} ;
typedef struct struct_SFC_TYP SFC_TYP; /* Typ fuer die Surfaces einer Subdomain */
#define SFC_SURF(p)                             (p->surface)
#define SFC_NEXT(p)                             (p->next)

struct struct_SD_TYP {
  struct struct_SD_TYP *next;
  SFC_TYP *sfces;
  INT nmbofsfces;
  INT sbdname;
} ;
typedef struct struct_SD_TYP SD_TYP; /* Typ fuer die Subdomains */
#define SD_NAME(p)                              (p->sbdname)
#define SD_NEXT(p)                              (p->next)
#define SD_SFCS(p)                              (p->sfces)
#define SD_NMB_OF_SFCS(p)               (p->nmbofsfces)



typedef struct {
  int Knoten_i;
  int Knoten_j;
  int Knoten_k;
  int VierterKnoten;
  double SurfaceIdentifier;
} CAD_SFE_TYP;
#define CAD_SFE_ND_I(p)                         (p->Knoten_i)
#define CAD_SFE_ND_J(p)                         (p->Knoten_j)
#define CAD_SFE_ND_K(p)                         (p->Knoten_k)
#define CAD_SFE_ND_4(p)                         (p->VierterKnoten)
#define CAD_SFE_SFE_IDF(p)                      (p->SurfaceIdentifier)

typedef struct {
  int Element_i;
  int Side_j;
  double SurfaceIdentifier;
} BND_SFE_TYP;
#define BND_SFE_ELEMENTID(p)                            (p->Element_i)
#define BND_SFE_SIDEID(p)                               (p->Side_j)
#define BND_SFE_SFCIDF(p)                               (p->SurfaceIdentifier)


typedef struct {
  INT nmb_of_SFEs;
  INT nmb_of_BndNds;
  CAD_SFE_TYP *SFE_Array;
  DOUBLE *n_koord_array;
  DOUBLE radius;
  DOUBLE MidPoint[3];
} EXCHNG_TYP1; /* Typ fuer den Austausch von Daten zw. Cadconvert()
                  und Ansys2LGM */
#define EXCHNG_TYP1_NMB_OF_SFES(p)                      (p->nmb_of_SFEs)
#define EXCHNG_TYP1_NMB_OF_BNDNDS(p)                    (p->nmb_of_BndNds)
#define EXCHNG_TYP1_SFE_ARRAY(p)                        (p->SFE_Array)
#define EXCHNG_TYP1_KOORDS(p)                           (p->n_koord_array)
#define EXCHNG_TYP1_RADIUS(p)                           (p->radius)
#define EXCHNG_TYP1_MIDPOINT_K(p,k)                     (p->MidPoint[k])

typedef struct {
  SF_TYP *rootsfc;       /*Root Pointer auf Surfaces*/
  SD_TYP *rootsd;       /*Root Pointer auf Subdomains*/
  PL_TYP *rootpl;       /*Root Pointer auf Polylines*/
  SFE_KNOTEN_TYP **SFE_HashTable;
  LI_KNOTEN_TYP **LI_HashTable;
} EXCHNG_TYP2; /* Typ fuer den Austausch von Daten zw. Cadconvert()
                  und Ansys2LGM */
#define EXCHNG_TYP2_ROOT_SFC(p)                 (p->rootsfc)
#define EXCHNG_TYP2_ROOT_SBD(p)                 (p->rootsd)
#define EXCHNG_TYP2_ROOT_PLY(p)                 (p->rootpl)
#define EXCHNG_TYP2_SFE_HASHTAB(p)              (p->SFE_HashTable)
#define EXCHNG_TYP2_LI_HASHTAB(p)               (p->LI_HashTable)







#define T 1
#define F 0
#define FERTIG 3 /*Aller guten Dinge sind ...*/
#define SEC_SFC_IDF_DEFAULT_VAL -1
#define SFE_KNID_4_DEFAULT_VAL -1
#define SF_RL_SBD_NOT_SET_YET -1
#define SEC_SFC_NAME_DEFAULT_VAL 0.0 /*geaendert in 0.0, da dann
                                        "Aussen in EvalLeftRightOfSfcs
                                        automatisch zu "0" werden*/
#define PI_HALBE 1.570796327
#define LI_IFS_ERROR 2
#define SORTED 0
#define REACHED 3 /*Aller guten Dinge sind ...*/

#define SPLIT_SURFACE_DISTINGUISH       0.0001 /*z.B. aus der in 5 Teile gesplitteten Surface 15.78
                                                                                        entstehen die 5 neuen SUrfaces
                                                                                        15.7801, 15.7802, 15.7803, 15.7804 und 15.7805*/

#define MAXLINE                         100
#define NMBOFCNDS                       30

#define CORNERS_OF_BND_SEG      3
#define CORNERS_OF_BND_SIDE     3
#define CORNERS_OF_ELEMENT      4

#define SEG_TRI                         1


#define MAX_NUB_OF_SBDMS        101

#define USE_SFC_DETECTOR        1

#define DIMENSION           3

/*#define STATISTICAL_INFORMATIONS      1*/
/* if defined  ===> characteristic values of the hashtable
    are written in ug_shell
    specified informations see
    http://dom.ica3.uni-stuttgart.de/~dirk/
 */
/*explanation of char*/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/
int LGM_ANSYS_ReadDomain                        (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info, INT MarkKey);
int LGM_ANSYS_ReadSizes                         (LGM_SIZES *lgm_sizes);
int LGM_ANSYS_ReadSubDomain         (int i, LGM_SUBDOMAIN_INFO *subdom_info);
int LGM_ANSYS_ReadLines                         (int i, LGM_LINE_INFO *line_info);
int LGM_ANSYS_ReadPoints                        (LGM_POINT_INFO *lgm_point_info);
int LGM_ANSYS_ReadSurface                       (int i, LGM_SURFACE_INFO *surface_info);
int LGM_ANSYS_ReadMesh              (HEAP *Heappointer, LGM_MESH_INFO *theMesh, int MarkKey);
