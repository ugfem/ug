// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ansys2lgm.c                                                   */
/*                                                                          */
/* Purpose:   convert Ansys-files to the UG-Domain-structure		   */
/*                                                                          */
/* Author:    Dirk Feuchter                                                                     */
/*			  Institut fuer Computeranwendungen III	                                        */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27					*/
/*			  70569 Stuttgart, Germany				*/
/*			  email: ug@ica3.uni-stuttgart.de			*/
/*										*/
/* History:   27.5.97 begin, ug version 3.7	                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/


#include <config.h>
#include <assert.h>

#include "defaults.h"
#include "domain.h"
#include "lgm_domain.h"
#include "lgm_transfer.h"
#include "heaps.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#include "ugtypes.h"
#include "fileopen.h"
#include "general.h"
#include "heaps.h"
#include "gm.h"
#include "misc.h"

#include "ugdevices.h"
#include "cmdline.h"
/*#include "problem.h"*/

/*#include "readcadfile.h"*/
#include "heaps.h"
/* #include "construct.h" */


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#include "ansys2lgm.h" /*%%%%%*/
/*#include "cadconvert.h"*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "namespace.h"

USING_UG_NAMESPACE
  USING_UGDIM_NAMESPACE

#define NUOFCLMS 70
#define NU_SFCES_BNDP 9
/*LINE_NEU*/
#define NU_LINES_BNDP 7


static EXCHNG_TYP1 ExchangeVar_1;
static EXCHNG_TYP1 *ExchangeVar_1_Pointer;
static EXCHNG_TYP2 ExchangeVar_2;
static EXCHNG_TYP2 *ExchangeVar_2_Pointer;
static INT SFE_p, LI_p  ;
static INT nmb_of_trias_of_sf; /*wird verwendet in Ansys2lgmCreateTriaOrientations*/

static INT nbofnds;

/*just for debugging*/
SD_TYP *sd_global;


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
/* see paramlist of Ansys2lgm and see Ansys2lgmInit */
/*static SFE_KNOTEN_TYP *SFE_HashTable;*/
/*static LI_KNOTEN_TYP *LI_HashTable;*/

/*neu EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer) */
/*static SF_TYP *rootsfc;*/ /*Root Pointer auf Surfaces*/


/*neu EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer) */
/*static SD_TYP *rootsd;*//*Root Pointer auf Subdomains*/

/*neu EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer) */
/*static PL_TYP *rootpl; *//*Root Pointer auf Polylines*/


/*sowie deren Subfunktionen als Vergleichszahl*/

static INT zaehler;   /*wird verwendet in Ansys2lgmCreateTriaOrientations
                         sowie deren Subfunktionen als Zaehler, bzw. Abbruchkriterium
                         dafuer, ob alle Triangles einer Surface schon einmal durchlaufen worden
                                sind.*/
static DOMAIN_INFO_TYP DomainInfo;
static DOMAIN_INFO_TYP *DomainInfo_Pointer;

static char ProblemName[31];

double ZoomFactorX;
double ZoomFactorY;
double ZoomFactorZ;

INT komponentenzaehler;
INT *KomponentenIndexArray; /* diese Array beinhaltet die IDs der einzelnen
                               Komponenten aus dem CAD File,
                               die am Ende der Datei in der Form
                                  K,1,salzgebiet
                               ergaenzt wurden.
                               Das array wir von 1 bis AnzahlAngaben gefuellt
                               und hat die selbe Reihenfol�ge wie in der CAD-DAtei*/
/* DIRKS NEU : static bisherige_ID_array */
static INT *bisherige_ID_array;    /* diese Array beinhaltet die Subdomain/Komponenten IDs
                                      aus dem CAD-Array ==> gefuellt sin die Plaetze 1
                                      bis AnzahlSbdms. Die Reihenfolge bezieht sich auf die
                                      SubdomainIDs des UGs !!!*/
char *KomponentenNamenArray;

static HEAP *theHeap;
static INT ANS_MarkKey;

static INT ansysfilepathes_set;


static INT triangle_found; /*used in Create_RealSurfaces and Find_SFE_Triangle*/

static char *TmpMemArray;

static TRIANGLE_TYP *New_Triangle_List;
INT nmb_of_triangles;

/* DIRKS NEU*/
static INT NuClmsMINUS1;

/* DIRKS NEU ... */
static INT *node_element_matrix;
/*      possesses for each node a list of all adjacent elements                 */
static INT *el_array;
/*      possesses all elements; each element is represented by 8 numbers:       */
/*      NodeID_i, NodeID_j, NodeID_k, NodeID_l,                                 */
/*      NbourElem_ijk, NbourElem_ijl, NbourElemt_jkl and NbourElem_ilk          */
static DOUBLE *n_koord_array_UG;
static INT *point_array;
static INT *point_array_UG_CAD;
static int nmbOfTetrhdrOfThisSbd;
static int nmbOfSidesOfThisSbd;
static INT *el_besucht_array;
static INT *elemflag_array;
/*Bei der Bestimmung der Tetraeder-Subdomain-Zugehoerigkeit
   wird el_besucht_array als Kennzeichnungsfeld verwendet
   Jedes Tetraederelement erhaelt dort seine zugehoerige SubdomainID.
   Initialisierung mit 0
   In nmbOfTetrhdrOfThisSbd wird die Anzahl der Tetraeder je Subdomain gezaehlt.*/

/* DIRKS NEU statistik[6]*/
static INT statistik[7];
/*
        statistic informations about the CAD-file :
        ===========================================

        statistik[0], used for number of inner nodes
        statistik[1], used for number of boundary nodes
        statistik[2], used for number of inner elements
        statistik[3], used for number of boundary elements
        statistik[4], used for number of boundary segments
        statistik[5], used for number of different boundary conditions
        statistik[6], used for number of elements == ( statistik[2] + statistik[3] )
 */
/* ... DIRKS NEU */

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);



#ifdef STATISTICAL_INFORMATIONS
static int ST_INF_Anz_SFE;
static int ST_INF_Anz_ds;
static int ST_INF_Anz_SFE_real;
static int ST_INF_m;
static int ST_INF_gef;
static double ST_INF_fg;
static double ST_INF_Kollishf;
static int ST_INF_2er;
static double ST_INF_2er_P;
static int ST_INF_3er;
static double ST_INF_3er_P;
static int ST_INF_4er;
static double ST_INF_4er_P;
static int ST_INF_5er;
static double ST_INF_5er_P;
static int ST_INF_maxK;
static int ST_INF_Klsstellen;
static int ST_INF_Kollis;
static int ST_INF_Anz_LI;
/*static int ST_INF_Anz_dl;*/
static int ST_INF_Anz_LI_real;
static double ST_INF_LI_real_durch_SFE_real;
static int ST_INF_mL;
static int ST_INF_gefL;
static double ST_INF_fgL;
static double ST_INF_KhL;
static int ST_INF_2erL;
static double ST_INF_2erL_P;
static int ST_INF_3erL;
static double ST_INF_3erL_P;
static int ST_INF_4erL;
static double ST_INF_4erL_P;
static int ST_INF_5erL;
static double ST_INF_5erL_P;
static int ST_INF_maxKL;
static int ST_INF_KlsstellenL;
static int ST_INF_KollisL;
#endif



/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/


/****************************/
/*aus ehemals readcadfile.c:*/




/****************************************************************************/
/**
*\brief reads the CAD file and fills the arrays od the intermediate format.
*
* This function reads the essential and necessary datas from the CAD OutputFile
* and fills different arrays. These arrays describe an intermediate format
* which is used by ConvertCADGrid to build the UG-Grid with its boundary segments
* and boundary descriptions.
*
* @param CADOutputFileName  name of the CAD Outputfile.
* @param node  pointer to an array used for the intermediate format.
* @tetrahedron  pointer to an array used for the intermediate forma.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
/************************************************************************************/
/*	INPUT  : actual file pointer concerning the CADOutputFile                                               */
/*	OUTPUT : 0 if oK																*/
/*		        complete line of the CADOutputFile from **actfp to  next '\n' given in	*/
/*			refenrence parameter buffer												*/
/************************************************************************************/

/*/ / ReadLine already exists in ANSI-C twice : getline() or fgets() !!!*/

INT ReadLine(char *linebuffer, INT mxln, FILE *filePtr)
{
  INT ind = 0;

  do
  {
    linebuffer[ind] = fgetc(filePtr);
    ind++;
  }
  while((linebuffer[ind-1]!=10)&&(linebuffer[ind-1]!=13)&&(ind < MAXLINE));

  return(ind);
}

/*/ / ReadLine already exists in ANSI-C twice : getline() or fgets() !!!*/


/*/ / ReadLine already exists in ANSI-C twice : getline() or fgets() !!!*/


/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/* purpose: function reads zoomfactor for z direction				*/
/************************************************************************************/
                                                int ZoomFct(char *linebuffer,char axis)
{
  INT index,elem_surf_cond;
  char *endp,*s;

  index = 1;
  elem_surf_cond = 0;


  endp = &(linebuffer[index]);

  s = &(endp[1]);

  switch(axis)
  {
  case 'x' :
    ZoomFactorX=(double)(strtod(s,&endp));
    if(ZoomFactorX <0.0)
    {
      PrintErrorMessage('E',"ZoomFct","ZoomFactorX ivalid use ZF_X<0.0");
      return(1);
    }
    break;
  case 'y' :
    ZoomFactorY=(double)(strtod(s,&endp));
    if(ZoomFactorY <0.0)
    {
      PrintErrorMessage('E',"ZoomFct","ZoomFactorY ivalid use ZF_Y<0.0");
      return(1);
    }
    break;
  case 'z' :
    ZoomFactorZ=(double)(strtod(s,&endp));
    if(ZoomFactorZ <0.0)
    {
      PrintErrorMessage('E',"ZoomFct","ZoomFactorZ ivalid use ZF_Z<0.0");
      return(1);
    }
    break;
  default :
  {
    PrintErrorMessage('E',"ZoomFct","axis must be x, y or z");
    return(1);
  }

  }


  return(0);
}





/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/* purpose: function reads names for SUbdomains resp. CAD-components				*/
/************************************************************************************/
int KomponentFct(char *linebuffer)
{
  INT index,elem_surf_cond;
  char *endp,*s;

  INT i,offset,stop;

  index = 1;
  elem_surf_cond = 0;


  endp = &(linebuffer[index]);

  s = &(endp[1]);
  komponentenzaehler++;
  if(komponentenzaehler == MAX_NUB_OF_SBDMS)
  {
    PrintErrorMessage('E',"cadconvert"," Komponentenzaehler bigger than MAX_NUB_OF_SBDMS");
    return(1);
  }
  KomponentenIndexArray[komponentenzaehler]=(INT)(strtol(s,&endp,10));

  s = &(endp[1]);
  offset = komponentenzaehler * 31;

  /*KomponentenNamenArray[komponentenzaehler * 31]=String einlesen;*/

  i=0;
  stop =0;

  while((s[i] != '\n') && (stop == 0))
  {
    if (i == 30)
    {
      PrintErrorMessage('W',"cadconvert"," KomponentName in ansFile is too long=> use first 30 bytes");
      KomponentenNamenArray[offset]  = '\0';
      i++;
      stop = 1;
    }
    else
    {
      KomponentenNamenArray[offset]  = (s[i]);
      offset++;
      i++;
    }
  }
  if(stop==0)
  {
    KomponentenNamenArray[offset]  = '\0';
  }
  return(0);
}



/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/*	changed Parameter:                                                                                              */
/* purpose: function reads problemname                                                                                  */
/************************************************************************************/
INT ProbNameFct(char *linebuffer)
{
  INT index,xyz,i,stop;

  index = 0;
  xyz = 0;

  do
  {
    index++;
  }
  while(linebuffer[index] != ',');
  index++;

  i=0;
  stop =0;

  while((linebuffer[index] != '\n') && (stop == 0))
  {
    if (i == 30)
    {
      PrintErrorMessage('W',"cadconvert"," Problemname in ansFile is too long=> use first 30 bytes");
      ProblemName[i] = '\0';
      i++;
      stop = 1;
    }
    else
    {
      ProblemName[i] = (linebuffer[index]);
      index++;
      i++;
    }
  }
  if(stop==0)
  {
    ProblemName[i] = '\0';
  }


  return(0);
}



/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/*	changed Parameter: n_koord_arrays												*/
/* purpose: function reads nodcoordinates                                                                                       */
/************************************************************************************/
INT NodeLineFct(INT nz,DOUBLE *n_koord_array,char *linebuffer)
{
  INT index,xyz;
  char *endp,*s;


  index = 2;
  xyz = 0;

  do
  {
    index++;
  }
  while(linebuffer[index] != ',');

  endp = &(linebuffer[index]);

  while((endp[0] != 10)&&(endp[0] != 13))
  {
    s = &(endp[1]);
    n_koord_array[nz*3+xyz]=(DOUBLE)(strtod(s,&endp));
    xyz++;
  }
  return(0);
}



/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/*	changed Parameters: el_array and  node_element_matrix							*/
/* purpose: function reads element node context and fills node element matrix           */
/************************************************************************************/
int ElementLineFct(INT *ez,INT *el_array,INT *node_element_matrix,char *linebuffer)
{
  INT index,ijkl,clm;
  char *endp,*s;
  INT help1;
  INT help2;

  index = 3;
  ijkl = 0;

  /*only if linebuffer isnt the "ET-line"*/

  if (linebuffer[1] != 'T')
  {

    do
    {
      index++;
    }
    while(linebuffer[index] != ',');

    endp = &(linebuffer[index]);

    while((endp[0] != 10)&&(endp[0] != 13))
    {
      s = &(endp[1]);
      help1 = (*ez)*8+ijkl;
      help2 = el_array[help1] = (INT)(strtol(s,&endp,0));
      help2 *= NUOFCLMS;
      clm = 0;
      while (node_element_matrix[help2+clm] != 0)
      {
        clm ++;
        /* DIRKS NEU ... :*/
        if (clm == NuClmsMINUS1)
        {
          PrintErrorMessage('E',"ElementLineFct","more than NUOFCLMS-1  elements corresponding to one node");
          return(1);
        }
        /* ... DIRKS NEU */
      }
      node_element_matrix[help2+clm] = *ez;
      ijkl++;
    }

  }
  else
    (*ez)--;


  return(0);
}



/************************************************************************************/
/*	INPUT  : linebuffer,  current line  of the CADOutputFile                                                */
/*	changed Parameters: el_array,bndseg_array and  nodeflag_array					*/
/*           besides statistik[5] and condition_array								*/
/*  st...[5]: number of different boundary conditions                                                           */
/*  condition_array: the different bndcndnumbers of the CAD-user which                          */
/*					characterisize the surface of the 3Dgeometry					*/
/* purpose: function reads SFEs later on used for boundary segments,				*/
/*                      selects bondary nodes and  boundary elements                                                    */
/************************************************************************************/
int SurfaceLoadFct(INT sz,INT *statistik, INT *condition_array,BND_SFE_TYP *bndseg_array,INT *nodeflag_array,INT *elemflag_array,INT *el_array,char *linebuffer)
{
  INT w,u;
  INT index,elem_surf_cond;
  char *endp,*s;

  INT help;
  DOUBLE aid;
  help = sz;

  index = 3;
  elem_surf_cond = 0;


  endp = &(linebuffer[index]);

  s = &(endp[1]);
  BND_SFE_ELEMENTID((&(bndseg_array[help])))=(INT)(strtol(s,&endp,10));

  s = &(endp[1]);
  BND_SFE_SIDEID((&(bndseg_array[help])))=(INT)(strtol(s,&endp,10));

  /* mark element in elemflag_array :  */
  elemflag_array[BND_SFE_ELEMENTID((&(bndseg_array[help])))] = TRUE;

  do
    endp++;
  while(endp[0] != ',');      /*Komma after PRESS*/
  do
    endp++;
  while(endp[0] != ',');      /*Komma after number after PRESS*/

  s = &(endp[1]);
  aid = BND_SFE_SFCIDF((&(bndseg_array[help]))) = strtod(s,&endp);


  /*/ / concerning number of bnd conditions and the condition_array : ....*/
  /*NOT YET
          ind = 0;
          already = FALS;
          while((condition_array[ind] != 0) && (already == FALS))
          {
                  if (condition_array[ind] == aid)
                          already = TRU;
                  ind++;
          }
     NOT YET */
  /*NOT YETif (already == FALS)NOT YET */               /*// still !*/
  /*NOT YET
          {
                  condition_array[ind] = aid;
                  statistik[5]++;
     NOT YET */
  /*NOT YET
                  bndcndflag_array[sz-1] = ind;

          }
          else
     NOT YET */
  /*NOT YET
                  bndcndflag_array[sz-1] = ind-1;
     NOT YET */
  /*/ /  .... concerning number of bnd conditions and the condition_array* /
  /*NOT YET Achtung das geht jetzt nicht mehr! muesste aber auch einfach mit
            Zuweisung von -1 gehen, ich habe zumindest keine Stelle gefunden wo wirklich auf
            den negierten SFE-IdentifierValue zugewiesen wird.
          el_array[8*bndseg_array[help]+3+bndseg_array[help+1]]=(-1)*bndseg_array[help+2];
     NOT YET */

  /* DIRKS NEU : */
  /*	el_array[8*BND_SFE_ELEMENTID((&(bndseg_array[help])))+3+BND_SFE_SIDEID((&(bndseg_array[help])))]= -1;*/
  /* ... DIRKS NEU schon wieder alt : siehe in ReadANsysFile */


  /* BND_SFE_ELEMENTID((&(bndseg_array[help])))  = zugehoerige ElementID */
  /* BND_SFE_SIDEID((&(bndseg_array[help]))) = zugehoerige SideID */
  /* BND_SFE_SFCIDF((&(bndseg_array[help]))) = zugehoerige Lastzahl */
  /*help = sz = sfe-zaehler*/


  w = BND_SFE_SIDEID((&(bndseg_array[sz])));
  u = BND_SFE_ELEMENTID(8*(&(bndseg_array[sz])));
  switch(w)
  {
  case 1 :              /*ijk012*/
    nodeflag_array[el_array[u]]=TRUE;
    nodeflag_array[el_array[1+u]]=TRUE;
    nodeflag_array[el_array[2+u]]=TRUE;
    break;
  case 2 :              /*ijl013*/
    nodeflag_array[el_array[u]]=TRUE;
    nodeflag_array[el_array[1+u]]=TRUE;
    nodeflag_array[el_array[3+u]]=TRUE;
    break;
  case 3 :              /*jkl123*/
    nodeflag_array[el_array[1+u]]=TRUE;
    nodeflag_array[el_array[2+u]]=TRUE;
    nodeflag_array[el_array[3+u]]=TRUE;
    break;
  case 4 :              /*ikl023*/
    nodeflag_array[el_array[u]]=TRUE;
    nodeflag_array[el_array[2+u]]=TRUE;
    nodeflag_array[el_array[3+u]]=TRUE;
  }

  return(0);
}




INT ReadCADFile(char *CADOutputFileName, INT *statistik, DOUBLE *abx,INT *condition_array, INT *nodeflag_array,INT *elemflag_array,INT *el_array,DOUBLE *n_koord_array,BND_SFE_TYP *bndseg_array,INT *node_element_matrix)
{
  INT nz,ez,sz;

  INT z,zm,rgw;

  INT help1,help2;

  char linebuffer[80];
  INT linebufferlength;

  FILE *filePtr;
  /*
          ansysfilepathes_set = 0;

          if (ReadSearchingPaths(DEFAULTSFILENAME,"ansysfilepathes")==0)
                  ansysfilepathes_set = 1;
   */
  /* reading */
  nz=ez=sz=0;

  if (ansysfilepathes_set)
    filePtr = FileOpenUsingSearchPaths(CADOutputFileName,"r","ansysfilepathes");
  else
    filePtr = fileopen(CADOutputFileName,"r");
  if (filePtr==NULL)
  {
    UserWriteF("cannot open file %s\n",CADOutputFileName);
    return(1);
  }

  /*/ /    linebufferadrs = fgets(linebuffer,MAXLINE,filePtr);*/
  linebufferlength = ReadLine(linebuffer,MAXLINE,filePtr);

  /*	if (ReadLine(&actfp,linebuffer) != 0)
          {
                  PrintErrorMessage('E',"ReadCADFile/ReadLine","execution failed");
                  return (CMDERRORCODE);
          }
   */

  /*DIRKS NEU :	*/
  NuClmsMINUS1 = NUOFCLMS -1 ;

  while(linebuffer[0] != 'F')
  {
    switch(linebuffer[0])
    {
    case 'X' :
      if ((rgw =ZoomFct(linebuffer,'x')) != 0)
      {
        PrintErrorMessage('E',"ZoomFct","execution failed");
        return(1);
      }
      break;
    case 'Y' :
      if ((rgw =ZoomFct(linebuffer,'y')) != 0)
      {
        PrintErrorMessage('E',"ZoomFct","execution failed");
        return(1);
      }
      break;
    case 'Z' :
      if ((rgw =ZoomFct(linebuffer,'z')) != 0)
      {
        PrintErrorMessage('E',"ZoomFct","execution failed");
        return(1);
      }
      break;
    case 'K' :
      if ((rgw =KomponentFct(linebuffer)) != 0)
      {
        PrintErrorMessage('E',"KomponentFct","execution failed");
        return(1);
      }
      break;
    case 'P' :
      if ((rgw =ProbNameFct(linebuffer)) != 0)
      {
        PrintErrorMessage('E',"ProbNameFct","execution failed");
        return(1);
      }
      break;
    case 'N' :       nz++;
      if ((rgw =NodeLineFct(nz,n_koord_array,linebuffer)) != 0)
      {
        PrintErrorMessage('E',"NodeLineFct","execution failed");
        return(1);
      }
      /* evaluations for the bounding box*/
      help1 = nz*3;
      for(z=0; z<3; z++)
      {
        help2 = help1+z;
        zm = z+3;
        if(n_koord_array[help2] < abx[z])
          abx[z] = n_koord_array[help2];
        if(n_koord_array[help2] > abx[zm])
          abx[zm] = n_koord_array[help2];
      }
      /* evaluations for the bounding box*/
      break;
    case 'E' :       ez++;
      if((rgw = ElementLineFct(&ez,el_array,node_element_matrix,linebuffer))!=0)
      {
        PrintErrorMessage('E',"ElementLineFct","execution failed");
        return(1);
      }
      break;
    case 'S' :       sz++;
      if((rgw = SurfaceLoadFct(sz,statistik,condition_array,bndseg_array,nodeflag_array,elemflag_array,el_array,linebuffer)) != 0)
      {
        PrintErrorMessage('E',"SurfaceLoadFct","execution failed");
        return(1);
      }
    }

    /*/ /            linebufferadrs = fgets(linebuffer,MAXLINE,filePtr);*/
    linebufferlength = ReadLine(linebuffer,MAXLINE,filePtr);
    /*		if (ReadLine(&actfp,linebuffer) != 0)
                    {
                            PrintErrorMessage('E',"ReadCADFile/ReadLine","execution failed");
                            return (CMDERRORCODE);
                    }
     */
  }

  /* DIRKS NEU : An dieser Stelle waere eine statistische Fuellgradfunktion fuer die
     node_element_matrix ganz interessant um zu entscheiden ob ein Umstieg auf Listen sinvoll waere,
     bzw um NUOFCLMS ideal einzustellen.*/

  fclose(filePtr);




  /* sort element array */

  /*complete statistics*/
  for(z=1; z<=statistik[0]; z++)
  {
    if(nodeflag_array[z] ==1)             /* = boundary node */
      statistik[1]++;                             /* increment number of boundary nodes*/
  }
  statistik[0]-=statistik[1];      /* evaluate number of inner nodes: total number of nodes minus number of boundary nodes */
  for(z=1; z<=statistik[2]; z++)
  {
    if(elemflag_array[z] ==1)             /* = boundary element */
      statistik[3]++;                             /* increment number of boundary elements*/
  }
  /*DIRKS NEU ... */
  statistik[6]=statistik[2];       /* merke die GEsamtAnzahl der Elemente*/
  statistik[2]-=statistik[3];      /* evaluate number of inner elements: total number of elements minus number of boundary elements */

  return(0);
}





/***************************/
/*aus ehemals cadconvert.c:*/
/****************************************************************************/
/**
*\brief builds a complete UG-Grid out of a CAD-Output-File (ANSYS-format).
*
*  This function builds a complete UG-Grid out of a the "IntermediateFormat",
*  which was created by the function ReadCADFile(...);
*  Thereby the multigrid-structure is filled completely.
*  BoundarySegments/Patches are generated automatically.
*  BoundaryConditions are created automatically.
*  Multigrid, BndSegments/Patches and all BndConditions are created out of an
*  ANSYS-file. At the moment we use the ANSYS-output from the CAD-Software
*  "Pro/ENGINEER" (from Parametrics), available at the "StutCAD" of the university
*  Stuttgart. The grids are created with "Pro/Mesh" ,the tetrahedron mesh genarator
*  of "Pro/ENGINEER". Besides in the CAD-program each "geometrical patch"
*  (not each surface triangle!) gets a "significant number".
*  The OUTPUT-ANSYS-file possesses all the (3D-tetrahedron)grid informations and boundary
*  informations.
*  For each surface triangle of the (3D-tetrahedron)grid "ConvertCADGrid"
*  creates a BoundarySegment/Patch.
*  Not for each surface triangle but for each "geometrical patch" (= each "significant number")
*  "ConvertCADGrid" creates a BoundaryCondition (default:GenBndConditionDc = general Dirichlet-
*  Boundary Condition with x=0.0, y=0.0 and z=0.0  ===>  defined in construct.c). Other
*  BoundaryConditions can be easily  added in the same way as GenBndConditionDc.
*  BoundaryConditions and BoundarySegments/Patches are matched.
*  (more about the "IntermediateFormat" ==> ask Dirk Feuchter)
*
* @param Var_1  Variable of ExchangeTyp  must be filled ( ==> INT nmb_of_SFEs;
*                                       CAD_SFE_TYP *SFE_Array;
*                                       DOUBLE *n_koord_array;).
* @param filename name of the CAD-Output-File.
*
* @return  MULTIGRID *theMG if ok,  NULL if error occured.
   */
/****************************************************************************/
/*MULTIGRID *ConvertCADGrid  (char *theFormat,char *CADOutputFileName,unsigned long heapSize)*/
INT ReadAnsysFile(char *filename)
{
  FILE *filePtr;

  char linebuffer[80];
  char s[5];

  INT point[CORNERS_OF_BND_SEG];
  INT linebufferlength;
  INT totalels;
  INT offsSGI;
  INT offsHP;
  INT nbofelms;
  INT cntbnodes,cntinodes;

  INT ii,rv,iminus1;

  INT ppp;
  DOUBLE nnnn;

  /* abutmentbox */
  DOUBLE abx[6];
  DOUBLE radius,MidPoint[3];

  unsigned long memforif;



  /********************************************************************************/
  /*      arrays for the "IntermediateFormat" ===> . . .                          */
  /*                                                                              */
  /*      more about the "IntermediateFormat" ==> ask Dirk Feuchter               */
  /*                                                                              */
  INT *nodeflag_array;
  /*      distinguishes between InnerNodes and BoundaryNodes                      */
  /*                                                                              */
  /*	INT *point_array;  DIRKS NEU static*/
  /*      possesses the correct node IDs for UG                                   */
  NODE **UGID_NdPtrarray;
  /*  possesses the correct node adresses for UG */
  /*                                                                              */
  /*INT *elemflag_array;	DIRKS NEU: jetzt statisch !								*/
  /*      distinguishes between InnerElements and BoundaryElements                */
  /*                                                                              */
  /*                                                                              */
  BND_SFE_TYP *bndseg_array;
  /*      = "Surface-Part" of the ANSYS-File                                      */
  /*                                                                              */
  CAD_SFE_TYP *CAD_SFE_Part;
  /*      = "CAD Surface-Part" of the ANSYS-File mit jeweils 4.Idf als zusaetliche*/
  /*      Angabe, Angabe in UG-IDs                                                                              */
  /*                                                                              */
  INT condition_array[NMBOFCNDS];
  /*      stores the "significant numbers" (see Description), which the CAD-User  */
  /*      added to the different "geometrical patches"                            */
  /*                                                                              */
  DOUBLE *bndcndflag_array;
  /*      matches each Patch(resp. SurfaceTriangle/BoundarySegment) with the      */
  /*      the corresponding BoundaryCondition                                                                   */
  /*                                                                              */
  /*                                                                              */
  DOUBLE *n_koord_array;
  /*      coordinates in the sequence of the CAD-file                             */
  /*                                                                              */
  /* DIRKS NEU static   DOUBLE *n_koord_array_UG; */
  /*      UGID - coordinates od the CAD-file                                                    */
  /*                                                                              */
  DOUBLE *koord_array;
  /*      used for the nine koordinates necessary for one boundary segment!       */
  /*      reused in GeneralBoundary                                               */
  /*                                                                              */
  /*      . . . <=== arrays for the IntermediateFormat                            */
  /********************************************************************************/


  INT help1,ofs0,ofs1,ofs2,ofs3,i,dummy, nmofnds,ug_side_offset;

  if(ansysfilepathes_set != 1)       /*wenn noch nicht gesetzt*/
  {
    if (ReadSearchingPaths(DEFAULTSFILENAME,"ansysfilepathes")==0)
      ansysfilepathes_set = 1;
  }


  dummy=0;


  /* Init of statistik[i] and abx */
  abx[0]=1.0E+38;                       /*used for x-coordinate of minimum corner of the bounding box*/
  abx[1]=1.0E+38;                       /*used for y-coordinate of minimum corner of the bounding box*/
  abx[2]=1.0E+38;                       /*used for z-coordinate of minimum corner of the bounding box*/
  abx[3]=-1.0E+38;                      /*used for x-coordinate of maximum corner of the bounding box*/
  abx[4]=-1.0E+38;                      /*used for y-coordinate of maximum corner of the bounding box*/
  abx[5]=-1.0E+38;                      /*used for z-coordinate of maximum corner of the bounding box*/
  statistik[0]=0;                       /*used for number of inner nodes*/
  statistik[1]=0;                       /*used for number of boundary nodes*/
  statistik[2]=0;                       /*used for number of inner elements*/
  statistik[3]=0;                       /*used for number of boundary elements*/
  statistik[4]=0;                       /*used for number of boundary segments*/
  statistik[5]=0;                       /*used for number of different boundary conditions*/
  /* DIRKS NEU */
  statistik[6]=0;                       /*used for number of elements == statistik[2] + statistik[3]*/



  /* opening CADFile */

  if (ansysfilepathes_set)
    filePtr = FileOpenUsingSearchPaths(filename,"r","ansysfilepathes");
  else
    filePtr = fileopen(filename,"r");
  if (filePtr==NULL)
  {
    UserWriteF("cannot open file %s\n",filename);
    return(1);
  }


  /*	filePtr = fileopen(filename,"r");*/


  linebufferlength = ReadLine(linebuffer,MAXLINE,filePtr);
  /* reading file only to get statistic informations */
  while(linebuffer[0] != 'F')
  {
    switch(linebuffer[0])
    {
    case 'N' :       statistik[0]++;break;
    case 'E' :       if(linebuffer[1] != 'T') statistik[2]++;break;
    case 'S' :       statistik[4]++;
    }
    linebufferlength = ReadLine(linebuffer,MAXLINE,filePtr);
  }
  fclose(filePtr);
  /* closing CADFile */


  memset(condition_array,0,NMBOFCNDS*sizeof(INT));


  /* necessary memory :*/
  /*offsMac   = 184;*/
  offsHP   = 256;
  /*offsSGI = 1250;*/
  memforif =  (statistik[0]+1)*((2+NUOFCLMS)*sizeof(INT) + sizeof(NODE*) + 3*sizeof(DOUBLE)) +
             (statistik[2]+1)*(9*sizeof(INT)) +
             (statistik[4]+1)*(9*sizeof(DOUBLE)+3*sizeof(INT)+sizeof(CAD_SFE_TYP)+sizeof(BND_SFE_TYP)) +
             (statistik[4])*(sizeof(INT)) +
             offsSGI +
             sizeof(BLOCK);



  /************************************************************************/
  /* Get memory for different  arrays for the "IntermediateFormat ... --> */
  /************************************************************************/
  /* nodeflag_array */
  nodeflag_array = (INT *) GetTmpMem(theHeap,(statistik[0]+1)*sizeof(INT),ANS_MarkKey);
  if ( nodeflag_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(nodeflag_array,0,(statistik[0]+1)*sizeof(INT));
  /***************/
  /* point_array */
  nbofnds = statistik[0];
  point_array = (INT*)GetTmpMem(theHeap,(statistik[0]+1)*sizeof(INT),ANS_MarkKey);
  if ( point_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(point_array,0,(statistik[0]+1)*sizeof(INT));
  /***************/
  /* point_array_UG_CAD */
  point_array_UG_CAD = (INT*)GetTmpMem(theHeap,(statistik[0])*sizeof(INT),ANS_MarkKey);
  if ( point_array_UG_CAD == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory for point_array_UG_CAD");
    return(1);
  }
  memset(point_array,0,(statistik[0])*sizeof(INT));
  /*******************/
  /* UGID_NdPtrarray */
  UGID_NdPtrarray = (NODE **)GetTmpMem(theHeap,(statistik[0]+1)*sizeof(NODE*),ANS_MarkKey);
  if ( UGID_NdPtrarray == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory for UGID_NdPtrarray");
    return(1);
  }
  memset(UGID_NdPtrarray,0,(statistik[0]+1)*sizeof(NODE*));
  /******************/
  /* elemflag_array */
  totalels = statistik[2];
  elemflag_array = (INT*)GetTmpMem(theHeap,(totalels+1)*sizeof(INT),ANS_MarkKey);
  if ( elemflag_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(elemflag_array,0,(totalels+1)*sizeof(INT));
  /************/
  /* el_array */
  nbofelms = statistik[2];
  el_array = (INT*)GetTmpMem(theHeap,(statistik[2]+1)*sizeof(INT)*8,ANS_MarkKey);
  if ( el_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(el_array,0,(statistik[2]+1)*sizeof(INT)*8);
  /****************/
  /* bndseg_array */
  bndseg_array = (BND_SFE_TYP *)GetTmpMem(theHeap,(statistik[4]+1)*sizeof(BND_SFE_TYP),ANS_MarkKey);
  if ( bndseg_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  for(i=1; i<=statistik[4]; i++)
  {
    BND_SFE_ELEMENTID((&(bndseg_array[i]))) = -1;
    BND_SFE_SIDEID((&(bndseg_array[i]))) = -1;
    BND_SFE_SFCIDF((&(bndseg_array[i]))) = -1.0;
  }
  /****************/
  /* CAD_SFE_Part */
  CAD_SFE_Part = ( CAD_SFE_TYP *)GetTmpMem(theHeap,(statistik[4])*sizeof(CAD_SFE_TYP),ANS_MarkKey);
  if ( CAD_SFE_Part == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  for(i=0; i<statistik[4]; i++)
  {
    CAD_SFE_ND_I((&(CAD_SFE_Part[i]))) = -1;
    CAD_SFE_ND_J((&(CAD_SFE_Part[i]))) = -1;
    CAD_SFE_ND_K((&(CAD_SFE_Part[i]))) = -1;
    CAD_SFE_ND_4((&(CAD_SFE_Part[i]))) = -1;
    CAD_SFE_SFE_IDF((&(CAD_SFE_Part[i]))) = -1.0;
  }

  /********************/
  /* bndcndflag_array */
  bndcndflag_array = (DOUBLE *)GetTmpMem(theHeap,(statistik[4])*sizeof(DOUBLE),ANS_MarkKey);
  if (statistik[4] != 0)
    if ( bndcndflag_array == NULL )
    {
      PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
      return(1);
    }
  memset(bndcndflag_array,-1,(statistik[4]+1)*sizeof(DOUBLE));
  /***********************/
  /* node_element_matrix */
  node_element_matrix = (INT*)GetTmpMem(theHeap,(statistik[0]+1)*sizeof(INT)*NUOFCLMS,ANS_MarkKey);
  if ( node_element_matrix == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(node_element_matrix,0,(statistik[0]+1)*sizeof(INT)*NUOFCLMS);
  /*****************/
  /* n_koord_array */
  n_koord_array = (DOUBLE *)GetTmpMem(theHeap,(statistik[0]+1)*sizeof(DOUBLE)*3,ANS_MarkKey);
  if ( n_koord_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(n_koord_array,-1,(statistik[0]+1)*sizeof(DOUBLE)*3);
  /*****************/
  /* n_koord_array_UG */
  n_koord_array_UG = (DOUBLE *)GetTmpMem(theHeap,(statistik[0])*sizeof(DOUBLE)*3,ANS_MarkKey);
  if ( n_koord_array_UG == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(n_koord_array_UG,-1,(statistik[0])*sizeof(DOUBLE)*3);
  /***************/
  /* koord_array */
  koord_array = (DOUBLE *)GetTmpMem(theHeap,(statistik[4]+1)*sizeof(DOUBLE)*9,ANS_MarkKey);
  if ( koord_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(koord_array,-1,(statistik[4]+1)*sizeof(DOUBLE)*9);
  /*              */
  /* ... <--- Get memory for different  arrays for the "IntermediateFormat*/
  /************************************************************************/
  /*Neu Initialisierung der Variablen fuer die SubdomainIdentifizierung mit Namen*/
  komponentenzaehler = 0;
  KomponentenIndexArray = (INT *)GetTmpMem(theHeap,(MAX_NUB_OF_SBDMS)*sizeof(INT),ANS_MarkKey);
  if ( KomponentenIndexArray == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(KomponentenIndexArray,-1,(MAX_NUB_OF_SBDMS)*sizeof(INT));

  bisherige_ID_array  = (INT *)GetTmpMem(theHeap,(MAX_NUB_OF_SBDMS)*sizeof(INT),ANS_MarkKey);
  if ( bisherige_ID_array == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  memset(bisherige_ID_array,-1,(MAX_NUB_OF_SBDMS)*sizeof(INT));

  KomponentenNamenArray = (char *)GetTmpMem(theHeap,(MAX_NUB_OF_SBDMS)*31*sizeof(char),ANS_MarkKey);
  if ( KomponentenNamenArray == NULL )
  {
    PrintErrorMessage('E',"cadconvert"," ERROR: No memory !!! error in cadconvertfunction <ConvertCADGrid>");
    return(1);
  }
  /*An der Stelle 0 des KomponentenNamenArray den NotSet-String eintragen :*/
  strcpy(KomponentenNamenArray, "< NOT SET IN ANSYS-FILE >");



  /*********************************************************************************************************/
  /*********************************************************************************************************/

  /*function "ReadCADFile" will now read the ANSY-FIle and fill different arrays of the "Int.Format"*/

  /*********************************************************************************************************/
  /*********************************************************************************************************/
  /*GOGON TODO HERE Achtung bndseg_array und bndcndflag_array haben neue Daten typen
     ab hier weiter ueberarbeiten zuerst die Fuktion ,dann ReadCADFile*/
  /*	OLD : if (ReadCADFile(CADOutputFileName,statistik,abx,condition_array,nodeflag_array,elemflag_array,el_array,
                          n_koord_array,bndseg_array,node_element_matrix,bndcndflag_array) != 0)
   */
  if ((rv = ReadCADFile(filename,statistik,abx,condition_array,nodeflag_array,elemflag_array,el_array, n_koord_array,bndseg_array,node_element_matrix)) != 0)
  {
    PrintErrorMessage('E',"ReadCADFile","execution failed");
    return (1);
  }

  /*********************************************************************************************************/
  /*********************************************************************************************************/




  /* calculating MidPoint and radius of abutment ball*/
  MidPoint[0] = (abx[0]+abx[3])/2;
  MidPoint[1] = (abx[1]+abx[4])/2;
  MidPoint[2] = (abx[2]+abx[5])/2;
  radius = sqrt( (abx[3]-abx[0])*(abx[3]-abx[0])/4 +
                 (abx[4]-abx[1])*(abx[4]-abx[1])/4 +
                 (abx[5]-abx[2])*(abx[5]-abx[2])/4
                 );
  EXCHNG_TYP1_MIDPOINT_K(ExchangeVar_1_Pointer,0) = MidPoint[0];
  EXCHNG_TYP1_MIDPOINT_K(ExchangeVar_1_Pointer,1) = MidPoint[1];
  EXCHNG_TYP1_MIDPOINT_K(ExchangeVar_1_Pointer,2) = MidPoint[2];
  EXCHNG_TYP1_RADIUS(ExchangeVar_1_Pointer) = radius;



  cntbnodes = 0; cntinodes = statistik[1];
  for(i=1; i<=nbofnds; i++)
  {
    if(nodeflag_array[i] == 0)                    /* inner node !!! */
    {
      point_array[i] = cntinodes;
      point_array_UG_CAD[cntinodes] = i;
      cntinodes++;
    }
    else                  /* bnd node !!! */
    {
      point_array[i] = cntbnodes;
      point_array_UG_CAD[cntbnodes] = i;
      cntbnodes++;
    }
  }
  /*Probe just for debugging :*/
  for(i=0; i<nbofnds; i++)
  {
    if(point_array_UG_CAD[i] == 0)
    {
      PrintErrorMessage('E',"ReadCADFile","point_array_UG_CAD contains 0!");
      return (1);
    }
  }


  for(i=1; i<=statistik[4]; i++)

  {

    /* DIRK NEU Setzen der -1 Eintraege bei den Elementnachbarn !*/
    /*bzw. genauer -1 * Vorkammastelle der zugewiesenen CAD-Last.*/
    switch(BND_SFE_SIDEID((&(bndseg_array[i]))))
    {
    case 1 : ug_side_offset = 4; break;
    case 2 : ug_side_offset = 7; break;
    case 3 : ug_side_offset = 5; break;
    case 4 : ug_side_offset = 6;
    }
    el_array[ug_side_offset + (8 * (BND_SFE_ELEMENTID((&(bndseg_array[i])))))] = (int)(-1 * (floor(BND_SFE_SFCIDF((&(bndseg_array[i]))))));
    /* . . . ENDE DIRKS NEU */

    iminus1 = i-1;
    s[3] = '\0';

    switch(BND_SFE_SIDEID((&(bndseg_array[i]))))
    {
    case 1 : ofs0 = 1; ofs1 = 0; ofs2 = 2; ofs3 = 3; break;
    case 2 : ofs0 = 0; ofs1 = 1; ofs2 = 3; ofs3 = 2; break;
    case 3 : ofs0 = 1; ofs1 = 2; ofs2 = 3; ofs3 = 0; break;
    case 4 : ofs0 = 2; ofs1 = 0; ofs2 = 3; ofs3 = 1;
    }

    help1 = BND_SFE_ELEMENTID((&(bndseg_array[i])))*8 ;


    point[0] = el_array[help1+ofs0];             /* node index of CAD !!!*/
    point[1] = el_array[help1+ofs1];             /* node index of CAD !!!*/
    point[2] = el_array[help1+ofs2];             /* node index of CAD !!!*/

    /*Concerning das fuer ansys2UG neue Feld  CAD_SFE_Part*/
    CAD_SFE_ND_I((&(CAD_SFE_Part[iminus1]))) = point_array[point[0]];
    CAD_SFE_ND_J((&(CAD_SFE_Part[iminus1]))) = point_array[point[1]];
    CAD_SFE_ND_K((&(CAD_SFE_Part[iminus1]))) = point_array[point[2]];
    CAD_SFE_ND_4((&(CAD_SFE_Part[iminus1]))) = point_array[el_array[help1+ofs3]];             /*der vierte Knoten, den ansys2Ug
                                                                                                 benoetigt um innen und aussen einer Surface zu unterscheiden !!!*/
    CAD_SFE_SFE_IDF((&(CAD_SFE_Part[iminus1]))) = BND_SFE_SFCIDF((&(bndseg_array[i])));

  }
  /**********************************************************************************************/

  /*n_koord_array_UG setzen :*/
  nmofnds = statistik[0]+statistik[1];
  for(ii =1; ii <= nmofnds ; ii++)
  {
    ppp = 3*point_array[ii];
    nnnn =  n_koord_array[ii*3];
    n_koord_array_UG[ppp] = nnnn;
    /*	n_koord_array_UG[point_array[ii]] = n_koord_array[ii*3];*/

    ppp = 3*point_array[ii]+1;
    nnnn = n_koord_array[ii*3+1];
    n_koord_array_UG[ppp] = nnnn;
    /*	n_koord_array_UG[point_array[ii]+1] = n_koord_array[ii*3+1];*/

    ppp = 3*point_array[ii]+2;
    nnnn = n_koord_array[ii*3+2];
    n_koord_array_UG[ppp] = nnnn;
    /*	n_koord_array_UG[point_array[ii]+2] = n_koord_array[ii*3+2];*/
  }
  EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer) = n_koord_array_UG; /*TODO UG-Indexes !!!*/
  EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer) = CAD_SFE_Part;
  EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer) = cntbnodes;
  EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer) = statistik[4]; /*TODO TODO TODO ??????*/

  return(0);

}






/****************************************************************************/
/**
* \brief hashfunctions for the SFE-Knoten-Eintraege.
*
*  berechnet aus den IDs der 3 Knoten eines Oberfl�chendreiecks
*  mit einer "Hashfunktion" einen Hashwert, der wiederum als Index
*  f�r die Abspeicherung verwendet wird.
*  Somit wird eine quadratische Alg.komplexit�t durch eine lineaere
*  optimiert bei der ldgl. der konstante Aufwand der Hashfunktion
*  hinzukommt.
*
* @param val1 niederwertigste KnotenID desOberfl�chendreiecks.
* @param val2 "mittlel"wertigste KnotenID desOberfl�chendreiecks.
* @param val3 hoechstwertigste KnotenID desOberfl�chendreiecks.
*
* @return the evaluated Hashvalue.
   */
/****************************************************************************/
INT The_SFE_hashfunction(INT val1, INT val2, INT val3)
{
  unsigned int hash_value, help;
  /*  beste Hasfunktion :*/
  /*  beste Hashfunktion: Fuellgrad bei m = 2 * #SFEs == 30 bis 40 Prozent
      je mehr CAD-Flaechen aufeinander liegen, desto niedriger ist der
      Fuellgrad .
      Die Kollisionshaeufigkeit betraegt 16 bis 21 Prozent . Der hoechste Wert
      wird bei einem EinzeolKomponentenGeometrie erreicht*/
  /* Die untersuchten Geometrien waren
     a2l_01_fein, a2l_100, a2l_80, a2l_72 und a2l_200 */

  hash_value = val1;
  hash_value = hash_value << 16;
  hash_value = hash_value ^ val2;
  help = val3;
  help = help << 8;
  hash_value = hash_value ^ help;
  hash_value = hash_value % SFE_p;

  /* andere Hashfunktionen:  auskommentiert !!!
         hash_value = val1 + val2 + val3 + val1 * val1 + val2 * val2 + val3 * val3;
          hash_value = hash_value % SFE_p;
   */


  /* noch eine andere Hashfunktion
          hash_value = val1;
          help = val2;
          bits = 1;
          while(help != 1)
          {
                  help = help >>1;
                  bits ++;
          }
          hash_value = hash_value << bits;
          hash_value = hash_value ^ val2;
          help = val3;
          bits = 1;
          while(help != 1)
          {
                  help = help >>1;
                  bits ++;
          }
          hash_value = hash_value << bits;
          hash_value = hash_value ^ val3;
          hash_value = hash_value % SFE_p;
   */

  return (INT)(hash_value);
  /*die eigentliche Hashfunktion - TODO Gibt es noch eine bessere ???*/
}


/****************************************************************************/
/**
*\brief hashfunctions for the LI-Knoten-Eintraege.
*
* berechnet aus den IDs der 2 Knoten einer "Oberfl�chenline"
* mit einer "Hashfunktion" einen Hashwert, der wiederum als Index
* f�r die Abspeicherung verwendet wird.
* Somit wird eine quadratische Alg.komplexit�t durch eine lineaere
* optimiert bei der ldgl. der konstante Aufwand der Hashfunktion
* hinzukommt.
*
* @param val1 niederwertige KnotenID derOberfl�chenLine.
* @param val2 hoeherwertige KnotenID derOberfl�chenLine.
*
* @return the evaluated Hashvalue.
   */
/****************************************************************************/
INT the_LI_hashfunction(INT val1, INT val2)
{
  unsigned int hash_value;
  /*unsigned int help, bits;*/
  /*  beste Hashfunktion: Fuellgrad bei m = 3 * #LIs == 30 bis 40 Prozent
      je mehr CAD-Flaechen aufeinander liegen, desto niedriger ist der
      Fuellgrad .
      Die Kollisionshaeufigkeit betraegt 13 bis 19 Prozent . Der hoechste Wert
      wird bei einem EinzeolKomponentenGeometrie erreicht */
  /*    Die untersuchten Geometrien waren
      a2l_01_fein, a2l_100, a2l_80, a2l_72 und a2l_200 */
  hash_value = val1;
  hash_value = hash_value << 16;
  hash_value = hash_value ^ val2;
  hash_value = hash_value % LI_p;


  /* andere Hash funktion
     hash_value = val1 + val2 + val1 * val1 + val2 * val2;
     hash_value = hash_value % LI_p;
   */
  /* weitere Hashfunktion
          hash_value = val1;
          help = val2;
          bits = 1;
          while(help != 1)
          {
                  help = help >>1;
                  bits ++;
          }
          hash_value = hash_value << bits;
          hash_value = hash_value ^ val2;
          hash_value = hash_value % LI_p;
   */


  return (INT)(hash_value);
  /*die eigentliche Hashfunktion - TODO Gibt es noch eine bessere
     oder TODO Untersuchungen/BenchmaarkTests/Statistiken ???*/
}


/****************************************************************************/
/**
*\brief initializations for Ansys2lgm.
*
* Die Sortierung ist notwendig, damit die SFE_Hashtabelle eindeutig bleibt und
* doppelt aufgef�hrte (im AnsysFile) Oberfl�chendreiecke erkannt werden !!!.
*
* @param bndsegids_array Das BoundarySegmentArray aus cadconvert
*                        wird hier so sortiert, das die 3 Nodeids eines SFEs aufsteigend sortiert
*                        sind. Achtung die 4-te ID steght fuer die KnotenID des zugehoerigen
*                        vierten innerern Knoten des jeweiligen Tetraeders: Diese CornerID
*                        wird bei der Sortierung natuerlich nicht beruecksichtigt.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT SortBndSegArray()
{
  INT lf,a,b,c,h;
  CAD_SFE_TYP *help1;
  CAD_SFE_TYP help2;
  for(lf=0; lf<(EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer)); lf++)
  {
    help1 = EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer);
    help2 = help1[lf];
    help2 = (EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf];
    help1 = &help2;
    a= CAD_SFE_ND_I(help1);
    a = CAD_SFE_ND_I((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    b = CAD_SFE_ND_J((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    c = CAD_SFE_ND_K((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));

    if ( (a<0) || (b<0) || (c<0) )
    {
      PrintErrorMessage('E',"SortBndSegArray","There are one or several ids with negative values !!");
      return(1);
    }

    if ( (a==b) || (a==c) || (b==c) )
    {
      PrintErrorMessage('E',"SortBndSegArray","There are twoids with the same value !!");
      return(1);
    }

    /*der Gr��e nach sortieren*/
    if (a>b)
    {
      h = a;
      a = b;
      b = h;
    }
    if (b>c)
    {
      h = b;
      b = c;
      c = h;
    }
    if (a>b)
    {
      h = a;
      a = b;
      b = h;
    }

    CAD_SFE_ND_I((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf]))) = a;
    CAD_SFE_ND_J((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf]))) = b;
    CAD_SFE_ND_K((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf]))) = c;

  }


  return(0);
}






/****************************************************************************/
/**
*\brief fetchs a good prime number.
*
* changes TheNumber to the next good prime number for a good Hashing.
*
* @param TheNumber The reference of the input number.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT NextGoodPrimeNumber(INT *TheNumber)
{


  /*vermutlich schnellere Alternative : SAtz des Eratostenes
     siehe DTV-Atlas der Mathematik von C Wrobel
     Prinzip : Man beginnt mit einem boolschen Feld, das bei allen Zahlen
     auf TRUE- gesetzt ist . Anschliessend streicht mal (beginnend von klein nach
     gross die Vielfachen von Primzahel aus (False)*/


  int nmbofoprimzahls,i;
  int *primzahlarray;
  int primzahl,index;
  int index_used;
  int wurzel,wurzel2;
  int newindex;

  int rv;

  wurzel2 = (int)floor(2*sqrt(*TheNumber));


  if ((primzahlarray = (int *)GetTmpMem(theHeap,wurzel2*sizeof(int),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"NextGoodPrimeNumber","  ERROR: No memory for primzahlarray");
    return(1);
  }


  primzahlarray[0] = 2;
  primzahlarray[1] = 3;
  newindex = 2;
  for (i = 2; i<wurzel2; i++)
  {
    primzahlarray[i] = 0;
  }
  /*	for (i = 0; i<1000; i++)
          {
                  usedprimenumbers[i] = 0;
          }
   */
  nmbofoprimzahls = 1;

  index_used = 0;
  rv = 1;
  /*for (i = 4; i<TheNumber; i++)*/
  i=4;
  while(rv == 1)
  {

    primzahl = T;             /*bis jetzt ist es noch eine Primzahl.*/
    index =0;
    wurzel = (int)floor(sqrt(i));
    while((primzahl == T)&&(primzahlarray[index] != 0)&&(primzahlarray[index] <= wurzel))
    {
      if (i % primzahlarray[index] == 0)
      {
        primzahl = F;                         /*dann war's also wohl doch keine Primzahl*/
      }
      index++;
    }


    if(primzahl == T)
    {
      primzahlarray[newindex] = i;

      if(newindex == wurzel2)
      {
        i = *TheNumber;
      }

      nmbofoprimzahls ++;

      if (i > *TheNumber)
      {
        /*Liegt eine neue Primzahl fuer ansys2ug vor ?*/
        /*neuer 500er Schritt !!*/
        if (
          (abs(primzahlarray[newindex] - 128) > 15) &&
          (abs(primzahlarray[newindex] - 256) > 15) &&
          (abs(primzahlarray[newindex] - 512) > 15) &&
          (abs(primzahlarray[newindex] - 1024) > 15) &&
          (abs(primzahlarray[newindex] - 2048) > 15) &&
          (abs(primzahlarray[newindex] - 4096) > 15) &&
          (abs(primzahlarray[newindex] - 8192) > 15) &&
          (abs(primzahlarray[newindex] - 16384) > 15) &&
          (abs(primzahlarray[newindex] - 32768) > 15) &&
          (abs(primzahlarray[newindex] - 65536) > 15) &&
          (abs(primzahlarray[newindex] - 131072) > 15) &&
          (abs(primzahlarray[newindex] - 262144) > 15) &&
          (abs(primzahlarray[newindex] - 524288) > 15) &&
          (abs(primzahlarray[newindex] - 1048576) > 15) &&
          (abs(primzahlarray[newindex] - 100) > 15) &&
          (abs(primzahlarray[newindex] - 1000) > 15) &&
          (abs(primzahlarray[newindex] - 10000) > 15) &&
          (abs(primzahlarray[newindex] - 100000) > 15) &&
          (abs(primzahlarray[newindex] - 1000000) > 15)
          )
        {
          /* *TheNumber = primzahlarray[newindex]; */
          rv = primzahlarray[newindex];
        }
      }

      newindex++;
    }

    i++;
  }

  if(rv ==1)
  {
    PrintErrorMessage('E',"NextGoodPrimeNumber","calculation of the next good prime number did not succeed !");
    return(1);
  }
  else
  {
    *TheNumber = rv;
    return(0);
  }

}





/****************************************************************************/
/**
*
* @param nmber_of_SFEs  The number of nodes, evaluated in readcadfile/cadconvert,
*                       now used for the hashfunctions
*                       In Ansys2lgmInit its value is used to set the static variable
*                       "number_of_nodes" of the file ansys2ug.c
* @param bndsgmntndids_array Das BoundarySegmentArray aus cadconvert
*                            wird hier an die SortierroutineSortBndSegArray(..)= weitergegeben.
*                            Hat Quadrupeleintraege: die 3 IDS der Ecken des Oberflaechen-
*                            dreieck plus zusaetzlich die ID des 4-ten innener Knotens des
*                            zugehoerigen Tetraeders.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT Ansys2lgmInit()
{
  INT rv,lf;


  /*Root Pointer auf Surfaces und Root Pointer auf Subdomains mit NULL initialisieren */
  EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer) = NULL;       /*Root Pointer auf Surfaces*/
  EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer) = NULL;       /*Root Pointer auf Subdomains*/

  SFE_p = EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer) * 2;
  LI_p = EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer) * 3;
  /*das waeren schlechte Zahlen :SFE_p = 1024; fuehrten zu zahlreichen Kollisionen !!! in LI-HsTab
     LI_p = 1024;*/

  if ((rv = NextGoodPrimeNumber(&SFE_p)) == 1)
  {
    PrintErrorMessage('E',"Ansys2lgmInit","got ERROR from function NextGoodPrimeNumber");
    return(1);
  }



  if ((rv = NextGoodPrimeNumber(&LI_p)) == 1)
  {
    PrintErrorMessage('E',"Ansys2lgmInit","got ERROR from function NextGoodPrimeNumber");
    return(1);
  }




  if ((rv = SortBndSegArray()) == 1)
  {
    PrintErrorMessage('E',"Ansys2lgmInit","got ERROR Response from function SortBndSegArray");
    return(1);
  }

  /*SFE_HashTable initialisieren!!! */
  /* OLD SFE_HashTable =	malloc(SFE_p * sizeof(SFE_KNOTEN_TYP));*/
  if ((EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer) = (SFE_KNOTEN_TYP **)GetTmpMem(theHeap,SFE_p * sizeof(SFE_KNOTEN_TYP*),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"Ansys2lgmInit","  ERROR: No memory for SFE-Hashtable");
    return(1);
  }
  for(lf=0; lf<SFE_p; lf++)
  {
    /*SFE_HashTable[lf] = NULL;*/
    (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lf] = NULL;
  }

  /*LI_HashTable initialisiern!!! */
  /*LI_HashTable =	malloc(LI_p * sizeof(LI_KNOTEN_TYP));*/
  if ((EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer) = (LI_KNOTEN_TYP **)GetTmpMem(theHeap,LI_p * sizeof(LI_KNOTEN_TYP*),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"Ansys2lgmInit","  ERROR: No memory for LI-Hashtable");
    return(1);
  }
  for(lf=0; lf<LI_p; lf++)
  {
    /*LI_HashTable[lf] = NULL;*/
    (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[lf] = NULL;
  }

  /*Initialisations of DomainInformation*/
  NMB_OF_SBDMS(DomainInfo_Pointer) = 0;
  NMB_OF_SFCES(DomainInfo_Pointer) = 0;
  NMB_OF_PLINES(DomainInfo_Pointer) = 0;
  NMB_OF_POINTS(DomainInfo_Pointer) = EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer);
  /* = Anzahl BoundaryPoints aus cadconvert */


  return(0);       /*alles in Ordnung*/
}



/****************************************************************************/
/**
* holt Speicher f�r neuen SFE-Eintrag und f�llt diesen mit Hilfe der Parameter.
*
* @param ii, jj, kk  NodeIDs of the  new SFE-Knoten.
* @param ss  der Identifiervalue des neuen SFE-Knotens.
* @param id_4  die Id des vierten Knotens des neuen SFEs.
*
* @return normally returns pointer to the new SFE_KNOTEN_TYP / Triangle, NULL if no memory is available.
   */
/****************************************************************************/
SFE_KNOTEN_TYP *GetMemAndFillNewSFE(INT ii, INT jj, INT kk, INT id_4, DOUBLE ss)
{
  SFE_KNOTEN_TYP *newsfemem;

  if ((newsfemem = (SFE_KNOTEN_TYP *)GetTmpMem(theHeap,sizeof(SFE_KNOTEN_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemAndFillNewSFE","  ERROR: No memory for a SFE_Hashtab_Entry, see ansys2lgm.c");
    return(NULL);
  }
  SFE_NDID1(newsfemem) = ii;
  SFE_NDID2(newsfemem) = jj;
  SFE_NDID3(newsfemem) = kk;
  SFE_NEXT(newsfemem)  = NULL;
  SFE_NGHB1(newsfemem) = NULL;
  SFE_NGHB2(newsfemem) = NULL;
  SFE_NGHB3(newsfemem) = NULL;
  SFE_IDF_0(newsfemem) = ss;
  SFE_IDF_1(newsfemem) = SEC_SFC_NAME_DEFAULT_VAL;
  SFE_4ND_0(newsfemem) = id_4;
  SFE_4ND_1(newsfemem) = SFE_KNID_4_DEFAULT_VAL;
  SFE_ORIENTATION_FLAG(newsfemem) = F;       /*benoetigt fuer Orientierung der Dreiecke / SFEs
                                                                T bedeutet Orientierung der SFE/DreiecksIDS
                                                                wurde ueberprueft.
                                                                siehe Ansys2lgmCreateTriaOrientations
                                                                und Subfunctions*/

  return(newsfemem);
}



/****************************************************************************/
/**
* Diese Funktion �berpr�ft, ob zwei aufsteigend nach Gr��e sortierte ID-Tripel
* identisch (return 1) sind oder nicht (return 0).
*
* @param iii, jjj, kkk  NodeIDs of the possible new SFE.
* @param sv_ids  pointer to the nodeid[3] of already existing and possibly same SFE.
*
* @return value 0:diferent. 1:IDs are the same.
   */
/****************************************************************************/
INT SameSFE(INT iii, INT jjj, INT kkk, INT *sv_ids)
{
  if ( ( iii == sv_ids[0] ) && ( jjj == sv_ids[1] ) && ( kkk == sv_ids[2] ) )
  {
    return(1);
  }
  else
  {
    return(0);
  }
}


/****************************************************************************/
/**
* Diese Funktion tr�gt ein SFE in die Hashtabelle ordnungsgem�� ein.
* Als Input hat die Funktion die 3 IDs der Ecken ijk des SFEs sowie den Identifier s.
* Ferner wird auch der vierte Knoten des Tetraeders, zu dem das SFE gehoert
* mitabgespeichert. Seine Info benoertig man spaeter um festzustelle was rechts und
* was links von einer Surface liegt ...
* Hash_SFE sucht in der SFE_HashTable, ob es schon mal einen Eintrag mit den
* selbern Ecken ijk gibt.
* Wenn ja ist zu diesem SFE-Eintrag ein Zweiter gefunden
* und der zweite Identifiers des SKE-Knotens wird von -1.0 auf s gesetzt.
* Wenn nein, dann wird Speicher fuer einen neuen SFE-Knoten angefordert und dierser
* gef�llt.
*
* @param i, j, k  NodeIDs of the SFE
* @param id4  Identifier of the fourth  "quasi inner" node
* @param s  Identifier of the SFE
*
* @return pointer to the new/actual SFE/triangle if ok, NULL if error occured.
   */
/****************************************************************************/
SFE_KNOTEN_TYP *Hash_SFE(INT i, INT j, INT k, INT id4, DOUBLE s)
{
  INT hw, schonvorhanden;
  SFE_KNOTEN_TYP *hp,*merke, *thenewmem;

  hw = The_SFE_hashfunction(i,j,k);
  hp = (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[hw];

  if (hp == NULL)       /*wenn es an dieser Stelle der Hashtabelle noch gar keinen Eintrag gibt*/
  {
    if ((thenewmem = GetMemAndFillNewSFE(i,j,k,id4,s)) == NULL)
    {
      PrintErrorMessage('E',"Hash_SFE","did receive nilpointer from GetMemAndFillNewSFE");
      return (NULL);
    }
    (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[hw] = thenewmem;
    return(thenewmem);
  }
  else
  {
    /*suchen, ob es schon einmal einen Hashtabeintrag mit den IDs i,j und k gegeben hat*/
    schonvorhanden = 0;
    while ((hp != NULL) && (schonvorhanden == 0))
    {
      schonvorhanden = SameSFE(i,j,k,SFE_NDID_ARRAY(hp));
      merke = hp;
      hp = SFE_NEXT(hp);
    }

    if (schonvorhanden == 0)             /*es gibt also noch keinen SFE-Eintrag der Art "i_j_k" !*/
    {
      if ((thenewmem = GetMemAndFillNewSFE(i,j,k,id4,s)) == NULL)
      {
        PrintErrorMessage('E',"Hash_SFE","did receive nilpointer from GetMemAndFillNewSFE");
        return (NULL);
      }
      SFE_NEXT(merke) = thenewmem;                   /*eingef�gt wird an das Ende der SFE-HW-Liste*/
      return(thenewmem);
    }
    else             /* es gibt also schon einen SFE-Eintrag der Art "i_j_k" !*/
    {
      if( (SFE_IDF_1(merke) == SEC_SFC_NAME_DEFAULT_VAL) &&
          (SFE_4ND_1(merke) == SFE_KNID_4_DEFAULT_VAL)       )
      {
        /*dieses SFE erhaelt somit seinen zweiten Identifier*/
        /*Die beiden Identifier muessen folglich der Groesse nach aufsteigend*/
        /*in den SFE-Knoten eingetragen werden.*/
        if(s > SFE_IDF_0(merke))
        {
          SFE_IDF_1(merke) = s;                              /*den Neuen an die zweite Stelle eintragen*/
          SFE_4ND_1(merke) = id4;                              /*die 4te KnotenID an die zweite Stelle eintragen*/
        }
        else
        {
          /*fuer Identifiers und 4teKnotenIDs ist Folgendes zu tun :*/
          /*den Ersten an die zweite Stelle eintragen
             und den Neuen an die erste Stelle eintragen*/
          SFE_IDF_1(merke) = SFE_IDF_0(merke);
          SFE_4ND_1(merke) = SFE_4ND_0(merke);
          SFE_IDF_0(merke) = s;
          SFE_4ND_0(merke) = id4;
        }
        return(merke);
      }
      else
      {
        PrintErrorMessage('E',"Hash_SFE","could not insert SFE_ijk for the second time because \nsecond value of IDFIis no more SEC_SFC_NAME_DEFAULT_VAL !!! or\n second value of 4ID is no more SFE_KNID_4_DEFAULT_VAL");
        return (NULL);
      }
    }
  }
}


/****************************************************************************/
/**
* erzeugt und f�llt eine neue Struktur vom IDF_TYP.
* Das ist eines der Listenelemente der Line-internen identifiers Liste.
*
* @param k  dritte NodeID des SFEs, zu dem die Line in diesem Fall geh�rt.
* @param s  Identifier of the SFE, zu dem die Line in diesem Fall geh�rt.
* @param act_tria  Zeiger auf zugeh�riges SFE_KNOTEN_TYP-Element.
*
* @return pointer to the new Identifier, NULL if error occured.
   */
/****************************************************************************/
IDF_TYP *GetMemandFillNewIDF(INT k, DOUBLE s,SFE_KNOTEN_TYP *act_tria)
{
  IDF_TYP *idf_new;

  if ((idf_new = (IDF_TYP *)GetTmpMem(theHeap,sizeof(IDF_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemandFillNewIDF","  ERROR: No memory for a IDF_TYP_Entry, see ansys2lgm.c");
    return(NULL);
  }

  IDF_VAL(idf_new) = s;
  IDF_SFE_TRIA(idf_new) = act_tria;
  IDF_ID3(idf_new) = k;
  LI_NEXT(idf_new) =NULL;

  return(idf_new);
}


/****************************************************************************/
/**
* f�gt den neuen Identifier idf_new in die Identifierliste der Line merke ein.
*
* @param idf_new dritte NodeID des SFEs, zu dem die Line in diesem Fall geh�rt.
* @param merke Identifier of the SFE, zu dem die Line in diesem Fall geh�rt.
* @param s SurfaceIdentifier(SFE-Wert aus AnsysDatei) of the SFE, zu dem die Line in diesem Fall geh�rt.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT InsertNewIdfIIntoIdfslist(IDF_TYP *idf_new, LI_KNOTEN_TYP *merke, DOUBLE s)
{
  IDF_TYP *idf, *pred_idf;
  INT Einfuegestelle;

  idf = LI_IDFS(merke);
  pred_idf = NULL;
  Einfuegestelle = 0;
  while((idf != NULL) && (Einfuegestelle == 0))
  {
    if(s <= IDF_VAL(idf))
    {
      Einfuegestelle = 1;
      if (pred_idf == NULL)                   /*d.h. ganz vorne mu� eingef�gt werden*/
      {
        LI_IDFS(merke) = idf_new;
        LI_NEXT(idf_new) =idf;
      }
      else
      {
        LI_NEXT(pred_idf) =idf_new;
        LI_NEXT(idf_new) =idf;
      }
    }
    pred_idf = idf;
    idf = IDF_NXT(idf);
  }
  if(Einfuegestelle == 0)       /*noch immer*/
  {
    /*d.h. ganz hinten mu� eingef�gt werden*/
    LI_NEXT(pred_idf) =idf_new;
    /*LI_NEXT(idf_new) =NULL; */ /* ist ja schon in GetMemandFillNewIDF() passiert. */
  }
  return(0);
}


/****************************************************************************/
/**
* erzeugt und f�llt eine neue Struktur vom LI_KNOTEN_TYP.
* Das ist eines der LIKNoten der LI-Hashtabelle
* f�r die interne Identifierliste, die zu jeder Line angelegt wird,
* wird die Funktion GetMemandFillNewIDF aufgerufen.
*
* @param  i  niederwertige NodeID der Line.
* @param  j  h�herwertige NodeID der Line.
* @param  k  dritte NodeID des SFEs, zu dem die Line in diesem Fall geh�rt.
* @param  s  Identifier of the SFE, zu dem die Line in diesem Fall geh�rt.
* @param  act_tria  Zeiger auf zugeh�riges SFE_KNOTEN_TYP-Element.
*
* @return pointer to the new LIne, NULL if error occured (no memory for example available).
   */
/****************************************************************************/
LI_KNOTEN_TYP *GetMemAndFillNewLI(INT i, INT j, INT k, DOUBLE s, SFE_KNOTEN_TYP *act_tria)
{
  LI_KNOTEN_TYP *line_new;

  if ((line_new = (LI_KNOTEN_TYP *)GetTmpMem(theHeap,sizeof(LI_KNOTEN_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemAndFillNewLI","  ERROR: No memory for a LI_Hashtab_Entry, see ansys2lgm.c");
    return(NULL);
  }


  LI_NDID1(line_new) = i;
  LI_NDID2(line_new) = j;
  LI_NEXT(line_new) = NULL;
  if ((LI_IDFS(line_new) = GetMemandFillNewIDF(k,s,act_tria)) == NULL)
  {
    PrintErrorMessage('E',"GetMemAndFillNewLI","did receive a nil ptr out of GetMemandFillNewIDF\n instead of a pointer to a new Identifier");
    return (NULL);
  }

  return(line_new);
}


/****************************************************************************/
/**
* Diese Funktion �berpr�ft, ob zwei aufsteigend nach Gr��e sortierte ID-Paare
* identisch (return 1) sind oder nicht (return 0).
*
* @param i, j  NodeIDs of the possible new LIne.
* @param lv_ids  pointer to the nodeid[2] of already existing and possibly same LIne.
*
* @return value   0:diiferent. 1:IDs are the same.
   */
/****************************************************************************/
INT SameLI(INT i, INT j, INT *lv_ids)
{
  if ( ( i == lv_ids[0] ) && ( j == lv_ids[1] ) )
  {
    return(1);
  }
  else
  {
    return(0);
  }
}


/****************************************************************************/
/**
* Diese Funktion tr�gt ein LI in die Hashtabelle ordnungsgem�� ein.
* Als Input hat die Funktion die 2 IDs der Ecken ij der Line sowie den Identifier s.
* Ferner steht in k die dritte Knotenid des zugeh. SFEs.
* Hash_LI sucht in der LI_HashTable, ob es schon mal einen Eintrag mit den
* selbern Ecken ij gibt.
* Wenn ja, dann ist zu diesem LI-Eintrag ein Weiterer gefunden,der dann in
* die identifiers-Liste der Gr��e nach aufsteigend einsortiert wird.
* Wenn nein, dann wird Speicher fuer einen neuen LI-Knoten angefordert und dieser
* gef�llt.
*
* @param i, j  aufsteigende NodeIDs of the new Line.
* @param k  dritte NodeID des SFEs, zu dem die Line geh�rt.
* @param s  Identifier of the SFE, zu dem die Line geh�rt:
* @param act_tria  Zeiger auf zugeh�riges SFE_KNOTEN_TYP-Element.
*
* @return pointer to the new/actual SFE/triangle if ok, NULL if error occured.
   */
/****************************************************************************/
LI_KNOTEN_TYP *Hash_LI(INT i, INT j, INT k, DOUBLE s, SFE_KNOTEN_TYP *act_tria)
{
  INT hw, schonvorhanden, rv;
  LI_KNOTEN_TYP *hp,*merke, *thenewmem;
  IDF_TYP *idf_new;

  hw = the_LI_hashfunction(i,j);
  hp = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[hw];

  if (hp == NULL)       /*wenn esa dieser Stelle der Hashtabelle noch gar keinen Eintrag gibt*/
  {
    if ((thenewmem = GetMemAndFillNewLI(i,j,k,s,act_tria)) == NULL)
    {
      PrintErrorMessage('E',"Hash_LI","did receive nilpointer from GetMemAndFillNewLI");
      return (NULL);
    }
    (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[hw] = thenewmem;
    return(thenewmem);
  }
  else
  {
    /*suchen, ob es schon einmal einen Hashtabeintrag mit den IDs i,j gegeben hat*/
    schonvorhanden = 0;
    while ((hp != NULL) && (schonvorhanden == 0))
    {
      schonvorhanden = SameLI(i,j,LI_NDID(hp));
      merke = hp;
      hp = LI_NEXT(hp);
    }

    if (schonvorhanden == 0)             /*es gibt also noch keinen LI-Eintrag der Art "i_j" !*/
    {
      if ((thenewmem = GetMemAndFillNewLI(i,j,k,s,act_tria)) == NULL)
      {
        PrintErrorMessage('E',"Hash_LI","did receive nilpointer from GetMemAndFillNewLI");
        return (NULL);
      }
      LI_NEXT(merke) = thenewmem;                   /*eingef�gt wird an das Ende der SFE-HW-Liste*/
      return(thenewmem);
    }
    else             /*es gibt also  bereits einen LI-Eintrag der Art "i_j" !*/
    {
      if((idf_new = GetMemandFillNewIDF(k,s,act_tria)) == NULL)
      {
        PrintErrorMessage('E',"Hash_LI","did receive nilpointer from GetMemandFillNewIDF");
        return (NULL);
      }

      if((rv = InsertNewIdfIIntoIdfslist(idf_new,merke,s)) == 1)
      {
        PrintErrorMessage('E',"InsertNewIdfIIntoIdfslist","did receive returnvalue = 1  ==> ERROR !");
        return (NULL);
      }

      return(merke);
    }
  }
}



/****************************************************************************/
/**
* l�uft �ber alle Oberfl�chendreiecke des ANSYS-Files und f�llt dabei
* die beiden Hashtabellen SFE-Hashtable und LI-HashTable.

* @param  bndsegids_array  Pointer to an array created (TODO) in cadconvert
*                                which includes 3erKnotentripel, die den ANSYS-SFE-Part
*                                wiederspiegeln.
* @param  bndsegidentifier_array  Pointer to an array created (TODO) in cadconvert
*                                which includes Identifiers, die den ANSYS-SFE-Part
*                                wiederspiegeln und zu dem zugeh�rigen 3er-KnotenID-Triple plus
*                                auch der vierte innere Knoten geh�ren, der spaeter fuer die
*                                Orientierung der Surfaces bzgl Innen/Aussen wichtig ist.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT Ansys2lgmCreateHashTables()
{
  INT lf;
  SFE_KNOTEN_TYP *actual_triangle_ptr;
  LI_KNOTEN_TYP *actual_line_ptr;
  DOUBLE sfcidf;
  INT id_i,id_j,id_k, id_4;

        #ifdef STATISTICAL_INFORMATIONS
  SFE_KNOTEN_TYP *sfeptr;
  INT lff,zaehler,zaehlerL;
  INT LI_HT_Index;
  LI_KNOTEN_TYP *li_ptr;
        #endif


  /*only for debugging*/
  SFE_KNOTEN_TYP **hilfusSFE;
  LI_KNOTEN_TYP **hilfusLI;

        #ifdef STATISTICAL_INFORMATIONS
  /*Initialisierungen*/
  ST_INF_Anz_SFE = EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer);
  ST_INF_Anz_LI = ST_INF_Anz_SFE *3;
  ST_INF_Anz_ds = 0;
  /*ST_INF_Anz_dl = 0;*/
  ST_INF_2er = 0;
  ST_INF_3er = 0;
  ST_INF_4er = 0;
  ST_INF_5er = 0;
  ST_INF_maxK = 1;
  ST_INF_maxKL = 1;
  ST_INF_Klsstellen = 0;
  ST_INF_Kollis = 0;
  ST_INF_2erL = 0;
  ST_INF_3erL = 0;
  ST_INF_4erL = 0;
  ST_INF_5erL = 0;
  ST_INF_KlsstellenL = 0;
  ST_INF_KollisL = 0;

  ST_INF_m = SFE_p;
  ST_INF_gef = 0;
  ST_INF_mL = LI_p;
  ST_INF_gefL = 0;

  ST_INF_Anz_LI_real = 0;

        #endif





  for(lf=0; lf<(EXCHNG_TYP1_NMB_OF_SFES(ExchangeVar_1_Pointer)); lf++)       /*laufe �ber alle SFE's*/
  {
    /* init help variable for nodeids and surfaceidentifiers : */
    sfcidf = CAD_SFE_SFE_IDF((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    id_i = CAD_SFE_ND_I((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    id_j = CAD_SFE_ND_J((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    id_k = CAD_SFE_ND_K((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));
    id_4 = CAD_SFE_ND_4((&((EXCHNG_TYP1_SFE_ARRAY(ExchangeVar_1_Pointer))[lf])));

    /* Hashtabelleneintrag in die SFE-HashTabelle : */
    /*GOON HERE :*/
    if ((actual_triangle_ptr = Hash_SFE(id_i,id_j,id_k,id_4,sfcidf)) == NULL)
    {
      PrintErrorMessage('E',"Ansys2lgmCreateHashTables","got nil-ptr out of hashSFE()");
      return(1);
    }

    /* Hashtabelleneintr�ge in die LI-HashTabelle : */
    /* (1.) Hashtabelleneintrag der "Line vom ersten zum zweiten node (i-->j)"  in die LI-HashTabelle : */
    if ((actual_line_ptr = Hash_LI(id_i,id_j,id_k,sfcidf,actual_triangle_ptr)) == NULL)
    {
      PrintErrorMessage('E',"Ansys2lgmCreateHashTables","got NULL Response from fct hash_LI (i-->j) ");
      return(1);
    }
    /* (2.) Hashtabelleneintrag der "Line vom ersten zum dritten node (i-->k)" in die LI-HashTabelle : */
    if ((actual_line_ptr = Hash_LI(id_i,id_k,id_j,sfcidf,actual_triangle_ptr)) == NULL)
    {
      PrintErrorMessage('E',"Ansys2lgmCreateHashTables","got NULL Response from fct hash_LI (i-->k)");
      return(1);
    }
    /* (3.) Hashtabelleneintrag der "Line vom zweiten zum dritten node (j-->k)" in die LI-HashTabelle : */
    if ((actual_line_ptr = Hash_LI(id_j,id_k,id_i,sfcidf,actual_triangle_ptr)) == NULL)
    {
      PrintErrorMessage('E',"Ansys2lgmCreateHashTables","got NULL Response from fct hash_LI (j-->k)");
      return(1);
    }
  }

  /*for debugging*/
  if ((hilfusSFE = (SFE_KNOTEN_TYP **)GetTmpMem(theHeap,SFE_p * sizeof(SFE_KNOTEN_TYP*),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"Ansys2lgmCreateHashTables","  ERROR: No memory for hilfusSFE");
    return(1);
  }
  if ((hilfusLI = (LI_KNOTEN_TYP **)GetTmpMem(theHeap,LI_p * sizeof(LI_KNOTEN_TYP*),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"Ansys2lgmCreateHashTables","  ERROR: No memory for hilfusLI");
    return(1);
  }
  for(lf=0; lf < SFE_p; lf++)
  {
    hilfusSFE[lf] = (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lf];
  }
  for(lf=0; lf < LI_p; lf++)
  {
    hilfusLI[lf] = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[lf];
  }


        #ifdef STATISTICAL_INFORMATIONS
  /*Berechnungen*/

  /*laufe ueber die gesamte SFE-Hashtabelle*/
  /*laeuft ueber die gesamte SFE-Hashtabelle*/
  for(lff=0; lff < SFE_p; lff++)
  {
    /*wenn an dieser Stelle ueberhaupt ein Eintrag erfolgte*/
    if ((EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lff] != NULL)
    {
      /*an dieser Stelle der SFE-Hashtabelle steht etwas*/
      ST_INF_gef ++;

      sfeptr = (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lff];

      /*wenn es hier mehr als diesen einen Eintrag gibt . . . */
      if(SFE_NEXT(sfeptr) != NULL)
      {
        /*dann liegt hier eine Kollisionsstelle vor.*/
        ST_INF_Klsstellen++;
        /*laufe ueber alle Eintraege an dieser Stelle der Hashtabelle*/
        zaehler = 0;
        while(sfeptr != NULL)
        {
          zaehler++;
          /*wenn hier eine Doppelsurface vorliegt*/
          if(SFE_IDF_1(sfeptr) != SEC_SFC_NAME_DEFAULT_VAL)
          {
            ST_INF_Anz_ds++;
          }


          /*weiter gehts mit dem naechsten SFE Eintrag an dieser Stelle der SFE-Hashtabelle*/
          sfeptr = SFE_NEXT(sfeptr);
        }

        if(zaehler > ST_INF_maxK)
        {
          ST_INF_maxK = zaehler;
        }

        /*Update der Kollisionsanzahl:*/
        ST_INF_Kollis = ST_INF_Kollis + zaehler -1;

        switch(zaehler)
        {
        case 0 : PrintErrorMessage('E',"Ansys2lgmCreateHashTables","STATISTICAL_INFORMATIONS: 0 not possible");
          return (1);
          break;
        case 1 : PrintErrorMessage('E',"Ansys2lgmCreateHashTables","STATISTICAL_INFORMATIONS: 1 not possible");
          return (1);
          break;
        case 2 : ST_INF_2er++;
          break;
        case 3 : ST_INF_3er++;
          break;
        case 4 : ST_INF_4er++;
          break;
        default :                                       /*d.h. >= 5*/
          ST_INF_5er++;
        }

      }
      else
      {
        /*wenn hier eine Doppelsurface vorliegt*/
        if(SFE_IDF_1(sfeptr) != SEC_SFC_NAME_DEFAULT_VAL)
        {
          ST_INF_Anz_ds++;
        }
      }
    }
  }

  /*weitere statistische Berechnungen*/
  ST_INF_Anz_SFE_real = ST_INF_Anz_SFE - ST_INF_Anz_ds;
  ST_INF_fg = ((double)ST_INF_gef) / ((double)SFE_p);
  ST_INF_Kollishf = ((double)ST_INF_Kollis) / ((double)ST_INF_Anz_SFE_real);
  ST_INF_2er_P = ((double)ST_INF_2er) / ((double)ST_INF_Klsstellen);
  ST_INF_3er_P = ((double)ST_INF_3er) / ((double)ST_INF_Klsstellen);
  ST_INF_4er_P = ((double)ST_INF_4er) / ((double)ST_INF_Klsstellen);
  ST_INF_5er_P = ((double)ST_INF_5er) / ((double)ST_INF_Klsstellen);




  /*here begins the LI-Part : . . .*/
  /*laufe ueber die gesamte LI-Hashtabelle*/
  for(LI_HT_Index = 0; LI_HT_Index < LI_p; LI_HT_Index++)
  {
    /*Existiert hier ueberhaupt ein Eintrag ?*/
    li_ptr = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[LI_HT_Index];
    if(li_ptr != NULL)
    {
      ST_INF_gefL ++;
      ST_INF_Anz_LI_real++;                           /*mind eine Line*/
      /*wenn es hier mehr als diesen einen Eintrag gibt . . . */
      if(LI_NEXT(li_ptr) != NULL)
      {
        /*dann liegt hier eine Kollisionsstelle vor.*/
        ST_INF_KlsstellenL++;
        /*laufe ueber alle Eintraege an dieser Stelle der Hashtabelle*/
        /*laufe ueber die lineare Liste, die mit diesem Hashpointer li_ptr beginnt*/
        zaehlerL = 0;
        while(li_ptr != NULL)
        {
          zaehlerL++;

          li_ptr = LI_NEXT(li_ptr);
        }
        ST_INF_Anz_LI_real = ST_INF_Anz_LI_real + zaehlerL - 1;
        if(zaehlerL > ST_INF_maxKL)
        {
          ST_INF_maxKL = zaehlerL;
        }

        /*Update der Kollisionsanzahl:*/
        ST_INF_KollisL = ST_INF_KollisL + zaehlerL -1;

        switch(zaehlerL)
        {
        case 0 : PrintErrorMessage('E',"Ansys2lgmCreateHashTables","STATISTICAL_INFORMATIONS: 0 LI not possible");
          return (1);
          break;
        case 1 : PrintErrorMessage('E',"Ansys2lgmCreateHashTables","STATISTICAL_INFORMATIONS: 1 LI not possible");
          return (1);
          break;
        case 2 : ST_INF_2erL++;
          break;
        case 3 : ST_INF_3erL++;
          break;
        case 4 : ST_INF_4erL++;
          break;
        default :                                       /*d.h. >= 5*/
          ST_INF_5erL++;
        }

      }
    }
  }

  /*weitere statistische Berechnungen*/
  ST_INF_fgL = ((double)ST_INF_gefL) / ((double)LI_p);
  ST_INF_KhL = ((double)ST_INF_KollisL) / ((double)ST_INF_Anz_LI_real);
  ST_INF_2erL_P = ((double)ST_INF_2erL) / ((double)ST_INF_KlsstellenL);
  ST_INF_3erL_P = ((double)ST_INF_3erL) / ((double)ST_INF_KlsstellenL);
  ST_INF_4erL_P = ((double)ST_INF_4erL) / ((double)ST_INF_KlsstellenL);
  ST_INF_5erL_P = ((double)ST_INF_5erL) / ((double)ST_INF_KlsstellenL);

  /*jetzt moeglich:*/
  ST_INF_LI_real_durch_SFE_real = ((double)ST_INF_Anz_LI_real) / ((double)ST_INF_Anz_SFE_real);



        #endif




        #ifdef STATISTICAL_INFORMATIONS
  UserWrite("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  UserWrite("Statistical Informations about the two hashtables are following :   \n");
  UserWrite("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
  UserWriteF("Anz_SFE                 %d \n",(int)ST_INF_Anz_SFE);
  UserWriteF("Anz_ds                  %d \n",(int)ST_INF_Anz_ds);
  UserWriteF("Anz_SFE_real            %d \n",(int)ST_INF_Anz_SFE_real);
  UserWriteF("m                       %d \n",(int)ST_INF_m);
  ST_INF_fg = (floor((ST_INF_fg * 1000) + 0.5)/10);
  UserWriteF("fg                      %f \n",(double)ST_INF_fg);
  ST_INF_Kollishf = (floor((ST_INF_Kollishf * 1000) + 0.5)/10);
  UserWriteF("Kollishf                %f \n",(double)ST_INF_Kollishf);
  UserWriteF("2er                     %d \n",(int)ST_INF_2er);
  ST_INF_2er_P = (floor((ST_INF_2er_P * 1000) + 0.5)/10);
  UserWriteF("2erP                     %f \n",(double)ST_INF_2er_P);
  UserWriteF("3er                     %d \n",(int)ST_INF_3er);
  ST_INF_3er_P = (floor((ST_INF_3er_P * 1000) + 0.5)/10);
  UserWriteF("3erP                     %f \n",(double)ST_INF_3er_P);
  UserWriteF("4er                     %d \n",(int)ST_INF_4er);
  ST_INF_4er_P = (floor((ST_INF_4er_P * 1000) + 0.5)/10);
  UserWriteF("4erP                     %f \n",(double)ST_INF_4er_P);
  UserWriteF("5er                     %d \n",(int)ST_INF_5er);
  ST_INF_5er_P = (floor((ST_INF_5er_P * 1000) + 0.5)/10);
  UserWriteF("5erP                     %f \n",(double)ST_INF_5er_P);
  UserWriteF("maxK                    %d \n",(int)ST_INF_maxK);
  UserWriteF("Klsstellen              %d \n",(int)ST_INF_Klsstellen);
  UserWriteF("Anz_LI                  %d \n",(int)ST_INF_Anz_LI);
  /*UserWriteF("Anz_dl                  %d \n",(int)ST_INF_Anz_dl);*/
  UserWriteF("LI_real                 %d \n",(int)ST_INF_Anz_LI_real);
  UserWriteF("mL                       %d \n",(int)ST_INF_mL);
  UserWriteF("LI_real_durch_SFE_real  %f \n",(double)ST_INF_LI_real_durch_SFE_real);
  ST_INF_fgL = (floor((ST_INF_fgL * 1000) + 0.5)/10);
  UserWriteF("fgL                     %f \n",(double)ST_INF_fgL);
  ST_INF_KhL = (floor((ST_INF_KhL * 1000) + 0.5)/10);
  UserWriteF("KhL                     %f \n",(double)ST_INF_KhL);
  UserWriteF("2erL                    %d \n",(int)ST_INF_2erL);
  ST_INF_2erL_P = (floor((ST_INF_2erL_P * 1000) + 0.5)/10);
  UserWriteF("2erLP                    %f \n",(double)ST_INF_2erL_P);
  UserWriteF("3erL                    %d \n",(int)ST_INF_3erL);
  ST_INF_3erL_P = (floor((ST_INF_3erL_P * 1000) + 0.5)/10);
  UserWriteF("3erLP                    %f \n",(double)ST_INF_3erL_P);
  UserWriteF("4erL                    %d \n",(int)ST_INF_4erL);
  ST_INF_4erL_P = (floor((ST_INF_4erL_P * 1000) + 0.5)/10);
  UserWriteF("4erLP                    %f \n",(double)ST_INF_4erL_P);
  UserWriteF("5erL                    %d \n",(int)ST_INF_5erL);
  ST_INF_5erL_P = (floor((ST_INF_5erL_P * 1000) + 0.5)/10);
  UserWriteF("5erLP                    %f \n",(double)ST_INF_5erL_P);
  UserWriteF("maxKL                   %d \n",(int)ST_INF_maxKL);
  UserWriteF("KlsstellenL             %d \n",(int)ST_INF_KlsstellenL);
        #endif




  return(0);       /*alles glatt durchgelaufen*/
}



/****************************************************************************/
/**
* diese FUnktion allokiert den Speicher fuer eine neue Surface und fuehrt die anfaenglichen
* Initialisierungen durch.
*
* @param surfacename Pointer auf DOUBLE-Array, das die eine oder die beiden
*                    Surfacezahlen beinhaltet
*
* @return pointer to new respectively already existing Surface, NULL if error occured.
   */
/****************************************************************************/
SF_TYP *GetMemandFillNewSF(DOUBLE *surfacename)
{
  SF_TYP *Surface;


  if((Surface = (SF_TYP *)GetTmpMem(theHeap,sizeof(SF_TYP),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"GetMemandFillNewSF","got  no memory  for a new Surface !?!");
    return(NULL);
  }
  SF_NEXT(Surface) = NULL;
  SF_TRIAS(Surface) = NULL;
  SF_NAME1(Surface) = surfacename[0];
  SF_NAME2(Surface) = surfacename[1];
  SF_NMB_OF_TRIAS(Surface) = 0;
  SF_NMB_OF_POINTS(Surface) = 0;
  SF_RIGHT_SBD(Surface)  = SF_RL_SBD_NOT_SET_YET;
  SF_LEFT_SBD(Surface)   = SF_RL_SBD_NOT_SET_YET;
  SF_NMB_OF_POLYLINES(Surface) = 0;
  SF_POLYLINES(Surface) = NULL;
  /* SURFACEDETECTOR . . .  */
  SF_NMB_OF_POLYLI_ZYK(Surface) = 0;
  SF_NMB_OF_REALSFCS(Surface) =   0;
  SF_POLYLI_ZYK(Surface) = NULL;
  SF_REALSFCS(Surface) = NULL;
  /* . . . SURFACEDETECTOR */

  /*Update of statistical Domain Info*/
  NMB_OF_SFCES(DomainInfo_Pointer) = NMB_OF_SFCES(DomainInfo_Pointer) + 1;


  return(Surface);
}


/****************************************************************************/
/**
* Funktion laeuft ueber die Surfaces und sucht, ob es bereits eine Surface gibt,
* die die beiden DOUBLE-Zahlen als Merkmale hat besitzt.
*
* @param surfacename  Pointer auf DOUBLE-Array, das die eine oder die beiden
*                     Surfacezahlen beinhaltet
*
* @return pointer to new respectively already existing Surface, NULL if error occured.
   */
/****************************************************************************/
SF_TYP *CreateOrFetchSurface(DOUBLE *surfacename)
{
  SF_TYP *merk_sfc,*lf_sfc;
  int gibtsschon;

  /*wenn es noch keine einzige Surface in der SF-Liste gibt ...*/
  if (EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer) == NULL)
  {
    if((EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer) = GetMemandFillNewSF(surfacename)) == NULL)
    {
      PrintErrorMessage('E',"CreateSF","got nil-ptr out of GetMemandFillNewSF() no memory ?!?");
      return(NULL);
    }
    return(EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer));
  }

  /* wenn es jedoch bereits eine oder mehrere Surfaces in der SF-Liste gibt...*/
  else
  {
    lf_sfc = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
    gibtsschon = F;

    /*laufe ueber die Surfaces und suche ob es die Surface mglweise schon gibt ...*/
    while((lf_sfc != NULL) && (gibtsschon == F))
    {
      if((surfacename[0] == SF_NAME1(lf_sfc)) && (surfacename[1] == SF_NAME2(lf_sfc)))
      {
        gibtsschon = T;
      }
      merk_sfc = lf_sfc;
      lf_sfc = SF_NEXT(lf_sfc);
    }

    /* wenn es diese Surface wirklich noch nicht gibt ...*/
    if(gibtsschon == F)
    {
      if((SF_NEXT(merk_sfc) = GetMemandFillNewSF(surfacename)) == NULL)
      {
        PrintErrorMessage('E',"CreateSF","got nil-ptr out of GetMemandFillNewSF() no memory ?!?");
        return(NULL);
      }
      return(SF_NEXT(merk_sfc));
    }
    else
    {
      return(merk_sfc);
    }
  }       /*von else ... es gibt bereits eine oder mehrere SDs*/
}



/****************************************************************************/
/**
* holt tatsaechlich den Speicher fuer einen neuen SFC-Entry ,
* der weiter oben bei der zugehoerigen SFC-Entry-Liste einer Subdomain
* angesiedelt ist.
*
* @param theSurface  Gibt Zeiger auf theSurface an.
*
* @return pointer to new SFC_TYP-Entry bzw. den schon bestehendne SFC_TYP-Entry, NULL if error occured.
   */
/****************************************************************************/
SFC_TYP *GetMemandFillNewSFC(SF_TYP *theSurface)
{
  SFC_TYP *SfceEntry;

  if ((SfceEntry = (SFC_TYP *)GetTmpMem(theHeap,sizeof(SFC_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemandFillNewSFC","  got no Memory out of GetTmpMem..., see ansys2lgm.c");
    return(NULL);
  }


  SFC_NEXT(SfceEntry) = NULL;
  SFC_SURF(SfceEntry) = theSurface;


  return(SfceEntry);
}



/****************************************************************************/
/**
* Diese Funktion verknuepft eine Subdomain mit einer Surface
* bzw. untersucht ob es die Verknuepfung vielleicht sogar schon gibt.
*
* @param sbdm Pointer auf Subdomain, die auf die "neue" Surface theSurface zeigen soll.
* @param theSurface Gibt Zeige auf theSurface an.
*
* @return pointer to new SFC_TYP-Entry bzw. den schon bestehenden SFC_TYP-Entry, NULL if error occured.
   */
/****************************************************************************/
SFC_TYP *CreateAndConnectSfceEntryWithSbd(SD_TYP *sbdm, SF_TYP *theSurface)
{
  SFC_TYP *lf_sfc, *merke_sfc;
  int gibtsschon;

  /*wenn es noch gar keine Surfaceeintraege gibt zu dieser Subdomain ...*/
  if(SD_SFCS(sbdm) == NULL)
  {
    if((SD_SFCS(sbdm) = GetMemandFillNewSFC(theSurface)) == NULL)
    {
      PrintErrorMessage('E',"ConnectSdWithSfce","got no SFC-Ptr out of GetMemandFillNewSFC");
      return(NULL);
    }
    /* Subdomain sbdm hat nun eine Surface mehr ! Inkrementation !*/
    SD_NMB_OF_SFCS(sbdm) = SD_NMB_OF_SFCS(sbdm) + 1;
    return(SD_SFCS(sbdm));
  }
  /*wenn diese Subdomain aber schon Surfaceeintraege besitzt ...*/
  else
  {
    lf_sfc = SD_SFCS(sbdm);
    gibtsschon = F;
    while((lf_sfc != NULL) && (gibtsschon == F))
    {
      if( SFC_SURF(lf_sfc) == theSurface)
      {
        gibtsschon = T;
      }
      merke_sfc = lf_sfc;
      lf_sfc = SFC_NEXT(lf_sfc);
    }
    /*wenn es den SFC-Eintrag wirklich noch nicht gibt ...*/
    if(gibtsschon == F)
    {
      /*wenn die Surface eine einfache Surface ist ==> hinten in Liste einfuegen*/
      if(SF_NAME2(theSurface) == SEC_SFC_NAME_DEFAULT_VAL)
      {
        if ((SFC_NEXT(merke_sfc) = GetMemandFillNewSFC(theSurface)) == NULL )
        {
          PrintErrorMessage('E',"ConnectSdWithSfce","got no SFC-Ptr out of GetMemandFillNewSFC");
          return(NULL);
        }
        /* Subdomain sbdm hat nun eine Surface mehr ! Inkrementation !*/
        SD_NMB_OF_SFCS(sbdm) = SD_NMB_OF_SFCS(sbdm) + 1;
        return(SFC_NEXT(merke_sfc));
      }
      else                   /*es handelt sich um eine doppelte Surface ==> diese muss an den Listenanfang
                                eingefuegt werden, da sonst ConnectPolylineWithSurfaces Probleme machen kann.*/
      {
        merke_sfc = SD_SFCS(sbdm);                         /*bisherigen Listenanfang merken*/
        if ((SD_SFCS(sbdm) = GetMemandFillNewSFC(theSurface)) == NULL )
        {
          PrintErrorMessage('E',"ConnectSdWithSfce","got no SFC-Ptr out of GetMemandFillNewSFC");
          return(NULL);
        }
        SFC_NEXT(SD_SFCS(sbdm)) = merke_sfc;
        /* Subdomain sbdm hat nun eine Surface mehr ! Inkrementation !*/
        SD_NMB_OF_SFCS(sbdm) = SD_NMB_OF_SFCS(sbdm) + 1;
        /*Achtung in diesem Fall muss natuerlich der Listenanfang zurueckgegeben
           werden, da doppelteSurfaces ja am Listenkopf eingefuegt werden*/
        return(SD_SFCS(sbdm));
      }

    }
    /*wenn es den SFC-Eintrag aber bereits schon gibt ...*/
    else
    {
      /*muss auch nicht neu angelegt werden*/
      return(merke_sfc);
    }
  }
}



/****************************************************************************/
/**
* holt oder erzeugt sich mit weiteren Subfunktionen  die
* aktuelle Surface und laesst sie mit den Subdomains verbinden.
*
* @param sfe_pter Pointer auf SFE-HTab-Eintrag, dessen Identifierinformation verwendet wird.
* @param sbdm0 Gibt das zugehoerige Subdomain an.
* @param sbdm1 Gibt ggf. ein zweites Subdomain an, zu dem die neue Surface auch gehoert.
*
* @return pointer to new respectively already existing Surface, NULL if error occured.
   */
/****************************************************************************/
SF_TYP *ConnectSdWithSfce(SFE_KNOTEN_TYP *sfe_ptr, SD_TYP *sbdm0, SD_TYP *sbdm1)
{
  DOUBLE surfacename[2];
  SF_TYP *theSurface;
  SFC_TYP *Sfce_Ret_Val;

  /*wenn es keine gemeinsame Surface ist ...*/
  if(sbdm1 == NULL)
  {
    surfacename[0] = SFE_IDF_0(sfe_ptr);
    surfacename[1] = SEC_SFC_NAME_DEFAULT_VAL;
  }
  /*wenn es aber eine gemeinsame Surface ist ...*/
  else
  {
    /*Achtung hier wird aufsteigend sortiert,
       damit man spaeter in CreateOrFetchSurface auch was eindeutig
       wieder finden kann*/
    if (SFE_IDF_0(sfe_ptr) < SFE_IDF_1(sfe_ptr))             /*richtig sortieren !!!*/
    {
      surfacename[0] = SFE_IDF_0(sfe_ptr);
      surfacename[1] = SFE_IDF_1(sfe_ptr);
    }
    else
    {
      surfacename[1] = SFE_IDF_0(sfe_ptr);
      surfacename[0] = SFE_IDF_1(sfe_ptr);
    }
  }

  /*only for debugging*/
  sd_global =  sbdm0;

  /*jetzt wird die Surface erzeugt bzw. geholt, falls es sie schon gibt*/
  if((theSurface = CreateOrFetchSurface(surfacename)) == NULL)
  {
    PrintErrorMessage('E',"ConnectSdWithSfce","got no surface out of CreateOrFetchSurface");
    return(NULL);
  }

  /*sbdm0 mit der neuen Surface u.U. verbinden*/
  if((Sfce_Ret_Val = CreateAndConnectSfceEntryWithSbd(sbdm0,theSurface)) == NULL)
  {
    PrintErrorMessage('E',"ConnectSdWithSfce","got no SbdSfceEntry for sbdm0 out of CreateAndConnectSfceEntryWithSbd");
    return(NULL);
  }

  /*sbdm1 mit der neuen Surface u.U. verbinden*/
  if (sbdm1 != NULL)       /*nur dann notwendig ...*/
  {
    if((Sfce_Ret_Val = CreateAndConnectSfceEntryWithSbd(sbdm1,theSurface)) == NULL)
    {
      PrintErrorMessage('E',"ConnectSdWithSfce","got no SbdSfceEntry for sbdm1 out of CreateAndConnectSfceEntryWithSbd");
      return(NULL);
    }
  }

  return(theSurface);

}



/****************************************************************************/
/**
* allocates memory for a new subdomain and initializes the components.
*
* @param SubdomName INT-Value of the Subdomain to be created.
*
* @return pointer to new respectively already existing subdomain, NULL if error occured.
   */
/****************************************************************************/
SD_TYP *GetMemandFillNewSD(INT SubdomName)
{
  SD_TYP *Subdomain;

  /* neuen Speicher allokieren */

  if ((Subdomain = (SD_TYP *)GetTmpMem(theHeap,sizeof(SD_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemandFillNewSD","  got no MEM for the new subdomain, see ansys2lgm.c");
    return(NULL);
  }

  /*Initialisierungen : */
  SD_NEXT(Subdomain) = NULL;
  SD_SFCS(Subdomain) = NULL;
  SD_NAME(Subdomain) = SubdomName;
  SD_NMB_OF_SFCS(Subdomain) = 0;

  /*Update of statistical Domain Info*/
  NMB_OF_SBDMS(DomainInfo_Pointer) = NMB_OF_SBDMS(DomainInfo_Pointer) + 1;


  return(Subdomain);
}




/****************************************************************************/
/**
* Erzeugt einen neuen Eintrag in der Subdomainliste, �berprueft aber zuvor, ob es diesen
* bereits in derr Liste gibt !
*
* @param sfe_pter Pointer auf SFE-HTab-Eintrag, dessen Identifierinformation verwendet wird.
* @param f Gibt an ob das Subdomain mit SFE-Identifier[0] oder [1] erzeugt werden soll.
*
* @return pointer to new respectively already existing subdomain, NULL if error occured.
   */
/****************************************************************************/
SD_TYP *CreateSD(SFE_KNOTEN_TYP *sfe_pter, INT f)
{
  INT gibtsschon;
  SD_TYP *lf_sd, *merk_sd, *new_sd;
  INT subdomainname;

  /*name of the possible new subdomain*/
  subdomainname = (int) floor(SFE_IDF_N(sfe_pter,f));

  /*wenn es noch keine einzige Subdomain in der SD-Liste gibt ...*/
  if (EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer) == NULL)
  {
    if((EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer) = GetMemandFillNewSD(subdomainname)) == NULL)
    {
      PrintErrorMessage('E',"CreateSD","got nil-ptr out of GetMemandFillSD() no memory ?!?");
      return(NULL);
    }
    return(EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer));
  }

  /* wenn es jedoch bereits ein oder mehrere Subdomains in der SD-Liste gibt...*/
  else
  {
    lf_sd = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
    gibtsschon = F;

    /*laufe ueber die Subdomains und suche ob es die Subdomain mglweise schon gibt ...*/
    while((lf_sd != NULL) && (gibtsschon == F))
    {
      if(subdomainname == SD_NAME(lf_sd))
      {
        gibtsschon = T;
      }
      merk_sd = lf_sd;
      lf_sd = SD_NEXT(lf_sd);
    }

    /* wenn es diese Subdomain wirklich noch nicht gibt ...*/
    if(gibtsschon == F)
    {
      if((new_sd = GetMemandFillNewSD(subdomainname)) == NULL)
      {
        PrintErrorMessage('E',"CreateSD","got nil-ptr out of GetMemandFillSD() no memory ?!?");
        return(NULL);
      }
      SD_NEXT(merk_sd) = new_sd;
      return(SD_NEXT(merk_sd));
    }
    else
    {
      return(merk_sd);
    }
  }       /*von else ... es gibt bereits eine oder mehrere SDs*/

}



/****************************************************************************/
/**
* This function connects a Surface with a triangle , d.h.
* the triangle is part of this surface.
* Thereby the number of triangles is incremented.s.
*
* @param sf Pointer to the surface.
* @param sfeptr Pointer to the triangle to be connected with the surface.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT ConnectSfcTria(SF_TYP *sf, SFE_KNOTEN_TYP *sfeptr)
{
  TRIANGLE_TYP *merketria, *newtria;

  merketria = SF_TRIAS(sf);

  if ((newtria = (TRIANGLE_TYP *)GetTmpMem(theHeap,sizeof(TRIANGLE_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"ConnectSfcTria","  got no MEM for a new triangle, see ansys2lgm.c");
    return(1);
  }

  TRIA_SFE_KN(newtria) = sfeptr;
  TRIA_NEXT(newtria) = merketria;
  SF_TRIAS(sf) = newtria;

  /*increment number of triangles of this surface:*/
  SF_NMB_OF_TRIAS(sf) = SF_NMB_OF_TRIAS(sf) + 1;

  return(0);
}




/****************************************************************************/
/**
* sucht in der LI-Hashtabelle ob es einen Nachbarn fuer die Eingangsparameter gibt.
* Wenn ja wird die Verknuepfung in beide Richtungen durchgefuehrt.
*
* @param ii niederwertige KnotenID der Kante des SFEs, zu der ein Nachbar gesucht werden soll.
* @param jj hoeherwertige KnotenID der Kante des SFEs, zu der ein Nachbar gesucht werden soll.
* @param kte  ID der Kante des SFEs, zu der ein Nachbar gesucht werden soll.
* @param sfep  SFE, zu dem ein Nachbar gesucht werden soll.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT Neighbourhood(INT ii, INT jj, INT kte, SFE_KNOTEN_TYP *sfep)
{
  INT hw;
  LI_KNOTEN_TYP *hp;
  IDF_TYP *lfptr,*merke_lfptr;
  INT gefunden;

  /*Hashwert/Schluessel berechnen*/
  hw = the_LI_hashfunction(ii,jj);
  hp = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[hw];

  /*Beim Hashwert den richtigen Listeneintrag suchen*/
  if (hp != NULL)
  {
    while( (LI_NDID1(hp) != ii) || (LI_NDID2(hp) != jj) )
    {
      hp = LI_NEXT(hp);
    }
  }

  if(hp == NULL)
  {
    PrintErrorMessage('E',"Neighbourhood","could not find the Line in the LI-Hashtable");
    return (1);
  }

  lfptr = LI_IDFS(hp);
  if(lfptr == NULL)
  {
    PrintErrorMessage('E',"Neighbourhood","the found LI-HashTable-Entry has no(!) IDF-Pointer!");
    return (1);
  }

  gefunden = F;
  /* Laufe ueber alle Eintraege der Identifierliste der Line hp */
  while( (lfptr != NULL) && (gefunden == F) )
  {

    /*wenn die Identifiers des zugehoerigen Triangles mit denen von sfep zusammenpassen*/
    if ( (SFE_IDF_0(IDF_SFE_TRIA(lfptr)) == SFE_IDF_0(sfep)) &&
         (SFE_IDF_1(IDF_SFE_TRIA(lfptr)) == SFE_IDF_1(sfep))    )
    {
      /* das zugehoerige Triangles des IdflistEintrags darf natuerlich nicht sfep selbst sein */
      if(sfep != IDF_SFE_TRIA(lfptr))
      {
        gefunden = T;
      }
    }

    merke_lfptr = lfptr;
    lfptr = LI_NEXT(lfptr);
  }

  /*wenn der Nachbar gefunden wurde ...*/
  if(gefunden == T)
  {
    /*... dann wird jetzt wirklich mit einander verknuepft/vernachbart*/
    SFE_NGHB(sfep,kte) = IDF_SFE_TRIA(merke_lfptr);

    /* und auch umgekehrt :*/
    if( IDF_ID3(merke_lfptr) < ii )
    {
      /* d.h. Nachbar von zweiter Kante "jk"*/
      SFE_NGHB2(IDF_SFE_TRIA(merke_lfptr)) = sfep;
    }
    else if( IDF_ID3(merke_lfptr) > jj )
    {
      /* d.h. Nachbar von erster Kante "ij"*/
      SFE_NGHB1(IDF_SFE_TRIA(merke_lfptr)) = sfep;
    }
    else
    {
      /* d.h. Nachbar von dritter Kante "ik"*/
      SFE_NGHB3(IDF_SFE_TRIA(merke_lfptr)) = sfep;
    }
  }
  /*else bleibt der Nachbarschftspointer von sfep auf NULL ===>  kein NAchbar*/

  return(0);
}



/****************************************************************************/
/**
* ruft fuer die 3 Kanten vom OFD "sfep" die Nachbarschaftsfindefunktion
* "Neighbourhood(...)" auf
*
* @param sfep SFE, dessen Nachbarn gesucht und eingetragen werden.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT TriaNeighbourhood(SFE_KNOTEN_TYP *sfep)
{
  INT i,j,k;
  INT rw;

  /* die 3 KnotenIDs des eingegangenen SFES "sfep"*/
  i = SFE_NDID1(sfep);
  j = SFE_NDID2(sfep);
  k = SFE_NDID3(sfep);

  /*wenn bzgl Kante "ij" noch keinen Nachbar existiert*/
  if(SFE_NGHB1(sfep) == NULL)
  {
    if( (rw = Neighbourhood(i,j,0,sfep)) == 1 )
    {
      PrintErrorMessage('E',"TriaNeighbourhood","got ERROR from calling Neighbourhood");
      return (1);
    }
  }

  /*wenn bzgl Kante "jk" noch keinen Nachbar existiert*/
  if(SFE_NGHB2(sfep) == NULL)
  {
    if( (rw = Neighbourhood(j,k,1,sfep)) == 1 )
    {
      PrintErrorMessage('E',"TriaNeighbourhood","got ERROR from calling Neighbourhood");
      return (1);
    }
  }

  /*wenn bzgl Kante "ki" noch keinen Nachbar existiert*/
  if(SFE_NGHB3(sfep) == NULL)
  {
    if( (rw = Neighbourhood(i,k,2,sfep)) == 1 )
    {
      PrintErrorMessage('E',"TriaNeighbourhood","got ERROR from calling Neighbourhood");
      return (1);
    }
  }

  return(0);
}





/****************************************************************************/
/**
* laeuft ueber die gesamte SFE-Hashtabelle und erzeugt alle notwendigen
* Subdomains und auch Surfaces sowie deren Verbindung bzw. Referenzierung.
* Ferner wird bei diesem Durchlauf auch die Nachbarschaftsbeziehungen der Dreiecke erzeugt.
* sowie ebenso die Surface-Triangle-Beziehung.
*
* @return value 0:ok. 1:error.
   */
/****************************************************************************/
INT Ansys2lgmCreateSbdsSfcsTriaRelations()
{
  SFE_KNOTEN_TYP *sfeptr;
  SF_TYP *sf;
  SD_TYP *sd0, *sd1;
  INT lff,rv;
  SF_TYP *sf_lfv;
  TRIANGLE_TYP *triangle;

  /*laeuft ueber die gesamte SFE-Hashtabelle*/
  for(lff=0; lff < SFE_p; lff++)
  {
    /*wenn an dieser Stelle ueberhaupt ein Eintrag erfolgte*/
    if ((EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lff] != NULL)
    {
      sfeptr = (EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lff];

      /*laufe ueber alle Eintraege an dieser Stelle der Hashtabelle*/
      while(sfeptr != NULL)
      {

        /* erzeuge Subdomain bzgl. des ersten Identifiers des SFEs*/
        if((sd0 = CreateSD(sfeptr,0)) == NULL)
        {
          PrintErrorMessage('E',"Ansys2lgmCreateSbdsSfcsTriaRelations"," Returnvalue from CreateSD was nil instead of subdomain pointer");
          return (1);
        }

        /* wenn das SFE gar keinen zweiten Identifier besitzt*/
        if(SFE_IDF_1(sfeptr) == SEC_SFC_NAME_DEFAULT_VAL)
        {
          /* erzeuge Surface aus dem ersten Identifier*/
          if ((sf = ConnectSdWithSfce(sfeptr,sd0,NULL)) == NULL)
          {
            PrintErrorMessage('E',"Ansys2lgmCreateSbdsSfcsTriaRelations"," Returnvalue from CreateSF was NULL instead of a surface pointer");
            return (1);
          }

        }

        /*wenn es aber einen zweiten Identifier gibt ...*/
        else
        {
          /* erzeuge  auch  noch Subdomain bzgl. des zweiten Identifiers des SFEs*/
          if((sd1 = CreateSD(sfeptr, 1)) == NULL)
          {
            PrintErrorMessage('E',"Ansys2lgmCreateSbdsSfcsTriaRelations"," Returnvalue from CreateSD was NULL instead of subdomain pointer");
            return (1);
          }

          /* erzeuge eine Surface aus beiden Identifiers*/
          if ((sf = ConnectSdWithSfce(sfeptr,sd0,sd1)) == NULL)                               /*zweiter und dritter Parameter sind belegt.*/
          {
            PrintErrorMessage('E',"Ansys2lgmCreateSbdsSfcsTriaRelations"," Returnvalue from CreateSF was NULL instead of a surface pointer");
            return (1);
          }

        }

        /* f�r beide F�lle (ein resp. zwei Identifier) mu� Sfce-Tria sowie Neighbourhoodbez. aufgebaut werden:*/

        /*Verbinde Surface sf mit SFE bzw. Triangle sfeptr*/
        if( (rv = ConnectSfcTria(sf,sfeptr)) == 1)
        {
          PrintErrorMessage('E',"Ansys2lgmCreateSbdsSfcsTriaRelations"," Returnvalue of ConnectSfcTria was 1 Could not connect surface with SFE");
          return (1);
        }

        /*weiter gehts mit dem naechsten SFE Eintrag an dieser Stelle der SFE-Hashtabelle*/
        sfeptr = SFE_NEXT(sfeptr);
      }

    }             /* von "if ((EXCHNG_TYP2_SFE_HASHTAB(ExchangeVar_2_Pointer))[lff] != NULL)" */

  }       /*von for */



  /****************************************
     Erzeugung der Dreiecksnachbarschaften
  ****************************************/
  /*laufe ueber alle Surfaces */
  sf_lfv = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  while(sf_lfv != NULL)
  {
    /*laufe ueber alle Triangles dieser Surface*/
    triangle = SF_TRIAS(sf_lfv);
    while(triangle != NULL)
    {
      sfeptr = TRIA_SFE_KN(triangle);
      /*Erzeuge Nachbarschaftsbeziehungen von dem SFE/Triangle sfeptr*/
      if( (rv = TriaNeighbourhood(sfeptr)) == 1)
      {
        PrintErrorMessage('E',"TriaNeighbourhood"," Returnvalue of TriaNeighbourhood was 1 Could not create neighbourhood");
        return (1);
      }

      /*weiter gehts mit dem naechsten SFE Eintrag an dieser Stelle der SFE-Hashtabelle*/
      triangle = TRIA_NEXT(triangle);
    }

    /*naechste Surface*/
    sf_lfv = SF_NEXT(sf_lfv);
  }      /*von while*/
  /****************************************/



  return (0);
}


/****************************************************************************/
/**
* Funtion ueberprueft anhand mehrerer DOUBLE-Identifier, die zu einer ganz
* bestimmten Line (siehe aufrufende Funktion Ansys2lgmCreatePloylines gehoert),
* ob diese Line zu einer Polyline gehoert oder nicht. Dies ist dann der Fall,
* wenn mind. zwei unterschiedliche Identifier-Zahlen vorkommen und nicht der
* Spezialfall AABB vorliegt.(Ausfuehrliche Aufschluesselung alle moeglichen
* Faellebzw. Identifierkombinationen, die eine Polyline besitzen kann siehe
* Extrablatt im ansys2UG-Konzept von Dirk).
*
* @param identifiers_pointer Pointer first Identifier of the observed Line.
*
* @return T if LineIdentifers indicate that their line is Part of a Polyline,
*         F if LineIdentifers don't indicate that their line is Part of a Polyline,
*         LI_IFS_ERROR if identifiers_pointer is a nil pointer.
   */
/****************************************************************************/
INT Check_If_Line_On_Polyline(IDF_TYP *identifiers_pointer)
{
  IDF_TYP *idfovernext;

  if(identifiers_pointer == NULL)
  {
    PrintErrorMessage('E',"Check_If_Line_On_Polyline","The input parameter identifiers_pointer is nil ==> ERROR");
    return (2);
  }

  if(IDF_NXT(identifiers_pointer) == NULL)
  {
    PrintErrorMessage('E',"Check_If_Line_On_Polyline","Input identifiers_pointer has only one identifier.");
    return (2);
  }

  idfovernext = IDF_NXT(IDF_NXT(identifiers_pointer));       /*der dritte Identifier bzw. der �bern�chste*/

  /*wenn der erste Identifer gleich dem Zweiten ist ...*/
  if(IDF_VAL(identifiers_pointer) == IDF_VAL(IDF_NXT(identifiers_pointer)))
  {
    /*wenn es nur diese zwei Identifier gibt ...*/
    if(idfovernext == NULL)
    {
      return(0);                   /*Die Line liegt nicht auf einer Polyline !!!*/
    }
    else
    {
      /*wenn es mehr als 3 Identifier gibt...*/
      if(IDF_NXT(idfovernext) != NULL)
      {
        /*wenn auch der dritte Identifier gleich dem vierten Identifier ist ...*/
        if(IDF_VAL(idfovernext) == IDF_VAL(IDF_NXT(idfovernext)))
        {
          /*wenn es ausser diesen beiden Paaren (aabb) keinen weiteren Identifier gibt ...*/
          if( IDF_NXT(IDF_NXT(idfovernext)) == NULL)
          {
            return(0);                                     /*Die Line liegt nicht auf einer Polyline !!!*/
          }
          /*
             else ... d.h. die Identifierliste beginnt zwar mit zwei Paaren
                   es folgen jedoch weitere Identifier so, dass diese Line
                   auf einer Polyline liegen muss.
           */
        }
      }
      /*
         else,...d.h. es gibt genau 3 Identifier wobei die ersten zwei gleich sind !
       */
    }
  }       /*von wenn der erste Identifer gleich dem Zweiten ist ...*/
          /*wenn die unmittelbar obige if-Schleife durchlaufen wurde ohne dass ein return(F)
             erfolgte, dann liegt die Line mit Sicherheit auf einer Polyline und muss weiter
             beruecksichtigt werden ---> Rueckgabe also an die aufrufende Funktion Ansys2lgmCreatePloylines()
             "JA!!!" Die eingegangenen Iderntifiers der zugehoerigen Line best�tigen,
             das diese Line auf einer Polyline liegen mu� deshalb nun : ...*/
  return(1);
}



/****************************************************************************/
/**
* Diese Funktion untersucht, ob es zur eingehenden Line bereits eine Polyline gibt.
* Dazu laeuft sie ueber alle bereits bestehenden Polylines und vergleicht ob eine
* dieser Polylines die gleichen charkt. IDFs besitzt. Wenn ja ist die bereits
* existierende Polyline gefunden.
*
* TODO: Moeglichkeit zur optimierung: alle Polylines koennten anhand der Charakteris-
* tischen IdfsListe in einem ausgeglichenen Binaerbaum abgespeichert werden
* um eine logarithmische Komplexitaet zu erhalten.
*
* @param identifiers_pointer Pointer to the observed Line.
*
* @return Pointer to already existing Polyline, NULL if no Polyline exists resp. was found
   */
/****************************************************************************/
PL_TYP *Exist_Polyline(LI_KNOTEN_TYP *identifiers_pointer)
{
  PL_TYP *pl_lauf;
  IDF_TYP *idf_laufINPUT, *idf_laufPLINE;

  INT identischeIDs;


  if(identifiers_pointer == NULL)
  {
    PrintErrorMessage('E',"Exist_Polyline","Input-IDFsList of the function is NULL ==> ERROR !");
    return (NULL);
  }


  pl_lauf = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);

  /*laufe ueber alle Polylines und suche ob es bereits eine gibt mit den obigen Identifiers*/
  while(pl_lauf != NULL)
  {
    idf_laufINPUT = LI_IDFS(identifiers_pointer);
    idf_laufPLINE = PL_IDFS(pl_lauf);

    /*Annahme*/
    identischeIDs = T;

    if(idf_laufPLINE == NULL)
    {
      PrintErrorMessage('E',"Exist_Polyline","IDFsList of a Polyline is NIL ==> ERROR !");
      return (NULL);
    }
    else
    {
      /*Solange keines der beiden Identifiers-Listen erreicht ist bzw.
         solange keiner der beiden IDF-Laufpointer idf_laufINPUT oder idf_laufPLINE NULL ist.
         UND
         die jeweiligen IDF_Wert/Values identisch sind !!!*/
      while( (idf_laufPLINE != NULL) && (idf_laufINPUT != NULL) && (identischeIDs == T) )
      {
        /*Sind die Values aus den beiden IdentifiersListen an dieser Stelle unterschiedlich?*/
        if (IDF_VAL(idf_laufPLINE) != IDF_VAL(idf_laufINPUT))
        {
          identischeIDs = F;
        }
        idf_laufPLINE = IDF_NXT(idf_laufPLINE);
        idf_laufINPUT = IDF_NXT(idf_laufINPUT);
      }

      /*Wenn unterschiedliche ID-Values an einer Stelle gefunden wurden*/
      if(identischeIDs == T)
      {
        /*... und wenn  beide LaufPointer das ListenEnde erreicht habe / bzw. die beiden Listen
           auch wirklich die selbe Laenge besitzten ...*/
        if((idf_laufPLINE == NULL) && (idf_laufINPUT == NULL))
        {
          /*Die gesuchte Polyline gibt es tatsaechlich schon*/
          return(pl_lauf);
        }
      }
      /*
         else... diese Polyline ist nicht die Gesuchte !
       */
    }
    /*noch nicht gefunden ===> weiter mit naechster Polyline*/
    pl_lauf = PL_NXT(pl_lauf);
  }
  /*Die gesucht Polyline gibt es noch nicht bzw. wurde nicht gefunden*/
  return(NULL);
}



/****************************************************************************/
/**
* allokiert Speicher fuer eine neue PolylineLine, initialisiert deren
* Komponenten und haengt die neue PolylineLine in die Liste der zu-
* geh�rigen Polyline und inkrementiert auch den PolylineLineZ�hler.
*
* @param  linepointer Pointer on Line, which is part of.
* @param  polylinepointer this polyline.
*
* @return Pointer on new Polyline_Line if ok, NULL if error occured.
   */
/****************************************************************************/
PL_LINE_TYP *GetMemFillAddNewPolylineLine(LI_KNOTEN_TYP *linepointer, PL_TYP *polylinepointer)
{
  PL_LINE_TYP *merkefirstpl_line, *newpl_line;

  merkefirstpl_line = PL_LINES(polylinepointer);

  if ((newpl_line = (PL_LINE_TYP *)GetTmpMem(theHeap,sizeof(PL_LINE_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemFillAddNewPolylineLine","did not receive  memory for the new polyline_Line");
    return (NULL);
  }

  PL_LINES_NXT(newpl_line) = merkefirstpl_line;
  PL_LINES_LINE(newpl_line) = linepointer;
  PL_LINES(polylinepointer) = newpl_line;
  /* Die Polyline "polylinepointer" hat nun eine Line mehr! ==> Inkrementation:"*/
  PL_NMB_OF_POINTS(polylinepointer) ++;
  return(newpl_line);
}




/****************************************************************************/
/**
* allokiert Speicher fuer eine neue Ployline und initialisiert deren Komponenten.
*
* @param linepointer Pointer on Line, which is characteristic resp. representative
*                    for the New Polyline and which will be a part of the new Polyline.
*
* @return Pointer on new Polyline if ok, NULL if error occured.
   */
/****************************************************************************/
PL_TYP *GetMemFillAddNewPolyline(LI_KNOTEN_TYP *linepointer)
{
  PL_TYP *merkeplptr, *new_pl;
  IDF_TYP *li_idf_lf;

  merkeplptr = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);

  if ((new_pl = (PL_TYP *)GetTmpMem(theHeap,sizeof(PL_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"GetMemFillAddNewPolyline","did not receive  memory for the new polyline");
    return (NULL);
  }

  PL_NMB_OF_CH_IDFS(new_pl) = 0;
  /*charakteristische Identifiers zaehlen*/
  li_idf_lf = LI_IDFS(linepointer);
  while(li_idf_lf != NULL)
  {
    PL_NMB_OF_CH_IDFS(new_pl) ++;
    li_idf_lf = IDF_NXT(li_idf_lf);
  }
  PL_IDFS(new_pl) = LI_IDFS(linepointer);
  PL_NXT(new_pl) = merkeplptr;
  PL_LINES(new_pl) = NULL;
  PL_NMB_OF_POINTS(new_pl) = 1;       /*mu� bei 1 beginnen da Anzahl der Points einer Polyline immer eins mehr ist als
                                         Anzahl der Lines einer Polyline !!! */
  if((PL_LINES(new_pl) = GetMemFillAddNewPolylineLine(linepointer,new_pl)) == NULL)
  {
    PrintErrorMessage('E',"GetMemFillAddNewPolyline","did receive nilpointer from GetMemFillAddNewPolylineLine");
    return (NULL);
  }

  /*Es gibt eine Polyline mehr !*/
  NMB_OF_PLINES(DomainInfo_Pointer) ++;

  EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer) = new_pl;

  return(new_pl);
}



/****************************************************************************/
/**
* kopiert eine lineare Liste vom IDF_TYP in eine neue Liste vom Typ IDF_SHORT_TYP
* In die neue Liste werden natuerlich nur die Identifiers gelegt entsprechend dem
* abgespeckten Typ IDF_SHORT_TYP. Achtung die Reihenfolge (d.h. sortiert von klein
* nach gross) bleibt dabei erhalten !!!
*
* @param charact_identifiers Pointer to the Identifierlist of the observed Polyline.
*
* @return Pointer to copied IdentifierList if ok, NULL if error occured resp. no copied List could be created.
   */
/****************************************************************************/
IDF_SHORT_TYP *CopyCharacteristicList2HelpList(IDF_TYP *charact_identifiers)
{
  /*only for debugging*/
  DOUBLE helpvar;
  IDF_SHORT_TYP *helplf;


  IDF_SHORT_TYP *thedoubledlist, *merke, *listenstart;

  thedoubledlist = NULL;

  if(charact_identifiers == NULL)
  {
    PrintErrorMessage('E',"CopyCharacteristicList2HelpList","The value of the InputParameter was NULL");
    return (NULL);
  }

  /*das erste Listenelement*/
  if ((thedoubledlist = (IDF_SHORT_TYP *)GetTmpMem(theHeap,sizeof(IDF_SHORT_TYP),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"CopyCharacteristicList2HelpList","no memory obtained from GetMem(sizeof(IDF_SHORT_TYP))");
    return (NULL);
  }
  IDF_SHORT_NXT(thedoubledlist) = NULL;
  IDF_SHORT_VAL(thedoubledlist) = IDF_VAL(charact_identifiers);
  charact_identifiers = IDF_NXT(charact_identifiers);

  listenstart = thedoubledlist;


  while(charact_identifiers != NULL)
  {
    merke = thedoubledlist;             /*merkt sich das Listenende*/

    if ((thedoubledlist = (IDF_SHORT_TYP *)GetTmpMem(theHeap,sizeof(IDF_SHORT_TYP),ANS_MarkKey))==NULL)
    {
      PrintErrorMessage('E',"CopyCharacteristicList2HelpList","no memory obtained from GetMem(sizeof(IDF_SHORT_TYP))");
      return (NULL);
    }

    IDF_SHORT_NXT(thedoubledlist) = NULL;             /*das neue Listenende*/
    IDF_SHORT_NXT(merke) = thedoubledlist;
    IDF_SHORT_VAL(thedoubledlist) = IDF_VAL(charact_identifiers);

    charact_identifiers = IDF_NXT(charact_identifiers);
  }

  /*only for debugging*/
  helplf = listenstart;
  while (helplf != NULL)
  {
    helpvar = IDF_SHORT_VAL(helplf);
    helplf = IDF_SHORT_NXT(helplf);
  }

  return(listenstart);
}



/****************************************************************************/
/**
* Laeuft ueber alle Subdomains und sucht dabei die mit der ID sbdmid.
*
* @param sbdmid The Id of the searched Subdomain.
*
* @return Pointer to the found subdomain with the Id sbdmid, NULL if error occured resp. no such subdomain could be found
   */
/****************************************************************************/
SD_TYP *FindSubdomain(INT sbdmid)
{
  SD_TYP *lauf_sd;

  lauf_sd = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);


  while(lauf_sd != NULL)
  {
    if(SD_NAME(lauf_sd) == sbdmid)
    {
      /*gefunden!*/
      return(lauf_sd);
    }
    else
    {
      lauf_sd = SD_NEXT(lauf_sd);
    }
  }

  /*wenn er hier her kommt, dann wurde keine solche Subdomain gefunden*/
  PrintErrorMessage('E',"FindSubdomain","Did not find a subdomain with th ID sbdmid");
  return (NULL);
}



/****************************************************************************/
/**
* Die Surface sfce erhaelt einen neuen PolylineEintrag plptr !!!
* Trage die Adresse der Polyline plptr in die Surface sfce ein und inkremetiere
* die AnzahlderPolylinesVon der Surface um 1.
*
* @param  plptr Pointer to polyline which belongs to surface.
* @param  sfce Pointer to surface.
*
* @return Pointer to new SurfacePolylineEntry if OK,
*         NULL if error occured resp. no SurfacePolylineEntry could be created
   */
/****************************************************************************/
SFPL_TYP *MakeNewSfcPlEntry(PL_TYP *plptr, SF_TYP *sfce)
{
  SFPL_TYP *mp;

  /*merke Dir den Anfang der Polylineeintraege von der Surface sfce*/
  mp = SF_POLYLINES(sfce);

  if((SF_POLYLINES(sfce) = (SFPL_TYP *)GetTmpMem(theHeap,sizeof(SFPL_TYP),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"MakeNewSfcPlEntry","no memory obtained from GetMem(sizeof(SFPL_TYP))");
    return (NULL);
  }

  SFPL_NEXT(SF_POLYLINES(sfce)) = mp;
  SFPL_PL(SF_POLYLINES(sfce)) = plptr;

  /*Die Surface sfce hat eine Polylione mehr:*/
  SF_NMB_OF_POLYLINES(sfce) ++;

  return(SF_POLYLINES(sfce));
}



/****************************************************************************/
/**
* Funktion such in einer Identifier liste einen Eintrag mit Wert sfn.
* Der Listenanfang kommt mit Idfi in die Funktion
* Im Referenzparameter pre_Idfi wird der Vorgaenger des gesuchten/gefundenen Identifiers
* nach aussen zurueckgegeben.
*
* @param Idfi Pointer to the first Identifier of the IdfList.
* @param pre_Idfi Pointer to predescessor of Idfi == !!reference parameter, which is changed
*                 within the function resp. updated to the predescessor of the Identifier with
*                 SurfaceValue sfn.
* @param sfn DOUBLE-Value of the searched Identifier.
*
* @return Pointer to the Identifier with the Value sfn,
*         NULL if no such Identifier could be found
   */
/****************************************************************************/
IDF_SHORT_TYP *SearchPartner(IDF_SHORT_TYP *Idfi, IDF_SHORT_TYP **pre_Idfi, DOUBLE sfn)
{
  DOUBLE sfcevalue;

  /*Solange das Listenende der HilfsIdentifierliste nicht erreicht ist (das kann uU gleich zu Beginn sein)
     Suche nach dem Eintrag mit Wert sfn ...*/
  while(Idfi != NULL)
  {
    sfcevalue = IDF_SHORT_VAL(Idfi);

    if(sfcevalue == sfn)
    {
      return(Idfi);
    }
    else             /*weiter mit naechstem Identifier der linearen Liste*/
    {
      *pre_Idfi = Idfi;                   /*VorgaengerUpdate*/
      Idfi = IDF_SHORT_NXT(Idfi);
    }
  }

  /*wenn die while-Schleife vollstaendig durchlaufen wurde, so bedeutet das, dass kein Identifer
     mit dem Wert sfn bzw. kein Partner gefunden wurde ...*/
  return(NULL);
}



/****************************************************************************/
/**
* abhaengig von den charakteristischen Identifiers die Polyline bei den Surfaces eintragen
* verwendet wird dabei eine temporaere Hilfsliste helplist, die die charakteristische IDFs
* einer Polyline in aufsteigender Reihenfolge beinhaltet. (z.B aabcd)
*
* Idee: Zu jedem Eintrag  werden die Surfaces durchgesucht. Findet man eine mit einem
* gleichen Identifier, dann wird diese Polyline mit der gefundenen Surface verknuepft.
* Danach kann  der Eintrag der Helpliste geloescht werden und mit dem nexten fortgefahren werden.
* Findet man gar eine Surface mit 2 Identifiern, von denen einer mit dem Helplisteneintrag
* uebereinstimmt, so laeuft man mit dem zweiten ueber den Rest der Helplist und sucht ob
* man diesen Zweiten vieleicht auch noch finden kann. Ist dem so, so hat man eine Polyline
* von einer DoppelIDSurface gefunden - man kann verknuepfen und danach aus der Helplist
* gleich zwei Eintraege loeschen.
*
* @param plptr Pointer to new Polyline.
*
* @return value 0: ok, 1:error.
   */
/****************************************************************************/
INT ConnectPolylineWithSurfaces(PL_TYP *plptr)
{
  IDF_SHORT_TYP *helplist;       /*dient dazu die charak. IDFs der Polyline plptr temporaer zu speichern*/
  IDF_SHORT_TYP *Idf;       /* = Variable fuer den ersten Identifier der Hilfsliste*/
  IDF_SHORT_TYP *mIdf;       /* = Merk-Variable fuer den ersten Identifier der Hilfsliste bevor die
                                  Funktion SearchPartner aufgerufen wird mit Idf->next  --
                                  mIdf wird auch an SearchPartner mituebergeben - allerdings als Referenz-
                                  parameter und kann so in SearchPartner so verandert werden, das nach aussen
                                  in die aufruf. Fkt ConnectPolylineWithSurfaces der Vorgaenger des gefundenen
                                  IDentifierPartners uebergeben werden kann.*/
  IDF_SHORT_TYP *partnerIdf;       /* = Variable fuer den zweiten Identifier der zusammen mit dem ersten
                                      fuer eine Doppelsurface steht bzw. diese repraesentiert */
  INT sbdmid;
  SD_TYP *sbd;
  SFC_TYP *sfce;
  INT gefunden;       /* F bzw. T, wenn zur Polyline eine Surface gefunden wurde*/
  SFPL_TYP *new_sfce_plyln;

  /* Kopiere die charakteristische List der eingegangenen Polyline plptr in eine Hilfsliste helplist
     Diese Hilfsliste wird dann nach und nach bearbeitet - Ist sie leer, so ist die Polyline
     ausreichend bearbeitet. In die Helplist werden allerdinngs nur die DOUBLE-Identifiers eingtragen,
     und nicht mehr ! Verwendet wird daher auch nur der abgespeckte Datentyp IDF_SHORT_TYP */
  if((helplist = CopyCharacteristicList2HelpList(PL_IDFS(plptr))) == NULL)
  {
    PrintErrorMessage('E',"ConnectPolylineWithSurfaces","did receive nilpointer from CopyCharacteristicList2HelpList");
    return (1);
  }
  else
  {
    /*solange die helplist != NULL bzw. solange es etwas zum Arbeiten gibt bzw.
       solange die temporaere Hilfsliste etwas beinhaltet*/
    while(helplist != NULL)
    {
      Idf = helplist;                   /*man nehme den ersten Identifier der Hilfsliste, welcher aufgrund
                                           der Sortierung auch der niederwertigste ist. */
      sbdmid = (int) (floor(IDF_SHORT_VAL(Idf)));

      /* Suche die zugehoerige Subdomain : */
      if((sbd = FindSubdomain(sbdmid)) == NULL)
      {
        PrintErrorMessage('E',"ConnectPolylineWithSurfaces","no subdomain found: NULL returnd by FindSubdomain");
        return (1);
      }
      else
      {
        gefunden = F;
        sfce = SD_SFCS(sbd);
        /*laufe ueber die surfaces und suche Surfaces, die zur Polyline passen*/
        while((sfce != NULL) && (gefunden == F))
        {
          /*wenn der Surfacename aus nur einer DOUBLE-Zahl besteht (d.h. also keine Doppelsurface)*/
          if(SF_NAME2(SFC_SURF(sfce)) == 0.0)
          {
            /*wenn der Identifier der Hilfsliste mit dem einen Surfaceidentifier uebereinstimmt...*/
            if(IDF_SHORT_VAL(Idf) == SF_NAME1(SFC_SURF(sfce)))
            {
              /* passende Surface gefunden !!!*/
              if((new_sfce_plyln = MakeNewSfcPlEntry(plptr,SFC_SURF(sfce))) == NULL)
              {
                PrintErrorMessage('E',"ConnectPolylineWithSurfaces","no new SingleSurface-Polyline created : MakeNewSfcPlEntry returnd NULL");
                return (1);
              }

              /*Update der temoraeren Hilfsliste - Identifier abgearbeitet :*/
              gefunden = T;
              helplist = IDF_SHORT_NXT(helplist);

              /*No ReleaseTmpMem here ...*/
              /*free(Idf, sizeof(IDF_SHORT_TYP));*/
            }
          }
          else                               /* es liegt eine Doppelsurface vor, die zwei Subdomains gemeinsam als Trennsurface besitzen */
          {
            mIdf = Idf;                                     /*MerkDir den Vorgaenger von Idf->next ...*/
            /*wenn der Identifier gleich dem niederwertigeren Surfaceidentifer*/
            if(IDF_SHORT_VAL(Idf) == SF_NAME1(SFC_SURF(sfce)))
            {
              /*... dann suche in der restlichen IdfsHilfsliste nach dem hoeherwertigeren Surfaceidentifer:*/
              if((partnerIdf = SearchPartner(IDF_SHORT_NXT(Idf), &mIdf, SF_NAME2(SFC_SURF(sfce)))) != NULL)
              {
                if((new_sfce_plyln = MakeNewSfcPlEntry(plptr,SFC_SURF(sfce))) == NULL)
                {
                  PrintErrorMessage('E',"ConnectPolylineWithSurfaces","no new DoubleSurface-Polyline created : MakeNewSfcPlEntry returnd NULL");
                  return (1);
                }

                /*Update der temoraeren Hilfsliste - 2 weitere Identifiers sind abgearbeitet :*/
                gefunden = T;
                IDF_SHORT_NXT(mIdf) = IDF_SHORT_NXT(partnerIdf);
                /*free(partnerIdf, sizeof(IDF_SHORT_TYP)); *//* No ReleaseTmpMem here ! free ??? */
                helplist = IDF_SHORT_NXT(helplist);
                /*free(Idf, sizeof(IDF_SHORT_TYP));*/ /* No ReleaseTmpMem here !  free ??? */
              }
              /*wenn kein Partner gefunden wurde ==> weiter mit naechster Surface*/
            }
            /* TODO TOTHINK : der umgekehrte Fall, d.h IDF_SHORT_VAL(Idf) == SF_NAME2(SFC_SURF(sfce)) ist
               meiner Meinung nach nicht zu betrachten, da man dabei sowieso nichts finden wuerde
               Begruendung: Sowohl die Idfs der Hilfsliste als auch die jew. 2 Idfs einer Surface
               sind  aufsteigend sortiert. ===> wenn also Hilfsliste->firstIDF !=  Surface->firstIDF
               aber Hilfsliste->firstIDF == Surface->secondIDF, dann hat man sowieso keine Chance
               in der Hilfsliste den Partner Surface->firstIDF zu finden da gilt:
               Surface->firstIDF < Surface->secondIDF == Hilfsliste->firstIDF < Hilfsliste->allenachfolger
               siehe Quelltext im Kopnzept - hier in der Implementierung wurde er zunaechst einmal
               weggelassen*/
          }
          sfce = SFC_NEXT(sfce);                               /* weiter mit naechster Surface*/
        }                        /*while sfce ...*/
      }                  /*else*/
    }            /*while*/
  }      /*else*/
  return(0);
}






/****************************************************************************/
/**
* This functions splits a polyline, which coud not be sorted by SortPolyline() !
* The function will split the already sorted part of the polyline and create a new polyline
* out of it. The remaining part wil be the old  Polyline with change features.
*
*  Idee:
*  -----
*  der bereits sortierte Teile der Polyline wird zu einer neuen Polyline gemacht.
*  Der verbleibende Teil wird upgedatet.
*
* @param anfang refeerence parameter - will be changed in fct.
* @param rechtesMuster refeerence parameter - will be changed in fct.
* @param Polyline the corresp. Polyline - pointer will be used to change polyline features
* @param  idl2 value parameter
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT PolylineSplit(PL_LINE_TYP **anfang, PL_LINE_TYP **rechtesMuster, PL_TYP *Polyline, PL_LINE_TYP *idl2)
{
  PL_TYP *merkeplptr,*neuePolyline;
  PL_LINE_TYP *lauf_PLL,*endeabgespPL;
  INT anzbersorPLLPoints;
  INT rv;


  anzbersorPLLPoints = 2;

  lauf_PLL = idl2;
  while(lauf_PLL != (*rechtesMuster))
  {
    lauf_PLL = PL_LINES_NXT(lauf_PLL);
    anzbersorPLLPoints++;
  }


  endeabgespPL = *rechtesMuster;


  /*alte Polyline*/
  /*=> *anfang, *rechtesMuster := ganz vorne usw. ... UPDATEN */
  (*rechtesMuster) = PL_LINES_NXT((*rechtesMuster));      /*ganz vorne*/
  if((*rechtesMuster) == NULL)
  {
    PrintErrorMessage('E',"PolylineSplit","PolylineSpliiting makes no sense - no remaining Polyline");
    return (1);
  }
  else
  {
    (*anfang) = PL_LINES_NXT((*rechtesMuster));
  }

  PL_NMB_OF_POINTS(Polyline) = PL_NMB_OF_POINTS(Polyline) - anzbersorPLLPoints + 1;
  PL_LINES(Polyline) = *rechtesMuster;
  /*PL_IDFS, PL_NMB_OF_CH_IDFS und PL_NXT aendern sich nicht !!!*/



  /*abgespaltete Polyline Teil2 ..*/
  PL_LINES_NXT(endeabgespPL) = NULL;
  /*=> bereits sortierte PlLinelist aufabrbeiten (z.B.:NullPointer ans ende)*/

  /*mit der abgespaltenen eine neue anlegen ...*/
  /*laeuft von idl2 bis endeabgespPL.*/
  merkeplptr = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  if((neuePolyline = (PL_TYP *)GetTmpMem(theHeap,sizeof(PL_TYP),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"PolylineSplit","got no mem for the new polyline, which split");
    return (1);
  }
  /*dass sie vom SortPolyline-DURCHLAUF nicht mehr erfasst wird also am besten an den
     Listenanfang*/
  PL_NXT(neuePolyline) = merkeplptr;
  PL_NMB_OF_CH_IDFS(neuePolyline) = PL_NMB_OF_CH_IDFS(Polyline);
  PL_IDFS(neuePolyline) = PL_IDFS(Polyline);
  PL_LINES(neuePolyline) = idl2;
  PL_NMB_OF_POINTS(neuePolyline) = anzbersorPLLPoints;
  NMB_OF_PLINES(DomainInfo_Pointer) ++;
  EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer) = neuePolyline;
  if((rv = ConnectPolylineWithSurfaces(neuePolyline)) == 1)
  {
    PrintErrorMessage('E',"PolylineSplit","Error occured calling ConnectPolylineWithSurfaces");
    return (1);
  }


  return (0);
}






/****************************************************************************/
/**
*  This functions sorts the the sequence of points of the polyline Polyline
*  in the way : "a,b" "b,g", "g,j", "j,k" usw. und aendert dabei u.U. auch
*  die Reihenfolge der Identifiers der Lines.*
*
*  Idee
*  -----
*  => erste Line als Anfangs-Muster und Ende-Muster nehmen
*  => durchlaufe alle anderen und suche eine Partnerline, d.h eine Line, die
*     mit der Anfangs-Muster-Line oder aber der Ende-Muster-Line einen gemeinsamen Punkt besitzt.
*  => wenn Partnerline gefunden, d.f. Partnerline links oder rechts des bereits
*     sortierten Listenteils einfuegen und Liste an Entnahmestelle updaten.
*  => die gefundene Partnerline wird zum neuen Anfangs-Muster oder Ende-Muster
*  => von vorne anfangen bis man einmal durch die Liste durch ist
*
* @param Polyline The Polyline to be sorted
*
* @return SORTED if ok, 1 if error occured.
   D*/
/****************************************************************************/
INT SortPolyline(PL_TYP *Polyline)
{
  PL_LINE_TYP *idl, *idl2;       /*Laufvariable durch die Lineliste*/
  PL_LINE_TYP *idl_pred;       /*Variable fuer den Vorgaenger von idl*/
  PL_LINE_TYP *anfang;       /*zeigt auf die erste Line der noch zu sortierenden Lines*/
                             /*zeigt auf den Anfang der noch zu sortierenden LineListe (Schwanz)*/
                             /*wenn anfang NULL wird, ist der Sortiervorgang beendet*/
  PL_LINE_TYP *new_partnerline;       /*Variable fuer die gefundene neue Partnerline, die an den
                                                                        Listenkopf eingefuegt werden soll und die dann als neues
                                                                        Muster dienen soll*/
  PL_LINE_TYP *new_partnerline_pred;       /*Variable fuer den Vorgaenger der neuen Partnerline*/
  PL_LINE_TYP *mpll;       /*MerkVariable fuer PolylineLines*/
  PL_LINE_TYP *rechtesMuster;       /*PolylineLine, die vor dem noch zu sortierenden Listenschwanz steht
                                       also vor anfang ! Sie ist das Muster das nach rechts  verwendet wird -
                                       und bestimmt somit die Reihenfolge der Polyline mit
                                       rechtesMuster steht immer vor anfang*/

  INT predlinefound;       /*Flag, das angibt, ob eine Partner- bzw. Vorgaengerline gefunden wurde.*/
  INT succlinefound;       /*Flag, das angibt, ob eine Partner- bzw. Nachfolgerline RECHTS gefunden wurde.*/
  INT merkeID;
  INT rv;



  if((idl = PL_LINES(Polyline)) == NULL)
  {
    PrintErrorMessage('E',"SortPolyline","Polyline has no LineEntries !!!");
    return (1);
  }
  /*wenn die Polyline mind. 2 Lines besitzt ...*/
  else if(PL_LINES_NXT(idl) != NULL)       /*otherwise not necessary to sort !!! only one line!*/
  {
    rechtesMuster = idl;            /*LineListenelement vor anfang (Schwanz) merken, fuer den Fall
                                       dass anfang nach vorne muss und new_partnerline_pred noch Null ist*/

    /*der erste LineEintrag Idl der Polyline ist das Muster
       Nun wird eine Partner-/NachbarLIne gesucht beginnend mit dem Anfang Idl-next*/
    anfang = PL_LINES_NXT(idl);

    /*Solange es noch etwas zu sortieren gibt ...*/
    while(anfang != NULL)
    {
      idl = anfang;
      new_partnerline_pred = NULL;
      idl_pred = NULL;
      predlinefound = F;
      succlinefound = F;

      /*Solange noch keine Nachbarline zum Muster gefunden ist UND
         falls das Listenende noch nicht erreicht ist : VorgaengerLineSuchen*/
      while((idl != NULL)&&(predlinefound == F)&&(succlinefound == F))
      {
        /*Ist die Line idl der Vorgaenger der LINKKENLine
           bzw haben diese beiden Lines einen gemeinsamen Punkt*/
        idl2 = PL_LINES(Polyline);
        if(  (LI_NDID1(PL_LINES_LINE(idl)) == LI_NDID1(PL_LINES_LINE(idl2))) ||
             (LI_NDID1(PL_LINES_LINE(idl)) == LI_NDID2(PL_LINES_LINE(idl2))) ||
             (LI_NDID2(PL_LINES_LINE(idl)) == LI_NDID1(PL_LINES_LINE(idl2))) ||
             (LI_NDID2(PL_LINES_LINE(idl)) == LI_NDID2(PL_LINES_LINE(idl2)))    )
        {
          /*Vorgaengerline gefunden !*/
          predlinefound = T;
          new_partnerline = idl;
          new_partnerline_pred = idl_pred;
        }
        /*... oder aber ein Nachfolger der RECHTEN Line*/
        else if(  (LI_NDID1(PL_LINES_LINE(rechtesMuster)) == LI_NDID1(PL_LINES_LINE(idl))) ||
                  (LI_NDID1(PL_LINES_LINE(rechtesMuster)) == LI_NDID2(PL_LINES_LINE(idl))) ||
                  (LI_NDID2(PL_LINES_LINE(rechtesMuster)) == LI_NDID1(PL_LINES_LINE(idl))) ||
                  (LI_NDID2(PL_LINES_LINE(rechtesMuster)) == LI_NDID2(PL_LINES_LINE(idl)))    )
        {
          /*Nachfolgerline gefunden !*/
          succlinefound = T;
          new_partnerline = idl;
          new_partnerline_pred = idl_pred;
        }
        idl_pred = idl;
        idl = PL_LINES_NXT(idl);
      }                   /*while((idl != NULL)&&(predlinefound == F))*/

      /*wenn keine Vorgaengerline gefunden wurde ...*/
      if(predlinefound == F)
      {
        if(succlinefound == F)
        {
          if((rv = PolylineSplit(&anfang,&rechtesMuster,Polyline,idl2)) == 1)
          {
            PrintErrorMessage('E',"PolylineSplit","returned ERROR");
            return(1);
          }
        }
        else                         /*ein RECHTER Nachfolger wurde gefunden*/
        {

          /* Fuer den FAll  dass (new_partnerline_pred == NULL) ist gilt:
             Reihenfolge der Lines kann in diesem Fall gleich bleiben
             Die gefundene neue Partnerline liegt bereits an der richtigen
             Stelle: Sie ist der unmittelbare Nachfolger des RECHTEN Musters!

             Wenn nicht ...==>UPDATE DER LISTE ...*/
          if(new_partnerline_pred != NULL)
          {
            /*gefundene Partnerline an gef. Stelle entnehmen und RECHTS in Liste einhaengen
               sowie Update an der Entnahmestelle und Update am RECHTEN Listenende durchfuehren :*/
            mpll = PL_LINES_NXT(rechtesMuster);                                     /*Nachfolger des rechten Musters merken*/
            PL_LINES_NXT(rechtesMuster) = new_partnerline;                                     /*neuer Nachfolger des rechten Musters ist
                                                                                                  die gefundene Line new_partnerline*/
            PL_LINES_NXT(new_partnerline_pred) = PL_LINES_NXT(new_partnerline);                                    /*Nachfolger der gefundenen Line merken*/
            PL_LINES_NXT(new_partnerline) = mpll;
          }

          /*Update der IDs*/
          /*Nur wenn die ID-Reihenfolge  des soeben gefundenen neuen PartnerLine nicht (!)
             zu der ID-Reihenfolge des bisherigen RECHTENListenendes passen
             dann muessen die beiden IDS der neuen Partnerline vertauscht werden.*/
          /*Der Fall, dass es sich um den ersten Nachbar zum ersten Muster handelt landet
             nie hier sondern geht stets in die LINKS-Schleife, da zuerst auf links abgefragft wird
             und indiesem Falle ja links und rechts einfuegbar waere.*/
          if(LI_NDID2(PL_LINES_LINE(rechtesMuster)) != LI_NDID1(PL_LINES_LINE(new_partnerline)))
          {
            merkeID = LI_NDID2(PL_LINES_LINE(new_partnerline));
            LI_NDID2(PL_LINES_LINE(new_partnerline)) = LI_NDID1(PL_LINES_LINE(new_partnerline));
            LI_NDID1(PL_LINES_LINE(new_partnerline)) = merkeID;


          }


          /*Im Fall eines neuen RECHTEN Nachbars :*/
          /* Es muessen stets anfang und rechtes Muster jeweils um 1 weiter geschaltet werden*/
          rechtesMuster = PL_LINES_NXT(rechtesMuster);
          /*ferner aendert sich in diesem Fall der Anfang des noch zu bearbeitenden Listenteils :*/
          anfang = PL_LINES_NXT(rechtesMuster);


        }                         /* else . . . ein RECHTER Nachfolger wurde gefunden*/
      }
      else                   /*ein LINKER Vorgaenger wurde gefunden*/
      {
        /*wenn die gefundene Line gleich dem Anfang des noch zu sortierenden
           Teils (Schwanz) der Lineliste ist ...*/
        if(new_partnerline_pred == NULL)
        {
          /*Da new_partnerline_pred NOCH NULL IST bzw die gefundene Line == anfang ist,
             wir als new_partnerline_pred das erste Muster, d.h. der Listenkopf der unsortierten
             Liste verwendet.*/
          new_partnerline_pred = rechtesMuster;

          /*ferner aendert sich in diesem Fall der Anfang des noch zu bearbeitenden Listenteils :*/
          anfang = PL_LINES_NXT(anfang);

        }

        /*gefundene Partnerline an gef. Stelle entnehmen und vorne in Liste einhaengen
           sowie Update an der Entnahmestelle und Update am Listenkopf durchfuehren :*/
        mpll = PL_LINES(Polyline);                         /*Listenanfang merken*/
        PL_LINES(Polyline) = new_partnerline;                         /*neuer Listenkopf == gefundene Line*/
        PL_LINES_NXT(new_partnerline_pred) = PL_LINES_NXT(new_partnerline);                         /*Entnahmestelle ueberbruecken*/
        PL_LINES_NXT(new_partnerline) = mpll;                        /*der bisherige Listenkopf wird als Nachfolger vom neuen Listen-
                                                                        kopf eingetragen*/
        /*anfang ist richtig gesetzt und wurde nur im obigen Falle
           if(new_partnerline_pred == NULL)
           weitergeschaltet*/

        /*Update der IDs*/
        /*Nur wenn die ID-Reihenfolge  des soeben gefundenen neuen PartnerLine nicht (!)
           zu der ID-Reihenfolge des bisherigen Listenkopfes passen
           dann muessen die beiden IDS der neuen Partnerline vertauscht werden.*/
        if(LI_NDID2(PL_LINES_LINE(new_partnerline)) != LI_NDID1(PL_LINES_LINE(PL_LINES_NXT(new_partnerline))))
        {
          merkeID = LI_NDID2(PL_LINES_LINE(new_partnerline));
          LI_NDID2(PL_LINES_LINE(new_partnerline)) = LI_NDID1(PL_LINES_LINE(new_partnerline));
          LI_NDID1(PL_LINES_LINE(new_partnerline)) = merkeID;

          /*Nur Nach Ausrichtung der beiden ersten Lines kann sein das die KnotenIDs
             immer noch nicht stimmen - die erste Line wird immer linksvorne eingetragen*/
          if(LI_NDID2(PL_LINES_LINE(new_partnerline)) != LI_NDID1(PL_LINES_LINE(PL_LINES_NXT(new_partnerline))))
          {
            merkeID = LI_NDID2(PL_LINES_LINE(PL_LINES_NXT(new_partnerline)));
            LI_NDID2(PL_LINES_LINE(PL_LINES_NXT(new_partnerline))) = LI_NDID1(PL_LINES_LINE(PL_LINES_NXT(new_partnerline)));
            LI_NDID1(PL_LINES_LINE(PL_LINES_NXT(new_partnerline))) = merkeID;

            /*Fuer die Faelle ab_ac, ab_ca und ba_ac stimmt nun alles
               Fuer den Fall ba_ca gilt aber:
               Aus ba_ca wurde zunaechst ab_ca und dann soeben ab_ac ===>
               Es stimmt alls immer noch nicht und es mus eine dritte Vertauschung
               zu ba_ac erfolgen . Dann stimmts aber*/
            if(LI_NDID2(PL_LINES_LINE(new_partnerline)) != LI_NDID1(PL_LINES_LINE(PL_LINES_NXT(new_partnerline))))
            {
              merkeID = LI_NDID2(PL_LINES_LINE(new_partnerline));
              LI_NDID2(PL_LINES_LINE(new_partnerline)) = LI_NDID1(PL_LINES_LINE(new_partnerline));
              LI_NDID1(PL_LINES_LINE(new_partnerline)) = merkeID;
            }

          }
        }

      }                   /* else ein LINKER Vorgaenger wurde gefunden*/
    }
  }
  return(0);
}



/****************************************************************************/
/**
*  Erzeugt alle Polylines durch kompletten Lauf ueber die LI-HAshtabelle
*  und sortieret sie anschliessend
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgmCreatePloylines()
{
  INT rv, LI_HT_Index;
  LI_KNOTEN_TYP *li_ptr;
  PL_TYP *plptr, *newpl, *lf_plptr;
  PL_LINE_TYP *newpllineptr;

  INT justfordebugging;
  justfordebugging =0;
  /*laufe ueber die gesamte LI-Hashtabelle*/
  for(LI_HT_Index = 0; LI_HT_Index < LI_p; LI_HT_Index++)
  {
    /*Existiert hier ueberhaupt ein Eintrag ?*/
    li_ptr = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[LI_HT_Index];
    if(li_ptr != NULL)
    {
      /*laufe ueber die lineare Liste, die mit diesem Hashpointer li_ptr beginnt*/
      while(li_ptr != NULL)
      {


        /*Liegt diese Line ueberhaupt auf einer Polyline?*/
        if((rv = Check_If_Line_On_Polyline(LI_IDFS(li_ptr))) == T)
        {
          /*Existiert die Polyline zu dieser Line vielleicht bereits schon
             oder besteht sie noch nicht ???*/
          if((plptr = Exist_Polyline(li_ptr)) == NULL)
          {
            /*neue Polyline anlegen*/
            if((newpl = GetMemFillAddNewPolyline(li_ptr)) == NULL)
            {
              PrintErrorMessage('E',"Ansys2lgmCreatePloylines","did receive nilpointer from GetMemAndFillNewPolyline");
              return (1);
            }

            /*diese neue Polyline in zugehoerigen Surfaces eintragen*/
            if((rv = ConnectPolylineWithSurfaces(newpl)) == 1)
            {
              PrintErrorMessage('E',"Ansys2lgmCreatePloylines","did receive nilpointer from GetMemAndFillNewPolyline");
              return (1);
            }
          }
          else                               /*die Polyline existiert also bereits*/
          {
            if((newpllineptr = GetMemFillAddNewPolylineLine(li_ptr,plptr)) == NULL)
            {
              PrintErrorMessage('E',"Ansys2lgmCreatePloylines","did receive nilpointer from GetMemFillAddNewPolylineLine");
              return (1);
            }
          }
        }
        else if (rv == LI_IFS_ERROR)
        {
          PrintErrorMessage('E',"Ansys2lgmCreatePloylines","did receive ERRORVALUE from fct Check_If_Line_On_Polyline");
          return (1);
        }
        /*else rv = 0, d.h. diese Line liegt auf keiner Polyline!*/
        li_ptr = LI_NEXT(li_ptr);
      }
    }
  }      /*for*/


  /* laufe ueber alle Polylines und sortiere die einzelnen Lines (=Kantenzuege) so,
     dass eine Sortierung der Art ij, jk, ko, oc,cv,vf, usw. usw., bzgl der KnotenIDs vorliegt.*/
  lf_plptr = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  while(lf_plptr != NULL)
  {
    if((rv = SortPolyline(lf_plptr)) == SORTED)
    {
      lf_plptr = PL_NXT(lf_plptr);
    }
    else
    {
      PrintErrorMessage('E',"Ansys2lgmCreatePloylines","did not receive SORTED-Message from SortPolyline");
      return (1);
    }
  }

  return(0);
}




/****************************************************************************/
/**
* this function create a new Polylinecycle for surface theSurface.
*
* @param  anfang reference parameter -> will be changed in this Fct
* @param  rechtesMuster reference parameter -> will be changed in this fct
* @param  theSurface Pointer to the correspomdinmg surface - itf features will be
*                    changed in this fct.
* @param  idl2 value parameter
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT GetMemAndFillNewPlz(SFPL_TYP **anfang,SFPL_TYP **rechtesMuster,SF_TYP *theSurface,SFPL_TYP *idl2)
{
  /* vergleiche PolylineSplit()*/
  PLZ_TYP *merkeplzptr,*neuerPolylinezyk;
  SFPL_TYP *lauf_PLZ,*endeabgespSF_PLs;
  INT anzbersorPLs;
  PL_LINE_TYP *last_Line_of_idl2, *last_Line_of_endeabgespSF_PLs;

  anzbersorPLs = 1;

  lauf_PLZ = idl2;
  while(lauf_PLZ != (*rechtesMuster))
  {
    lauf_PLZ = SFPL_NEXT(lauf_PLZ);
    anzbersorPLs++;
  }


  endeabgespSF_PLs = *rechtesMuster;


  /*alte Surface->Polylines*/
  /*=> *anfang, *rechtesMuster := ganz vorne usw. ... UPDATEN */
  (*rechtesMuster) = SFPL_NEXT((*rechtesMuster));      /*ganz vorne*/
  if((*rechtesMuster) == NULL)
  {
    /*In diesen Teil laeuft das Programm nur im Falle, dass in Create_PLZN  mit dem verbleibenden Rest
       auch noch ein PLZ erzeugt wird.*/
    if((*anfang) != NULL)
    {
      PrintErrorMessage('E',"GetMemAndFillNewPlz","anfang == NULL is not possible");
      return (1);
    }
  }
  else
  {
    (*anfang) = SFPL_NEXT((*rechtesMuster));
  }
  /* die Anzahl der Polyline lassen, da die Surface ja immer noch alle Polylines hat! */
  /*SF_NMB_OF_POLYLINES(theSurface) = SF_NMB_OF_POLYLINES(theSurface) - anzbersorPLs;*/
  SF_POLYLINES(theSurface) = *rechtesMuster;



  /*abgespaltete Surface->Polylines Teil2 ..*/
  SFPL_NEXT(endeabgespSF_PLs) = NULL;
  /*=> bereits sortierte SfcePlinelist aufabrbeiten (d.h.:NullPointer ans ende)*/

  /*mit der abgespaltenen eine neue anlegen ...*/
  /*laeuft von idl2 bis endeabgespSF_PLs.*/
  merkeplzptr = SF_POLYLI_ZYK(theSurface);
  if((neuerPolylinezyk = (PLZ_TYP *)GetTmpMem(theHeap,sizeof(PLZ_TYP),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"GetMemAndFillNewPlz","got no mem for the new polylinecycle");
    return (1);
  }
  /*eingefuegt wird am Listenanfang*/
  PLZ_NEXT(neuerPolylinezyk) = merkeplzptr;
  PLZ_POLYLINES(neuerPolylinezyk) = idl2;
  PLZ_NMB_OF_POLYLINES(neuerPolylinezyk) = anzbersorPLs;
  SF_NMB_OF_POLYLI_ZYK(theSurface)++;
  SF_POLYLI_ZYK(theSurface) = neuerPolylinezyk;


  /*Probe: Ist der neue Polylinezyklus wirklich zyklisch ?*/
  /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
  /*erste oder zweite LI_NDID der ersten Line muss einen Gemeinsamen besitzen mit
     der ersten oder zweiten LI_NDID der letzten Line sein . . .*/
  last_Line_of_idl2 = PL_LINES(SFPL_PL(idl2));
  while(PL_LINES_NXT(last_Line_of_idl2) != NULL)
  {
    last_Line_of_idl2 = PL_LINES_NXT(last_Line_of_idl2);
  }
  last_Line_of_endeabgespSF_PLs = PL_LINES(SFPL_PL(endeabgespSF_PLs));
  while(PL_LINES_NXT(last_Line_of_endeabgespSF_PLs) != NULL)
  {
    last_Line_of_endeabgespSF_PLs = PL_LINES_NXT(last_Line_of_endeabgespSF_PLs);
  }
  if(
    (LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(idl2)))) != LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(endeabgespSF_PLs)))))
    && (LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(idl2)))) != LI_NDID2(PL_LINES_LINE(last_Line_of_endeabgespSF_PLs)))
    && (LI_NDID2(PL_LINES_LINE(last_Line_of_idl2)) != LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(endeabgespSF_PLs)))))
    && (LI_NDID2(PL_LINES_LINE(last_Line_of_idl2)) != LI_NDID2(PL_LINES_LINE(last_Line_of_endeabgespSF_PLs)))
    )
  {
    PrintErrorMessage('E',"Create_PLZN","Surface has got a PolylineZyklus which is not cyclic !");
    return (1);
  }


  return (0);
}


/****************************************************************************/
/**
* this function detects the polylinecycles of a surface
* a polylinecycle is a set of polylines, which can be combined to
* a close polygon.
* Each surface has at leats one Polylinecycle.
* The polylines itself are already sorted.
*
* @param theSurface - the Surface to detect its polylinecycles.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Create_PLZN(SF_TYP *theSurface)
{

  SFPL_TYP *idl, *idl2;       /*Laufvariable durch die PolyLineliste*/
  SFPL_TYP *idl_pred;       /*Variable fuer den Vorgaenger von idl*/
  SFPL_TYP *anfang;       /*zeigt auf die erste Polyline der noch zu untersuchenden PolyLines*/
                          /*zeigt auf den Anfang der noch zu sortierenden PolyLineListe (Schwanz)*/
                          /*wenn anfang NULL wird, ist der Sortiervorgang beendet*/
  SFPL_TYP *new_partnerline;       /*Variable fuer die gefundene neue PartnerPolyline, die an den
                                                                        Listenkopf eingefuegt werden soll und die dann als neues
                                                                        Muster dienen soll*/
  SFPL_TYP *new_partnerline_pred;       /*Variable fuer den Vorgaenger der neuen PartnerPolyline*/
  SFPL_TYP *mpll;       /*MerkVariable fuer Polylines*/
  SFPL_TYP *rechtesMuster;       /*Polyline, die vor dem noch zu sortierenden Listenschwanz steht
                                     also vor anfang ! Sie ist das Muster das nach rechts  verwendet wird -
                                     und bestimmt somit die Reihenfolge der Polylinezyklen
                                      rechtesMuster steht immer vor anfang*/

  INT predlinefound;       /*Flag, das angibt, ob eine Partner- bzw. VorgaengerPolyline gefunden wurde.*/
  INT succlinefound;       /*Flag, das angibt, ob eine Partner- bzw. NachfolgerPolyline RECHTS gefunden wurde.*/
  INT rv;
  INT lff;
  PL_LINE_TYP *polylinelinesofidl;

  INT start_idl,start_idl2,start_rechtesMuster,end_idl,end_idl2,end_rechtesMuster;



  if((idl = SF_POLYLINES(theSurface)) == NULL)
  {
    PrintErrorMessage('E',"Create_PLZN","Surface has no PolyLineEntries !!!");
    return (1);
  }
  /*wenn die Surface mind. 2 PolyLines besitzt ...*/
  else if(SFPL_NEXT(idl) != NULL)       /*otherwise not necessary to sort !!! only one Polyline!*/
  {
    rechtesMuster = idl;            /*PolyLineListenelement vor anfang (Schwanz) merken, fuer den Fall
                                       dass anfang nach vorne muss und new_partnerline_pred noch Null ist*/

    /*der erste PolyLineEintrag Idl der Surface ist das Muster
       Nun wird eine Partner-/NachbarLIne gesucht beginnend mit dem Anfang Idl-next*/
    anfang = SFPL_NEXT(idl);

    /*Solange es noch etwas zu untersuchen gibt ...*/
    while(anfang != NULL)
    {
      idl = anfang;
      new_partnerline_pred = NULL;
      idl_pred = NULL;
      predlinefound = F;
      succlinefound = F;

      /*Solange noch keine NachbarPolyline zum Muster gefunden ist UND
         falls das Listenende noch nicht erreicht ist : VorgaengerPolyLineSuchen*/
      while((idl != NULL)&&(predlinefound == F)&&(succlinefound == F))
      {
        idl2 = SF_POLYLINES(theSurface);

        /*Im folgenden werden die Start und Endpunkte der 3 betrachteten Polylines bestimmt*/

        start_idl = LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(idl))));

        start_idl2 = LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(idl2))));

        start_rechtesMuster = LI_NDID1(PL_LINES_LINE(PL_LINES(SFPL_PL(rechtesMuster))));

        polylinelinesofidl = PL_LINES(SFPL_PL(idl));
        for (lff = 2; lff < PL_NMB_OF_POINTS(SFPL_PL(idl)); lff++)                         /*Lauf zur letzten Line dieser Polyline*/
        {
          polylinelinesofidl = PL_LINES_NXT(polylinelinesofidl);
        }
        end_idl = LI_NDID2(PL_LINES_LINE(polylinelinesofidl));


        polylinelinesofidl = PL_LINES(SFPL_PL(idl2));
        for (lff = 2; lff < PL_NMB_OF_POINTS(SFPL_PL(idl2)); lff++)                         /*Lauf zur letzten Line dieser Polyline*/
        {
          polylinelinesofidl = PL_LINES_NXT(polylinelinesofidl);
        }
        end_idl2 = LI_NDID2(PL_LINES_LINE(polylinelinesofidl));

        polylinelinesofidl = PL_LINES(SFPL_PL(rechtesMuster));
        for (lff = 2; lff < PL_NMB_OF_POINTS(SFPL_PL(rechtesMuster)); lff++)                         /*Lauf zur letzten Line dieser Polyline*/
        {
          polylinelinesofidl = PL_LINES_NXT(polylinelinesofidl);
        }
        end_rechtesMuster = LI_NDID2(PL_LINES_LINE(polylinelinesofidl));


        /*Ist die SfcePOlyLine idl der Vorgaenger der LINKKEN_SfcePolyLine
           bzw. ist einer der vier STart-/Endpunkte bei beiden SFCPLYLNS vertreten ?*/
        if(  (start_idl == start_idl2) ||
             (start_idl == end_idl2) ||
             (end_idl == end_idl2) ||
             (end_idl == start_idl2)    )
        {
          /*VorgaengerPolyline gefunden !*/
          predlinefound = T;
          new_partnerline = idl;
          new_partnerline_pred = idl_pred;
        }
        /*... oder aber ein Nachfolger der RECHTEN PolyLine*/
        else if(  (start_idl == start_rechtesMuster) ||
                  (start_idl == end_rechtesMuster) ||
                  (end_idl == end_rechtesMuster) ||
                  (end_idl == start_rechtesMuster)    )
        {
          /*NachfolgerPolyline gefunden !*/
          succlinefound = T;
          new_partnerline = idl;
          new_partnerline_pred = idl_pred;
        }
        idl_pred = idl;
        idl = SFPL_NEXT(idl);
      }                   /*while((idl != NULL)&&(predlinefound == F))*/

      /*wenn keine VorgaengerPolyline gefunden wurde ...*/
      if(predlinefound == F)
      {
        /*wenn auch keine NachfolgerPolyline gefunden wurde ...*/
        if(succlinefound == F)
        {
          if((rv = GetMemAndFillNewPlz(&anfang,&rechtesMuster,theSurface,idl2)) == 1)
          {
            PrintErrorMessage('E',"GetMemAndFillNewPlz","returned ERROR");
            return(1);
          }
        }
        else                         /*ein RECHTER Nachfolger wurde gefunden*/
        {

          /* Fuer den FAll  dass (new_partnerline_pred == NULL) ist gilt:
             Reihenfolge der PolyLines kann in diesem Fall gleich bleiben
             Die gefundene neue PartnerPolyline liegt bereits an der richtigen
             Stelle: Sie ist der unmittelbare Nachfolger des RECHTEN Musters!

             Wenn nicht ...==>UPDATE DER LISTE ...*/
          if(new_partnerline_pred != NULL)
          {
            /*gefundene PartnerPolyline an gef. Stelle entnehmen und RECHTS in Liste einhaengen
               sowie Update an der Entnahmestelle und Update am RECHTEN Listenende durchfuehren :*/
            mpll = SFPL_NEXT(rechtesMuster);                                     /*Nachfolger des rechten Musters merken*/
            SFPL_NEXT(rechtesMuster) = new_partnerline;                                     /*neuer Nachfolger des rechten Musters ist
                                                                                                die gefundene PolyLine new_partnerline*/
            SFPL_NEXT(new_partnerline_pred) = SFPL_NEXT(new_partnerline);                                    /*Nachfolger der gefundenen PolyLine merken*/
            SFPL_NEXT(new_partnerline) = mpll;
          }


          /*Im Fall eines neuen RECHTEN Nachbars :*/
          /* Es muessen stets anfang und rechtes Muster jeweils um 1 weiter geschaltet werden*/
          rechtesMuster = SFPL_NEXT(rechtesMuster);
          /*ferner aendert sich in diesem Fall der Anfang des noch zu bearbeitenden Listenteils :*/
          anfang = SFPL_NEXT(rechtesMuster);


        }                         /* else . . . ein RECHTER Nachfolger wurde gefunden*/
      }
      else                   /*ein LINKER Vorgaenger wurde gefunden*/
      {
        /*wenn die gefundene PolyLine gleich dem Anfang des noch zu sortierenden
           Teils (Schwanz) der PolyLineliste ist ...*/
        if(new_partnerline_pred == NULL)
        {
          /*Da new_partnerline_pred NOCH NULL IST bzw die gefundene PolyLine == anfang ist,
             wird als new_partnerline_pred das erste Muster, d.h. der Listenkopf der unsortierten
             Liste verwendet.*/
          new_partnerline_pred = rechtesMuster;

          /*ferner aendert sich in diesem Fall der Anfang des noch zu bearbeitenden Listenteils :*/
          anfang = SFPL_NEXT(anfang);

        }

        /*gefundene Partnerpolyline an gef. Stelle entnehmen und vorne in Liste einhaengen
           sowie Update an der Entnahmestelle und Update am Listenkopf durchfuehren :*/
        mpll = SF_POLYLINES(theSurface);                         /*Listenanfang merken*/
        SF_POLYLINES(theSurface) = new_partnerline;                         /*neuer Listenkopf == gefundene PolyLine*/
        SFPL_NEXT(new_partnerline_pred) = SFPL_NEXT(new_partnerline);                         /*Entnahmestelle ueberbruecken*/
        SFPL_NEXT(new_partnerline) = mpll;                        /*der bisherige Listenkopf wird als Nachfolger vom neuen Listen-
                                                                       kopf eingetragen*/
        /*anfang ist richtig gesetzt und wurde nur im obigen Falle
           if(new_partnerline_pred == NULL)
           weitergeschaltet*/


      }                   /* else ein LINKER Vorgaenger wurde gefunden*/
    }            /* von while */
  }      /*von ... wenn die Surface mind. 2 PolyLines besitzt ...*/

  /*Erzeuge mit dem Rest auch noch einen PLZyklus - aber nur wenn es nicht der einzige ist ...*/
  if(SF_NMB_OF_POLYLI_ZYK(theSurface) >0)
  {
    /* Achtung: idl2 benoetigt ein Update*/
    idl2 = SF_POLYLINES(theSurface);
    if((rv = GetMemAndFillNewPlz(&anfang,&rechtesMuster,theSurface,idl2)) == 1)
    {
      PrintErrorMessage('E',"GetMemAndFillNewPlz","returned ERROR");
      return(1);
    }
  }
  return(0);
}




/****************************************************************************/
/**
* This fct searches the (only!) correspondig SFE_Triangle of the Polyline-Line
* line with the same surface identifiers like theSurface.

* @param line - pointer to the Line of a polyline
* @param theSurface - pointer to the Surface
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
SFE_KNOTEN_TYP *Find_SFE_Triangle(LI_KNOTEN_TYP *line,SF_TYP *theSurface)
{
  INT hw;
  LI_KNOTEN_TYP *hp;
  IDF_TYP *lfptr;
  INT gefunden;
  SFE_KNOTEN_TYP *returnSFETriangle;

  /*Hashwert/Schluessel berechnen*/
  if(LI_NDID1(line) < LI_NDID2(line))
  {
    hw = the_LI_hashfunction(LI_NDID1(line),LI_NDID2(line));
  }
  else
  {
    hw = the_LI_hashfunction(LI_NDID2(line),LI_NDID1(line));
  }
  hp = (EXCHNG_TYP2_LI_HASHTAB(ExchangeVar_2_Pointer))[hw];

  /*Beim Hashwert den richtigen Listeneintrag suchen*/
  if (hp != NULL)
  {
    while( (LI_NDID1(hp) != LI_NDID1(line)) && (LI_NDID2(hp) != LI_NDID2(line)) )
    {
      hp = LI_NEXT(hp);
    }
  }

  if(hp == NULL)
  {
    PrintErrorMessage('E',"Find_SFE_Triangle","could not find the Line in the LI-Hashtable");
    return (NULL);
  }

  lfptr = LI_IDFS(hp);
  if(lfptr == NULL)
  {
    PrintErrorMessage('E',"Find_SFE_Triangle","the found LI-HashTable-Entry has no(!) IDF-Pointer!");
    return (NULL);
  }

  gefunden = F;
  /* Laufe ueber alle Eintraege der Identifierliste der Line hp */
  while( (lfptr != NULL) )
  {

    /*wenn die Identifiers des zugehoerigen Triangles mit denen der Surface zusammenpassen*/
    if ( (SFE_IDF_0(IDF_SFE_TRIA(lfptr)) == SF_NAME1(theSurface)) &&
         (SFE_IDF_1(IDF_SFE_TRIA(lfptr)) == SF_NAME2(theSurface))    )
    {
      if(gefunden == F)
      {
        gefunden = T;
        returnSFETriangle = IDF_SFE_TRIA(lfptr);
      }
      else
      {
        /*wenn ein zweites Dreieck gefunden wurde...*/
        if( (IDF_SFE_TRIA(lfptr)) != returnSFETriangle)
        {
          PrintErrorMessage('E',"Find_SFE_Triangle","es wurden zwei(!!!) moegliche SFE_Triangles gefunden");
          return (NULL);
        }
        /*else adressengleiches Dreieck wurde gefunden ->
           Dies ist unbedenklich, da im Spezialfalle
           einer zu splittenden Doppelsurface passieren kann
           vergleiche ersten discovered Bug im Oktober*/
      }
    }

    lfptr = LI_NEXT(lfptr);
  }

  if(gefunden == F)
  {
    PrintErrorMessage('E',"Find_SFE_Triangle","did not find the SFE_Triangle");
    return(NULL);
  }
  else
  {
    /* SFE_Triangle wurde gefunden :*/
    return(returnSFETriangle);
  }
}



/****************************************************************************/
/**
*
* @param SFE_search
* @param SFE_destination
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT TriangleNeighbourSearcher(SFE_KNOTEN_TYP *SFE_search,SFE_KNOTEN_TYP *SFE_destination)
{
  SFE_KNOTEN_TYP *Nachbar_SFE;
  INT kante;       /*Laufvariable ueber die 3 Kanten eines Dreiecks*/
  INT neubesetzt[3];
  INT rgbwrt;


  neubesetzt[0] = F; neubesetzt[1] = F; neubesetzt[2] = F;

  /* Laufe ueber die 3 Nachbarn von SFE_search */
  for(kante = 0; kante < 3; kante++)
  {
    Nachbar_SFE =  SFE_NGHB(SFE_search,kante);
    if(Nachbar_SFE != NULL)
    {
      /*natuerlich nur wenn dieser Nachbar nicht schon betrachtet wurde ...*/
      if(SFE_ORIENTATION_FLAG(Nachbar_SFE) == F)
      {

        SFE_ORIENTATION_FLAG(Nachbar_SFE) = T;
        neubesetzt[kante] = T;
        if(Nachbar_SFE == SFE_destination)
        {
          triangle_found = T;
          return(FERTIG);
        }
        /*************************************************************/
      }
    }
  }       /* von for */

  /*weiterer Lauf ueber die 3 Nachbarn, um die RekursionsHierarchie etwas kleiner zu halten*/
  for(kante = 0; kante < 3; kante++)
  {
    /*nur die soeben neu betrachteten Dreiecke sind fuer einen Rekursionschritt interessant*/
    if(neubesetzt[kante] == T)
    {
      rgbwrt = TriangleNeighbourSearcher(SFE_NGHB(SFE_search,kante),SFE_destination);
      if(triangle_found == T)
      {
        return(FERTIG);
      }
    }
  }       /* von for */

  return(FERTIG);

}




/****************************************************************************/
/**
* this function creates a new RealSurface with one or different PLZs for surface theSurface.
*
* @param  anfang reference parameter -> will be changed in this Fct
* @param  rechtesMuste reference parameter -> will be changed in this fct
* @param  theSurface Pointer to the correspomdinmg surface - itf features will be
*                    changed in this fct.
* @param  idl2 value parameter.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT GetMemAndFillNewRlSfc(PLZ_TYP **anfang,PLZ_TYP **rechtesMuster,SF_TYP *theSurface,PLZ_TYP *idl2)
{
  /* vergleiche PolylineSplit()*/
  RS_TYP *merkeplzptr,*neueRealSurface;
  PLZ_TYP *lauf_PLZ,*endeabgespSF_PLZs;
  INT anzbersorPLZs;

  anzbersorPLZs = 1;

  lauf_PLZ = idl2;
  while(lauf_PLZ != (*rechtesMuster))
  {
    lauf_PLZ = SFPL_NEXT(lauf_PLZ);
    anzbersorPLZs++;
  }


  endeabgespSF_PLZs = *rechtesMuster;


  /*alte Surface->Polylinecycles*/
  /*=> *anfang, *rechtesMuster := ganz vorne usw. ... UPDATEN */
  (*rechtesMuster) = PLZ_NEXT((*rechtesMuster));      /*ganz vorne*/
  if((*rechtesMuster) == NULL)
  {
    /*In diesen Teil laeuft das Programm nur im Falle, dass in Create_RealSurfaces mit dem verbleibenden Rest
       auch noch eine RS erzeugt wird.*/
    if((*anfang) != NULL)
    {
      PrintErrorMessage('E',"GetMemAndFillNewRlSfc","anfang == NULL is not possible");
      return (1);
    }
  }
  else
  {
    (*anfang) = PLZ_NEXT((*rechtesMuster));
  }
  /* die Anzahl der Polylinezyklen belassen, da die Surface ja immer noch alle Polylinezyklen hat! */
  /*SF_NMB_OF_POLYLI_ZYK(theSurface) = SF_NMB_OF_POLYLI_ZYK(theSurface) - anzbersorPLZs;*/
  SF_POLYLI_ZYK(theSurface) = *rechtesMuster;



  /*abgespaltete Surface->Polylinezyklen Teil2 ..*/
  PLZ_NEXT(endeabgespSF_PLZs) = NULL;
  /*=> bereits sortierte SfcePlineCyclelist aufarbeiten (d.h.:NullPointer ans ende)*/

  /*mit der abgespaltenen eine neue anlegen ...*/
  /*laeuft von idl2 bis endeabgespSF_PLZs.*/
  merkeplzptr = SF_REALSFCS(theSurface);
  if((neueRealSurface = (RS_TYP *)GetTmpMem(theHeap,sizeof(RS_TYP),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"GetMemAndFillNewRlSfc","got no mem for the new realsurface");
    return (1);
  }
  /*eingefuegt wird am Listenanfang*/
  RS_NEXT(neueRealSurface) = merkeplzptr;
  RS_PL_ZKLN(neueRealSurface) = idl2;
  RS_NMB_OF_PL_ZKLN(neueRealSurface) = anzbersorPLZs;
  SF_NMB_OF_REALSFCS(theSurface)++;
  SF_REALSFCS(theSurface) = neueRealSurface;

  return (0);
}



/****************************************************************************/
/**
*  This function reconstructs the Surface->Polyline-Entries in the case of
*  no realsurfaces but different created polylinecycles
*  Fct. tests whether all polylines have been put back.
*  Fct. tests whether all polylinecycles do have all their polylines.
*
* @param  theSurface pointer to the original surface which perhaps has to be split.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT ReconstructSurfacePolylines(SF_TYP *theSurface)
{
  PLZ_TYP *laufvar_plz;
  SFPL_TYP *lauf_sfce_polyline,*merkeSfceplanfang,*anfang, *ende;
  INT nmbofplz;
  INT i,numberofpolylinesofpolylinecycle,numberofpolylinesofsurface;

  /*laufe ueber alle Polylinezyklen von theSurface*/
  nmbofplz = SF_NMB_OF_POLYLI_ZYK(theSurface);
  if(nmbofplz < 2)
  {
    PrintErrorMessage('E',"ReconstructSurfacePolylines","Surface schoud have at least 2 PLZs");
    return (1);
  }
  laufvar_plz = SF_POLYLI_ZYK(theSurface);
  if(laufvar_plz == NULL)
  {
    PrintErrorMessage('E',"ReconstructSurfacePolylines","Surface has no PLZ  at all");
    return (1);
  }
  numberofpolylinesofsurface =0;
  for(i=1; i<=nmbofplz; i++)
  {
    if(laufvar_plz == NULL)
    {
      PrintErrorMessage('E',"ReconstructSurfacePolylines","Surface has not enough PLZs");
      return (1);
    }

    /*Wiedereinbauen der Polylines von laufvar_plz in Surface->Polylines*/
    lauf_sfce_polyline = PLZ_POLYLINES(laufvar_plz);
    anfang = lauf_sfce_polyline;
    numberofpolylinesofpolylinecycle = 0;
    while(lauf_sfce_polyline != NULL)
    {
      numberofpolylinesofpolylinecycle++;
      ende = lauf_sfce_polyline;
      lauf_sfce_polyline = SFPL_NEXT(lauf_sfce_polyline);

    }
    /*Probe*/
    if(numberofpolylinesofpolylinecycle != PLZ_NMB_OF_POLYLINES(laufvar_plz))
    {
      PrintErrorMessage('E',"ReconstructSurfacePolylines","A PLZ has too much or too less polylines");
      return (1);
    }

    numberofpolylinesofsurface += numberofpolylinesofpolylinecycle;

    merkeSfceplanfang = SF_POLYLINES(theSurface);
    SF_POLYLINES(theSurface) = anfang;
    SFPL_NEXT(ende) = merkeSfceplanfang;
    /*durch diese 3 Zeilen wurden alle Polylines des betrachteten Poylinezyklus laufvar_plz
       an die alte Stelle Surface->Polyline zurueckgehaengt.*/

    laufvar_plz = PLZ_NEXT(laufvar_plz);
  }

  /*Probe*/
  if(numberofpolylinesofsurface != SF_NMB_OF_POLYLINES(theSurface))
  {
    PrintErrorMessage('E',"ReconstructSurfacePolylines","Surface has reconstructed too much or too less polylines  with PLZs");
    return (1);
  }

  if(laufvar_plz != NULL)
  {
    PrintErrorMessage('E',"ReconstructSurfacePolylines","Surface has too much PLZs");
    return (1);
  }
  /*... und haenge die Polylinelsten in die Surface->Polylines*/
  /* ueberpruefe dabei die Anzahlen der Polylinelisten pointermaessig*/
  return(0);
}





/****************************************************************************/
/**
*  This fct splits the polylinecycles of a surface in one or several group(s)
*  of polylinecycles which describe a "real" surface.
*  idea: Any triangle T1 of a "real" surface can be reached from any other triangle T2
*        of the same "real" surface using triangleneighbourhoodrelations.
*
* @param  theOrigSfce pointer to the original surface which perhaps has to be split
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Create_RealSurfaces(SF_TYP *theOrigSfce)
{

  PLZ_TYP *idl, *idl2;       /*Laufvariable durch die PolyLineZyklusliste*/
  PLZ_TYP *idl_pred;       /*Variable fuer den Vorgaenger von idl*/
  PLZ_TYP *anfang;       /*zeigt auf den ersten PolylineZyklus der noch zu untersuchenden PolyLineZyklen*/
  /*zeigt auf den Anfang der noch zu sortierenden PolyLineZyklenListe (Schwanz)*/
  /*wenn anfang NULL wird, ist der Sortiervorgang beendet*/
  PLZ_TYP *new_partnerplz;       /*Variable fuer den gefundenen neuen PartnerPolylineZyklus, der stets an den
                                                                        Listenkopf eingefuegt wird  und der dann als neues
                                                                        Muster dient*/
  PLZ_TYP *new_partnerplz_pred;       /*Variable fuer den Vorgaenger der neuen PartnerPolylineZyklus*/
  PLZ_TYP *mpll;       /*MerkVariable fuer PolylineZyklen*/
  PLZ_TYP *rechtesMuster;       /*Polylinezyklus, der vor dem noch zu sortierenden Listenschwanz steht
                                     also vor anfang ! Er ist das Muster das nach rechts wird aber nie verwendet
                                     da stets links am Listenanfang eingefuegt wird im Falle der Erzeugeung von RSs -
                                      rechtesMuster steht immer vor anfang*/

  INT predplzfound;       /*Flag, das angibt, ob ein Partner- bzw. VorgaengerPolylineZyklus gefunden wurde.*/
  INT rv,rw;

  LI_KNOTEN_TYP *erste_line_von_idl,*erste_line_von_idl2;
  SFE_KNOTEN_TYP *SFE_destination,*SFE_search;
  TRIANGLE_TYP *lauf_tria;




  if((idl = SF_POLYLI_ZYK(theOrigSfce)) == NULL)
  {
    PrintErrorMessage('E',"Create_RealSurfaces","Surface has no PolyLineCycles !!!");
    return (1);
  }
  /*wenn die Surface mind. 2 PolyLineZyklen besitzt ...*/
  else if(PLZ_NEXT(idl) != NULL)       /*otherwise not necessary to sort !!! only one Polylinecycle!*/
  {
    rechtesMuster = idl;            /*PolyLineZyklusListenelement vor anfang (Schwanz) merken, fuer den Fall
                                       dass anfang nach vorne muss und new_partnerplz_pred noch Null ist*/

    /*der erste PolyLineZyklusEintrag Idl der Surface ist das Muster
       Nun wird eine Partner-/NachbarLIne gesucht beginnend mit dem Anfang Idl-next*/
    anfang = PLZ_NEXT(idl);

    /*Solange es noch etwas zu untersuchen gibt ...*/
    while(anfang != NULL)
    {
      idl = anfang;
      new_partnerplz_pred = NULL;
      idl_pred = NULL;
      predplzfound = F;

      /*Solange noch kein NachbarPolylineZyklus zum Muster gefunden ist UND
         falls das Listenende noch nicht erreicht ist : VorgaengerPolyLineZyklusSuchen*/
      while((idl != NULL)&&(predplzfound == F))
      {
        idl2 = SF_POLYLI_ZYK(theOrigSfce);

        /*Ist der SfcePOlyLineZyklus idl der Vorgaenger ders LINKKEN_SfcePolyLinezyklus idl2 ?*/

        /* d.h. nehme erste Line von der ersten Polyline von idl*/
        erste_line_von_idl = PL_LINES_LINE(PL_LINES(SFPL_PL(PLZ_POLYLINES(idl))));

        /* und ermiitle zugehoeriges SFE-Dreieck SFE_destination*/
        if( (SFE_destination = Find_SFE_Triangle(erste_line_von_idl,theOrigSfce)) == NULL)
        {
          PrintErrorMessage('E',"Create_RealSurfaces","could not find SFE_destination with fct. Find_SFE_Triangle");
          return (1);
        }

        /* nehme erste Line von der ersten Polyline von idl2 */
        erste_line_von_idl2 = PL_LINES_LINE(PL_LINES(SFPL_PL(PLZ_POLYLINES(idl2))));

        /* und ermiitle zugehoeriges SFE-Dreieck SFE_search*/
        if( (SFE_search = Find_SFE_Triangle(erste_line_von_idl2,theOrigSfce)) == NULL)
        {
          PrintErrorMessage('E',"Create_RealSurfaces","could not find SFE_search with fct. Find_SFE_Triangle");
          return (1);
        }

        SFE_ORIENTATION_FLAG(SFE_search) = T;
        triangle_found = F;
        /*Spez.fall: Stimmen womoeglich gleich die ersten beiden ueberein?*/
        if(SFE_search == SFE_destination)
        {
          triangle_found = T;
          SFE_ORIENTATION_FLAG(SFE_search) = F;
        }

        /*Normalfall:*/
        else if ( (rw = TriangleNeighbourSearcher(SFE_search,SFE_destination)) != FERTIG)
        {
          PrintErrorMessage('E',"Create_RealSurfaces"," Returnvalue of TriangleNeighbourSearcher was not FERTIG - Problems with searching triangle");
          return (1);
        }
        else
        {
          /*Achtung ! Update der orientationflags aller Dreiecke dieser Surface ist notwendig*/
          /*SFE_ORIENTATION_FLAG*/
          /*laufe ueber alle Dreiecke dieser Surface und setze das OrientationFlag wieder auf F*/
          /* damit es fuer TriangleIDOrientations wieder stimmt*/
          lauf_tria = SF_TRIAS(theOrigSfce);
          while(lauf_tria != NULL)
          {
            SFE_ORIENTATION_FLAG(TRIA_SFE_KN(lauf_tria)) = F;
            lauf_tria = TRIA_NEXT(lauf_tria);
          }
        }

        if(triangle_found == T)
        {
          /*VorgaengerPolylineZyklus gefunden !*/
          predplzfound = T;
          new_partnerplz = idl;
          new_partnerplz_pred = idl_pred;
        }

        idl_pred = idl;
        idl = SFPL_NEXT(idl);


      }                   /*while((idl != NULL)&&(predplzfound == F))*/


      /*wenn kein VorgaengerPolylineZyklus gefunden wurde ...*/
      if(predplzfound == F)
      {
        /*hier laeuft das Programm ldgl. bei der letzten RS nicht 'rein.*/
        if((rv = GetMemAndFillNewRlSfc(&anfang,&rechtesMuster,theOrigSfce,idl2)) == 1)
        {
          PrintErrorMessage('E',"GetMemAndFillNewRlSfc","returned ERROR");
          return(1);
        }
      }
      else                   /*ein LINKER Vorgaenger wurde gefunden*/
      {
        /*wenn der gefundene PolyLineZyklus gleich dem Anfang des noch zu sortierenden
           Teils (Schwanz) der PolyLineZyklusliste ist ...*/
        if(new_partnerplz_pred == NULL)
        {
          /*Da new_partnerplz_pred NOCH NULL IST bzw der gefundene PolyLineZyklus == anfang ist,
             wird als new_partnerplz_pred das erste Muster, d.h. der Listenkopf der unsortierten
             Liste verwendet.*/
          new_partnerplz_pred = rechtesMuster;

          /*ferner aendert sich in diesem Fall der Anfang des noch zu bearbeitenden Listenteils :*/
          anfang = PLZ_NEXT(anfang);

        }

        /*gefundener PartnerpolylineZyklus an gef. Stelle entnehmen und vorne in Liste einhaengen
           sowie Update an der Entnahmestelle und Update am Listenkopf durchfuehren :*/
        mpll = SF_POLYLI_ZYK(theOrigSfce);                         /*Listenanfang merken*/
        SF_POLYLI_ZYK(theOrigSfce) = new_partnerplz;                         /*neuer Listenkopf == gefundener PolyLineZyklus*/
        PLZ_NEXT(new_partnerplz_pred) = PLZ_NEXT(new_partnerplz);                         /*Entnahmestelle ueberbruecken*/
        PLZ_NEXT(new_partnerplz) = mpll;                        /*der bisherige Listenkopf wird als Nachfolger vom neuen Listen-
                                                                       kopf eingetragen*/
        /*anfang ist richtig gesetzt und wurde nur im obigen Falle
           if(new_partnerplz_pred == NULL)
           weitergeschaltet*/

      }                   /* else ein LINKER Vorgaenger wurde gefunden*/
    }            /* von while */
  }      /*von ... wenn die Surface mind. 2 PolyLineZyklen besitzt ...*/

  /*Erzeuge mit dem Rest auch noch eine RealSurface - aber nur wenn es nicht die einzige ist ...*/
  if(SF_NMB_OF_REALSFCS(theOrigSfce) > 0)
  {
    /* Achtung: idl2 benoetigt ein Update*/
    idl2 = SF_POLYLI_ZYK(theOrigSfce);

    if((rv = GetMemAndFillNewRlSfc(&anfang,&rechtesMuster,theOrigSfce,idl2)) == 1)
    {
      PrintErrorMessage('E',"GetMemAndFillNewPlz","returned ERROR");
      return(1);
    }
  }
  else
  {
    /*Achtung, wenn nur eine einzige RealSurface entstehen wuerde dann muessen die Surface-Polylines
       wieder angelegt werden einfach aus den PLZ heraus.*/
    if((rv = ReconstructSurfacePolylines(theOrigSfce)) == 1)
    {
      PrintErrorMessage('E',"Create_RealSurfaces","ReconstructSurfacePolylines returned ERROR");
      return(1);
    }
  }
  return(0);
}




/****************************************************************************/
/**
*  Diese rekursive Funktion durchlaeuft alle Dreiecke einer Surface
*  mit HIlfe der schon bestehenden Nachbarschaftsbeziehungen zwischen
*  den Dreiecken. Dabei wird nach jedem ueberprueften resp. einge-
*  tragenen Dreieck die statische fuer die Funktion global verwendbare
*  Variable "nmb_of_triangles" inkrementiert. Im Gegensatz zur Fkt TRIANGLEIDORIENTATIONS
*  wird bei dieser Fkt. die gesamte Rekursionshierarchie durchlaufen.
*  Besuchte Dreiecke erhalten das Flag und werden dabei nich nochmals besucht.
*
* @param dasDreieck  SFE, dessen drie Nachbarn durchlaufen werden und in die SurfaceListe
*                    eingetragen werden.
*
* @return FERTIG if alle Dreiecke einmal besucht bzw.
*         static INT zaehler == static INT nmb_of_trias_of_sf,
*         1 if error occured.
   */
/****************************************************************************/
INT FetchAllTriangles(SFE_KNOTEN_TYP *dasDreieck)
{
  SFE_KNOTEN_TYP *Nachbar_SFE;
  INT kante;       /*Laufvariable ueber die 3 Kanten eines Dreiecks*/
  INT neubesetzt[3];
  INT rgbwrt;
  TRIANGLE_TYP *merke_triangle_pointer;

  neubesetzt[0] = F; neubesetzt[1] = F; neubesetzt[2] = F;

  /* Laufe ueber die 3 Nachbarn von dasDreieck */
  for(kante = 0; kante < 3; kante++)
  {
    Nachbar_SFE =  SFE_NGHB(dasDreieck,kante);
    if(Nachbar_SFE != NULL)
    {
      /*natuerlich nur wenn dieser Nachbar nicht schon aufgenommen wurde ...*/
      if(SFE_ORIENTATION_FLAG(Nachbar_SFE) == F)
      {
        /* Dreieck vorne in Liste einfuegen ... */
        merke_triangle_pointer = New_Triangle_List;
        if((New_Triangle_List = (TRIANGLE_TYP *)GetTmpMem(theHeap,sizeof(TRIANGLE_TYP),ANS_MarkKey))== NULL)
        {
          PrintErrorMessage('E',"SplitSurface","got  no memory  for  New_Triangle_List !?!");
          return(1);
        }
        TRIA_SFE_KN(New_Triangle_List) = Nachbar_SFE;
        TRIA_NEXT(New_Triangle_List) = merke_triangle_pointer;
        nmb_of_triangles++;




        /*************************************************************/
        /* Dieser Nachbar ist nun in die New_Triangle_List eingefuegt, folglich . . . */
        SFE_ORIENTATION_FLAG(Nachbar_SFE) = T;
        neubesetzt[kante] = T;
        /*************************************************************/
      }
    }
  }       /* von for */

  /*weiterer Lauf ueber die 3 Nachbarn, um die RekursionsHierarchie etwas kleiner zu halten*/
  for(kante = 0; kante < 3; kante++)
  {
    /*nur die soeben neu besetzten Dreiecke sind fuer einen Rekursionschritt interessant*/
    if(neubesetzt[kante] == T)
    {
      rgbwrt = FetchAllTriangles(SFE_NGHB(dasDreieck,kante));
    }
  }       /* von for */

  return(FERTIG);

}




/****************************************************************************/
/**
*  fct splits surface theSurface using its realSurfaces
*  and delets the old surface
*  all necessary aupdate are executed.
*
*
* @param theSurface the Surface which must be spliited
* @param thePredSurface  the predecessor of the surface
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT SplitSurface(SF_TYP *theSurface, SF_TYP *thePredSurface)
{
  SF_TYP *sf_lfv2;
  SF_TYP *newSurface,*merkeSurface;
  SD_TYP *the_sbd, *the_sbd_Double;
  SFC_TYP *lauf_sd_sfces, *lauf_sd_sfces_Double, *pred_sfc_OF_the_sbd, *pred_sfc_OF_the_sbd_Double,  *sf_of_sd_lauf,  *sf_of_sd_lauf_Double, *neue_sfc_of_sbd, *neue_sfc_of_sbd_Double;
  INT sbdmid,sbdmid_Double,lauf_rs,lauf_plz,numberofpolylinesofpolylinecycle,rv;
  RS_TYP *theRealSurface;
  PLZ_TYP *thePolylineCycle;
  SFPL_TYP *lauf_sfce_polyline,*merkeSfceplanfang,*anfang, *ende;
  LI_KNOTEN_TYP *ersteline;
  SFE_KNOTEN_TYP *erstesdreieck;
  INT dreiecksanzahl;
  TRIANGLE_TYP *lauf_tria;
  DOUBLE new_sfc_name,new_sfc_name2;
  TRIANGLE_TYP *lauf_sftria;
  INT DoubleSurfaceCase;      /*Flag ob Surface Beruehrungssurface ist oder nicht*/

  DoubleSurfaceCase = F;

  if(SF_NAME2(theSurface) != SEC_SFC_NAME_DEFAULT_VAL)
  {
    DoubleSurfaceCase = T;             /*Es liegt eine Beruehrungssurface vor*/
  }


  /*suche Subdomain the_sbd der Surface*/
  sbdmid = (int) (floor(SF_NAME1(theSurface)));
  if(DoubleSurfaceCase == T)
  {
    sbdmid_Double = (int) (floor(SF_NAME2(theSurface)));
  }
  /* Suche die zugehoerige Subdomain : */
  if((the_sbd = FindSubdomain(sbdmid)) == NULL)
  {
    PrintErrorMessage('E',"SplitSurface","no subdomain found: NULL returnd by FindSubdomain");
    return (1);
  }
  if(DoubleSurfaceCase == T)
  {
    if((the_sbd_Double = FindSubdomain(sbdmid_Double)) == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","no subdomain found:sbdmid_Double  NULL returnd by FindSubdomain");
      return (1);
    }
  }


  /*suche den Surfaceeintrag von the_sbd und ermittle dessen Vorgaenger pred_sfc_OF_the_sbd*/
  pred_sfc_OF_the_sbd = NULL;
  lauf_sd_sfces = SD_SFCS(the_sbd);
  if(lauf_sd_sfces == NULL)
  {
    PrintErrorMessage('E',"SplitSurface","lauf_sd_sfces is NULL at begin and theSurface was not found yet ");

    return (1);
  }
  while(SFC_SURF(lauf_sd_sfces) != theSurface)
  {
    pred_sfc_OF_the_sbd = lauf_sd_sfces;
    lauf_sd_sfces = SFC_NEXT(lauf_sd_sfces);
    if(lauf_sd_sfces == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","lauf_sd_sfces has just become NULL and theSurface was not found yet ");
      return (1);
    }
  }

  if(DoubleSurfaceCase == T)
  {
    /*suche den Surfaceeintrag von the_sbd_Double und ermittle dessen Vorgaenger pred_sfc_OF_the_sbd_Double*/
    pred_sfc_OF_the_sbd_Double = NULL;
    lauf_sd_sfces_Double = SD_SFCS(the_sbd_Double);
    if(lauf_sd_sfces_Double == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","lauf_sd_sfces_Double is NULL at begin and theSurface was not found yet ");
      return (1);
    }
    while(SFC_SURF(lauf_sd_sfces_Double) != theSurface)
    {
      pred_sfc_OF_the_sbd_Double = lauf_sd_sfces_Double;
      lauf_sd_sfces_Double = SFC_NEXT(lauf_sd_sfces_Double);
      if(lauf_sd_sfces_Double == NULL)
      {
        PrintErrorMessage('E',"SplitSurface","lauf_sd_sfces_Double has just become NULL and theSurface was not found yet ");
        return (1);
      }
    }
  }


  lauf_rs=0;
  theRealSurface = SF_REALSFCS(theSurface);

  dreiecksanzahl = 0;

  /*laufe ueber die real surfaces und erzeuge neue Surfaces*/
  while(theRealSurface != NULL)
  {


    lauf_rs++;

    /*erzeuge neue surface und trage sie vorne in SF-Liste ein.*/
    merkeSurface = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
    if((newSurface = (SF_TYP *)GetTmpMem(theHeap,sizeof(SF_TYP),ANS_MarkKey))== NULL)
    {
      PrintErrorMessage('E',"SplitSurface","got  no memory  for  newSurface !?!");
      return(1);
    }
    SF_NEXT(newSurface) = merkeSurface;
    EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer) = newSurface;

    /*Attention! hier wird noch der alte Surfacename verwendet, da diese noc in Find_SFE_Triangle
       benoetigt wird und ausserdem alle zugehoerigen Dreiecke nach FetchAllTriangles den neuen
       SUrfacenamen erhalten sollten
       TODO ?? Sollten auch die Lines d.h. die LI-Hashtabelle upgedatet werden ?
       Meiner Meinung <25.8.97> nach NICHT NOETIG, da es im UG Lines in meinem Verstaendnis
       nicht mehr gibt Triangles sehr wohl
       Die Triangle->SUrfaceNamen werden allerdings spaeter auch nicht mehr benoetigt
       oder ?? In jedewm Fall besteht in dieser Funktion <siehe unten> ein gute Moeglichkeit
       zum Update ...*/
    SF_NAME1(newSurface) = SF_NAME1(theSurface);
    SF_NAME2(newSurface) = SF_NAME2(theSurface);

    SF_RIGHT_SBD(newSurface)  = SF_RIGHT_SBD(theSurface);
    SF_LEFT_SBD(newSurface)   = SF_LEFT_SBD(theSurface);
    /*Setzen der PolylineInformationen:*/
    SF_NMB_OF_POLYLINES(newSurface) = 0;
    SF_POLYLINES(newSurface) = NULL;

    thePolylineCycle = RS_PL_ZKLN(theRealSurface);
    if(thePolylineCycle == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","theRealSurface has no polylinecycle");
      return(1);
    }
    lauf_plz = 0;
    while(thePolylineCycle != NULL)
    {
      lauf_plz++;
      SF_NMB_OF_POLYLINES(newSurface) += PLZ_NMB_OF_POLYLINES(thePolylineCycle);

      /*Einbauen der Polylines von thePolylineCycle in newSurface->Polylines*/
      /* qwert */
      lauf_sfce_polyline = PLZ_POLYLINES(thePolylineCycle);
      anfang = lauf_sfce_polyline;
      numberofpolylinesofpolylinecycle = 0;
      while(lauf_sfce_polyline != NULL)
      {
        numberofpolylinesofpolylinecycle++;
        ende = lauf_sfce_polyline;
        lauf_sfce_polyline = SFPL_NEXT(lauf_sfce_polyline);

      }
      /*Probe*/
      if(numberofpolylinesofpolylinecycle != PLZ_NMB_OF_POLYLINES(thePolylineCycle))
      {
        PrintErrorMessage('E',"SplitSurface","A PLZ has too much or too less polylines");
        return (1);
      }

      merkeSfceplanfang = SF_POLYLINES(newSurface);
      SF_POLYLINES(newSurface) = anfang;
      SFPL_NEXT(ende) = merkeSfceplanfang;
      /*durch diese 3 Zeilen wurden alle Polylines des betrachteten Poylinezyklus laufvar_plz
         an die  Stelle newSurface->Polylines gehaengt.*/


      thePolylineCycle = PLZ_NEXT(thePolylineCycle);
    }
    if(lauf_plz != RS_NMB_OF_PL_ZKLN(theRealSurface))
    {
      PrintErrorMessage('E',"SplitSurface","theRealSurface has too much or too less polylinecycles");
      return (1);
    }

    /* SURFACEDETECTOR . . .  werden hier wieder mit Nullen initialisiert.*/
    SF_NMB_OF_POLYLI_ZYK(newSurface) = 0;
    SF_NMB_OF_REALSFCS(newSurface) =        0;
    SF_POLYLI_ZYK(newSurface) = NULL;
    SF_REALSFCS(newSurface) = NULL;
    /* . . . SURFACEDETECTOR */


    /* erzeuge dabei Dreiecke der neuen Surface mit GetMem TRIANGLE_TYP*/
    SF_TRIAS(newSurface) = NULL;
    SF_NMB_OF_TRIAS(newSurface) = 0;
    SF_NMB_OF_POINTS(newSurface) = 0;

    /*erstes Dreieck der ersten Line*/
    ersteline = PL_LINES_LINE(PL_LINES(SFPL_PL(SF_POLYLINES(newSurface))));
    if((erstesdreieck = Find_SFE_Triangle(ersteline,newSurface)) == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","did not get a triangle out of Find_SFE_Triangle");
      return (1);
    }

    nmb_of_triangles = 0;
    /*New_Triangle_List ist global und wird in der Funktion FetchAllTriangles weiter gefuellt*/
    if((New_Triangle_List = (TRIANGLE_TYP *)GetTmpMem(theHeap,sizeof(TRIANGLE_TYP),ANS_MarkKey))== NULL)
    {
      PrintErrorMessage('E',"SplitSurface","got  no memory  for  New_Triangle_List !?!");
      return(1);
    }
    TRIA_SFE_KN(New_Triangle_List) = erstesdreieck;
    TRIA_NEXT(New_Triangle_List) = NULL;
    nmb_of_triangles++;

    SFE_ORIENTATION_FLAG(erstesdreieck) = T;
    /*In der folg. Fkt werden die statics New_Triangle_List und nmb_of_triangles upgedatet*/
    if((rv = FetchAllTriangles(erstesdreieck)) != FERTIG)
    {
      PrintErrorMessage('E',"SplitSurface","ret-value of FetchAllTriangles was not FERTIG");
      return(1);
    }

    SF_TRIAS(newSurface) = New_Triangle_List;
    SF_NMB_OF_TRIAS(newSurface) = nmb_of_triangles;

    /*Update des Surfacenamens von newSurface und dessen Dreiecken - bei den Lines wird auf ein Update verzichtet*/
    new_sfc_name = SF_NAME1(newSurface) + (SPLIT_SURFACE_DISTINGUISH*lauf_rs);
    SF_NAME1(newSurface) = new_sfc_name;

    /*wenn hier eine Doppelsurface gesplittet wurde, so muss der erste und der zweite Identifier einen Splitoffset bzw.
       ein Update erhalten.*/
    if(SF_NAME2(newSurface) != 0.0)
    {
      new_sfc_name2 = SF_NAME2(newSurface) + (SPLIT_SURFACE_DISTINGUISH*lauf_rs);
      SF_NAME2(newSurface) = new_sfc_name2;

      /*Update der Tria->Surfacezugehoerigkeit im Falle einer Doppelsurface bzgl dem ersten und zweiten Identifier:*/
      lauf_sftria = SF_TRIAS(newSurface);
      while(lauf_sftria != NULL)
      {
        SFE_IDF_0(TRIA_SFE_KN(lauf_sftria)) = new_sfc_name;
        SFE_IDF_1(TRIA_SFE_KN(lauf_sftria)) = new_sfc_name2;
        lauf_sftria = TRIA_NEXT(lauf_sftria);
      }
    }
    else
    {
      /*Update der Tria->Surfacezugehoerigkeit im Falle einer einfachen Surface*/
      lauf_sftria = SF_TRIAS(newSurface);
      while(lauf_sftria != NULL)
      {
        SFE_IDF_0(TRIA_SFE_KN(lauf_sftria)) = new_sfc_name;
        lauf_sftria = TRIA_NEXT(lauf_sftria);
      }
    }


    dreiecksanzahl = dreiecksanzahl + nmb_of_triangles;

    /* NOT TODO :	SF_NMB_OF_POINTS(newSurface) WIRD in EvalNmbOfPointsOfSfcs ermittelt*/

    /*Update of statistical Domain Info es gibt eine Surface mehr*/
    NMB_OF_SFCES(DomainInfo_Pointer) = NMB_OF_SFCES(DomainInfo_Pointer) + 1;


    /*  das noch prgrammieren*/
    /*trage neue surfaces auch bei the_sbd ein*/
    /* the_sbd->NMBsurfaces ++*/
    if((neue_sfc_of_sbd = (SFC_TYP *)GetTmpMem(theHeap,sizeof(SFC_TYP),ANS_MarkKey))== NULL)
    {
      PrintErrorMessage('E',"SplitSurface","got  no SFC_TYP memory  for  neue_sfc_of_sbd !?!");
      return(1);
    }
    SFC_SURF(neue_sfc_of_sbd) = newSurface;
    SFC_NEXT(neue_sfc_of_sbd) = SD_SFCS(the_sbd);                    /*vorne einfuegen*/
    SD_SFCS(the_sbd) = neue_sfc_of_sbd;                     /*neue Wurzel*/
    SD_NMB_OF_SFCS(the_sbd) = SD_NMB_OF_SFCS(the_sbd) + 1;


    if(DoubleSurfaceCase == T)
    {
      /*  das noch prgrammieren*/
      /*trage neue surfaces auch bei the_sbd_Double ein*/
      /* the_sbd_Double->NMBsurfaces ++*/
      if((neue_sfc_of_sbd_Double = (SFC_TYP *)GetTmpMem(theHeap,sizeof(SFC_TYP),ANS_MarkKey)) == NULL)
      {
        PrintErrorMessage('E',"SplitSurface","got  no SFC_TYP memory  for  neue_sfc_of_sbd_Double !?!");
        return(1);
      }
      SFC_SURF(neue_sfc_of_sbd_Double) = newSurface;
      SFC_NEXT(neue_sfc_of_sbd_Double) = SD_SFCS(the_sbd_Double);                          /*vorne einfuegen*/
      SD_SFCS(the_sbd_Double) = neue_sfc_of_sbd_Double;                           /*neue Wurzel*/
      SD_NMB_OF_SFCS(the_sbd_Double) = SD_NMB_OF_SFCS(the_sbd_Double) + 1;
    }


    theRealSurface = RS_NEXT(theRealSurface);
  }      /* von laufe ueber die real surfaces und erzeuge neue Surfaces*/



  /*Dreiecksflags wieder zuruecksetzen.!!!
     alle Dreiecke der alten Surface muessen gesetzt sein
     wenn es noch ein unbesetztes gibt ==>  ERROR*/
  lauf_tria = SF_TRIAS(theSurface);
  while(lauf_tria != NULL)
  {
    if(SFE_ORIENTATION_FLAG(TRIA_SFE_KN(lauf_tria)) == F)             /*wenn dieses Dreieck nicht verwendet wurde*/
    {
      PrintErrorMessage('E',"SplitSurface","settting back flags : tria found, which was not used");
      return(1);
    }
    SFE_ORIENTATION_FLAG(TRIA_SFE_KN(lauf_tria)) = F;
    lauf_tria = TRIA_NEXT(lauf_tria);
  }



  /* Probe ANzahl der Dreieck = ANzahl der ehemealigen Gesamtdereeckke ?*/
  if(dreiecksanzahl != SF_NMB_OF_TRIAS(theSurface))
  {
    PrintErrorMessage('E',"SplitSurface","total number of trias of real surfaces does not match SF_NMB_OF_TRIAS of theSurface");
    return(1);
  }

  if(lauf_rs != SF_NMB_OF_REALSFCS(theSurface))
  {
    PrintErrorMessage('E',"SplitSurface","number of real surfaces does not match SF_NMB_OF_REALSFCS");
    return(1);
  }

  /*loesche die gesplittete Surface aus der liste von  the_sbd:*/
  SD_NMB_OF_SFCS(the_sbd) = SD_NMB_OF_SFCS(the_sbd) - 1;
  if(pred_sfc_OF_the_sbd == NULL)       /*d.h. die gesplittete Surface ist die allererste der Sbd->SurfaceListe*/
  {
    /* vorne sind aber bereits mind. 2 neue Surfaces in der Liste durch das Split
       oder aber schon viel mehr ...*/
    sf_of_sd_lauf = SD_SFCS(the_sbd);
    if(sf_of_sd_lauf == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","sf_of_sd_lauf became NULL but theSurface was not found yet");
      return(1);
    }
    /*laufe von der allerersten SUrface bis zur unmittelbar vor theSurface liegenden Surface*/
    while(SFC_SURF(SFC_NEXT(sf_of_sd_lauf)) != theSurface)
    {
      sf_of_sd_lauf = SFC_NEXT(sf_of_sd_lauf);
      if(sf_of_sd_lauf == NULL)
      {
        PrintErrorMessage('E',"SplitSurface","sf_of_sd_lauf became NULL but theSurface was not found yet");
        return(1);
      }
    }
    SFC_NEXT(sf_of_sd_lauf) = SFC_NEXT(SFC_NEXT(sf_of_sd_lauf));            /*loeschen der alten gesplitteten sf_lfv aus der Surfaceliste*/
  }
  else
  {
    SFC_NEXT(pred_sfc_OF_the_sbd) = SFC_NEXT(SFC_NEXT(pred_sfc_OF_the_sbd));
  }


  if(DoubleSurfaceCase == T)
  {
    /*loesche die gesplittete Surface aus der liste von  the_sbd_Double:*/
    SD_NMB_OF_SFCS(the_sbd_Double) = SD_NMB_OF_SFCS(the_sbd_Double) - 1;
    if(pred_sfc_OF_the_sbd_Double == NULL)             /*d.h. die gesplittete Surface ist die allererste der Sbd_Double->SurfaceListe*/
    {
      /* vorne sind aber bereits mind. 2 neue Surfaces in der Liste durch das Split
         oder aber schon viel mehr ...*/
      sf_of_sd_lauf_Double = SD_SFCS(the_sbd_Double);
      if(sf_of_sd_lauf_Double == NULL)
      {
        PrintErrorMessage('E',"SplitSurface","sf_of_sd_lauf_Double became NULL but theSurface was not found yet");
        return(1);
      }
      /*laufe von der allerersten Surface bis zur unmittelbar vor theSurface liegenden Surface*/
      while(SFC_SURF(SFC_NEXT(sf_of_sd_lauf_Double)) != theSurface)
      {
        sf_of_sd_lauf_Double = SFC_NEXT(sf_of_sd_lauf_Double);
        if(sf_of_sd_lauf_Double == NULL)
        {
          PrintErrorMessage('E',"SplitSurface","sf_of_sd_lauf_Double became NULL but theSurface was not found yet");
          return(1);
        }
      }
      SFC_NEXT(sf_of_sd_lauf_Double) = SFC_NEXT(SFC_NEXT(sf_of_sd_lauf_Double));                  /*loeschen der alten gesplitteten sf_lfv aus der Surfaceliste*/
    }
    else
    {
      SFC_NEXT(pred_sfc_OF_the_sbd_Double) = SFC_NEXT(SFC_NEXT(pred_sfc_OF_the_sbd_Double));
    }
  }



  /*loesche die gesplittete Surface aus der Gesamtliste :*/
  NMB_OF_SFCES(DomainInfo_Pointer) = NMB_OF_SFCES(DomainInfo_Pointer) - 1;
  if(thePredSurface == NULL)       /*d.h. die gesplittete Surface ist die allererste der alten SurfaceListe*/
  {
    /*vorne sind aber bereits mind 2 neue Surfaces in der Liste durch das Split eingetragen worden
       oder aber schon viel mehr ...*/
    sf_lfv2 = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
    if(sf_lfv2 == NULL)
    {
      PrintErrorMessage('E',"SplitSurface","sf_lfv2 became NULL but theSurface was not found yet");
      return(1);
    }
    /*laufe von der allerersten SUrface bis zur unmittelbar vor theSurface liegenden Surface*/
    while(SF_NEXT(sf_lfv2) != theSurface)
    {
      sf_lfv2 = SF_NEXT(sf_lfv2);
      if(sf_lfv2 == NULL)
      {
        PrintErrorMessage('E',"SplitSurface","sf_lfv2 became NULL but theSurface was not found yet");
        return(1);
      }
    }
    SF_NEXT(sf_lfv2) = SF_NEXT(theSurface);            /*loeschen der alten gesplitteten sf_lfv aus der Surfaceliste*/
  }
  else
  {
    SF_NEXT(thePredSurface) = SF_NEXT(theSurface);
  }


  return(0);
}



/****************************************************************************/
/**
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgmSurfaceDetecting()
{
  SF_TYP *sf_lfv, *pred_sf_lfv;
  INT rw, nmb_of_plzs_polylines, llf;
  PLZ_TYP *lauf_plz;

  pred_sf_lfv = NULL;
  sf_lfv = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  while(sf_lfv != NULL)
  {
    /*Doppelsurfaces werden nicht gesplittet, folglich ...*/
    /* FALSCH !! 7.10.97 */
    /*Auch Doppelsurfaces muesen u.U. gesplittet werden:
       siehe a2l_80.ans*/
    /* egal ob Surface nur einen oder aber zwei SfcIdentifier besitzt ...*/
    /*erzeuge die Polylinezyklen zu dieser Surface ...*/
    if((rw = Create_PLZN(sf_lfv)) == 1)
    {
      PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","did receive ERROR from Create_PLZN");
      return (1);
    }

    /*Ist genau  ein Polylinezyklus entstanden ?*/
    /*Das kann nicht sein, da im Falle eines Zyklus die Anzahl 0 vorliegt. siehe Create_PLZN*/
    if(SF_NMB_OF_POLYLI_ZYK(sf_lfv) == 1)
    {
      PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","did receive exactly 1 PLZ from Create_PLZN but sfce must have at least 2 PLZs or none of it");
      return (1);
    }

    /*Ist mehr als ein Polylinezyklus entstanden ?*/
    if(SF_NMB_OF_POLYLI_ZYK(sf_lfv) >1)
    {
      /*Probe1: Die Surface darf keine Polylines mehr uber dirketen Zugriff besitzen
               sondern muss alle ihre Polylines in den Polylinzyklen haben !*/
      if(SF_POLYLINES(sf_lfv) != NULL)
      {
        PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","Surface->Polylines ist not NULL after successfull Create_PLZN");
        return (1);
      }

      /*Probe 2: Die Polylines existieren aber noch in den Polylinezyklen ==>
                 Nun erfolgt ein quantitative Ueberpruefung, ob noch alle da sind.*/
      nmb_of_plzs_polylines = 0;
      lauf_plz = SF_POLYLI_ZYK(sf_lfv);
      if(lauf_plz == NULL)
      {
        PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","Surface should have Polylinecycle<s>");
        return (1);
      }
      nmb_of_plzs_polylines += PLZ_NMB_OF_POLYLINES(lauf_plz);
      for(llf = 2; llf<=SF_NMB_OF_POLYLI_ZYK(sf_lfv); llf++)
      {
        lauf_plz = PLZ_NEXT(lauf_plz);
        if(lauf_plz == NULL)
        {
          PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","Surface doesnt have all Polylinecycle<s>");
          return (1);
        }
        nmb_of_plzs_polylines += PLZ_NMB_OF_POLYLINES(lauf_plz);
      }
      if(nmb_of_plzs_polylines != SF_NMB_OF_POLYLINES(sf_lfv))
      {
        PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","Surface doesnt have as much Polylines as all its PLZs together");
        return (1);
      }

      /* GOON HERE */
      /*erzeuge die RealSurfaces aus den Polylinezyklen*/
      if((rw = Create_RealSurfaces(sf_lfv)) == 1)
      {
        PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","did receive ERROR from Create_RealSurfaces");
        return (1);
      }

      /*Ist mehr als eine RealSurface entstanden ?*/
      if(SF_NMB_OF_REALSFCS(sf_lfv) >1)
      {
        /*die wahren Surfaces koennen nun erzeugt werden*/
        if((rw = SplitSurface(sf_lfv,pred_sf_lfv)) == 1)
        {
          PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","did receive ERROR from SplitSurface");
          return (1);
        }
        /*in diesem Fall pred_sf_lfv belassen und nur sf_lfv fortschalten,
           da <Splitten> bedeutet, dass sf_lfv ra�sfaellt --> also muss man den
           Vorgaenger pred_sf_lfv als Predescessor belassen !*/
      }
      else
      {
        if (SF_NMB_OF_REALSFCS(sf_lfv) ==1)
        {
          PrintErrorMessage('E',"Ansys2lgmSurfaceDetecting","SF_NMB_OF_REALSFCS(sf_lfv) == 1 ->impossible");
          return (1);
        }
        /* else OK:  es wurden eben keine Realsurfaces erstellt : Surface-Polylines wurden in
           Create_RealSurfaces ugdedatet.*/
        pred_sf_lfv = sf_lfv;
      }
    }
    else
    {
      pred_sf_lfv = sf_lfv;
    }

    sf_lfv = SF_NEXT(sf_lfv);            /*wird im Gegensatz zu pred_sf_lfv stets fortgeschaltet
                                                                        pred_sf_lfv im Falle SplitSurface nicht !!!*/
  }      /*von while*/


  return(0);
}




/****************************************************************************/
/**
*  dreht die Reihenfolge der NodeIDs vom OFD "dasSFE" um:
*  Dazu werden erster mit zweitem Node vertauscht sowie
*  zweiter mit drittem Node vertauscht
*  So ist alles umgedreht (siehe auch theoretische Ueberlegungen im Konzept)
*
* @param dasSFE SFE, dessen Knotenreihenfolge korrigiert bzw. umgedreht werden muss
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT ChangeOrientation(SFE_KNOTEN_TYP *dasSFE)
{
  INT merk_id;
  SFE_KNOTEN_TYP *merk_SFEknoten;

  /*Vertauschung des ersten mit dem zweiten Knoten*/
  merk_id = SFE_NDID1(dasSFE);
  SFE_NDID1(dasSFE) = SFE_NDID2(dasSFE);
  SFE_NDID2(dasSFE) = merk_id;

  /*Vertauschung des zweiten mit dem dritten Nachbarn*/
  merk_SFEknoten = SFE_NGHB2(dasSFE);
  SFE_NGHB2(dasSFE) = SFE_NGHB3(dasSFE);
  SFE_NGHB3(dasSFE) = merk_SFEknoten;

  return(0);
}



/****************************************************************************/
/**
*  Passt die KnotenIDReihenfolge des NachbarSFEs Nachbar_SFE an die KnReihenfolge
*  von Muster_SFE an.
*
*  Die beiden Reihenfolgen muessen so beschaffen sein,
*  dass sie an dem gemeinsamen Streckenzug entgegengesetzte Richtung besitzen.
*
*  Die Funktion sucht zunaechst die beiden Knoten beim Nachbarelement
*  und prueft dann, ob sie dort in gleicher oder umgekehrter Reihenfolge auftreten.
*  Bei gleicher Reihenfolge ist es noetig diese zu korrigieren.
*  Dazu wird die Subfunktion "ChangeOrientation(..)" aufgerufen.
*
*  Fuer die Kantenreihenfolge eines Dreiecks mit den KnotenIDs i,j und k gilt
*  erste Kante von i nach j, zweite Kante von j nach k sowie dritte Kante von k nach i
*
* @param Muster_SFE  SFE, deren Knotenreihenfolge bereits stimmt
* @param Nachbar_SFE  SFE, deren Knotenreihenfolge an die von Muster_SFE angepasst werden soll
* @param kante  KantenID von Muster_SFE, bzgl. der die Nachbarschaft zu Nachbar_SFE besteht
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ausrichtung(SFE_KNOTEN_TYP *Muster_SFE, SFE_KNOTEN_TYP *Nachbar_SFE, INT kante)
{
  INT ID_firstnode, ID_secndnode;
  INT merki1, merki2, i;
  INT rw;

  switch(kante)
  {
  case 0 :         ID_firstnode = SFE_NDID1(Muster_SFE);
    ID_secndnode = SFE_NDID2(Muster_SFE);
    break;
  case 1 :         ID_firstnode = SFE_NDID2(Muster_SFE);
    ID_secndnode = SFE_NDID3(Muster_SFE);
    break;
  case 2 :         ID_firstnode = SFE_NDID3(Muster_SFE);
    ID_secndnode = SFE_NDID1(Muster_SFE);
    break;
  default :        PrintErrorMessage('E',"Ausrichtung","got wrong Input-Value: kante != {0|1|2}");
    return (1);
  }

  /*die beiden KnotenIDs beim Nachbarn suchen*/
  merki2 = -1;
  merki1 = -1;
  for(i=0; i<3; i++)
  {
    if(SFE_NDID(Nachbar_SFE,i) == ID_secndnode)
    {
      merki2 = i;
    }
    else if(SFE_NDID(Nachbar_SFE,i) == ID_firstnode)
    {
      merki1 = i;
    }
  }
  if((merki1 != -1) && (merki2 != -1))
  {
    /* Ist die NachfolgeKnotenID von ID_firstnode == ID_secndnode */
    /* bzw.: Ist die Reihenfolge wieder genau identisch ?*/
    if( ((merki1+1)%3) == merki2)
    {
      /*Die Reihenfolge muss aber immer andersherum sein ==> Aendern ...*/
      if( (rw = ChangeOrientation(Nachbar_SFE)) == 1 )
      {
        PrintErrorMessage('E',"Ausrichtung","got ERROR from calling ChangeOrientation");
        return (1);
      }
    }
  }
  else
  {
    PrintErrorMessage('E',"Ausrichtung","die beiden SFEs sind ja gar keine Nachbarn");
    return (1);
  }

  return(0);
}



/****************************************************************************/
/**
*  Diese rekursive Funktion durchlaeuft alle Dreiecke einer Surface
*  mit HIlfe der schon bestehenden Nachbarschaftsbeziehungen zwischen
*  den Dreiecken. Dabei wird nach jedem ueberprueften resp. ausgerich-
*  teten Dreieck die statische fuer die Funktion global verwendbare
*  Variable "zaehler" inkrementiert. Wenn dieser Zaehler den Wert der
*  GesamtDreiecksAnzahl der betrachteten Surface erreicht, so wird die
*  Rekursionshierarchie mit return(FERTIG) auf dem schnellsten Wege abgebaut.
*
* @param  Muster_SFE SFE, das als Muster verwendet wird, dessen
*                    Umlaufrichtung also bereits stimmt
*
* @return  FERTIG if alle Dreiecke einmal besucht bzw.
*          static INT zaehler == static INT nmb_of_trias_of_sf,
*          1 if error occured.
   D*/
/****************************************************************************/
INT TriangleIDOrientations(SFE_KNOTEN_TYP *Muster_SFE)
{
  SFE_KNOTEN_TYP *Nachbar_SFE;
  INT kante;       /*Laufvariable ueber die 3 Kanten eines Dreiecks*/
  INT neubesetzt[3];
  INT rv,rgbwrt;

  neubesetzt[0] = F; neubesetzt[1] = F; neubesetzt[2] = F;

  /* Laufe ueber die 3 Nachbarn von Muster_SFE */
  for(kante = 0; kante < 3; kante++)
  {
    Nachbar_SFE =  SFE_NGHB(Muster_SFE,kante);
    if(Nachbar_SFE != NULL)
    {
      /*natuerlich nur wenn dieser Nachbar nicht schon orientiert wurde ...*/
      if(SFE_ORIENTATION_FLAG(Nachbar_SFE) == F)
      {
        if((rv = Ausrichtung(Muster_SFE, Nachbar_SFE, kante)) == 1)
        {
          PrintErrorMessage('E',"TriangleIDOrientations"," Returnvalue of Ausrichtung was 1 ===> ERROR");
          return (1);
        }

        /*************************************************************/
        /* Dieser Nachbar ist nun richtig orientiert, folglich . . . */
        SFE_ORIENTATION_FLAG(Nachbar_SFE) = T;
        zaehler ++;
        neubesetzt[kante] = T;
        if(zaehler == nmb_of_trias_of_sf)
        {
          return(FERTIG);
        }
        /*************************************************************/
      }
    }
  }       /* von for */

  /*weiterer Lauf ueber die 3 Nachbarn, um die RekursionsHierarchie etwas kleiner zu halten*/
  for(kante = 0; kante < 3; kante++)
  {
    /*nur die soeben neu besetzten Dreiecke sind fuer einen Rekursionschritt interessant*/
    if(neubesetzt[kante] == T)
    {
      rgbwrt = TriangleIDOrientations(SFE_NGHB(Muster_SFE,kante));
      if(zaehler == nmb_of_trias_of_sf)
      {
        return(FERTIG);
      }
    }
  }       /* von for */

  return(FERTIG);

}



/****************************************************************************/
/**
*  erzeugt Dreiecksnachbarschaften :
*          d.h. jedes OFlaecDreieck einer Surface bzw jedem SFE werden die 3
*          Nachbar SFES ermittelt und zugewiesen.
*          Liegt ein solches ObFlDr am Rand einer Surface,
*          d.h es besitzt bzgl einer Kante kein Nachbarobfldreieck auf dieser
*          Surface, dann erfolgt ein "NULL"-Eintrag.
*  erzeugt ferner Dreiecksorientierungen :
*               d.h. die Reihenfolge der KnotenIds wird fuer alle Oberflaechendreiecke
*               so ausgerichtet, dass die "rechte Handregel" immer das selbe Ergebnis
*               liefert.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgmCreateTriaOrientations()
{
  SFE_KNOTEN_TYP *firstOFD;
  SF_TYP *sf_lfv;
  INT rv;

  /****************************************
     Pruefung der Dreiecksorientierungen
  ****************************************/
  /*laufe ueber alle Surfaces und erzeuge einheitliche Richtungen der
     Triangles (d.h. der TriangleNodeIDReihenfolgen */
  sf_lfv = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  while(sf_lfv != NULL)
  {
    /*Anzahl der Dreiecke von dieser Surface der statischen Variable nmb_of_trias_of_sf
       zuweisen : */

    nmb_of_trias_of_sf = SF_NMB_OF_TRIAS(sf_lfv);
    firstOFD = TRIA_SFE_KN(SF_TRIAS(sf_lfv));
    SFE_ORIENTATION_FLAG(firstOFD) = T;
    zaehler = 1;
    if((rv = TriangleIDOrientations(firstOFD)) != FERTIG)
    {
      PrintErrorMessage('E',"Ansys2lgmCreateTriaOrientations"," Returnvalue of TriangleIDOrientations was not FERTIG - Problems with checking ID-Orientations");
      return (1);
    }
    /*das erste der nmb_of_trias_of_sf Dreiecke wurde damit besucht -
       Deshalb wird jetzt zaehler auf 1 gesetzt
       Seine Richtung (Reihenfolge der KnotenIDs ist nun das Muster fuer alle anderen ...
       Jedes andere Dreieck muss nun einmal besucht werden und auf richtige Reihenfolge
       geprueft u. ggf. geaendert werden. Schluss ist wenn zaehler == nmb_of_trias_of_sf,
       da dann alle Dreiecke genau einmal gecheckt wurden. */


    sf_lfv = SF_NEXT(sf_lfv);
  }      /*von while*/
  /****************************************/
  return(0);
}



/****************************************************************************/
/**
*  laeuft ueber alle Surfaces und berechnet je Surface die Anzahl der Points,
*  dazu wird ein temporaeres Knotenfeld verwendet, bei dem fuer alle Knoten ein
*  Integer reserviert ist. Fuer jede Surface laeuft die Funktion ueber alle
*  zugehoerigen Triangles und markiert im temp. Feld die zugeh. Knoten.
*  Bei jeder neuen Markierung der NumberofPoints-Zaehler der Surface inkrementiert.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT EvalNmbOfPointsOfSfcs()
{
  SF_TYP *lauf_sf;
  TRIANGLE_TYP *lauf_tria;
  INT *TempNodeArray;
  INT index,i,n;

  lauf_sf = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);

  /*einmal fuer alle Surfaces ein temporaeres Hilfsfeld anlegen:*/
  /*erster Eintrag steht fuer keienn Node*/
  if ((TempNodeArray = (INT *)GetTmpMem(theHeap,EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer)*sizeof(INT),ANS_MarkKey))==NULL)
  {
    PrintErrorMessage('E',"EvalNmbOfPointsOfSfcs","  got no MEM for the TempNodeArray, see ansys2lgm.c");
    return(1);
  }


  /*laufe ueber alle Surfaces*/
  while(lauf_sf != NULL)
  {
    /*fuer jede Surface das temporaere Hilfsfeld mit F initialisieren*/
    for(i=0; i<EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer); i++)
    {
      TempNodeArray[i] = F;
    }

    /*laufe ueber alle Dreiecke*/
    lauf_tria = SF_TRIAS(lauf_sf);
    while(lauf_tria != NULL)
    {
      for(n=0; n<3; n++)
      {
        index = SFE_NDID(TRIA_SFE_KN(lauf_tria),n);
        if(TempNodeArray[index] == F)
        {
          TempNodeArray[index] = T;
          SF_NMB_OF_POINTS(lauf_sf) = SF_NMB_OF_POINTS(lauf_sf) + 1;
        }
      }
      lauf_tria = TRIA_NEXT(lauf_tria);
    }
    lauf_sf = SF_NEXT(lauf_sf);
  }
  return(0);
}



/****************************************************************************/
/**
*  berechnet mit Spat- und Kreuzprodukt, ob die RechteHandRegel mit i, j und k
*  nach "aussen" zeigt. Zur Feststellung dient v. Dieses "v" liegt bzgl der von ijk
*  definierten Ebene "innen"
*
* @param  i,j,k KnotenIDs in repraesentativer Reihenfolge fuer gesamte Surface
*               seit dem Durchlauf der FUnktion TriangleIDOrientation()
* @param  v  "v" wie vierte ID ==> Liegt dieser vierte Knoten auf der selben Seite,
*            wie die, auf die der Normalenvektor von i,j,k zeigt dann ist
*            "ijk"(genau in dieser Reihenfolge!) nicht nach aussen sondern nach
*            innen orientiert. Der Rueckgabewert ist also "F"
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT NachAussenOrientiert(INT i, INT j, INT k, INT v)
{
  DOUBLE II[3],JJ[3],KK[3],VV[3];
  DOUBLE A[3],B[3],C[3],V[3];
  DOUBLE cos_alpha,Laenge_C,Laenge_V;
  INT offs;
  INT imal3,jmal3,kmal3,vmal3;

  imal3 = i * 3;
  jmal3 = j * 3;
  kmal3 = k * 3;
  vmal3 = v * 3;

  /* II[3],JJ[3],KK[3],VV[3] fuellen mit cadconvert-Feldern und den Parametern i,j,k und v*/
  for(offs = 0; offs<3; offs++)
  {
    II[offs] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[imal3+offs];
    JJ[offs] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[jmal3+offs];
    KK[offs] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[kmal3+offs];
    VV[offs] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[vmal3+offs];
  }

  /*Berechnung des Vektors A von II nach JJ sowie  des Vektors B von JJ nach KK
     sowie  des Vektors V von II nach VV*/
  A[0] = JJ[0] - II[0]; A[1] = JJ[1] - II[1]; A[2] = JJ[2] - II[2];
  B[0] = KK[0] - JJ[0]; B[1] = KK[1] - JJ[1]; B[2] = KK[2] - JJ[2];
  V[0] = VV[0] - II[0]; V[1] = VV[1] - II[1]; V[2] = VV[2] - II[2];

  /*Berechnung des auf A,B senkrechten Vekotors C mit Hilfe des Kreuzprodukts*/
  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];

  /*Berechnung des Winkels zwischen C und V*/
  Laenge_C = sqrt( C[0]*C[0] + C[1]*C[1] + C[2]*C[2] );
  Laenge_V = sqrt( V[0]*V[0] + V[1]*V[1] + V[2]*V[2] );
  cos_alpha = ( ( C[0]*V[0] + C[1]*V[1] + C[2]*V[2] ) / (Laenge_C) / (Laenge_V) );

  if(cos_alpha > 0.0)
    return(F);
  else
    return(T);
}



/****************************************************************************/
/**
*  laeuft ueber alle Surfaces und weist diesen die SubdomainIds zu, die
*  "links" bzw. "rechts" liegen. Dazu verwendet es die vierteID,
*  die als Eingangsparameter kommt.
*  F�r die eigentliche Links-RechtsBerechnung wird "NachAussenOrientiert(...)"
*  aufgerufen.
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT EvalLeftRightOfSfcs()
{
  SF_TYP *lauf_sf;
  SFE_KNOTEN_TYP *sfe;
  INT i,j,k,v;

  lauf_sf = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);

  /*laufe ueber alle Surfaces*/
  while(lauf_sf != NULL)
  {
    /*man verwende den SFE-Pointer des ersten Triangles diese Surface lauf_sf*/
    sfe = TRIA_SFE_KN(SF_TRIAS(lauf_sf));

    /*die IDs der Eckknoten in einer Reihenfolge, die fuer die gesamte Surface
       repraesentativ ist.*/
    i = SFE_NDID1(sfe);
    j = SFE_NDID2(sfe);
    k = SFE_NDID3(sfe);
    v = SFE_4ND_0(sfe);            /*die vierte KnotenID bzgl. des ersten Identifiers*/

    /*Sind ijk so orientiert, dass die RechteHandregel aus der Subdomain raus zeigt,
       bzw. nicht (!)  in Richtung des vierten Knotens zeigt.*/
    if(NachAussenOrientiert(i,j,k,v) == T)
    {
      SF_LEFT_SBD(lauf_sf) = (int) (floor(SFE_IDF_1(sfe)));                  /*RechtHndReg mit ijk zeigt nach "links" ind das
                                                                                andere "linke-C.Tapp" Gebiet SFE_IDF_1(sfe)*/
      SF_RIGHT_SBD(lauf_sf) = (int) (floor(SFE_IDF_0(sfe)));                  /*LinkHndReg mit ijk zeigt nach "rechts" in das "(r)echte"
                                                                                                         Gebiet SFE_IDF_0(sfe) selbst*/
    }
    else
    {
      SF_RIGHT_SBD(lauf_sf) = (int) (floor(SFE_IDF_1(sfe)));
      SF_LEFT_SBD(lauf_sf) = (int) (floor(SFE_IDF_0(sfe)));
    }

    lauf_sf = SF_NEXT(lauf_sf);
  }
  return(0);
}


/****************************************************************************/
/**
* ruft die beiden Funktionen EvalNmbOfPointsOfSfc und EvalLeftRightOfSfcs auf
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgmEvalSurfaceInformations()
{
  INT rv;

  if ((rv = EvalNmbOfPointsOfSfcs()) != 0)
  {
    PrintErrorMessage('E',"Ansys2lgmEvalSurfaceInformations","ERR-Return-Val from EvalNmbOfPointsOfSfcs");
    return (1);
  }
  if ((rv = EvalLeftRightOfSfcs())!= 0)
  {
    PrintErrorMessage('E',"Ansys2lgmEvalSurfaceInformations","ERR-Return-Val from EvalLeftRightOfSfcs");
    return (1);
  }

  return(0);
}



char GetCharact(int input)
{
  switch(input)
  {
  case 0 : return('0');break;
  case 1 : return('1');break;
  case 2 : return('2');break;
  case 3 : return('3');break;
  case 4 : return('4');break;
  case 5 : return('5');break;
  case 6 : return('6');break;
  case 7 : return('7');break;
  case 8 : return('8');break;
  case 9 : return('9');
  }
  PrintErrorMessage('E',"GetCharact","no correct valu for var. input");
  return ('E');

}


/****************************************************************************/
/**
*  evalates an Output String with 19+1 chars ==> necessary for the user
*  of the CAD-Interface to distinguish between different LGM-Surfaces
*
* @param sfce_name1 Surface Identifer 1
* @param sfce_name1 Surface Identifer 2 extra number if double surface, else 0.0
* @param The_String reference Parameter, used for the evaluated String
* @param mft Pointer on flag : with or without first triangle
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT SurfaceNamer(double sfce_name1, double sfce_name2, char *The_String, int *mft)
{
  int hilf;


  The_String[19] = '\0';
  The_String[9] = '_';      /*Trenner zwischen erster und zweiter Surface Zahl*/
  The_String[4] = '.';      /*GleitkommazahlPunkt*/
  The_String[14] = '.';      /*GleitkommazahlPunkt*/

  /*erster Surfaceidentifier*/
  /*Vorkommastellen :*/
  /*Besetzung von The_String[0]*/
  hilf = (INT)( floor(sfce_name1)) / 1000;
  The_String[0] = GetCharact(hilf);

  /*Besetzung von The_String[1]*/
  sfce_name1 = sfce_name1 - (double)(hilf * 1000);
  hilf = (INT)(floor(sfce_name1)) / 100;
  The_String[1] = GetCharact(hilf);

  /*Besetzung von The_String[2]*/
  sfce_name1 = sfce_name1 - (double)(hilf * 100);
  hilf = (INT)(floor(sfce_name1)) / 10;
  The_String[2] = GetCharact(hilf);

  /*Besetzung von The_String[3]*/
  sfce_name1 = sfce_name1 - (double)(hilf * 10);
  hilf = (INT)(floor(sfce_name1));
  The_String[3] = GetCharact(hilf);

  /*Nachkommastellen :*/
  /*Besetzung von The_String[5]*/
  sfce_name1 = (floor(0.5 + 10000 *(sfce_name1 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name1));
  The_String[5] = GetCharact(hilf);

  /*Besetzung von The_String[6]*/
  sfce_name1 = (floor(0.5 + 10000 *(sfce_name1 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name1));
  The_String[6] = GetCharact(hilf);

  /*Besetzung von The_String[7]*/
  sfce_name1 = (floor(0.5 + 10000 *(sfce_name1 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name1));
  The_String[7] = GetCharact(hilf);

  /*Besetzung von The_String[8]*/
  sfce_name1 = (floor(0.5 + 10000 *(sfce_name1 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name1));
  if(hilf >0)
  {
    *mft = 1;
  }
  The_String[8] = GetCharact(hilf);


  /*zweiter Surfaceidentifier*/
  /*Vorkommastellen :*/
  /*Besetzung von The_String[10]*/
  hilf = (INT)(floor(sfce_name2)) / 1000;
  The_String[10] = GetCharact(hilf);

  /*Besetzung von The_String[11]*/
  sfce_name2 = sfce_name2 - (double)(hilf * 1000);
  hilf = (INT)(floor(sfce_name2)) / 100;
  The_String[11] = GetCharact(hilf);

  /*Besetzung von The_String[12]*/
  sfce_name2 = sfce_name2 - (double)(hilf * 100);
  hilf = (INT)(floor(sfce_name2)) / 10;
  The_String[12] = GetCharact(hilf);

  /*Besetzung von The_String[13]*/
  sfce_name2 = sfce_name2 - (double)(hilf * 10);
  hilf = (INT)(floor(sfce_name2));
  The_String[13] = GetCharact(hilf);

  /*Nachkommastellen :*/
  /*Besetzung von The_String[15]*/
  sfce_name2 = (floor(0.5 + 10000 *(sfce_name2 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name2));
  The_String[15] = GetCharact(hilf);

  /*Besetzung von The_String[16]*/
  sfce_name2 = (floor(0.5 + 10000 *(sfce_name2 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name2));
  The_String[16] = GetCharact(hilf);

  /*Besetzung von The_String[17]*/
  sfce_name2 = (floor(0.5 + 10000 *(sfce_name2 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name2));
  The_String[17] = GetCharact(hilf);

  /*Besetzung von The_String[18]*/
  sfce_name2 = (floor(0.5 + 10000 *(sfce_name2 - (double)hilf))/10000) * 10.0;
  hilf = (INT)(floor(sfce_name2));
  The_String[18] = GetCharact(hilf);


  return(0);
} /*end of SurfaceNamer*/




/****************************************************************************/
/**
*
* @param integer_zahl int to be converted int an string
* @param generated_string pointer to new string
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT int2string(int integer_zahl, char *generated_string )
{
  int i,j,k;
  char helpchar;

  if (integer_zahl > MAX_INT)
  {
    UserWrite("ERROR: in int2string: integer_zahl > MAX_INT");
    return (1);
  }
  if (integer_zahl < 0)
  {
    UserWrite("ERROR: in int2string: integer_zahl < 0");
    return (1);
  }


  i = 0;

  do
  {
    generated_string[i++] = integer_zahl % 10 + '0';
    integer_zahl = integer_zahl /10;
  }
  while(integer_zahl > 0);

  generated_string[i] = '\0';       /*Stringabschluss*/

  /*reverse :*/
  k = i-1 ;       /*Stelle des letzten Zeichens*/
  i = i/2;       /*Anzahl der Zeichen DIV 2 */
  for(j=0; j<i; j++)       /*bis zur MItte*/
  {
    helpchar = generated_string[j];
    generated_string[j] = generated_string[k];
    generated_string[k] = helpchar;
    k--;
  }



  return(0);
} /* of function int2string*/






/****************************************************************************/
/**
*  update the Subdomainidentifiers (=Names) to the UG numeration
*  beginning with 1, 2, 3, ... N
*  and updating of the surfaces wich reference these subdomains
*  with their left- and right-components
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgmUpdateSbdmIDs()
{
  SD_TYP *sd_pointer;
  SFC_TYP *scf_of_sbdm;
  int i, bisherige_ID,hlp,rv,mit_first_tria,a,b,c;
  SF_TYP *lauf_sfce;

  /*DIRKS NEU*/
  char TheString[20];

  int cad_id[2];
  /*ALT	int lf;
          int cad_id_gefunden[2];
          int anzahlgefundene;
     ALT*/

  FILE *stream;
  char *buffer;
  int necessary,rval;
  char generated_string[11];      /*2hoch31 = 2147483648 hat 10 Stellen*/

  /*Wieviel Speicher ist notwendig ?*/
  necessary = sizeof(char);

  /*Herleitung siehe File abschaetzung auf saentis unter Text*/
  necessary = 300 + 150 * NMB_OF_SBDMS(DomainInfo_Pointer) + 150 * NMB_OF_SFCES(DomainInfo_Pointer);


  if((buffer = (char *)GetTmpMem(theHeap,necessary * sizeof(char),ANS_MarkKey))== NULL)
  {
    PrintErrorMessage('E',"Ansys2lgmUpdateSbdmIDs","no memory obtained for buffer");
    return (1);
  }
  *buffer ='\0';       /*mit Leerstring beginnen*/


  UserWrite("\n");
  if((buffer = strcat(buffer,(char*)"\n")) == NULL) return (1);
  UserWrite("-------------------------------------------\n");
  if((buffer = strcat(buffer,(char*)"-------------------------------------------\n")) == NULL) return (1);

  UserWrite("CAD-Subdomain-ID =>LGM-Subdomain-ID => Identifierstring :\n");
  if((buffer = strcat(buffer,"CAD-Subdomain-ID =>LGM-Subdomain-ID => Identifierstring :\n")) == NULL) return (1);


  sd_pointer = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
  /*laufe ueber alle Subdomains*/
  /*merke Dir dabei die Beziehung alte Sbd ---> neue Sbd in bisherige_ID_array*/

  for(i=1; i<=NMB_OF_SBDMS(DomainInfo_Pointer); i++)
  {
    if(sd_pointer != NULL)
    {
      bisherige_ID = SD_NAME(sd_pointer);
      bisherige_ID_array[i] = bisherige_ID;
      SD_NAME(sd_pointer) = i;                   /*Zuweisung der neuen UG-verstaendlichen ID*/
      /*laufe ueber KomponentenIndexArray und suche die bisherige ID:*/
      hlp = 1;
      while((KomponentenIndexArray[hlp] != -1) && (KomponentenIndexArray[hlp] != bisherige_ID))
      {
        hlp++;
      }
      if(KomponentenIndexArray[hlp] == -1)
      {
        UserWriteF("  %d                 %d                       %s\n",bisherige_ID,i,&(KomponentenNamenArray[0]));
        if((buffer = strcat(buffer,(char*)"  ")) == NULL) return (1);
        if((rval = int2string(bisherige_ID,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                 ")) == NULL) return (1);
        if((rval = int2string(i,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                       ")) == NULL) return (1);
        if((buffer = strcat(buffer,&(KomponentenNamenArray[0]))) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"\n")) == NULL) return (1);
      }
      else
      {
        UserWriteF("  %d                 %d                       %s\n",bisherige_ID,i,&(KomponentenNamenArray[hlp*31]));
        if((buffer = strcat(buffer,(char*)"  ")) == NULL) return (1);
        if((rval = int2string(bisherige_ID,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                 ")) == NULL) return (1);
        if((rval = int2string(i,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                       ")) == NULL) return (1);
        if((buffer = strcat(buffer,&(KomponentenNamenArray[hlp*31]))) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"\n")) == NULL) return (1);
      }
    }
    else
    {
      UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: Subdoamin is missing !!");
      return (1);
    }
    scf_of_sbdm = SD_SFCS(sd_pointer);             /*die erste SubdomainSurfaceStruktur diese Sbd*/
    /*laufe ueber alle SubdomainSurfaceTeile */
    while(scf_of_sbdm != NULL)
    {
      if (SF_LEFT_SBD(SFC_SURF(scf_of_sbdm)) == bisherige_ID)
      {
        SF_LEFT_SBD(SFC_SURF(scf_of_sbdm)) = i*(-1);                        /*Zuweisung der neuen SbdID zunaechst mit neg. Vorzeichen*/
        /*negatives Vorzeichen, da sonst Surfaces entstehen koennen mit rechtsWert = linksWert*/
      }
      else if (SF_RIGHT_SBD(SFC_SURF(scf_of_sbdm)) == bisherige_ID)
      {
        SF_RIGHT_SBD(SFC_SURF(scf_of_sbdm)) = i*(-1);                        /*Zuweisung der neuen SbdID zunaechst mit neg. Vorzeichen*/
        /*negatives Vorzeichen, da sonst Surfaces entstehen koennen mit rechtsWert = linksWert*/
      }
      else
      {
        UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: wrong SurfaceLeftRightInformation");
        return (1);
      }

      scf_of_sbdm = SFC_NEXT(scf_of_sbdm);
    }
    /*weiter mit naechster Subdomain:*/
    sd_pointer = SD_NEXT(sd_pointer);
  }

  /*laufe ueber alle Surfaces und multipliziere die LinksRechts-Infosmit -1
     damit */
  lauf_sfce = EXCHNG_TYP2_ROOT_SFC (ExchangeVar_2_Pointer);
  /* DIRKS NEW : Ausgabe der Surfacinformationen . . .*/
  UserWrite("\n");
  UserWrite("\n");
  UserWrite("\n");
  if((buffer = strcat(buffer,(char*)"\n\n\n")) == NULL) return (1);
  UserWrite("-------------------------------------------\n");
  if((buffer = strcat(buffer,(char*)"-------------------------------------------\n")) == NULL) return (1);
  UserWrite("CAD-Surface-ID       =>   LGM-Surface-ID \n");
  if((buffer = strcat(buffer,(char*)"CAD-Surface-ID       =>   LGM-Surface-ID \n")) == NULL) return (1);
  /**********12345678901234567890*/
  for(i=0; i<NMB_OF_SFCES(DomainInfo_Pointer); i++)
  {
    if(lauf_sfce != NULL)
    {
      SF_RIGHT_SBD(lauf_sfce) = (-1) * SF_RIGHT_SBD(lauf_sfce);
      SF_LEFT_SBD(lauf_sfce) = (-1) * SF_LEFT_SBD(lauf_sfce);


      mit_first_tria = 0;
      if((rv = SurfaceNamer(SF_NAME1(lauf_sfce), SF_NAME2(lauf_sfce),TheString,&mit_first_tria)) == 1)
      {
        UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs : SurfaceNamer returns ERROR.");
        return (1);
      }

      if(mit_first_tria > 0)
      {
        a = SFE_NDID1(TRIA_SFE_KN(SF_TRIAS(lauf_sfce)));

        b = SFE_NDID2(TRIA_SFE_KN(SF_TRIAS(lauf_sfce)));
        c = SFE_NDID3(TRIA_SFE_KN(SF_TRIAS(lauf_sfce)));
        /*Achtung, a,b,c*/
        /* ALT ...	*/
        /*just for debugging*/
        /*Bestimmung des ersten Ausgabewertes der Fehlermeldung*/

        /*point_array durchlaufen und CAD-ID zu "i" suchen !!!!*/
        /*				lf = 1;
                                        cad_id_gefunden[0] = F;
                                        cad_id_gefunden[1] = F;
                                        cad_id_gefunden[2] = F;
                                        anzahlgefundene = 0;
         */				/*TODO u.U. kann die folgende do Schleife
                                          weggelassen werden ud durch Koordausgabe
                                          ersetzt werden*/
        /*				while( (lf <= nbofnds) && (anzahlgefundene < 3) )
                                        {
                                                if(point_array[lf] == a)
                                                {
                                                        cad_id[0] = lf;
                                                        cad_id_gefunden[0] = T;
                                                        anzahlgefundene++;
                                                }
                                                else if(point_array[lf] == b)
                                                {
                                                        cad_id[1] = lf;
                                                        cad_id_gefunden[1] = T;
                                                        anzahlgefundene++;
                                                }
                                                else if(point_array[lf] == c)
                                                {
                                                        cad_id[2] = lf;
                                                        cad_id_gefunden[2] = T;
                                                        anzahlgefundene++;
                                                }

                                                lf++;
                                        }

                                        if ((cad_id_gefunden[0] == F) || (cad_id_gefunden[1] == F) || (cad_id_gefunden[2] == F))
                                        {
                                                PrintErrorMessage('E',"ansys2lgm","Ansys2lgmUpdateSbdmIDs : cad_id_gefunden[i] == F nicht alle gefunden");
                                                return(1);
                                        }
         *//* ... ALT*/

        /* NEU ...*/
        cad_id[0] = point_array_UG_CAD[a];
        cad_id[1] = point_array_UG_CAD[b];
        cad_id[2] = point_array_UG_CAD[c];
        /* ... NEU */



        UserWriteF("%s                     %d          FirstTria:%d,%d,%d\n",&(TheString[0]),i,cad_id[0],cad_id[1],cad_id[2]);
        if((buffer = strcat(buffer,&(TheString[0]))) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                     ")) == NULL) return (1);
        if((rval = int2string(i,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"          FirstTria:")) == NULL) return (1);
        if((rval = int2string(cad_id[0],generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)",")) == NULL) return (1);
        if((rval = int2string(cad_id[1],generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)",")) == NULL) return (1);
        if((rval = int2string(cad_id[2],generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"\n")) == NULL) return (1);
      }
      else
      {
        UserWriteF("%s                     %d          \n",&(TheString[0]),i);
        if((buffer = strcat(buffer,&(TheString[0]))) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"                     ")) == NULL) return (1);
        if((rval = int2string(i,generated_string)) == 1)
        {
          UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: int2string failed");
          return (1);
        }
        if((buffer = strcat(buffer,generated_string)) == NULL) return (1);
        if((buffer = strcat(buffer,(char*)"          \n")) == NULL) return (1);
      }
      /**********12345678901234567890*/
      /* . . .DIRKS NEW : Ausgabe der Surfacinformationen */

    }
    else
    {
      UserWrite("ERROR: in Ansys2lgmUpdateSbdmIDs: Surface is missing !!");
      return (1);
    }
    lauf_sfce = SF_NEXT(lauf_sfce);
  }

  stream = LGM_WriteOpenFile("sid.geo");
  if(stream == NULL)
  {
    UserWrite("ERROR: in LGM_LoadDomain: could not open stream ...");
    return (1);
  }

  if (fprintf(stream,"%s\n",(char*)buffer)<0) return (1);
  /* close file */
  if (fclose(stream)==EOF) return (1);



  return(0);
}




/****************************************************************************/
/**
*\brief converts an ANSYS output file to an ug-understandable Domain.
*
* @param  ExchangeVar_1 beinmhaltet die folgenden 4 Komponenten
* @param nmb_of_SFEs Anzahl der Oberfl�chendreiecke, die im ANSYS-FIle mit
*                    dem Schl�sselwort "SFE" eingeleitet werden.
* @param  SFE_Array Pointer to an array which posesses components of the data structure
*                           CAD_SFE_TYP (the whole SFE-information out of the ANSYS-File)
* @param  n_koord_array Pointer to an array created () in cadconvert
*                               which includes the real coordiantes, Indizes sind geshiftete UG-Indizes
*                               beinhaltet alle Punkte, also auch die inneren !
* @param  ExchangeVar_2 beinmhaltet die folgenden fuer sich sprechenden 5 Komponenten
* @param  rootsfc Pointer auf erste Surface
* @param  rootsd Pointer auf erste Subdomain
* @param  rootpl Pointer auf erste Polyline
* @param  SFE_HashTable Pointer auf die SFE-Hashtabelle
* @param  LI_HashTable Pointer auf die LI-Hashtabelle
*
* @return value 0:ok, 1:error.
   */
/****************************************************************************/
INT Ansys2lgm  ()
{
  INT rv;



  if((rv = Ansys2lgmInit()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmInit returns ERROR.");
    return (1);
  }

  if((rv = Ansys2lgmCreateHashTables()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmCreateHashTables returns ERROR.");
    return (1);
  }

  if((rv = Ansys2lgmCreateSbdsSfcsTriaRelations()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmCreateSbdsSfcsTriaRelations returns ERROR.");
    return (1);
  }

  if((rv = Ansys2lgmCreatePloylines()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmCreatePloylines returns ERROR.");
    return (1);
  }

  if(USE_SFC_DETECTOR == 1)
  {
    if((rv = Ansys2lgmSurfaceDetecting()) == 1)
    {
      UserWrite("ERROR: in Ansys2lgm : Ansys2lgmSurfaceDetecting returns ERROR.");
      return (1);
    }
  }

  if((rv = Ansys2lgmCreateTriaOrientations()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmCreateTriaOrientations returns ERROR.");
    return (1);
  }

  if((rv = Ansys2lgmEvalSurfaceInformations()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmEvalSurfaceInformations returns ERROR.");
    return (1);
  }


  if((rv = Ansys2lgmUpdateSbdmIDs()) == 1)
  {
    UserWrite("ERROR: in Ansys2lgm : Ansys2lgmUpdateSbdmIDs returns ERROR.");
    return (1);
  }

  /*
     return(rv = Ansys2lgmUse());
     return(rv = Ansys2lgmDispose());
   */

  return(0);
}




/****************************************************************************/
/*
   LGM_ANSYS_ReadDomain - reads general domain information from file

   SYNOPSIS:
        int LGM_ANSYS_ReadDomain (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info)

   PARAMETERS:
   .  Heap - pointer to Heap
   .  filename - name of ANSYS-file to read
   .  domain_info - information

   DESCRIPTION:
   function reads an ANSYS-file and sets the domain informations

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
/* out of lgm_transfer3d.c ...
   static int nSubdomain, nSurface, nLine, nPoint;
   static fpos_t filepos;
   static fpos_t fileposline;
   static fpos_t filepossurface;*/

int NS_DIM_PREFIX LGM_ANSYS_ReadDomain (HEAP *Heap, char *filename, LGM_DOMAIN_INFO *domain_info, INT MarkKey)
{
  INT ret_val;
  char helpstring[50];
  INT i;

  /*Default ZoomFactorZ in z-Direction in Promille; */
  ZoomFactorX = 1.0;      /*so gibt es keine Aenderung*/
  ZoomFactorY = 1.0;      /*so gibt es keine Aenderung*/
  ZoomFactorZ = 1.0;      /*so gibt es keine Aenderung*/

  TmpMemArray = NULL;
  ExchangeVar_2_Pointer = &ExchangeVar_2;
  ExchangeVar_1_Pointer = &ExchangeVar_1;
  DomainInfo_Pointer   = &DomainInfo;

  /* store heapptr and MarkKey */
  if (Heap==NULL) return (1);
  theHeap = Heap;       /*static heap pointer setzen*/
  ANS_MarkKey = MarkKey;

  /*OLDif((ret_val =ReadAnsysFile(filename)) == 1)*/

  if((ret_val =ReadAnsysFile(filename)) == 1)
  {
    UserWrite("ERROR: in LGM_ANSYS_ReadDomain ReadAnsysFile returns ERROR.");
    return (1);
  }

  if((ret_val = Ansys2lgm()) == 1)
  {
    UserWrite("ERROR: in LGM_ANSYS_ReadDomain Ansys2lgm returns ERROR.");
    return (1);
  }



  /***************************************************************************/
  /*Setzen der Domaininfowerten analog zu LGM_ReadDomain in lgm_transfer3d.c*/
  i=0;
  while(filename[i] != '.' )
  {
    helpstring[i] = filename[i];
    i++;
  }

  helpstring[i] = '.';
  i++;
  helpstring[i] = 'l';
  i++;
  helpstring[i] = 'g';
  i++;
  helpstring[i] = 'm';
  i++;
  helpstring[i] = '\0';
  strcpy(domain_info->Name,filename);


  if (ProblemName[0] == '\0')
  {
    UserWrite("Warning: in LGM_ANSYS_ReadDomain no problemname defined in ANSYS-File\n");
    /*UserWrite("Warning: using _Configurable_ as default value\n");
       strcpy(domain_info->ProblemName,"_Configurable_"); */
    UserWrite("Warning: using elder_problem as default value\n");
    strcpy(domain_info->ProblemName,"elder_problem");
  }
  else
  {
    strcpy(domain_info->ProblemName,ProblemName);
  }




  domain_info->Dimension = 3;

  /*domain_info->Convex = 0; temoraer abgeaendert 200198*/
  domain_info->Convex = 1;

  domain_info->nSubDomain = NMB_OF_SBDMS(DomainInfo_Pointer);

  domain_info->nSurface = NMB_OF_SFCES(DomainInfo_Pointer);

  domain_info->nPolyline = NMB_OF_PLINES(DomainInfo_Pointer);

  domain_info->nPoint = NMB_OF_POINTS(DomainInfo_Pointer);


  return(0);
}


/****************************************************************************/
/*
   LGM_ANSYS_ReadSizes - reads general domain information from file

   SYNOPSIS:
   int LGM_ANSYS_ReadSizes (LGM_SIZES *lgm_sizes)

   PARAMETERS:
   .  lgm_sizes -

   DESCRIPTION:
   Zuweisung der Anzahl der Points je Line, Anzahl der Surfaces je Subdomain,
   sowie der Anzahl der Polylines, Points und Trianmgles je Surface,


   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int NS_DIM_PREFIX LGM_ANSYS_ReadSizes (LGM_SIZES *lgm_sizes)
{
  int i;
  PL_TYP *plyln;
  SD_TYP *sbdmn;
  SF_TYP *srfce;

  /*Anzahl der Points je Line : */
  plyln = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  for (i=0; i < NMB_OF_PLINES(DomainInfo_Pointer); i++)
  {
    if(plyln != NULL)
    {
      lgm_sizes->Polyline_nPoint[i] = PL_NMB_OF_POINTS(plyln);
      plyln = PL_NXT(plyln);
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSizes: Polyline is missing !!");
      return (1);
    }
  }


  /*Anzahl der Surfaces je Subdomain zuweisen*/
  sbdmn = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
  for (i=1; i <= NMB_OF_SBDMS(DomainInfo_Pointer); i++)
  {
    if(sbdmn != NULL)
    {
      lgm_sizes->Subdom_nSurf[i] = SD_NMB_OF_SFCS(sbdmn);
      sbdmn = SD_NEXT(sbdmn);
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSizes: Subdomain is missing !!");
      return (1);
    }
  }


  /*Anzahl der Polylines, Points und Trianmgles je Surface zuweisen*/
  srfce = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  for (i=0; i < NMB_OF_SFCES(DomainInfo_Pointer); i++)
  {
    if(srfce != NULL)
    {
      lgm_sizes->Surf_nPoint[i] = SF_NMB_OF_POINTS(srfce);
      lgm_sizes->Surf_nPolyline[i] = SF_NMB_OF_POLYLINES(srfce);
      lgm_sizes->Surf_nTriangle[i] = SF_NMB_OF_TRIAS(srfce);
      srfce = SF_NEXT(srfce);
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSizes: Surface is missing !!");
      return (1);
    }
  }

  return(0);
}


/****************************************************************************/
/*
   LGM_ANSYS_ReadSubDomain - reads general domain information from file

   SYNOPSIS:
        int LGM_ANSYS_ReadSubDomain (int subdom_i, LGM_SUBDOMAIN_INFO *subdom_info)

   PARAMETERS:
   .  subdom_i - die i-te Subdomain (1,2,3,....N).
   .  subdom_info - to fill .

   DESCRIPTION:
   ermittelt alle Surfaces, die zur Subdomain subdom_i gehoeren und traegt deren
   ID (0,...,NMBOFSFCESOFDOMAIN-1) in aufsteigender Reihenfolge in subdom_info ein.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int NS_DIM_PREFIX LGM_ANSYS_ReadSubDomain (int subdom_i, LGM_SUBDOMAIN_INFO *subdom_info)
{
  int s,n,neue_ID,bisherige_ID,i,hlp;
  SF_TYP *sfce;
  SD_TYP *sd_pointer;

  n = 0;



  /*Zuweisung des Subdomainnamens :*/
  /*laufe zur angegebenen Subdomain :*/
  sd_pointer = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
  if(sd_pointer == NULL)
  {
    UserWrite("ERROR: in LGM_ANSYS_ReadSubDomain: Subdoamin is missing !!");
    return (1);
  }

  for(i=1; i<subdom_i; i++)
  {
    sd_pointer = SD_NEXT(sd_pointer);
    if(sd_pointer == NULL)
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSubDomain: Subdoamin is missing !!");
      return (1);
    }
  }

  neue_ID = SD_NAME(sd_pointer);
  if ( (neue_ID >= MAX_NUB_OF_SBDMS) || (neue_ID<=0) )
  {
    UserWrite("ERROR: in LGM_ANSYS_ReadSubDomain: neue_ID is too big or too small");
    return (1);
  }
  bisherige_ID = bisherige_ID_array[neue_ID];
  if(bisherige_ID <= 0)
  {
    UserWrite("ERROR: in LGM_ANSYS_ReadSubDomain: bisherige_ID is <= 0 !!");
    return (1);
  }
  /*laufe ueber KomponentenIndexArray und suche die bisherige ID:*/
  hlp = 1;
  while((KomponentenIndexArray[hlp] != -1) && (KomponentenIndexArray[hlp] != bisherige_ID))
  {
    hlp++;
  }
  if(KomponentenIndexArray[hlp] == -1)
  {
    strcpy(subdom_info->Unit,&(KomponentenNamenArray[0]));             /*no name given for this Subdomain in ANSYS-FIle !!!*/
  }
  else
  {
    strcpy(subdom_info->Unit,&(KomponentenNamenArray[hlp*31]));
  }




  sfce = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  /*laufe ueber alle Surfaces*/
  for(s=0; s<NMB_OF_SFCES(DomainInfo_Pointer); s++)
  {
    if(sfce != NULL)
    {
      if(SF_RIGHT_SBD(sfce) == subdom_i)
      {
        subdom_info->SurfaceNumber[n] = s;
        n = n+1;
      }
      else if (SF_LEFT_SBD(sfce) == subdom_i)
      {
        subdom_info->SurfaceNumber[n] = s;
        n = n+1;
      }
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSubDomain: Surface is missing !!");
      return (1);
    }
    sfce = SF_NEXT(sfce);
  }


  return(0);
}


/****************************************************************************/
/*
   LGM_ANSYS_ReadSurface -

   SYNOPSIS:
        int LGM_ANSYS_ReadSurface (int sfcnumber, LGM_SURFACE_INFO *surface_info)

   PARAMETERS:
   .  sfcnumber - Number of the Surface.
   .  surface_info - bla.

   DESCRIPTION:
   bla bla

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int NS_DIM_PREFIX LGM_ANSYS_ReadSurface (int sfcnumber, LGM_SURFACE_INFO *surface_info)
{
  SF_TYP *sfce;
  TRIANGLE_TYP *erstes_Tria;
  SFE_KNOTEN_TYP *zugehoeriger_SFE_Knoten, *sfe_nachbar;
  SFPL_TYP *plli_sfce;
  PL_TYP *plli, *plli2;


  int s,p,tp,n,pl,t,found,hilfesvaria;

  /*die erste Surface :*/
  sfce = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);

  /*gehe zur Surface mit der ID sfcnumber */
  for(s=0; s<sfcnumber; s++)
  {
    if(sfce == NULL)
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSurface: Surface is missing !!");
      return (1);
    }
    else
    {
      sfce = SF_NEXT(sfce);
    }
  }

  /*Eintrag der Left-/Right-Information*/
  surface_info->left = SF_LEFT_SBD(sfce);
  surface_info->right = SF_RIGHT_SBD(sfce);

  /*alle Points der Surface => dazu Verwendung eines temporaeren Feldes:*/
  /*durch die folgende Abfrag wird nur einmal Speicher geholt und nicht fuer jede Surface nochmal ganz neu!*/
  if (TmpMemArray == NULL)       /*wenn noch NULL*/
  {
    /*... dann hole Speicher ...*/
    if((TmpMemArray = (char *)GetTmpMem(theHeap,EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer) * sizeof(char),ANS_MarkKey))== NULL)
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSurface: no memory obtained for TmpMemArray");
      return (1);
    }
  }

  for(p=0; p<EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer); p++)
  {
    TmpMemArray[p] = 0;
  }


  /* NEU : NEU : NEU: ...*/
  /*laufe ueber alle Dreieck der Surface*/
  erstes_Tria = SF_TRIAS(sfce);

  for(t=0; t<SF_NMB_OF_TRIAS(sfce); t++)
  {
    if(erstes_Tria != NULL)
    {
      zugehoeriger_SFE_Knoten = TRIA_SFE_KN(erstes_Tria);
      /* NEU : NEU : NEU: Erweiterung der SFE_Knotentyp_Struktur*/
      /*trage beim zugehoerigen SFE-Knoten die Surfacelokale TriangleID ein*/
      SFE_LOCAL_TRI_ID(zugehoeriger_SFE_Knoten) = t;
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSurface: Triangle is missing !!");
      return (1);
    }


    /*weiterschalten zum naechsten Dreieck:*/
    erstes_Tria = TRIA_NEXT(erstes_Tria);
  }
  /*  ... NEU : NEU : NEU:*/



  /*laufe ueber alle Dreieck der Surface*/
  erstes_Tria = SF_TRIAS(sfce);

  for(t=0; t<SF_NMB_OF_TRIAS(sfce); t++)
  {
    if(erstes_Tria != NULL)
    {
      zugehoeriger_SFE_Knoten = TRIA_SFE_KN(erstes_Tria);


      /*die 3 NodeIDs des Trias im Pointarray markieren*/
      for(tp=0; tp<3; tp++)
      {
        TmpMemArray[SFE_NDID(zugehoeriger_SFE_Knoten,tp)] = 1;

        /*bei dieser Gelegenheit : Trianglinfos zuweisen:*/
        surface_info->Triangle[t].corner[tp] = SFE_NDID(zugehoeriger_SFE_Knoten,tp);

        /*Achtung hier geht es um die TriangleID !!*/
        /*OLDsurface_info->Triangle[t].neighbor[tp] = SFE_LOCAL_TRI_ID(SFE_NGHB(zugehoeriger_SFE_Knoten,tp));*/
        /*Achtung richtige Nummerierung - Der i-te Nachbar bezieht sich auf die dem iten Knoten
           gegenueberliegende Kante*/
        hilfesvaria = (tp+1)%3;
        sfe_nachbar = (SFE_NGHB(zugehoeriger_SFE_Knoten,hilfesvaria));
        if(sfe_nachbar == NULL)
        {
          surface_info->Triangle[t].neighbor[tp] = -1;
        }
        else
        {
          surface_info->Triangle[t].neighbor[tp] = SFE_LOCAL_TRI_ID(sfe_nachbar);
        }
      }
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSurface: Triangle is missing !!");
      return (1);
    }


    /*weiterschalten zum naechsten Dreieck:*/
    erstes_Tria = TRIA_NEXT(erstes_Tria);
  }
  /*laufe ueber das TmpMemArray und trage die Indizes der markierten Stellen
     in das PointArray ein.*/
  n = 0;
  for(p=0; p<EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer); p++)
  {
    if(TmpMemArray[p] == 1)
    {
      surface_info->point[n++] = p;                  /*alle PointIDs - Reihenfolge stimmt hier sogar,
                                                        was gar nicht noetig waere.*/
    }
  }


  /*Einlesen der  PolylineIDs*/
  /*laufe ueber alle polylines der Domain: ...*/
  plli = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  n=0;
  for(pl = 0; pl<NMB_OF_PLINES(DomainInfo_Pointer); pl++)
  {
    if(plli != NULL)
    {
      /*Hat die betrachtete Surface diese Polyline ebenso ?*/

      plli_sfce= SF_POLYLINES(sfce);
      found = 0;
      /*durchsuche die Polylines der betrachteten Surface*/
      while((plli_sfce != NULL)&&(found == 0))
      {
        plli2 = SFPL_PL(plli_sfce);
        if(plli2 == plli)                         /*Adressengleichheit genuegt*/
        {
          found = 1;
          surface_info->line[n++] = pl;
        }
        /*Bug 26: plli2 = PL_NXT(plli2);*/
        plli_sfce = SFPL_NEXT(plli_sfce);
      }
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadSurface: Polyline is missing !!");
      return (1);
    }
    plli = PL_NXT(plli);
  }

  return(0);
}



/****************************************************************************/
/*
   LGM_ANSYS_ReadLines - reads general domain information from file

   SYNOPSIS:
   int LGM_ANSYS_ReadLines (int dummy, LGM_LINE_INFO *line_info)

   PARAMETERS:
   .  plline - diese Variable zeigt an welche Polyline gemeint ist.
   .  line_info - wird gefuellt mit den IDs der zugehoerigen Points
   in sequentielle Reihenfolge

   DESCRIPTION:
   blabla

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int NS_DIM_PREFIX LGM_ANSYS_ReadLines (int which_plline, LGM_LINE_INFO *line_info)
{
  PL_TYP *plyln;
  PL_LINE_TYP *pllyln;
  int i;


  plyln = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  /*bestimme gesuchte Polyline ...*/
  for(i=0; i<which_plline; i++)
  {
    if(plyln != NULL)
    {
      plyln = PL_NXT(plyln);
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadLines: Polyline is missing !!");
      return (1);
    }
  }
  pllyln = PL_LINES(plyln);
  /*die erste KnotenID der Polyline ...*/
  line_info->point[0] = LI_NDID1(PL_LINES_LINE(pllyln));
  /*weise KnotenIDs der gefundenen Polyline zu*/
  /*... laufe dazu ueber die Polylineteilstreckenzuege ...*/
  /*Der Fall der zyklischen Polyline ist von C.Tapp bisher nicht beruecksichtigt
     Ich uebergeb in diesm speziellen Fall einen Punkt zweimal
     d.h. am Anfang und am Ende,; ausserdem ist PL_NMB_OF_POINTS dann um 1 zu gross
     in Zukunft: sollte lgm-Domain so ueberarbeitet werden,
     dass die Pl-Struktur hiefuer ein Flag besitzt.*/
  /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
  for(i=1; i<PL_NMB_OF_POINTS(plyln); i++)
  {
    if(pllyln != NULL)
    {
      line_info->point[i] = LI_NDID2(PL_LINES_LINE(pllyln));
      /*jetzt immer den zweiten Knoten nehmen, da ja die
         Anzahl der Points einer Polyline um 1 groesser ist
         als die Anzahl der PL-LineSegmente derselbigen ...*/
      pllyln = PL_LINES_NXT(pllyln);
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadLines: PolylineLine is missing !!");
      return (1);
    }
  }

  return(0);
}



/****************************************************************************/
/*
   LGM_ANSYS_ReadPoints - reads general domain information from file

   SYNOPSIS:
   int LGM_ANSYS_ReadPoints (LGM_POINT_INFO *lgm_point_info)

   PARAMETERS:
   .  lgm_point_info - coordinates.

   DESCRIPTION:
   lreads all coordinates

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int NS_DIM_PREFIX LGM_ANSYS_ReadPoints (LGM_POINT_INFO *lgm_point_info)
{
  int p, hilf;

  /*laufe ueber alle Points <d.h. auch die inneren !>*/

  for(p=0; p<EXCHNG_TYP1_NMB_OF_BNDNDS(ExchangeVar_1_Pointer); p++)
  {
    hilf = p*3;

    if((ZoomFactorX != 1.0) || (ZoomFactorY != 1.0) || (ZoomFactorZ != 1.0))
    {
      lgm_point_info[p].position[0] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf] * (ZoomFactorX);
      hilf++;
      lgm_point_info[p].position[1] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf] * (ZoomFactorY);
      hilf++;
      /*lgm_point_info[p].position[2] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf] * ((float)ZoomFactorZ);*/
      lgm_point_info[p].position[2] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf] * (ZoomFactorZ);
    }
    else
    {
      lgm_point_info[p].position[0] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf];
      hilf++;
      lgm_point_info[p].position[1] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf];
      hilf++;
      lgm_point_info[p].position[2] = (EXCHNG_TYP1_KOORDS(ExchangeVar_1_Pointer))[hilf];
    }
  }



  return(0);
}





/****************************************************************************/
/**
*  reads coordinformations for boundary points as well as
*  innerpoint out of the array n_koord_array_UG
*
* @param  theMesh Pointer of the meshstructure to be filled
*
* @return value 1:error.
   */
/****************************************************************************/
int FillPositionInformations(LGM_MESH_INFO *theMesh)
{
  int bndpindex, innpindex, h;

  /*Anzahl der BoundaryPoints*/
  theMesh->nBndP = statistik[1];

  /*Anzahl der InnerPoints*/
  theMesh->nInnP = statistik[0];

  /*Feld anlegen fuer theMesh->BndPosition == Feld fuer Zeiger auf BndPointKoordinatenpaare-tripel*/
  if ((theMesh->BndPosition = (double **)GetTmpMem(theHeap,(statistik[1])*sizeof(DOUBLE*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillPositionInformations"," ERROR: No memory for theMesh->BndPosition");
    return(1);
  }

  h=0;      /*Hilfsvariable fuer Zugriff auf n_koord_array_UG --> wird bei den InnerPoints einfach
               ohne Neuinitialisierung weiterverwendet*/

  /* Uebergabe der Koordinatenwerte der BoundaryPoints mit den IDs 0,1,2,...,m */
  for(bndpindex=0; bndpindex<statistik[1]; bndpindex++)
  {
    if (((theMesh->BndPosition)[bndpindex]= (double *)GetTmpMem(theHeap,(DIMENSION)*sizeof(DOUBLE), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillPositionInformations"," ERROR: No memory for (theMesh->BndPosition)[bndpindex]");
      return(1);
    }
    if((ZoomFactorX != 1.0) || (ZoomFactorY != 1.0) || (ZoomFactorZ != 1.0))
    {
      ((theMesh->BndPosition)[bndpindex])[0] = n_koord_array_UG[h] * (ZoomFactorX) ;h++;
      ((theMesh->BndPosition)[bndpindex])[1] = n_koord_array_UG[h] * (ZoomFactorY) ;h++;
      /*((theMesh->BndPosition)[bndpindex])[2] = n_koord_array_UG[h] - ((float)ZoomFactorZ/1000.0) * n_koord_array_UG[h];h++; */
      ((theMesh->BndPosition)[bndpindex])[2] = n_koord_array_UG[h] * (ZoomFactorZ);h++;
    }
    else
    {
      ((theMesh->BndPosition)[bndpindex])[0] = n_koord_array_UG[h];h++;
      ((theMesh->BndPosition)[bndpindex])[1] = n_koord_array_UG[h];h++;
      ((theMesh->BndPosition)[bndpindex])[2] = n_koord_array_UG[h];h++;
    }
  }

  /*Feld anlegen fuer theMesh->InnPosition == Feld fuer Zeiger auf InnerPointKoordinatenpaare-tripel*/
  /* wenn ueberhaupt innere Knoten existieren ... */
  if(statistik[0] >0)
  {
    if ((theMesh->InnPosition = (double **)GetTmpMem(theHeap,(statistik[0])*sizeof(double*), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillPositionInformations"," ERROR: No memory for theMesh->InnPosition");
      return(1);
    }
  }

  /* Uebergabe der Koordinatenwerte der InnerPoints mit den IDs m,m+1,m+2,...n */
  for(innpindex=0; innpindex<statistik[0]; innpindex++)
  {
    if (((theMesh->InnPosition)[innpindex]= (double *)GetTmpMem(theHeap,(DIMENSION)*sizeof(double), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillPositionInformations"," ERROR: No memory for (theMesh->InnPosition)[innpindex]");
      return(1);
    }
    if((ZoomFactorX != 1.0) || (ZoomFactorY != 1.0) || (ZoomFactorZ != 1.0))
    {
      ((theMesh->InnPosition)[innpindex])[0] = n_koord_array_UG[h] * ZoomFactorX;h++;
      ((theMesh->InnPosition)[innpindex])[1] = n_koord_array_UG[h] * ZoomFactorY;h++;
      /*((theMesh->InnPosition)[innpindex])[2] = n_koord_array_UG[h] - ((float)ZoomFactor/1000.0) * n_koord_array_UG[h];h++;*/
      ((theMesh->InnPosition)[innpindex])[2] = n_koord_array_UG[h] * ZoomFactorZ;h++;
    }
    else
    {
      ((theMesh->InnPosition)[innpindex])[0] = n_koord_array_UG[h];h++;
      ((theMesh->InnPosition)[innpindex])[1] = n_koord_array_UG[h];h++;
      ((theMesh->InnPosition)[innpindex])[2] = n_koord_array_UG[h];h++;
    }
  }

  return (0);
}



/****************************************************************************/
/**
* sets neighbours of elements with a O(n)-algorithm using a node-element matrix.
*
* @param ne number of elements
*
* @return value 1:error.
   */
/****************************************************************************/
INT FindElNeighbours(INT ne)
{
  INT elind,nbijkl,gefunden,kna,knb,knc,ela,elb,elc,a,b,c,elaindex,elbindex,elcindex;
  INT realelind,realelemsurf,n[4],indexn,elgefstartindex,found,ofs;


  /*laufe ueber alle Elemente*/
  for(elind=1; elind<=ne; elind++)
  {
    /*laufe ueber die 4 Nachbarn ... 4 bis 7, da im elarray, erst nach den 4 nodeids 0bis3
                                     die Nachbarids in 4 bis 7 folgen*/
    for(nbijkl=4; nbijkl<8; nbijkl++)
    {
      realelind = elind * 8;
      realelemsurf =realelind + nbijkl;                  /* == der Nachbar, um den es nun geht .*/
      if(el_array[realelemsurf] == 0)                   /*no neighbour yet respectively no bndside <-1 * Vorkommanzahl der zugewiesenen CADLastzahl>*/
      /*concerning <-1 * Vorkommanzahl der zugewiesenen CADLastzahl>
         see ReadAnsysFile*/
      {
        gefunden = FALSE;                         /*noch keinen Nachbar gefunden*/
        switch(nbijkl)
        {
        /*TODO : Die nunfolgende Fallunterscheidung realisiert die UG-Reihenfolge fuer
                 Tetraeder aus elements.c bzgl. dem Referenzelement , d.h.
                 die Vernachbarung erfolgt im Sinne von ug nicht im SInne des ANS-Formats*/
        case 4 :
          kna = 0;
          knb = 2;
          knc = 1;
          break;
        case 5 :
          kna = 1;
          knb = 2;
          knc = 3;
          break;
        case 6 :
          kna = 0;
          knb = 3;
          knc = 2;
          break;
        case 7 :
          kna = 0;
          knb = 1;
          knc = 3;
        }                        /*"switch(nbijkl)"*/
        a = el_array[realelind + kna];
        elaindex = NUOFCLMS * a;
        ela = node_element_matrix[elaindex];                         /*first element of node a*/
        while ((ela != 0)  && (gefunden == FALSE))
        {
          if(ela != elind)
          {
            b = el_array[realelind + knb];
            elbindex = NUOFCLMS * b;
            elb = node_element_matrix[elbindex];                                     /*first element of node b*/
            while ((elb != 0) && (gefunden == FALSE))
            {
              if(ela == elb)
              {
                c = el_array[realelind + knc];
                elcindex = NUOFCLMS * c;
                elc = node_element_matrix[elcindex];                                                 /*first element of node c*/
                while ((elc != 0) && (gefunden == FALSE))
                {
                  if(elb == elc)
                  {
                    gefunden = TRUE;                                                             /*neighbour found*/
                    el_array[realelemsurf] = elc;                                                             /* forward connenction */
                    elgefstartindex = 8 * elc;
                    /*the four nodes of the found element: ...*/
                    n[0] = el_array[elgefstartindex];
                    n[1] = el_array[elgefstartindex+1];
                    n[2] = el_array[elgefstartindex+2];
                    n[3] = el_array[elgefstartindex+3];
                    indexn = -1; found = FALSE;
                    /*Welcher der 4 Knoten des gef. Tetraeders ist nicht betroffen*/
                    do
                    {
                      indexn++;
                      /* a,b,c : die 3 KnotenIDs der Sideflaeche, die hier betroffen ist*/
                      if(     (n[indexn] != a) &&
                              (n[indexn] != b) &&
                              (n[indexn] != c) ) found =TRUE;
                    } while(found == 0);
                    switch (indexn)                                                             /*indexn is misser!*/
                    {
                    /*die folgenden Offsets beziehen sich auf UG, nicht auf ANSYS !*/
                    case 0 : ofs = 1; break;
                    case 1 : ofs = 2; break;
                    case 2 : ofs = 3; break;
                    case 3 : ofs = 0;
                    }                                                            /*"switch (indexn)"*/
                    el_array[elgefstartindex+4+ofs] = elind;
                    /* ... = backward connenction */
                  }                                                      /*"if(elb == elc)"*/
                  else
                  {
                    elcindex++;
                    elc = node_element_matrix[elcindex];                                                             /*next element of node c*/
                  }
                }                                                /*"while ((elc != 0) && (gefunden == FALSE))"*/
              }                                          /*"if(ela == elb)"*/
              if(gefunden == FALSE)
              {
                elbindex++;
                elb = node_element_matrix[elbindex];                                                 /*next element of node a*/
              }
            }                                    /*"while ((elb != 0) && (gefunden == FALSE))"*/
          }                              /*"if(ela != elind)"*/
          if(gefunden == FALSE)
          {
            elaindex++;
            ela = node_element_matrix[elaindex];                                     /*next element of node a*/
          }
        }                        /*"while ((ela != 0)  && (gefunden == FALSE))"*/

        /*wenn kein Nachbar geunden wurde . . .*/
        if(found != TRUE)
        {
          PrintErrorMessage('E',"FindElNeighbours","no neighbour found");
          return(1);
        }
      }                  /*"if(el_array[realelemsurf] == 0)"*/
    }            /*"for(nbijkl=4;nbijkl<8;nbijkl++)"*/
  }      /* von "for(elind=1;elind<=ne;elind++)" */
  return (0);
} /*"INT FindElNeighbours(char *node_element_matrix,char *el_array,INT ne)"*/



/****************************************************************************/
/*
   FetchATetrahedronOfThisSbd - searches a tetrahedron element corresponding to Sbd

   SYNOPSIS:
   int FetchATetrahedronOfThisSbd(SD_TYP *sbd)

   PARAMETERS:
   .  sbd - Subdomainpointer

   DESCRIPTION:
   runs over el_array and searches the first entry with -1 *  SD_NAME(sbd)
   when succeeded the corresponding elementID is returned.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      -1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int FetchATetrahedronOfThisSbd(SD_TYP *sbd)
{
  /* Achtung Returnvalue -1 bei Fehler !! nicht 1 !!!*/
  int i,gefTetraederelementID,lff,sbd_name,gefunde,maxwert;

  /*laufe ueber el_array und suche ein Element mit einem -1 * SD_NAME(sbd) --- Eintrag*/
  sbd_name = (bisherige_ID_array[SD_NAME(sbd)]);
  /* Im el_array ist nach dem negierten Wert zu suchen*/
  sbd_name = -1 * sbd_name;

  lff = 12;       /*da kommt der erste Nachbareintrag des ersten Elements*/
  maxwert = (statistik[6] * 8) + 8;
  gefunde = FALSE;
  while((gefunde == FALSE)&&(lff < maxwert))
  {
    for (i=0; i<4; i++)
    {
      /* DIRKS NEU, sbd_name oben  mit bisherige_ID_array modifiziert da el_array CAD IDs hat*/
      if(sbd_name == el_array[lff])
      {
        gefunde = TRUE;
        gefTetraederelementID = lff / 8;                        /* DIV !!!*/
        return(gefTetraederelementID);
      }
      lff ++;
    }
    lff += 4;
  }

  PrintErrorMessage('E',"FetchATetrahedronOfThisSbd","did not find such a tetrahedron");
  return(-1);
}



/****************************************************************************/
/*
   SearchAllTetrahedronsOfThisSbd - searches all tetrahedron elements corresponding a special Sbd

   SYNOPSIS:
   int SearchAllTetrahedronsOfThisSbd(int tetra_el_ID, int sbdname)

   PARAMETERS:
   .  tetra_el_ID - TetraederelementID
   .  sbdname - Name of corresponding Subdomain

   DESCRIPTION:
   runs recursively over el_array and searches the neighbours of tetra_el_ID
   Thereby all elements of a special subdomain are found.
   The variables nmbOfTetrhdrOfThisSbd und el_besucht_array
   are updated.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if  error occured

   SEE ALSO:
 */
/****************************************************************************/
int SearchAllTetrahedronsOfThisSbd(int Muster_tetra_el, int sbdname)
{
  int Nachbar_Tetra, stelle;
  int tetraside;       /*Laufvariable ueber die 4 Tetraederseiten */
  int rgbwrt;
  int neubesetzt[4];

  neubesetzt[0] = F; neubesetzt[1] = F; neubesetzt[2] = F; neubesetzt[3] = F;

  /* Laufe ueber die 4 Tetraederseiten von Muster_tetra_el */
  stelle = Muster_tetra_el * 8 + 4;
  for(tetraside = 0; tetraside < 4; tetraside++)
  {
    Nachbar_Tetra =  el_array[stelle];
    if(Nachbar_Tetra > 0)            /*nur wenn ein Nachbartetr. existiert . . . */
    {
      /*natuerlich nur wenn dieser Nachbar nicht schon besucht wurde ...*/
      if(el_besucht_array[Nachbar_Tetra] == 0)
      {

        /*************************************************************/
        /* Dieser Nachbartetr. gehoert folglich zur selben Subdomain ...*/
        el_besucht_array[Nachbar_Tetra] = sbdname;
        nmbOfTetrhdrOfThisSbd++;
        neubesetzt[tetraside] = T;
        /*Update of nmbOfSidesOfThisSbd reicht auch noch in FillSubdomainInformations*/
        /*if(elemflag_array[Nachbar_Tetra] > 0)*/ /*wenn bndelement vorliegt...*/
        /*{
                stelle2 = Nachbar_Tetra * 8 + 4;
                for(l = 0; l < 4; l++)
                {
                        if(el_array[stelle2] < 0) *//*wenn das Element bzgl dieser Seite an eine Boundary angrenzt ... */
        /*	{
                        nmbOfSidesOfThisSbd ++;
                }
                stelle2++;
           }
           }*/
        /*************************************************************/
      }
      else
      {
        /*ldgl. eine Probe - kann spaeter wieder entfallen*/
        if(el_besucht_array[Nachbar_Tetra] != sbdname)
        {

          PrintErrorMessage('E',"SearchAllTetrahedronsOfThisSbd","tetr-element belongs to 2 diff sbds");
          return(1);
        }
      }
    }
    stelle ++;
  }       /* von for */

  /*weiterer Lauf ueber die 4 Nachbarn, um die RekursionsHierarchie etwas kleiner zu halten*/
  stelle = stelle - 4;
  for(tetraside = 0; tetraside < 4; tetraside++)
  {
    /*nur die soeben neu besuchten Tetraeder sind fuer einen Rekursionschritt interessant*/
    if(neubesetzt[tetraside] == T)
    {
      rgbwrt = SearchAllTetrahedronsOfThisSbd(el_array[stelle], sbdname);
      /*ldgl. Probe:*/
      if(rgbwrt == 1)
        return (1);
      /* . . . ldgl. Probe*/
    }
    stelle ++;
  }       /* von for */

  return(0);

}



/****************************************************************************/
/*
   FillSubdomainInformations - fills MeshStructure for UG with Subdomaininformations

   SYNOPSIS:
   int FillSubdomainInformations(LGM_MESH_INFO *theMesh, int SbdName, int ug_lgm_id)

   PARAMETERS:
   .  theMesh - Pointer to mesh-structure
   .  SbdName - ug_lgm_id-ID of the corresp. Sbd
   .  ug_lgm_id - LGM-ID of the corresp. Sbd

   DESCRIPTION:
   fills the informations for all subdomains:
   Number of Sides
   Number of Corners of each side
   SideCornerIDs
   Number of elements (tetrahedrons)
   Number of elementcorners (4)
   ElementCornerIDs
   ElementNeighbours


   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int FillSubdomainInformations(LGM_MESH_INFO *theMesh, int SbdName, int ug_lgm_id)
{
  SD_TYP *sbd;
  int nmbofsides,sides_zaehler,elems_zaehler, stelle, stelle2,help;
  SFC_TYP *sd_sfc;
  int lf,elem_lf,lfv;
  /*int ug_lgm_id_minus_1*/
  int ug_sd_offs[3];
  /*only for debugging */
  int hilfszaehler;
  /* bestimme Nmbofsides Info aus den bestehenden Austauschstrukturen*/
  nmbofsides = 0;sides_zaehler =0;elems_zaehler=0;
  /*ug_lgm_id_minus_1 = ug_lgm_id -1;*/
  sbd = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
  while (SD_NAME(sbd) != SbdName)
  {
    sbd = SD_NEXT(sbd);
  }
  sd_sfc = SD_SFCS(sbd);
  while (sd_sfc != NULL)
  {
    nmbofsides = nmbofsides + SF_NMB_OF_TRIAS(SFC_SURF(sd_sfc));
    sd_sfc = SFC_NEXT(sd_sfc);
  }
  /*Zuweisung Number of Sides*/
  nmbOfSidesOfThisSbd = nmbofsides;
  /* TO ASK KLAUS : eine Spalte freilassen, da keine SbdID 0 ? vermutlich nein ===> DOCH*/
  (theMesh->nSides)[ug_lgm_id] = nmbofsides;

  if(((theMesh->Side_corners)[ug_lgm_id] = (INT *)GetTmpMem(theHeap,(nmbofsides)*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for (theMesh->Side_corners)[ug_lgm_id]");
    return(1);
  }
  /*	memset(((theMesh->Side_corners)[ug_lgm_id]),CORNERS_OF_BND_SIDE,(nmbofsides)*sizeof(INT));hier immer 3 da Tetraederseiten gemeint sind*/
  for(lfv=0; lfv< nmbofsides; lfv++)
  {
    ((theMesh->Side_corners)[ug_lgm_id])[lfv] =  CORNERS_OF_BND_SIDE;
  }

  if(((theMesh->Side_corner_ids)[ug_lgm_id] = (INT **)GetTmpMem(theHeap,(nmbofsides)*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for (theMesh->Side_corner_ids)[ug_lgm_id]");
    return(1);
  }
  for(lf = 0; lf <nmbofsides; lf++)
  {
    if((((theMesh->Side_corner_ids)[ug_lgm_id])[lf] = (INT *)GetTmpMem(theHeap,(CORNERS_OF_BND_SIDE)*sizeof(INT), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for ((theMesh->Side_corner_ids)[ug_lgm_id])[lf]");
      return(1);
    }
  }


  /*elements*/
  if(((theMesh->Element_corners)[ug_lgm_id] = (INT *)GetTmpMem(theHeap,(nmbOfTetrhdrOfThisSbd)*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for (theMesh->Element_corners)[ug_lgm_id]");
    return(1);
  }
  /*	memset(((theMesh->Element_corners)[ug_lgm_id]),CORNERS_OF_ELEMENT,(nmbOfTetrhdrOfThisSbd)*sizeof(INT)); hier immer 4 da Tetraederseiten gemeint sind*/
  for (lfv=0; lfv < nmbOfTetrhdrOfThisSbd; lfv++)
  {
    ((theMesh->Element_corners)[ug_lgm_id])[lfv] = CORNERS_OF_ELEMENT;
  }
  if(((theMesh->Element_SideOnBnd)[ug_lgm_id] = (INT *)GetTmpMem(theHeap,(nmbOfTetrhdrOfThisSbd)*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for (theMesh->Element_SideOnBnd)[ug_lgm_id]");
    return(1);
  }
  memset(((theMesh->Element_SideOnBnd)[ug_lgm_id]),0,(nmbOfTetrhdrOfThisSbd)*sizeof(INT));      /*hier zunaechst 0; d.h. keine BndSide*/
  /*die Boundarysides muessen noch bitweise gesetzt werden; aber wo  ... TODO  GOON HERE ...*/

  if(((theMesh->Element_corner_ids)[ug_lgm_id] = (INT **)GetTmpMem(theHeap,(nmbOfTetrhdrOfThisSbd)*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for (theMesh->Element_corner_ids)[ug_lgm_id]");
    return(1);
  }
  for(lf = 0; lf <nmbOfTetrhdrOfThisSbd; lf++)
  {
    if((((theMesh->Element_corner_ids)[ug_lgm_id])[lf] = (INT *)GetTmpMem(theHeap,(CORNERS_OF_ELEMENT)*sizeof(INT), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillSubdomainInformations"," ERROR: No memory for ((theMesh->Element_corner_ids)[ug_lgm_id])[lf]");
      return(1);
    }
  }


  /*Lauf �ber el_besucht_array
     nehme Tetraeder, die zur Sbd gehoeren z�hle sie dabei zwecks Probe
          trage die 4 CornerIDs ein in der Reihenfolge wie sie kommen Achtung UGIDs-stattCAD-IDs verwenden
          trage ggf. die SideCornerIDs ein - Achtung in UG-Reihenfolge
          dabei Update of nbofsides zwecks Probe*/
  hilfszaehler = 0;       /*only for debugging*/
  for(elem_lf=1; elem_lf<=statistik[6]; elem_lf++)
  {
    if (el_besucht_array[elem_lf] == SbdName)
    {
      /*wenn es sich also um ein Element dieser Subdomain handelt ...*/
      stelle2 = stelle = elem_lf * 8;
      /*trage die 4 ElementNdIDs ein*/
      for (lf=0; lf<4; lf++)
      {
        (((theMesh->Element_corner_ids)[ug_lgm_id])[elems_zaehler])[lf] = point_array[el_array[stelle]];
        stelle++;
      }


      /*wenn �berhaupt ein Boundaryelement...*/
      if(elemflag_array[elem_lf] > 0)
      {

        hilfszaehler = 0;
        /*stelle stimmt bereits , s.o */
        for(lf = 0; lf < 4; lf++)
        {
          if(el_array[stelle] < 0)                               /*wenn das Element bzgl. dieser Seite an eine Boundary angrenzt ... */
          /* die Seitenindizes beziehen sich dabei auf die UG-TetraederSideNummerierung
             d.h. Side0:021; Side1:123; Side2:032; Side3:013
             Diese UG-Nummerierung wurde bei der Zuweisung in FindElNeighbours beruecksichtigt*/
          {
            /*SideCorners von dieser Side eintragen*/
            /* zusaetzlich an dieser Stelle: Update der bitweise gesetzten
               Komponente theMesh->Element_SideOnBnd . . .
               wird hier in der switch-Auswahl gesetzt.*/
            switch (lf)
            {
            /*die folgenden Zuw. beziehen sich auf UG, nicht auf ANSYS !*/
            case 0 : ug_sd_offs[0] = 0; ug_sd_offs[1] = 2; ug_sd_offs[2] = 1;
              /*estes Bit wird gesetzt fuer Side 0 . . . */
              ((theMesh->Element_SideOnBnd)[ug_lgm_id])[elems_zaehler] += 1;
              break;
            case 1 : ug_sd_offs[0] = 1; ug_sd_offs[1] = 2; ug_sd_offs[2] = 3;
              /*zweites  Bit wird gesetzt fuer Side 1 . . . */
              ((theMesh->Element_SideOnBnd)[ug_lgm_id])[elems_zaehler] += 2;
              break;
            case 2 : ug_sd_offs[0] = 0; ug_sd_offs[1] = 3; ug_sd_offs[2] = 2;
              /*drittes  Bit wird gesetzt fuer Side 2 . . . */
              ((theMesh->Element_SideOnBnd)[ug_lgm_id])[elems_zaehler] += 4;
              break;
            case 3 : ug_sd_offs[0] = 0; ug_sd_offs[1] = 1; ug_sd_offs[2] = 3;
              /*viertes  Bit wird gesetzt fuer Side 3 . . . */
              ((theMesh->Element_SideOnBnd)[ug_lgm_id])[elems_zaehler] += 8;
            }                                    /*"switch (lf)"*/
            for (help=0; help<3; help++)                                    /*ueber die 3 sidenodeIDs*/
            {
              (((theMesh->Side_corner_ids)[ug_lgm_id])[sides_zaehler])[help] = point_array[el_array[stelle2 + ug_sd_offs[help]]];
            }
            hilfszaehler ++;
            sides_zaehler ++;
          }
          stelle++;
        }
        if(hilfszaehler == 0)
        {
          UserWriteF("ERROR in FillSubdomainInformations Boundaryelement %d hat keine einzige BndSide\n",(int)elem_lf);
          return(1);
        }
      }
      elems_zaehler++;
    }
  }
  /*Probe:*/
  if(elems_zaehler != nmbOfTetrhdrOfThisSbd)
  {
    PrintErrorMessage('E',"FillSubdomainInformations","elems_zaehler != nmbOfTetrhdrOfThisSbd");
    return(1);
  }
  if(sides_zaehler != nmbOfSidesOfThisSbd)
  {
    PrintErrorMessage('E',"FillSubdomainInformations","sides_zaehler != nmbOfSidesOfThisSbd");
    return(1);
  }

  return(0);
} /*end of FillSubdomainInformations*/




/****************************************************************************/
/*
   FillBndPointInformations - fills the  BNDP-Part of the lgmExchangestructure

   SYNOPSIS:
   int FillBndPointInformations(LGM_MESH_INFO *theMesh, int *bnd_pnt_srfc, int *bnd_pnt_cntr, int *bnd_pnt_cor_TrID, int *bnd_pnt_case)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.
   .  bnd_pnt_srfc - includes all surfaces to each BNDP.
   .  bnd_pnt_cntr - includes numbers of surfaces of each BNDP.
   .  bnd_pnt_cor_TrID - includes triangleIDs of each BNDPsurfaceRelation.
   .  bnd_pnt_case - includes caseIDs of each BNDPsurfaceTriangleRelation.

   DESCRIPTION:
   Using the input data structures bnd_pnt_srfc, bnd_pnt_cntr, bnd_pnt_cor_TrID and bnd_pnt_case
   FillBndPointInformations fills the  BNDP-Part of the lgmExchangestructure

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int FillBndPointInformations(LGM_MESH_INFO *theMesh, int *bnd_pnt_srfc, int *bnd_pnt_cntr, int *bnd_pnt_cor_TrID, int *bnd_pnt_case)
{
  /* int nBndP;*/                       /* nb. of boundary points              */
  /* int *BndP_nSurf;*/                 /* nb. of surfaces per bound. point    */
  /* int **BndP_SurfID;*/               /* id of each surface                  */
  /* int **BndP_Cor_TriaID;*/           /* id of corresponding triangle of each surface*/
  /* float ***BndP_lcoord;*/            /* local coord of BndP on each surface */

  int b,s,stelle;
  /* Speicher holen f�r int *BndP_nSurf*/
  if((theMesh->BndP_nSurf = (INT *)GetTmpMem(theHeap,statistik[1]*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for theMesh->BndP_nSurf !!!");
    return(1);
  }

  /* Speicher holen f�r int **BndP_SurfID*/
  if((theMesh->BndP_SurfID = (INT **)GetTmpMem(theHeap,statistik[1]*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for theMesh->BndP_SurfID !!!");
    return(1);
  }

  /* Speicher holen f�r int **BndP_Cor_TriaID*/
  if((theMesh->BndP_Cor_TriaID = (INT **)GetTmpMem(theHeap,statistik[1]*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for theMesh->BndP_Cor_TriaID !!!");
    return(1);
  }

  /* Speicher holen f�r int ***BndP_lcoord*/
  if((theMesh->BndP_lcoord = (float ***)GetTmpMem(theHeap,statistik[1]*sizeof(float**), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for theMesh->BndP_lcoord !!!");
    return(1);
  }

  /* Lauf �ber die statistik[1] BoundaryPoints . . .*/
  for(b = 0; b < statistik[1]; b++)
  {
    /* Anzahl Surfaces je BndP setzen */
    (theMesh->BndP_nSurf)[b] = bnd_pnt_cntr[b];

    /* Speicher holen f�r int *BndP_SurfID */
    if(((theMesh->BndP_SurfID)[b] = (INT *)GetTmpMem(theHeap,bnd_pnt_cntr[b]*sizeof(INT), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for (theMesh->BndP_SurfID)[b] !!!");
      return(1);
    }

    /* Speicher holen f�r int *BndP_Cor_TriaID */
    if(((theMesh->BndP_Cor_TriaID)[b] = (INT *)GetTmpMem(theHeap,bnd_pnt_cntr[b]*sizeof(INT), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for (theMesh->BndP_Cor_TriaID)[b] !!!");
      return(1);
    }

    /* Speicher holen f�r int **BndP_lcoord */
    if(((theMesh->BndP_lcoord)[b] = (float **)GetTmpMem(theHeap,bnd_pnt_cntr[b]*sizeof(float*), ANS_MarkKey)) == NULL)
    {
      PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for (theMesh->BndP_lcoord)[b] !!!");
      return(1);
    }

    stelle = NU_SFCES_BNDP * b;
    for (s = 0; s < bnd_pnt_cntr[b]; s++)
    {
      /*  einen der bnd_pnt_cntr[b] Surface-IDs fuer den b-ten BndP setzen */
      ((theMesh->BndP_SurfID)[b])[s] = bnd_pnt_srfc[stelle];
      /*  die zugeh�rige TriangleID setzen */
      ((theMesh->BndP_Cor_TriaID)[b])[s] = bnd_pnt_cor_TrID[stelle];
      /*Speicher holen f�r die zugeh�rigen lokalen Koordinaten *BndP_lcoord */
      if((((theMesh->BndP_lcoord)[b])[s] = (float *)GetTmpMem(theHeap,NMBOFLOCLCOORDS*sizeof(float), ANS_MarkKey)) == NULL)
      {
        PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for ((theMesh->BndP_lcoord)[b])[s] !!!");
        return(1);
      }
      /*lokaleKoordinaten setzen*/
      switch (bnd_pnt_case[stelle])
      {
      case 0 :
        (((theMesh->BndP_lcoord)[b])[s])[0] = 1.0;
        (((theMesh->BndP_lcoord)[b])[s])[1] = 0.0;
        break;
      case 1 :
        (((theMesh->BndP_lcoord)[b])[s])[0] = 0.0;
        (((theMesh->BndP_lcoord)[b])[s])[1] = 1.0;
        break;
      case 2 :
        (((theMesh->BndP_lcoord)[b])[s])[0] = 0.0;
        (((theMesh->BndP_lcoord)[b])[s])[1] = 0.0;
        break;
      default :
        PrintErrorMessage('E',"FillBndPointInformations","kein Standardfall <0,1,2> bzgl.lok. Koords");
        return (1);
      }
      stelle++;
    }
  }     /* endof Lauf �ber die statistik[1] BoundaryPoints . . .*/
  return(0);
} /*end of FillBndPointInformations*/



/****************************************************************************/
/*
   EvalBndPointInformations - reads necessary informations for boundary points

   SYNOPSIS:
   int	EvalBndPointInformations(LGM_MESH_INFO *theMesh)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.

   DESCRIPTION:
   reads necessary informations for boundary points

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int     EvalBndPointInformations(LGM_MESH_INFO *theMesh)
{
  int *boundary_point_counter;      /*z�hlt Anzahl Surfaces pro BndP*/
  int *boundary_point_surface_array;      /*beinhaltet Surfaces f�r jeden BndPoint*/
  int *boundary_point_case_array;      /*beinhaltet die Fallzahl <1,2 oder3> f�r jede BNDPSFC-Rel.*/
  int *boundary_point_corresp_TriaID_array;      /*beinhaltet zu jeder BNDP-Sfce Relation die zugeh. TriaID*/
  int Sfc_ID, Tria_ID, i, BndPID_UG, stelle, maxstelle,einfuegestellegefunden, gibtsschon ;
  SF_TYP *Surface;
  TRIANGLE_TYP *Triangle;
  int cad_id;
  /*alt	int lf, cad_id_gefunden;*/


  /* Fetch Memory ...*/

  boundary_point_counter = (int *)GetTmpMem(theHeap,statistik[1]*sizeof(int),ANS_MarkKey);
  if ( boundary_point_counter == NULL )
  {
    PrintErrorMessage('E',"ansys2lgm"," ERROR: No memory for boundary_point_counter in EvalBndPointInformations ");
    return(1);
  }
  memset(boundary_point_counter,0,(statistik[1])*sizeof(int));

  boundary_point_surface_array = (int *)GetTmpMem(theHeap,statistik[1]*NU_SFCES_BNDP*sizeof(int),ANS_MarkKey);
  if ( boundary_point_surface_array == NULL )
  {
    PrintErrorMessage('E',"ansys2lgm"," ERROR: No memory for boundary_point_surface_array in EvalBndPointInformations ");
    return(1);
  }
  /*Achtung Init mit -1 da Surface "0" existiert*/
  memset(boundary_point_surface_array,-1,(statistik[1])*NU_SFCES_BNDP*sizeof(int));

  boundary_point_case_array = (int *)GetTmpMem(theHeap,statistik[1]*NU_SFCES_BNDP*sizeof(int),ANS_MarkKey);
  if ( boundary_point_case_array == NULL )
  {
    PrintErrorMessage('E',"ansys2lgm"," ERROR: No memory for boundary_point_case_array in EvalBndPointInformations ");
    return(1);
  }
  /*Achtung Init mit -1, da es den Fall "0" gibt*/
  memset(boundary_point_case_array,-1,(statistik[1])*NU_SFCES_BNDP*sizeof(int));

  boundary_point_corresp_TriaID_array = (int *)GetTmpMem(theHeap,statistik[1]*NU_SFCES_BNDP*sizeof(int),ANS_MarkKey);
  if ( boundary_point_corresp_TriaID_array == NULL )
  {
    PrintErrorMessage('E',"ansys2lgm"," ERROR: No memory for boundary_point_corresp_TriaID_array in EvalBndPointInformations ");
    return(1);
  }
  /*Achtung Init mit -1 da es Trias mit ID "0" gibt*/
  memset(boundary_point_corresp_TriaID_array,-1,(statistik[1])*NU_SFCES_BNDP*sizeof(int));


  /*Laufe �ber die Surfaces*/
  Surface = EXCHNG_TYP2_ROOT_SFC(ExchangeVar_2_Pointer);
  for (Sfc_ID = 0; Sfc_ID < NMB_OF_SFCES(DomainInfo_Pointer); Sfc_ID++)
  {
    /*just for debugging - kann sp�ter gel�scht werden*/
    if (Surface == NULL)
    {
      PrintErrorMessage('E',"EvalBndPointInformations","Surface-Laufpointer is NULL !!");
      return (1);
    }

    /*Laufe �ber die bereits orientierten Triangles jeder Surface*/
    Triangle = SF_TRIAS(Surface);
    for(Tria_ID = 0; Tria_ID < SF_NMB_OF_TRIAS(Surface); Tria_ID++)
    {
      /*just for debugging - kann sp�ter gel�scht werden*/
      if (Triangle == NULL)
      {
        PrintErrorMessage('E',"EvalBndPointInformations","Triangle-Laufpointer is NULL !!");
        return (1);
      }

      /*Laufe �ber die 3 BndPoints*/
      for(i=0; i<3; i++)
      {
        BndPID_UG = (TRIA_SFE_KN(Triangle))->nodeid[i];

        /*laufe beginnend bei boundary_point_surface_array[BndPID_UG*NU_SFCES_BNDP]
           bis zum ersten unbesetzten Feld <-1 !> und pr�fe dabei, ob es
           diese BndP-Surface-Relation nicht schon gibt ...*/
        stelle = NU_SFCES_BNDP * BndPID_UG;
        maxstelle = stelle + NU_SFCES_BNDP;
        gibtsschon = 0; einfuegestellegefunden = -1;
        do
        {
          if(boundary_point_surface_array[stelle] == -1)                              /*wenn noch kein SfcID-Eintrag*/
          {
            einfuegestellegefunden = stelle;
          }
          else if(boundary_point_surface_array[stelle] == Sfc_ID)                              /*wenn genau diese Surface bereits eingetragen wurde*/
          {
            gibtsschon = 1;
          }
          else
          {
            stelle++;
          }
        } while((einfuegestellegefunden == -1) && (gibtsschon == 0) && (stelle < maxstelle));

        if(stelle == maxstelle)
        {
          PrintErrorMessage('E',"ansys2lgm"," NU_SFCES_BNDP ist zu klein in EvalBndPointInformations");
          return(1);
        }
        else if(gibtsschon == 0)                         /*wenn es die Surface bei diesem BNDP noch nicht gibt*/
        {
          /*just for debug : dann muss aber auch eine Einfuegestelle gefunden worden sein*/
          if(einfuegestellegefunden == -1)
          {
            PrintErrorMessage('E',"ansys2lgm","<einfuegestellegefunden == -1> kann nicht sein in EvalBndPointInformations");
            return(1);
          }
          /*Es wird eine neue BndP-Surface-Relation eingetragen ... */
          boundary_point_surface_array[einfuegestellegefunden] = Sfc_ID;
          boundary_point_counter[BndPID_UG] ++;
          boundary_point_corresp_TriaID_array[einfuegestellegefunden] = Tria_ID;
          boundary_point_case_array[einfuegestellegefunden] = i;                               /*Fall 0,1 oder 2; je nachdem*/
        }                         /*von: wenn es die Surface bei diesem BNDP noch nicht gibt*/
      }                  /*vom Lauf �ber die 3 BoundaryPoints*/
      Triangle = TRIA_NEXT(Triangle);
    }            /* vom Lauf �ber die Triangles der Surface */
    Surface = SFC_NEXT(Surface);
  }      /* vom Lauf �ber die Surfaces*/

  /*just for debug Probe*/
  /*Probe: Ist ein BoundaryPoint unbesetzt geblieben ?*/
  for(i = 0; i< statistik[1]; i++)
  {
    if(boundary_point_counter[i] == 0)             /*d.h. Geh�rt dieser BoundaryPoint zu keiner Surface ?*/
    {
      /*just for debugging*/
      /*Bestimmung des ersten Ausgabewertes der Fehlermeldung*/

      /*point_array durchlaufen und CAD-ID zu "i" suchen !!!!*/
      /*	alt		lf = 1;
                              cad_id_gefunden = F;
                              do
                              {
                                      if(point_array[lf] == i)
                                      {
                                              cad_id = lf;
                                              cad_id_gefunden = T;
                                      }
                                      lf++;
                              }while((lf <= nbofnds)&&(cad_id_gefunden == F));

                              if (cad_id_gefunden == F)
                              {
                                      PrintErrorMessage('E',"ansys2lgm","EvalBndPointInformations : cad_id_gefunden == F");
                                      return(1);
                              }*/

      /*neu :*/
      cad_id = point_array_UG_CAD[i];


      UserWriteF("CAD_BndP%d bzw. UGBndP%d geh�rt zu keiner Surface ?!?! ERROR in EvalBndPointInformations\n",(int)cad_id,(int)i);
      return(1);
    }
  }

  /* Nun k�nnen die neuen Komponenten von lgm_mesh_info gef�llt werden: */
  if (FillBndPointInformations(theMesh,boundary_point_surface_array,boundary_point_counter,boundary_point_corresp_TriaID_array,boundary_point_case_array) != 0)
  {
    PrintErrorMessage('E',"EvalBndPointInformations->FillBndPointInformations","execution failed");
    return (1);
  }

  return(0);
} /*end of EvalBndPointInformations*/



/****************************************************************************/
/*
   BndPoint_Line_Alloc_Mem - allocates mem for the  BNDP-Line-Part

   SYNOPSIS:
   int BndPoint_Line_Alloc_Mem(LGM_MESH_INFO *theMesh, int *boundary_point_line_counter)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.
   .  boundary_point_line_counter - includes number of lines of each BNDP.

   DESCRIPTION:
   allokiert und initiailisierSpeicher fuer die BndP->linekomponenten
   int *BndP_nLine, int **BndP_LineID, float **BndP_lcoord_left und float **BndP_lcoord_right
   int *BndP_nLine wird dabei aus boundary_point_line_counter heraus gefuellt.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int BndPoint_Line_Alloc_Mem(LGM_MESH_INFO *theMesh, int *boundary_point_line_counter)
{
  int b;
  float float_help;
  int hlp;

  if((theMesh->BndP_nLine = (INT *)GetTmpMem(theHeap,statistik[1]*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"BndPoint_Line_Alloc_Mem"," ERROR: No memory for theMesh->BndP_nLine !!!");
    return(1);
  }

  if((theMesh->BndP_LineID = (INT **)GetTmpMem(theHeap,statistik[1]*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"BndPoint_Line_Alloc_Mem"," ERROR: No memory for theMesh->BndP_LineID !!!");
    return(1);
  }

  if((theMesh->BndP_lcoord_left = (float **)GetTmpMem(theHeap,statistik[1]*sizeof(float*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"BndPoint_Line_Alloc_Mem"," ERROR: No memory for theMesh->BndP_lcoord_left !!!");
    return(1);
  }

  if((theMesh->BndP_lcoord_right = (float **)GetTmpMem(theHeap,statistik[1]*sizeof(float*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"BndPoint_Line_Alloc_Mem"," ERROR: No memory for theMesh->BndP_lcoord_right !!!");
    return(1);
  }

  for(b = 0; b < statistik[1]; b++)
  {
    /* Line-Anzahl je BndP setzen */
    (theMesh->BndP_nLine)[b] = boundary_point_line_counter[b];

    /*Speicher anlegen fuer LineIDs und lokale Kordinaten*/
    /* TOASK: Gibt GetTmpMem NULL zurueck,wenn zweiter Parameter == 0 ?
       Antwort: JA ! siehe in low->heaps.c->GetMem*/
    /*****************************************************/
    /* Speicher holen f�r int *BndP_LineID */
    if (boundary_point_line_counter[b] != 0)
    {
      if(((theMesh->BndP_LineID)[b] = (INT *)GetTmpMem(theHeap,(boundary_point_line_counter[b])*sizeof(INT), ANS_MarkKey)) == NULL)
      {
        PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for <theMesh->BndP_LineID>[b] !!!");
        return(1);
      }
      memset((theMesh->BndP_LineID)[b],-1,(boundary_point_line_counter[b])*sizeof(INT));                  /*-1 da LineIDs von 0 bis .. gehen*/
    }
    else
    {
      (theMesh->BndP_LineID)[b] = NULL;
    }

    if (boundary_point_line_counter[b] != 0)
    {
      if(((theMesh->BndP_lcoord_left)[b] = (float *)GetTmpMem(theHeap,(boundary_point_line_counter[b])*sizeof(float), ANS_MarkKey)) == NULL)
      {
        PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for <theMesh->BndP_lcoord_left>[b] !!!");
        return(1);
      }
      /*Iint mit -2.0 da lok. LNkoords von 1.0, 0.0, 1.0, 2.0, ... gehen*/
      float_help = -2.0;
      for (hlp=0; hlp<boundary_point_line_counter[b]; hlp++)
      {
        ((theMesh->BndP_lcoord_left)[b])[hlp] = float_help;
      }
    }
    else
    {
      (theMesh->BndP_lcoord_left)[b] = NULL;
    }

    if (boundary_point_line_counter[b] != 0)
    {
      if(((theMesh->BndP_lcoord_right)[b] = (float *)GetTmpMem(theHeap,(boundary_point_line_counter[b])*sizeof(float), ANS_MarkKey)) == NULL)
      {
        PrintErrorMessage('E',"FillBndPointInformations"," ERROR: No memory for <theMesh->BndP_lcoord_right>[b] !!!");
        return(1);
      }
      /*Iint mit -2.0 da lok. LNkoords von 1.0, 0.0, 1.0, 2.0, ... gehen*/
      float_help = -2.0;
      for (hlp=0; hlp<boundary_point_line_counter[b]; hlp++)
      {
        ((theMesh->BndP_lcoord_right)[b])[hlp] = float_help;
      }
    }
    else
    {
      (theMesh->BndP_lcoord_right)[b] = NULL;
    }
  }

  return(0);
} /* of BndPoint_Line_Alloc_Mem*/




/****************************************************************************/
/*
   Put_BndPLineRelation_In_theMesh - fills a  BNDP-LINE-Relation of the lgmExchangestructure

   SYNOPSIS:
   int Put_BndPLineRelation_In_theMesh(LGM_MESH_INFO *theMesh, int BndP_UG_ID, int Line_ID, float linkelokkoord, float rechtelokkoord)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.
   .  BndP_UG_ID - the correct id of the oundaryPoint.
   .  Line_ID - the correspondug PoylylineID.
   .  linkelokkoord - left loc coord of this BndP-Polyline-Relation.
   .  rechtelokkoord - right loc coord of this BndP-Polyline-Relation.

   DESCRIPTION:
   traegt eine spezielle Bndp_Line_Relation in die Meshstruktur ein
   und prueft zuvor auf ausreichenden vorbereiteten Speicheplatz ab.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int Put_BndPLineRelation_In_theMesh(LGM_MESH_INFO *theMesh, int BndP_UG_ID, int Line_ID, float linkelokkoord, float rechtelokkoord)
{
  int freie_stelle;

  freie_stelle = 0;

  if((theMesh->BndP_LineID)[BndP_UG_ID] == NULL)
  {
    PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," ERROR: No memory prepared for <theMesh->BndP_LineID>[BndP_UG_ID]!!!");
    return(1);
  }

  /*suche erste freie Einfuegestelle*/
  while( ((theMesh->BndP_LineID)[BndP_UG_ID])[freie_stelle] != -1)       /*solange schon was eingefuegt wurde*/
  {
    freie_stelle++;
    if( freie_stelle == (theMesh->BndP_nLine)[BndP_UG_ID] )
    {
      PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," ERROR: <theMesh->BndP_LineID>[] already full!!!");
      return(1);
    }
  }

  /*zusaetzliche Checkfunktion kann spaeter evtl. geloescht werden:*/
  /*several tests . . .*/
  if((theMesh->BndP_lcoord_left)[BndP_UG_ID] == NULL)
  {
    PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," ERROR: No memory prepared for <theMesh->BndP_lcoord_left>[BndP_UG_ID]!!!");
    return(1);
  }
  if((theMesh->BndP_lcoord_right)[BndP_UG_ID] == NULL)
  {
    PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," ERROR: No memory prepared for <theMesh->BndP_lcoord_right>[BndP_UG_ID]!!!");
    return(1);
  }
  if(((theMesh->BndP_lcoord_left)[BndP_UG_ID])[freie_stelle] != -2.0)
  {
    PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," <<theMesh->BndP_lcoord_left>[BndP_UG_ID]>[freie_stelle] != -2.0>!!!");
    return(1);
  }
  if(freie_stelle > 0)
  {
    if(((theMesh->BndP_lcoord_left)[BndP_UG_ID])[freie_stelle-1] == -2.0)
    {
      PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," <<theMesh->BndP_lcoord_left>[BndP_UG_ID]>[freie_stelle-1] == -2.0>!!!");
      return(1);
    }
  }
  if(((theMesh->BndP_lcoord_right)[BndP_UG_ID])[freie_stelle] != -2.0)
  {
    PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," <<theMesh->BndP_lcoord_right>[BndP_UG_ID]>[freie_stelle] != -2.0>!!!");
    return(1);
  }
  if(freie_stelle > 0)
  {
    if(((theMesh->BndP_lcoord_right)[BndP_UG_ID])[freie_stelle-1] == -2.0)
    {
      PrintErrorMessage('E',"Put_BndPLineRelation_In_theMesh"," <<theMesh->BndP_lcoord_right>[BndP_UG_ID]>[freie_stelle-1] == -2.0>!!!");
      return(1);
    }
  }
  /* . . . end of several tests which may be deleted later on*/

  /*wenn bis hier her gekommen und alle OK ==> Einfuegestelle gefunden := freie_stelle*/
  /* !!! FILL IN ...*/
  ((theMesh->BndP_LineID)[BndP_UG_ID])[freie_stelle] = Line_ID;
  ((theMesh->BndP_lcoord_left)[BndP_UG_ID])[freie_stelle] = linkelokkoord;
  ((theMesh->BndP_lcoord_right)[BndP_UG_ID])[freie_stelle] = rechtelokkoord;


  return(0);
} /*von Put_BndPLineRelation_In_theMesh*/




/****************************************************************************/
/*
   Prepair_BndPointLineRelations_fortheMesh - fills the  BNDP-LINE-Part of the lgmExchangestructure

   SYNOPSIS:
   int Prepair_BndPointLineRelations_fortheMesh(LGM_MESH_INFO *theMesh)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.

   DESCRIPTION:
   laeuft ueber alle Lines bestimmt dabei die Bndp_Line_Relationen und uebergibt
   diese dann an Put_BndPLineRelation_In_theMesh

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int Prepair_BndPointLineRelations_fortheMesh(LGM_MESH_INFO *theMesh)
{
  PL_TYP *Line;
  int Plyln_ID;
  PL_LINE_TYP *pllyln, *pllyln_last;
  int i, der_erste_Point_der_Line, die_Point_ID;
  int BndP_UG_ID, Line_ID;
  float linkelokkoord;
  float rechtelokkoord;
  int Stop;

  /* TODO zweiter Lauf �ber die Lines
     dies'mal werden die Mesh->substuctures gesetzt*/
  Line = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  for (Plyln_ID = 0; Plyln_ID < NMB_OF_PLINES(DomainInfo_Pointer); Plyln_ID++)
  {
    /*just for debugging - kann sp�ter gel�scht werden*/
    if (Line == NULL)
    {
      PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh","Line-Laufpointer is NULL !!");
      return (1);
    }

    /*Laufe �ber die bereits sortierten Points jeder Line*/
    pllyln = PL_LINES(Line);
    /*die erste KnotenID der Polyline ...*/
    der_erste_Point_der_Line = LI_NDID1(PL_LINES_LINE(pllyln));
    Stop = PL_NMB_OF_POINTS(Line) - 2;
    /*der erste Point und die beiden letzten PL-Lines werden
       gesondert behandelt,
       da der zyklische Fall auftreten kann
       die letzten beiden PL-Lines ==> d.h. die letzten beiden Points
       da ja immer der zeite genommen wird
       Folgende Faelle muessen unterschieden werden koennen:
       <ab>,<bc>,<cd>,<de> von <ab>,<bc>,<cd>,<da>
       man benoetigt alo stets die letzten beiden PL-Lines gesondert.*/
    /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */

    /*weise KnotenIDs der gefundenen Polyline zu*/
    /*... laufe dazu ueber die Polylineteilstreckenzuege ...*/
    /*Der Fall der zyklischen Polyline ist von C.Tapp bisher nicht beruecksichtigt
       Ich uebergeb in diesm speziellen Fall einen Punkt zweimal
       d.h. am Anfang und am Ende,; ausserdem ist PL_NMB_OF_POINTS dann um 1 zu gross
       in Zukunft: sollte lgm-Domain so ueberarbeitet werden,
       dass die PL-Struktur hiefuer ein Flag besitzt.*/
    /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
    for(i=1; i<Stop; i++)
    {
      if(pllyln != NULL)
      {
        die_Point_ID = LI_NDID2(PL_LINES_LINE(pllyln));
        /*jetzt immer den zweiten Knoten nehmen, da ja die
           Anzahl der Points einer Polyline um 1 groesser ist
           als die Anzahl der PL-LineSegmente derselbigen ...*/

        /*BndPoint-Line-Relation fuer Eintrag in Mesh vorbereiten*/
        BndP_UG_ID = die_Point_ID;
        Line_ID = Plyln_ID;
        linkelokkoord = -1.0 + (float) i;
        rechtelokkoord = 1.0 + (float) i;
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }

        pllyln = PL_LINES_NXT(pllyln);
      }
      else
      {
        UserWrite("ERROR: in Prepair_BndPointLineRelations_fortheMesh: PolylineLine is missing Case1!!");
        return (1);
      }
    }            /*Ende von for<i=1; i<Stop; i++>*/

    if(pllyln == NULL)
    {
      UserWrite("ERROR: in Prepair_BndPointLineRelations_fortheMesh: PolylineLine is missing Case2!!");
      return (1);
    }


    /*gesonderte Betrachtung des ersten Points und der letzten beiden PL-Lines*/
    if(PL_NMB_OF_POINTS(Line) ==  2)
    {
      /*Fall1 Polyline besteht aus nur einer PL-Line*/

      /*Polyline muss mind. aus zwei versc. Points bestehen:*/
      if ( LI_NDID1(PL_LINES_LINE(pllyln)) == LI_NDID2(PL_LINES_LINE(pllyln)) )
      {
        UserWrite("ERROR: in Prepair_BndPointLineRelations_fortheMesh: only 1 PolylineLine with 2 identical nodes , cyclic ");
        return (1);
      }

      /*BndPoint-Line-Relation fuer Eintrag in Mesh vorbereiten*/

      /*der erste*/
      BndP_UG_ID = LI_NDID1(PL_LINES_LINE(pllyln));
      Line_ID = Plyln_ID;
      linkelokkoord = -1.0;
      rechtelokkoord = 1.0;
      /*Subfunktionsaufruf*/
      if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
      {
        PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
        return (1);
      }

      /*der zweite*/
      BndP_UG_ID = LI_NDID2(PL_LINES_LINE(pllyln));
      Line_ID = Plyln_ID;
      linkelokkoord = 0.0;
      rechtelokkoord = 1234567890.0;
      /*Subfunktionsaufruf*/
      if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
      {
        PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
        return (1);
      }

    }
    else
    {
      /* pllyln ist nun die vorletzte PL-Line*/

      /*die letzte PL-Line:*/
      pllyln_last = PL_LINES_NXT(pllyln);
      if(pllyln_last == NULL)
      {
        UserWrite("ERROR: in Prepair_BndPointLineRelations_fortheMesh: PolylineLine is missing Case3 pllyln_last!!");
        return (1);
      }

      /*TOCHECK: Liegt eine zyklische Polyline vor:*/
      /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
      if( LI_NDID2(PL_LINES_LINE(pllyln_last)) == der_erste_Point_der_Line)
      {
        /*Fall2: Spezialfall : <zyklische Polyline>*/
        /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */

        /*BndPoint-Line-Relation fuer Eintrag in Mesh vorbereiten*/

        /*der erste*/
        BndP_UG_ID = der_erste_Point_der_Line;
        Line_ID = Plyln_ID;
        linkelokkoord = -2.0 + (float) PL_NMB_OF_POINTS(Line);
        rechtelokkoord = 1.0;
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }

        /*der letzte*/
        BndP_UG_ID = LI_NDID2(PL_LINES_LINE(pllyln));
        Line_ID = Plyln_ID;
        linkelokkoord = -3.0 + (float) PL_NMB_OF_POINTS(Line);
        rechtelokkoord = 0.0;
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }

      }
      else
      {
        /*Fall3: Normalfall*/

        /*BndPoint-Line-Relation fuer Eintrag in Mesh vorbereiten*/

        /*der erste*/
        BndP_UG_ID = der_erste_Point_der_Line;
        Line_ID = Plyln_ID;
        linkelokkoord = -1.0;
        rechtelokkoord = 1.0;
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }

        /*der letzte*/
        BndP_UG_ID = LI_NDID2(PL_LINES_LINE(pllyln_last));
        Line_ID = Plyln_ID;
        linkelokkoord = -2.0 + (float) PL_NMB_OF_POINTS(Line);
        rechtelokkoord = 1234567890.0;
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }

        /*der vorletzte*/
        BndP_UG_ID = LI_NDID2(PL_LINES_LINE(pllyln));
        Line_ID = Plyln_ID;
        linkelokkoord = -3.0 + (float) PL_NMB_OF_POINTS(Line);
        rechtelokkoord = -1.0 + (float) PL_NMB_OF_POINTS(Line);
        /*Subfunktionsaufruf*/
        if (Put_BndPLineRelation_In_theMesh(theMesh, BndP_UG_ID, Line_ID, linkelokkoord, rechtelokkoord) != 0)
        {
          PrintErrorMessage('E',"Prepair_BndPointLineRelations_fortheMesh->Put_BndPLineRelation_In_theMesh","execution failed");
          return (1);
        }
      }

    }


    Line = PL_NXT(Line);
  }      /*vom zweiten Lauf �ber die Lines*/

  return(0);
} /*von Prepair_BndPointLineRelations_fortheMesh*/




/****************************************************************************/
/*
   EvalBndPoint_Line_Informations - reads necessary line-informations for boundary points

   SYNOPSIS:
   int	EvalBndPoint_Line_Informations(LGM_MESH_INFO *theMesh)

   PARAMETERS:
   .  theMesh - mesh structure which will be filled.

   DESCRIPTION:
   reads necessary line-informations for boundary points
   dazu werden zunaechst in einem Lauf ueber alle existierenden Lines
   die Anzahl der Lines fuer jeden einzelnen BndPoint berechnet
   Danach wird mit dieser Information, abgelegt im array boundary_point_line_counter,
   die Funktion BndPoint_Line_Alloc_Mem aufgerufen, die unter der Adresse theMesh
   den notwendigen Speicher sowie die notwendigen Stukturen bereitstellt.
   Abschliessend wird dieser allokierte Speicherbereich in der Funktion
   Prepair_BndPointLineRelations_fortheMesh, mit Hilfe eines zweiten Laufs ueber alle Lines
   gefuellt.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
int     EvalBndPoint_Line_Informations(LGM_MESH_INFO *theMesh)
{
  int *boundary_point_line_counter;      /*z�hlt Anzahl Lines pro BndP*/
  PL_TYP *Line;
  int Plyln_ID;
  PL_LINE_TYP *pllyln, *pllyln_last;
  int i, der_erste_Point_der_Line, die_Point_ID;
  int Stop;

  /* Fetch Memory ...*/
  boundary_point_line_counter = (int *)GetTmpMem(theHeap,statistik[1]*sizeof(int),ANS_MarkKey);
  if ( boundary_point_line_counter == NULL )
  {
    PrintErrorMessage('E',"ansys2lgm"," ERROR: No memory for boundary_point_line_counter in EvalBndPoint_Line_Informations ");
    return(1);
  }
  memset(boundary_point_line_counter,0,(statistik[1])*sizeof(int));


  /*erster Lauf �ber die Lines: boundary_point_line_counter inkrementieren*/
  Line = EXCHNG_TYP2_ROOT_PLY(ExchangeVar_2_Pointer);
  for (Plyln_ID = 0; Plyln_ID < NMB_OF_PLINES(DomainInfo_Pointer); Plyln_ID++)
  {
    /*just for debugging - kann sp�ter gel�scht werden*/
    if (Line == NULL)
    {
      PrintErrorMessage('E',"EvalBndPoint_Line_Informations","Line-Laufpointer is NULL !!");
      return (1);
    }

    /*Laufe �ber die bereits sortierten Points jeder Line*/
    pllyln = PL_LINES(Line);
    /*die erste KnotenID der Polyline ...*/
    der_erste_Point_der_Line = LI_NDID1(PL_LINES_LINE(pllyln));
    Stop = PL_NMB_OF_POINTS(Line) - 2;
    /*der erste Point und die beiden letzten PL-Lines werden
       gesondert behandelt,
       da der zyklische Fall auftreten kann
       die letzten beiden PL-Lines ==> d.h. die letzten beiden Points
       da ja immer der zweite genommen wird
       Folgende Faelle muessen unterschieden werden koennen:
       <ab>,<bc>,<cd>,<de> von <ab>,<bc>,<cd>,<da>
       man benoetigt alo stets die letzten beiden PL-Lines gesondert.*/
    /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */

    /*weise KnotenIDs der gefundenen Polyline zu*/
    /*... laufe dazu ueber die Polylineteilstreckenzuege ...*/
    /*Der Fall der zyklischen Polyline ist von C.Tapp bisher nicht beruecksichtigt
       Ich uebergeb in diesm speziellen Fall einen Punkt zweimal
       d.h. am Anfang und am Ende,; ausserdem ist PL_NMB_OF_POINTS dann um 1 zu gross
       in Zukunft: sollte lgm-Domain so ueberarbeitet werden,
       dass die PL-Struktur hiefuer ein Flag besitzt.*/
    /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
    for(i=1; i<Stop; i++)
    {
      if(pllyln != NULL)
      {
        die_Point_ID = LI_NDID2(PL_LINES_LINE(pllyln));
        /*jetzt immer den zweiten Knoten nehmen, da ja die
           Anzahl der Points einer Polyline um 1 groesser ist
           als die Anzahl der PL-LineSegmente derselbigen ...*/

        /*boundary_point_line_counter inkrementieren,*/
        boundary_point_line_counter[die_Point_ID] ++;

        pllyln = PL_LINES_NXT(pllyln);
      }
      else
      {
        UserWrite("ERROR: in EvalBndPoint_Line_Informations: PolylineLine is missing Case1!!");
        return (1);
      }
    }            /*Ende von for .. <Stop*/

    if(pllyln == NULL)
    {
      UserWrite("ERROR: in EvalBndPoint_Line_Informations: PolylineLine is missing Case2!!");
      return (1);
    }


    /*gesonderte Betrachtung des ersten Points und der letzten beiden PL-Lines*/
    if(PL_NMB_OF_POINTS(Line) ==  2)
    {
      /*Fall1 Polyline besteht aus nur einer PL-Line*/

      /*Polyline muss mind. aus zwei versc. Points bestehen:*/
      if ( LI_NDID1(PL_LINES_LINE(pllyln)) == LI_NDID2(PL_LINES_LINE(pllyln)) )
      {
        UserWrite("ERROR: in EvalBndPoint_Line_Informations: only 1 PolylineLine with 2 identical nodes , cyclic ");
        return (1);
      }

      /*boundary_point_line_counter inkrementieren,*/
      boundary_point_line_counter[LI_NDID1(PL_LINES_LINE(pllyln))] ++;                  /*erster*/
      boundary_point_line_counter[LI_NDID2(PL_LINES_LINE(pllyln))] ++;                  /*letzter*/

    }
    else
    {
      /* pllyln ist nun die vorletzte PL-Line*/

      /*die letzte PL-Line:*/
      pllyln_last = PL_LINES_NXT(pllyln);
      if(pllyln_last == NULL)
      {
        UserWrite("ERROR: in EvalBndPoint_Line_Informations: PolylineLine is missing Case3 pllyln_last!!");
        return (1);
      }

      /*TOCHECK: Liegt eine zyklische Polyline vor:*/
      /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */
      if( LI_NDID2(PL_LINES_LINE(pllyln_last)) == der_erste_Point_der_Line)
      {
        /*Fall2: Spezialfall : <zyklische Polyline>*/
        /* concerning cyclic Polylines see also file <zyklische_Polylines> on Dirk Mac */

        /*boundary_point_line_counter*/
        boundary_point_line_counter[der_erste_Point_der_Line] ++;                        /*erster*/
        boundary_point_line_counter[LI_NDID2(PL_LINES_LINE(pllyln))] ++;                        /*letzter*/

      }
      else
      {
        /*Fall3: Normalfall*/

        /*boundary_point_line_counter */

        boundary_point_line_counter[der_erste_Point_der_Line] ++;                        /*erster*/
        boundary_point_line_counter[LI_NDID2(PL_LINES_LINE(pllyln_last))] ++;                        /*letzter*/
        boundary_point_line_counter[LI_NDID2(PL_LINES_LINE(pllyln))] ++;                        /*vorletzter*/
      }

    }


    Line = PL_NXT(Line);
  }      /*vom ersten Lauf �ber die Lines*/


  /* Nun kann mit den Informationen in boundary_point_line_counter Speicher geholt werden*/
  if (BndPoint_Line_Alloc_Mem(theMesh, boundary_point_line_counter) != 0)
  {
    PrintErrorMessage('E',"EvalBndPoint_Line_Informations->BndPoint_Line_Alloc_Mem","execution failed");
    return (1);
  }

  /* Nun k�nnen die neuen Komponenten von Mesh bzw. lgm_mesh_info gef�llt werden: */
  if (Prepair_BndPointLineRelations_fortheMesh(theMesh) != 0)
  {
    PrintErrorMessage('E',"EvalBndPoint_Line_Informations->Prepair_BndPointLineRelations_fortheMesh","execution failed");
    return (1);
  }

  return(0);
} /*Ende von EvalBndPoint_Line_Informations*/




/****************************************************************************/
/*
   LGM_ANSYS_ReadMesh - reads mesh from ANSYS-Output and returns corret Mesh for UG

   SYNOPSIS:
   int LGM_ANSYS_ReadMesh (HEAP *theHeap, LGM_MESH_INFO *theMesh, int MarkKey)

   PARAMETERS:
   .  theHeap - heappointer.

   DESCRIPTION:
   lreads all coordinates

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/
/* DIRKS NEU : theHeap gibts doch schon seit LGM_ANSYS_ReadDomain
               und muss hier nicht nochmal uebergeben werden.*/
/*int LGM_ANSYS_ReadMesh (HEAP *theHeap, LGM_MESH_INFO *theMesh)*/
/* alte Version*/
int NS_DIM_PREFIX LGM_ANSYS_ReadMesh (char *name, HEAP *Heappointer, LGM_MESH_INFO *theMesh, int MarkKey) /* DIRKS NEU MarkKey*/
{
  SD_TYP *sbd;
  int i,TetraederelementID, SbdName,elem_lf;

  /* DIRKS NEU :*/
  theHeap = Heappointer;
  ANS_MarkKey = MarkKey;

  /* fill PositionInformations of Inner- and BndPoints of theMesh */
  if (FillPositionInformations(theMesh) != 0)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh/FillPositionInformations","execution failed");
    return (1);
  }

  /*In el_array benoetigen alle Tetraederseiten, die auf einer Boundary liegen, den Eintrag -1;
     bzw. genauer -1 * Vorkammastellle der zugehoerigen Lastzahl
     dies wird moeglich durch einen Lauf uber alle SFEs , dies geschieht in ReadAnsysFile, s.o.
     Wichtig, da sonst keine Unterscheidung  zwischen sich beruehrenden Subdomains
     in SearchAllTetrahedronsOfThisSbd moeglich. s.o. unter ug_side_offset */

  /* find neighbours */
  if (FindElNeighbours(statistik[6]) != 0)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh/FindElNeighbours","execution failed");
    return (1);
  }

  /*Initialsierung des notwendigen KontrollFlagFeldes zur Steuerung der Elem_Sbd_Zgh. */
  el_besucht_array = (INT *)GetTmpMem(theHeap,(statistik[6]+1)*sizeof(INT), ANS_MarkKey);
  if ( el_besucht_array == NULL )
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for el_besucht_array !!!");
    return(1);
  }
  memset(el_besucht_array,0,(statistik[6]+1)*sizeof(INT));


  /*the number of Subdomains:*/
  theMesh->nSubDomains = NMB_OF_SBDMS(DomainInfo_Pointer);
  /*allocate array for nmbofsides per sbd*/
  if((theMesh->nSides = (INT *)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->nSides !!!");
    return(1);
  }
  /*allocate array for nmbofelems per sbd*/
  if((theMesh->nElements = (INT *)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->nElements !!!");
    return(1);
  }
  /*allocate array for SideCorners*/
  if((theMesh->Side_corners = (INT **)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->nSides !!!");
    return(1);
  }
  /*allocate array for SideCornerIds*/
  if((theMesh->Side_corner_ids = (INT ***)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT**), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->Side_corner_ids !!!");
    return(1);
  }
  /*allocate array for Element_corners*/
  if((theMesh->Element_corners = (INT **)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->Element_corners !!!");
    return(1);
  }
  /*allocate array for Element_SideOnBnd*/
  if((theMesh->Element_SideOnBnd = (INT **)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT*), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->Element_SideOnBnd !!!");
    return(1);
  }
  /*allocate array for Element_corner_ids*/
  if((theMesh->Element_corner_ids = (INT ***)GetTmpMem(theHeap,(1 + NMB_OF_SBDMS(DomainInfo_Pointer))*sizeof(INT**), ANS_MarkKey)) == NULL)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR: No memory for theMesh->Element_corner_ids !!!");
    return(1);
  }
  /*allocate array for nbElements*/
  theMesh->nbElements = NULL;       /*not yet*/

  /* for all subdomains ...*/
  sbd = EXCHNG_TYP2_ROOT_SBD(ExchangeVar_2_Pointer);
  for (i=1; i <= NMB_OF_SBDMS(DomainInfo_Pointer); i++)
  {
    if(sbd != NULL)
    {
      SbdName = SD_NAME(sbd);
      nmbOfTetrhdrOfThisSbd = 0;                   /*Zaehler fuer Sbd->Tetraederanzahl*/
      nmbOfSidesOfThisSbd = 0;                   /*Zaehler fuer Sbd->Sidesanzahl*/


      if((TetraederelementID = FetchATetrahedronOfThisSbd(sbd)) == -1)                   /* -1 bei Fehler !!*/
      {
        UserWrite("ERROR: in LGM_ANSYS_ReadMesh: no tetrahedron out of FetchATetrahedronOfThisSbd");
        return (1);
      };

      /*erster Tetraeder der Subdomain sbd ist somit gefunden
         er wird nun im KontrollFlagFeld eingetragen */
      el_besucht_array[TetraederelementID] = SbdName;
      nmbOfTetrhdrOfThisSbd ++;
      /* Update of nmbOfSidesOfThisSbd des ersten Mustertetraeders - dieser liegt
         ja auf jeden Fall auf dem Rand*/
      /*reicht auch noch in FillSubdomainInformations */
      /*stelle = TetraederelementID * 8 + 4;
         for(ll = 0; ll < 4; ll++)
         {
              if(el_array[stelle] < 0) *//*wenn das Element bzgl dieser Seite an eine Boundary angrenzt ... */
      /*	{
                      nmbOfSidesOfThisSbd ++;
              }
              stelle++;
         }*/

      /* GOON HERE MONDAY die zwei Funktionen SearchAllTetrahedronsOfThisSbd und
         FillSubdomainInformations sind noch zu schreiben, ferner sollten Proben realisiert werden in der Art:
         Am schluss muss el_besucht_array ueberall ausser bei 0 genau einmal besucht worden sein
         In FillSubdomainInformations muss die AnzahlSidespersubdomain geprueft werden.
         In SearchAllTetrahedronsOfThisSbd muss die Sbd-SideAnzahl noch nicht gezaehlt werden
         Spaeter kann sie durch Addition der #Sbds->surfaces->'triangle gewaehlt werden*/

      /*Suche in einer rek. Funktion mit Hilfe der Nachbarschaftsbeziehungen
         alle Tetraeder die der aktuellen Subdomain angehoeren. Dabei muessen
         nmbOfTetrhdrOfThisSbd, nmbOfSidesOfThisSbd und el_besucht_array stets aktualisiert werden*/


      if (SearchAllTetrahedronsOfThisSbd(TetraederelementID,SbdName) == 1)
      {
        PrintErrorMessage('E',"LGM_ANSYS_ReadMesh"," ERROR out of SearchAllTetrahedronsOfThisSbd, = rekursive Funktion. !");
        return(1);
      }


      (theMesh->nElements)[i] = nmbOfTetrhdrOfThisSbd;

      if (FillSubdomainInformations(theMesh,SbdName,i) != 0)                   /*uses nmbOfTetrhdrOfThisSbd and el_besucht_array*/
      {
        PrintErrorMessage('E',"LGM_ANSYS_ReadMesh/FillSubdomainInformations","execution failed");
        return (1);
      }
    }
    else
    {
      UserWrite("ERROR: in LGM_ANSYS_ReadMesh: Subdomain is missing !!");
      return (1);
    }

    /*weiter mit naechster Subdomain*/
    sbd = SD_NEXT(sbd);
  }

  /*Probe:*/
  for(elem_lf=1; elem_lf<=statistik[6]; elem_lf++)
  {
    /*Probe : wenn ein Element keiner Subdomain zugeordnet werden konnte . . .*/
    if(el_besucht_array[elem_lf] == 0)
    {
      PrintErrorMessage('E',"LGM_ANSYS_ReadMesh","el_besucht_array nicht vollstaendig gefuellt");
      return(1);
    }
  }

  if (EvalBndPointInformations(theMesh) != 0)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh/EvalBndPointInformations","execution failed");
    return (1);
  }

  if (EvalBndPoint_Line_Informations(theMesh) != 0)
  {
    PrintErrorMessage('E',"LGM_ANSYS_ReadMesh/EvalBndPoint_Line_Informations","execution failed");
    return (1);
  }


  return(0);
  /*    UserWrite("ERROR: in LGM_ANSYS_ReadMesh: return1, da not finshed yet\n");
      return (1);*/
  /* TODO  BoundaryPoints !!! DIRKS NEU return 0  !!!*/
}
