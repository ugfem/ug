%{
/****************************************************************************/
/*                                                                          */
/* File:      grid.y                                                        */
/*                                                                          */
/* Purpose:   parser for gridfiles                                          */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen                             */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: ug@ica3.uni-stuttgart.de                            */
/*                                                                          */
/* History:   19.2.98 begin,                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "ng.h"

#define alloca(p)		malloc(p)
#define SP_COPY(d,s)    {(d)->surf_id=(s)->surf_id; \
						 (d)->tri_id=(s)->tri_id; \
						 (d)->local[0]=(s)->local[0]; \
						 (d)->local[1]=(s)->local[1];}
#define LP_COPY(d,s)	{(d)->line_id=(s)->line_id; (d)->local=(s)->local;}

static LINE_POSITION LinePos;
static SURFACE_POSITION SurfPos;
static BND_NODE BndNode;
static INNER_NODE InnerNode;
static ELEM_FACE ElemFace;
static ELEMENT Elem;



%}

%union 
{
	/* put RCS string here in order to get it into yacc-generated header file
	static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 1.4
	 1998/02/20 16:58:46 birken Exp $",UG_RCS_STRING);
	*/


	/* transfer lex->yacc */
    double dval;
    long ival;

	LINE_POSITION *lp;
	SURFACE_POSITION *sp;
	BND_NODE *bs;
	INNER_NODE *in;
	ELEM_FACE *ef;
	ELEMENT *el;
}

%token <dval> DOUBLE_VALUE
%token <ival> INT_VALUE
%token INODE BNODE SURFACE LINE ELEM FACE TEND

%type<el> ElemSpec InnerElemSpec
%type<ef> ElemFace
%type<in> InnerNode InnerNodeList
%type<bs> BndSpec BndNode BndNodeList Grid
%type<sp> SurfacePosition
%type<lp> LinePosition
%type<dval> Coord
%type<ival> Id

%%
Grid:
	BndNodeList InnerNodeList ElemList
	| BndNodeList ElemList
	;
ElemList:
    Elem
    | ElemList Elem
    ;
Elem:
	ELEM 											{
														Elem.n_c=Elem.n_f=0;
													}
	ElemSpec										{
														if (PutElement(&Elem)) YYABORT;
													} 
	TEND
	;
ElemSpec:
	InnerElemSpec
	| InnerElemSpec ElemFace
	| InnerElemSpec ElemFace ElemFace
	| InnerElemSpec ElemFace ElemFace ElemFace
	| InnerElemSpec ElemFace ElemFace ElemFace ElemFace
	| InnerElemSpec ElemFace ElemFace ElemFace ElemFace ElemFace
	| InnerElemSpec ElemFace ElemFace ElemFace ElemFace ElemFace ElemFace
	;
InnerElemSpec:
    Id Id Id Id Id									{
														Elem.subdom=(int)$1;
														Elem.c_id[0]=(int)$2;
														Elem.c_id[1]=(int)$3;
														Elem.c_id[2]=(int)$4;
														Elem.c_id[3]=(int)$5;
														Elem.n_c=4;
														$$=&Elem;
													}
	| Id Id Id Id Id Id								{
														Elem.subdom=(int)$1;
														Elem.c_id[0]=(int)$2;
														Elem.c_id[1]=(int)$3;
														Elem.c_id[2]=(int)$4;
														Elem.c_id[3]=(int)$5;
														Elem.c_id[4]=(int)$6;
														Elem.n_c=5;
														$$=&Elem;
													}
    | Id Id Id Id Id Id Id                          {
														Elem.subdom=(int)$1;
                                                        Elem.c_id[0]=(int)$2;
                                                        Elem.c_id[1]=(int)$3;
                                                        Elem.c_id[2]=(int)$4;
                                                        Elem.c_id[3]=(int)$5;
                                                        Elem.c_id[4]=(int)$6;
                                                        Elem.c_id[5]=(int)$7;
                                                        Elem.n_c=6;
                                                        $$=&Elem;
                                                    }
    | Id Id Id Id Id Id Id Id Id                    {
														Elem.subdom=(int)$1;
                                                        Elem.c_id[0]=(int)$2;
                                                        Elem.c_id[1]=(int)$3;
                                                        Elem.c_id[2]=(int)$4;
                                                        Elem.c_id[3]=(int)$5;
                                                        Elem.c_id[4]=(int)$6;
                                                        Elem.c_id[5]=(int)$7;
                                                        Elem.c_id[6]=(int)$8;
                                                        Elem.c_id[7]=(int)$9;
                                                        Elem.n_c=8;
                                                        $$=&Elem;
                                                    }
    ;
ElemFace:
    FACE Id Id Id									{
														Elem.face[Elem.n_f].c_id[0]=(int)$2;
														Elem.face[Elem.n_f].c_id[1]=(int)$3;
														Elem.face[Elem.n_f].c_id[2]=(int)$4;
														Elem.face[Elem.n_f].n_c=3;
														$$=&(Elem.face[Elem.n_f]);
														Elem.n_f++;
													}
    | FACE Id Id Id Id								{
														Elem.face[Elem.n_f].c_id[0]=(int)$2;
														Elem.face[Elem.n_f].c_id[1]=(int)$3;
														Elem.face[Elem.n_f].c_id[2]=(int)$4;
														Elem.face[Elem.n_f].c_id[3]=(int)$5;
														Elem.face[Elem.n_f].n_c=4;
														$$=&Elem.face[Elem.n_f];
														Elem.n_f++;
													}
    ;
InnerNodeList:
    InnerNode
    | InnerNodeList InnerNode
    ;
InnerNode:
    INODE Coord Coord Coord TEND                    {
														InnerNode.global[0]=$2;
														InnerNode.global[1]=$3;
														InnerNode.global[2]=$4;
														$$=&InnerNode;
														PutInnerNode(&InnerNode);
													}
    ;
BndNodeList:
    BndNode
    | BndNodeList BndNode
    ;
BndNode:
    BNODE Coord Coord Coord							{
														BndNode.n_lp=BndNode.n_sp=0;
														BndNode.global[0]=$2;
														BndNode.global[1]=$3;
														BndNode.global[2]=$4;
														$$=&BndNode;
													}
	BndSpec 										{PutBndNode(&BndNode);}
	TEND			
    ;
BndSpec:
    SurfacePosition									{
														SP_COPY(&(BndNode.sp[BndNode.n_sp]),$1);
														BndNode.n_sp++;
														$$=&BndNode;
													}
    | LinePosition									{
														LP_COPY(&(BndNode.lp[BndNode.n_lp]),$1);
														BndNode.n_lp++;
														$$=&BndNode;
													}
    | BndSpec SurfacePosition                       {
														SP_COPY(&(BndNode.sp[BndNode.n_sp]),$2);
														BndNode.n_sp++;
														$$=&BndNode;
													}
    | BndSpec LinePosition                          {
														LP_COPY(&(BndNode.lp[BndNode.n_lp]),$2);
														BndNode.n_lp++;
														$$=&BndNode;
													}
    ;
SurfacePosition:
	SURFACE Id Id Coord Coord                       {
														SurfPos.surf_id=(int)$2;
														SurfPos.tri_id=(int)$3;
														SurfPos.local[0]=(float)$4;
														SurfPos.local[1]=(float)$5;
														$$=&SurfPos;
													}
	;
LinePosition:
	LINE Id Coord 									{
														LinePos.line_id=(int)$2; 
														LinePos.local=(float)$3; 
														$$=&LinePos;
													}	
	;
Coord: 
	INT_VALUE										{$$=(double)$1;} 
	| DOUBLE_VALUE
	;
Id: 
	INT_VALUE										{$$=$1;}                                                            
	;


%%


yywrap (char *s)
{
    return (1);
}

yyerror (char *s)
{
	int line;
	char text[128];

	NP_Error(&line,text);
    NG_Print(ERROR_PREFIX "'%s', line %d\n",text,line);
    ngbreak();
}









