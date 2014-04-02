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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "ng2d.h"
#include "ngin-lex.h"

#ifndef alloca  
#define alloca(p)		malloc(p)
#endif

#define SP_COPY(d,s)    {(d)->surf_id=(s)->surf_id; \
						 (d)->tri_id=(s)->tri_id; \
						 (d)->local[0]=(s)->local[0]; \
						 (d)->local[1]=(s)->local[1];}
#define LP_COPY(d,s)	{(d)->line_id=(s)->line_id; (d)->local=(s)->local;}

#include "namespace.h"

USING_UG_NAMESPACE
USING_UGDIM_NAMESPACE
  
static LINE_POSITION LinePos;
static BND_NODE BndNode;
static INNER_NODE InnerNode;
static ELEM_SIDE ElemSide;
static NG_ELEMENT Elem;

 /* forward declare my own function (referenced by automatic parser) */ 
 int ng2derror (char *s);


%}

%union 
{
	/* put RCS string here in order to get it into yacc-generated header file
	static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 1.1
	 1998/02/20 16:58:46 birken Exp $",UG_RCS_STRING);
	*/


	/* transfer lex->yacc */
    double dval;
    long ival;

	NS_DIM_PREFIX LINE_POSITION *lp;
	NS_DIM_PREFIX BND_NODE *bs;
	NS_DIM_PREFIX INNER_NODE *in;
	NS_DIM_PREFIX ELEM_SIDE *es;
	NS_DIM_PREFIX NG_ELEMENT *el;
}

%token <dval> DOUBLE_VALUE
%token <ival> INT_VALUE
%token INODE BNODE LINE ELEM SIDE TEND

%type<el> ElemSpec InnerElemSpec
%type<es> ElemSide
%type<in> InnerNode InnerNodeList
%type<bs> BndSpec BndNode BndNodeList Grid
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
														Elem.n_c=Elem.n_s=0;
													}
	ElemSpec										{
														if (PutElement(&Elem)) YYABORT;
													} 
	TEND
	;
ElemSpec:
	InnerElemSpec
	| InnerElemSpec ElemSide
	| InnerElemSpec ElemSide ElemSide
	| InnerElemSpec ElemSide ElemSide ElemSide
	| InnerElemSpec ElemSide ElemSide ElemSide ElemSide
	;
InnerElemSpec:
    Id Id Id Id 									{
														Elem.subdom=(int)$1;
														Elem.c_id[0]=(int)$2;
														Elem.c_id[1]=(int)$3;
														Elem.c_id[2]=(int)$4;
														Elem.n_c=3;
														$$=&Elem;
													}
	| Id Id Id Id Id 								{
														Elem.subdom=(int)$1;
														Elem.c_id[0]=(int)$2;
														Elem.c_id[1]=(int)$3;
														Elem.c_id[2]=(int)$4;
														Elem.c_id[3]=(int)$5;
														Elem.n_c=4;
														$$=&Elem;
													}
    ;
ElemSide:
    SIDE Id Id 	    								{
														Elem.side[Elem.n_s].c_id[0]=(int)$2;
														Elem.side[Elem.n_s].c_id[1]=(int)$3;
														$$=&(Elem.side[Elem.n_s]);
														Elem.n_s++;
													}
    ;
InnerNodeList:
    InnerNode
    | InnerNodeList InnerNode
    ;
InnerNode:
    INODE Coord Coord TEND                         {
														InnerNode.global[0]=$2;
														InnerNode.global[1]=$3;
														$$=&InnerNode;
														PutInnerNode(&InnerNode);
													}
    ;
BndNodeList:
    BndNode
    | BndNodeList BndNode
    ;
BndNode:
    BNODE Coord Coord 								{
														BndNode.n_lp=0;
														BndNode.global[0]=$2;
														BndNode.global[1]=$3;
														$<bs>$=&BndNode;
													}
	BndSpec 										{PutBndNode(&BndNode);}
	TEND    										{}
    ;
BndSpec:
    LinePosition									{
														LP_COPY(&(BndNode.lp[BndNode.n_lp]),$1);
														BndNode.n_lp++;
														$$=&BndNode;
													}
    | BndSpec LinePosition                          {
														LP_COPY(&(BndNode.lp[BndNode.n_lp]),$2);
														BndNode.n_lp++;
														$$=&BndNode;
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

#ifdef __cplusplus
extern "C" 
{
#endif
int ng2dwrap (void )
{
    return (1);
}
#ifdef __cplusplus
}
#endif

int ng2derror (char *s)
{
	int line;
	char text[128];

	NP2d_Error(&line,text);
    NG_Print(ERROR_PREFIX "'%s', line %d\n",text,line);
    ngbreak();
	return 0;
}
