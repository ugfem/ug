%{
/****************************************************************************/
/*                                                                          */
/* File:      grid.l                                                        */
/*                                                                          */
/* Purpose:   lexer for gridfiles                                           */
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

#include <stdlib.h>
#include "ng2d.h"
#include "ngin-yacc.h"

static int noline=1;

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);


%}

COMMENT             #.*
KEY_DOUBLE 			-?(([0-9]+)|([0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?)
KEY_INT             [0-9]+

KEY_INODE           I
KEY_BNODE           B
KEY_LINE         	L
KEY_ELEM         	E
KEY_SIDE         	S
KEY_TEND        	;

%%

[ \t]+              ;
[\n]                {noline++;}
{COMMENT}           ;
{KEY_INT}			{yylval.ival=atol((const char *)yytext); return (INT_VALUE);}
{KEY_DOUBLE}		{yylval.dval=strtod((const char *)yytext,NULL); return (DOUBLE_VALUE);}
{KEY_INODE}			{return (INODE);}
{KEY_BNODE}			{return (BNODE);}
{KEY_LINE}			{return (LINE);}
{KEY_ELEM}			{return (ELEM);}
{KEY_SIDE}			{return (SIDE);}
{KEY_TEND}			{return (TEND);}
.                   {yyerror(NULL);}

%%

NP_Error (int *line, char *text)
{
	*line=noline;
	strcpy(text,yytext);
	return;
}


