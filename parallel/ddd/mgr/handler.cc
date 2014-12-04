/****************************************************************************/
/*                                                                          */
/* File:      handler.ct                                                    */
/*                                                                          */
/* Purpose:   template routines for handler definition functions            */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   970212 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* derived macros                                                           */
/*                                                                          */
/****************************************************************************/

#define _CAT(a,b)   a ## b
#define CAT(a,b)    _CAT(a,b)


#ifdef CPP_FRONTEND
#define HDLR_FUNCNAME CAT(DDD_Library::SetHandler,HDLR_NAME)
#else
#define HDLR_FUNCNAME CAT(DDD_SetHandler,HDLR_NAME)
#endif

#define HDLR_TYPENAME CAT(Handler,HDLR_NAME)
#define HDLR_VARNAME  CAT(handler,HDLR_NAME)


/****************************************************************************/
/*                                                                          */
/* auxiliary macros                                                         */
/*                                                                          */
/****************************************************************************/

#ifndef HDLR_AUX_MACROS
#define HDLR_AUX_MACROS

#define CONCAT(txt,fct)    txt FUNCNAME_STR(fct)
#define FUNCNAME_STR(f)   #f

#endif


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
void HDLR_FUNCNAME (DDD_TYPE type_id, HDLR_TYPENAME funcptr)
{
#endif

#ifdef F_FRONTEND
void HDLR_FUNCNAME (DDD_TYPE *fid, HDLR_TYPENAME funcptr)
{
	DDD_TYPE type_id = *fid;
#endif

	TYPE_DESC *desc = &(theTypeDefs[type_id]);

	if (desc->mode != DDD_TYPE_DEFINED)
	{
		sprintf(cBuffer, CONCAT("undefined DDD_TYPE %d in ",HDLR_FUNCNAME),
			type_id);
		DDD_PrintError('E', 9916, cBuffer);
		HARD_EXIT;
	}

	desc-> HDLR_VARNAME = funcptr;
}



/****************************************************************************/

/* undefs for all input macros and defines */

#undef HDLR_NAME
#undef HDLR_FUNCNAME
#undef HDLR_TYPENAME
#undef HDLR_VARNAME

#undef _CAT
#undef CAT


/* undefs for derived macros */



/****************************************************************************/

