// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ifcmds.c                                                      */
/*                                                                          */
/* Purpose:   routines concerning interfaces between processors             */
/*            part 2: usage of DDD interfaces                               */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/03/03 kb  complete rewrite                                 */
/*            94/09/12 kb  IFExchange & IFOneway rewrite, two bugs fixed    */
/*            94/09/21 kb  created from if.c                                */
/*            95/01/13 kb  added range functionality                        */
/*            95/07/26 kb  overlapping of gather/scatter and communication  */
/*            96/01/24 kb  added use of object shortcut tables              */
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

/* standard C library */
#include <stdlib.h>
#include <stdio.h>

#include "dddi.h"
#include "if.h"



/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/*
        description of input defines for ifcmd.ct

        IF_NAME:   name of interface function
        IF_CBR:    call by reference, otherwise call by address (CPP_FRONTEND)

        TODO: more doc!
 */


/* with normal arguments for gather/scatter */

#define IF_NAME Exchange
#define IF_EXCHANGE
#include "ifcmd.ct"

#ifdef CPP_FRONTEND
#define IF_NAME Exchange
#define IF_EXCHANGE
#define IF_CBR
#include "ifcmd.ct"
#endif



#define IF_NAME Oneway
#define IF_ONEWAY
#include "ifcmd.ct"

#ifdef CPP_FRONTEND
#define IF_NAME Oneway
#define IF_ONEWAY
#define IF_CBR
#include "ifcmd.ct"
#endif





#ifndef CPP_FRONTEND

#define IF_NAME ExecLocal
#define IF_EXECLOCAL
#include "ifcmd.ct"


#define IF_NAME AExchange
#define IF_EXCHANGE
#define IF_WITH_ATTR
#include "ifcmd.ct"

#define IF_NAME AOneway
#define IF_ONEWAY
#define IF_WITH_ATTR
#include "ifcmd.ct"

#define IF_NAME AExecLocal
#define IF_EXECLOCAL
#define IF_WITH_ATTR
#include "ifcmd.ct"


/* with extended arguments for gather/scatter */

#define IF_NAME ExchangeX
#define IF_EXCHANGE
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_NAME OnewayX
#define IF_ONEWAY
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_NAME ExecLocalX
#define IF_EXECLOCAL
#define IF_WITH_XARGS
#include "ifcmd.ct"


#define IF_NAME AExchangeX
#define IF_EXCHANGE
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_NAME AOnewayX
#define IF_ONEWAY
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_NAME AExecLocalX
#define IF_EXECLOCAL
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"

#endif

/****************************************************************************/


/*
        description of input defines for ifstd.ct

        IF_NAME:   name of interface function

        These includes generate internal functions for communication
        on the STD_INTERFACE (IF0). The Gather/Scatter-functions for
        the STD_INTERFACE will get the DDD_HDR as a parameter, not the
        DDD_OBJ.

        TODO: more doc!
 */


/* with normal arguments for gather/scatter */

#define IF_NAME Exchange
#define IF_EXCHANGE
#include "ifstd.ct"

#define IF_NAME ExecLocal
#define IF_EXECLOCAL
#include "ifstd.ct"


/* with extended arguments for gather/scatter */

#define IF_NAME ExchangeX
#define IF_EXCHANGE
#define IF_WITH_XARGS
#include "ifstd.ct"

#define IF_NAME ExecLocalX
#define IF_EXECLOCAL
#define IF_WITH_XARGS
#include "ifstd.ct"


/****************************************************************************/
