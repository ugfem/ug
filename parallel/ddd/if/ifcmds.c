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


/* with normal arguments for gather/scatter */

#define IF_FUNCNAME DDD_IFExchange
#define IF_EXCHANGE
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFOneway
#define IF_ONEWAY
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFExecLocal
#define IF_EXECLOCAL
#include "ifcmd.ct"


#define IF_FUNCNAME DDD_IFAExchange
#define IF_EXCHANGE
#define IF_WITH_ATTR
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFAOneway
#define IF_ONEWAY
#define IF_WITH_ATTR
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFAExecLocal
#define IF_EXECLOCAL
#define IF_WITH_ATTR
#include "ifcmd.ct"


/* with extended arguments for gather/scatter */

#define IF_FUNCNAME DDD_IFExchangeX
#define IF_EXCHANGE
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFOnewayX
#define IF_ONEWAY
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFExecLocalX
#define IF_EXECLOCAL
#define IF_WITH_XARGS
#include "ifcmd.ct"


#define IF_FUNCNAME DDD_IFAExchangeX
#define IF_EXCHANGE
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFAOnewayX
#define IF_ONEWAY
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"

#define IF_FUNCNAME DDD_IFAExecLocalX
#define IF_EXECLOCAL
#define IF_WITH_ATTR
#define IF_WITH_XARGS
#include "ifcmd.ct"


/****************************************************************************/
