// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* coupleddiscinit.c - initialization routine for the coupled discretization
 * numprocs and commands.
 * History:
 * Jun. 22, 2003 - created
 * D. Logashenko
 */

#include <config.h>
#include <stdio.h>
#include <math.h>

/* The UG namespaces: */
#include "namespace.h"
USING_UG_NAMESPACES

#include "gm.h"        /* for data structure               */
#include "evm.h"       /* for data structure               */
#include "shapes.h"    /* for data structure               */
#include "misc.h"      /* for MIN, MAX, PI, ...            */
#include "ugdevices.h" /* for UserWrite, PrintErrorMessage */
#include "general.h"
#include "rm.h"
#include "ugstruct.h"

/* The mudules to install: */
#include "globdisc.h"
#include "sparsepattern.h"
#include "multivd.h"

/* Own header: */
#include "coupleddiscinit.h"

INT NS_DIM_PREFIX Install_Coupled_Discretization ()
{
  if (Install_GlobalDisc ())
  {
    printf ("Error in main: Install_GlobalDisc returned error.\n");
    return 1;
  };

  if (Install_Sparsity_Descriptions ())
  {
    printf ("Error in main: Install_Sparsity_Descriptions returned error.\n");
    return 1;
  };

  if (Install_CombineVD ())
  {
    printf ("Error in main: Install_CombineVD returned error.\n");
    return 1;
  };

  return 0;
}

/* End of File */
