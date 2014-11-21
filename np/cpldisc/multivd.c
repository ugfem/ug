// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* multivd.c - implementation of a command that creates a vector descriptor
 * of several ones.
 * History:
 * Jun. 22, 2003 - created
 * D. Logashenko
 **
 * USAGE INSTRUCTIONS:
 * combinevd <name> $parts <vd1> <vd2> ...;
 * This command creates a vector data descriptor <name> that refers
 * to the components of vector descriptors <vd1>, <vd2>, ... in the
 * order in which these descriptors are listed in the command.
 */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* The UG namespaces: */
#include "namespace.h"
USING_UG_NAMESPACES

#include "gm.h"        /* for data structure               */
#include "misc.h"      /* for MIN, MAX, PI, ...            */
#include "ugdevices.h" /* for UserWrite, PrintErrorMessage */
#include "general.h"
#include "rm.h"
#include "ugstruct.h"
#include "sm.h"        /* for sparse block structures      */
#include "scan.h"
#include "cmdline.h"   /* for shell commands               */
#include "commands.h"
#include "udm.h"

/* The maximum number of the vector descriptors to combine: */
#define MAX_COMBINED 64

/* The name of the command: */
static char command_name [] = "combinevd";

/* read_vec_descs - reads the vector descriptors from the script.
 * Returns the number of the read vec. desc.'s
 * if OK, zero on an error:
 */
static INT read_vec_descs
(
  MULTIGRID * mg, /* the multigrid to consider */
  INT argc, char * * argv, /* the script specifications */
  VECDATA_DESC * vd_array [] /* to save the vec. descs. */
)
{
  char * str;
  char name [64];
  INT i;
  int n;

  /* Look for the specification: */
  i = 0;
  while (sscanf (argv [i], "%s %n", name, &n) != 1
         || strcmp (name, "parts") != 0)
    if (++i >= argc)
    {
      /* No specification found! */
      PrintErrorMessage ('E', "combinevd",
                         "No specification of the parts found!");
      return 0;
    };
  str = argv [i] + n;

  /* Parse the line: */
  for (i = 0; sscanf (str, "%s %n", name, &n) == 1; i++, str += n)
    if ((vd_array [i] = GetVecDataDescByName (mg, name)) == NULL)
    {
      PrintErrorMessageF ('E', "combinevd",
                          "Cannot get the vec. desc. '%s'", name);
      return 0;
    };

  /* Check whether the tail is empty: */
  if (str [0])
  {
    PrintErrorMessageF ('E', "combinevd", "Cannot parse '%s'", str);
    return 0;
  };

  return i;
}

/* create_combined_vd - creates the combined vd of the specified parts.
 * Besides assignes the names to the components of the vd. Returns the
 * pointer to the created vd if OK, or NULL if failed.
 */
static VECDATA_DESC * create_combined_vd
(
  MULTIGRID * mg, /* the multigrid to consider */
  char * name, /* the name for the new vd */
  INT n, /* number of the vec. desc.'s to combine */
  VECDATA_DESC * vd_array [] /* the specified parts */
)
{
  VECDATA_DESC * vd;
  char * cmps;
  INT i, k, len;

  /* Create the vec. desc.: */
  if ((vd = CombineVecDesc (mg, name, (const VECDATA_DESC * *) vd_array, n))
      == NULL)
  {
    PrintErrorMessageF ('E', "combinevd",
                        "Cannot create the combined vd '%s'", name);
    return NULL;
  };

  /* Assign the names: */
  k = 0;
  for (i = 0; i < n; i++)
  {
    len = strlen (cmps = VM_COMP_NAMEPTR (vd_array [i]));
    strcpy (VM_COMP_NAMEPTR (vd) + k, cmps);
    k += len;
  };

  /* Lock the vd: */
  if (LockVD (mg, vd))
    return NULL;

  /* Done: */
  return vd;
}

/* combinevd_command - the function that implements the command: */
static INT combinevd_command (INT argc, char * * argv)
{
  MULTIGRID * mg;
  char name [NAMESIZE];
  VECDATA_DESC * vd_array [MAX_COMBINED];
  INT n;

  if ((mg = GetCurrentMultigrid ()) == NULL)
  {
    PrintErrorMessage ('E', "combinevd", "No multigrid created!");
    return CMDERRORCODE;
  };

  if (ReadArgvChar (command_name, name, argc, argv) || name [0] == '\0')
  {
    PrintErrorMessage ('E', "combinevd", "No name for the vd specified");
    return CMDERRORCODE;
  };

  if ((n = read_vec_descs (mg, argc, argv, vd_array)) == 0)
    return CMDERRORCODE;

  if (create_combined_vd (mg, name, n, vd_array) == NULL)
    return CMDERRORCODE;

  return OKCODE;
}

/* The installer: */
INT NS_DIM_PREFIX Install_CombineVD ()
{
  return CreateCommand (command_name, combinevd_command) == NULL;
}

/* End of File */
