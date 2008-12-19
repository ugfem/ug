// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* sparsepattern.c - commands for creating a sparsity pattern from
 * a special description.
 * History:
 * May 30, 2003 - created
 * D. Logashenko
 **
 * Description of the sparsity pattern:
 * The description consists of variables places in the shell structure
 * ':sparsity_pattern'. These variables are:
 * size (an integer) - the size of the matrix blocks,
 * diag_line<i>, offd_line<i>  (1 <= i <= size) - descriptions of the
 *  sparsity patterns of the i-th line of diagonal and offdiagonal blocks,
 *  respectively. This description has the following syntax:
 *  {j | j1 - j2}
 *  (these specifications should be separated by space characters).
 *  j, j1 and j2 are positive integers. The specification "j" means a
 *  (potentially) nonzero entry in the position (i, j) in the matrix
 *  whereas "j1 - j2" means that positions (i, j) with j from j1 to j2
 *  are to be marked as nonzero.
 **
 * Hints:
 * 1. Create a shell structure (in :SparseFormats) for the variables
 * Dnn and Tnn that described the UG sparsity patterns. The structure can
 * be created using the UG command 'ms'. Go into this structure (using 'cs').
 * Then call the command 'create_sparsity_pattern' (supplied by this module)
 * to convert the sparsity descriptions described above to the UG
 * sparsity patterns.
 * 2. To set the variables diag_line<i> and offd_line<i> you can use the old
 * UG shell command 'set':
 * i = 4; j = 8;
 * set diag_line1 1 @i - @j;
 * sets diag_line1 to the string "1 4 - 8". The construction
 * diag_line1 = "1 @i - @j";
 * would set it to "1 @i - @j". The command 'set' does not fulfil any
 * arithmetics, so '-' in the description is not treated as an arithmetical
 * sign!
 *     The command 'set' does not accept arrays. for the automatic
 * pattern creation you can use something like
 * i = 2;
 * set var_name diag_line@i;
 * set @var_name @i;
 * This sets 'diag_line2' to "2".
 */
#include "config.h"
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
#include "cmdline.h"   /* for shell commands               */

/* Names of the variables: */
static char size_var [] = ":sparsity_pattern:size";
static char diag_line [] = ":sparsity_pattern:diag_line%d";
static char offd_line [] = ":sparsity_pattern:offd_line%d";

/* parse_line_desc - parses a description of a block line sparse
 * pattern. The function fills the entries of a given array with '*'
 * for potentially nonzero entries and with '0' for zero entries.
 * 'size' entries are initialized. The function returns 0 if OK,
 * nonzero on an error:
 */
static INT parse_line_desc
(
  INT size, /* the size of the block */
  char * desc, /* the description */
  char * result /* to fill with '*' and '0' */
)
{
  unsigned n, j, j1, j2;
  char * str;

  /* Fill the result with '0': */
  for (j = 0; j < size; j++)
    result [j] = '0';

  /* Skip the initial white space: */
  sscanf (desc, " %n", &n);
  str = desc + n;

  /* Interpret the description elements: */
  while (*str != '\0')
  {
    if (sscanf (str, "%u - %u %n", &j1, &j2, &n) == 2)
    {
      if (j1 < 1 || j2 > size || j1 > j2)
        return __LINE__;
      for (j = j1 - 1; j < j2; j++)
        result [j] = '*';
    }
    else if (sscanf (str, "%u %n", &j, &n) == 1)
    {
      if (j < 1 || j > size)
        return __LINE__;
      result [j - 1] = '*';
    };
    str += n;
  };

  return 0;
}

/* create_sparsity_pattern_command - implements a command that creates
 * the UG sparsity pattern structures from the sparsity pattern descriptions:
 */
static INT create_sparsity_pattern_command (INT argc, char * * argv)
{
  int size, i, j;
  char var_name [64];
  char * var_value;
  char * Tnn, * Dnn;

  if (argc != 1)
    PrintErrorMessage ('W', "create_sparsity_pattern",
                       "This command needs no arguments!\n");

  /* Get and compute the size, allocate the memory: */
  if (GetStringValueInt (size_var, &size))
  {
    PrintErrorMessageF ('E', "create_sparsity_pattern",
                        "Cannot get the value of the variables '%s'", size_var);
    return CMDERRORCODE;
  };

  i = (size + 1) * size;
  if ((Tnn = malloc (i)) == NULL)
  {
    PrintErrorMessageF ('E', "create_sparsity_pattern",
                        "Cannot allocate %d bites for the offdiagonal pattern", i);
    return CMDERRORCODE;
  };
  if ((Dnn = malloc (i)) == NULL)
  {
    PrintErrorMessageF ('E', "create_sparsity_pattern",
                        "Cannot allocate %d bites for the diagonal pattern", i);
    free (Tnn);
    return CMDERRORCODE;
  };

  /* Prepare the patterns: */
  i--;
  Tnn [i] = Dnn [i] = '\0';
  for (i = 1; i < size; i++)
  {
    j = i * (size + 1) - 1;
    Tnn [j] = Dnn [j] = ' ';
  };

  /* Parse the descriptions: */
  for (i = 1; i <= size; i++)
  {
    j = (i - 1) * (size + 1);
    /* The off-diagonal pattern: */
    sprintf (var_name, offd_line, i);
    if ((var_value = GetStringVar (var_name)) == NULL)
    {
      PrintErrorMessageF ('E', "create_sparsity_pattern",
                          "Cannot get variable '%s'", var_name);
      free (Dnn); free (Tnn);
      return CMDERRORCODE;
    };
    if (parse_line_desc (size, var_value, Tnn + j))
    {
      PrintErrorMessageF ('E', "create_sparsity_pattern",
                          "Cannot parse variable '%s': %s", var_name, var_value);
      free (Dnn); free (Tnn);
      return CMDERRORCODE;
    };
    /* The diagonal pattern: */
    sprintf (var_name, diag_line, i);
    if ((var_value = GetStringVar (var_name)) == NULL)
    {
      PrintErrorMessageF ('E', "create_sparsity_pattern",
                          "Cannot get variable '%s'", var_name);
      free (Dnn); free (Tnn);
      return CMDERRORCODE;
    };
    if (parse_line_desc (size, var_value, Dnn + j))
    {
      PrintErrorMessageF ('E', "create_sparsity_pattern",
                          "Cannot parse variable '%s': %s", var_name, var_value);
      free (Dnn); free (Tnn);
      return CMDERRORCODE;
    };
  };

  /* Set the UG structure variables: */
  if (SetStringVar ("Tnn", Tnn) || SetStringVar ("Dnn", Dnn))
  {
    PrintErrorMessage ('E', "create_sparsity_pattern",
                       "Cannot create variable 'Tnn' or 'Dnn'");
    free (Dnn); free (Tnn);
    return CMDERRORCODE;
  };

  /* Done: */
  free (Dnn); free (Tnn);
  return OKCODE;
}

/* The installer: */
INT NS_DIM_PREFIX Install_Sparsity_Descriptions ()
{
  return CreateCommand ("create_sparsity_pattern",
                        create_sparsity_pattern_command) == NULL;
}

/* End of File */
