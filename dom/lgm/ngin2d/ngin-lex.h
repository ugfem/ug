// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UG_NGIN_LEX_H
#define UG_NGIN_LEX_H

/* C++-compilers don't like this defined implicitly */
int ng2dlex();

/* error function used in YACC-parser as well */
int NP2d_Error (int *line, char *text);

#endif
