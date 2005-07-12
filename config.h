// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* README:

   This header is a dummy header so that both the classical UG build
   system and a automake/autoconf-approach can be used. arch/compiler.h
   switches it's functionality depending on symbols in config.h. As the
   compiler complains if no config.h is found we need this dummy.

   $Id$
 */

/* Define this if you want to use the full refinement rule set for tetrahedra */
#ifdef ModelP
#define TET_RULESET
#endif
