// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/** \file
    \brief Create a link time marker for the enable-dune build option.
    \author Oliver Sander

    If you want to use UG as part of DUNE, you have to build it with the
    build-option --enable-dune.  However, this option affects header
    files exclusively.  It is therefore impossible for the DUNE build
    system to check whether UG has been built with this opion.  People
    tend to forget it and then encounter strange errors.  To avoid this
    we add a data field here which gets compiled into the UG library
    if --enable-dune is set and can therefore be checked by the
    DUNE build system.
 */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

#include "namespace.h"
#include "dunemarker.h"

#ifdef FOR_DUNE

#ifndef __cplusplus
#error You must force the build system to compile as C++ when setting --enable-dune!
#endif

/** \brief Data field which is only there for the DUNE build system to check
    whether UG has been compiled with FOR_DUNE.
 */
int NS_PREFIX duneMarker = 0;
#endif
