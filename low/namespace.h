// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/** \file

   Defines macros to put symbols into namespaces if C++ compiler is
   used. Everything is put into the namespace UG, dimension-dependent
   functions into UG::D2 or UG::D3

   in a header-files use:

   #include "namespace.h"

    START_UG_NAMESPACE

    ...

    END_UG_NAMESPACE

   for stuff that is independent of the space dimension or

   #include "namespace.h"

    START_UGDIM_NAMESPACE

    ...

    END_UGDIM_NAMESPACE

   else.  In the implementation:

   #include "namespace.h"   // of course

    USING_UG_NAMESPACE       // for stuff from the namespace UG

   or
    USING_UGDIM_NAMESPACE    // for stuff from UG::D3 resp. UG::D2

   Write

    int NS_PREFIX function(...) { ... };

   if function is declared in UG and

    int NS_DIM_PREFIX function(...) { ... };

   if it is declared in a namespace with dimension.

 */

#ifndef UG_NAMESPACE_H
#define UG_NAMESPACE_H

#ifdef __cplusplus

#define START_UG_NAMESPACE namespace UG {
#define END_UG_NAMESPACE }
#define END_UGDIM_NAMESPACE }}
#define NS_PREFIX UG::
#define USING_UG_NAMESPACE using namespace UG;

#ifdef _3
#define START_UGDIM_NAMESPACE namespace UG { namespace D3 {
#define USING_UGDIM_NAMESPACE using namespace UG::D3;
#define USING_UG_NAMESPACES namespace UG {namespace D3 {} }; using namespace UG; using namespace UG::D3;
#define NS_DIM_PREFIX UG::D3::
#else
#define START_UGDIM_NAMESPACE namespace UG { namespace D2 {
#define USING_UGDIM_NAMESPACE using namespace UG::D2;
#define USING_UG_NAMESPACES namespace UG {namespace D2 {} }; using namespace UG; using namespace UG::D2;
#define NS_DIM_PREFIX UG::D2::
#endif

#else
/* normal C-compiler, no namespace-stuff */
# define START_UG_NAMESPACE
# define START_UGDIM_NAMESPACE
# define END_UG_NAMESPACE
# define END_UGDIM_NAMESPACE
# define NS_PREFIX
#define NS_DIM_PREFIX
# define USING_UG_NAMESPACE
# define USING_UGDIM_NAMESPACE
# define USING_UG_NAMESPACES
#endif

/* check if the required symbols exist */
#if !defined(NS_PREFIX) || !defined(START_UG_NAMESPACE) || !defined(END_UG_NAMESPACE)
# error missing symbol!
#endif

#endif
