// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/*

   defines macros to put symbols into namespaces if C++ compiler is
   used. Everything is put into the namespace UG, dimension-dependent
   functions into UG::D2 or UG::D3

   in a header-files use:

   #include "namespace.h"

    START_NAMESPACE

    ...

    END_NAMESPACE

   and in the implementation:

   #include "namespace.h"

    USING_UG_NAMESPACES

    int NS_PREFIX function(...) { ... };

 */

//#ifndef UG_NAMESPACE_H
//#define UG_NAMESPACE_H

#ifdef __cplusplus

// # if defined(__TWODIM__) || defined(__THREEDIM__)
//   /* namespace for dimension */
// #  ifdef __TWODIM__
// #   define NAMESPACE D2
// #   define NS_PREFIX D2::
// #  define USING_UG_NAMESPACES using namespace UG; using namespace UG::D2;
// #  else
// #   define NAMESPACE D3
// #   define NS_PREFIX D3::
// #  define USING_UG_NAMESPACES using namespace UG; using namespace UG::D3;
// #  endif
// #  define START_NAMESPACE namespace UG { namespace NAMESPACE {
// #  define END_NAMESPACE } }
// # else
//   /* no dimension set */
// #  define NS_PREFIX UG::
// #  define START_NAMESPACE namespace UG {
// #  define END_NAMESPACE }
// #  define USING_UG_NAMESPACES using namespace UG;
// # endif

#ifdef _3
#define START_NAMESPACE namespace UG3d {
#define END_NAMESPACE }
#define USING_UG_NAMESPACES using namespace UG3d;
#define NS_PREFIX UG3d::
#else
#define START_NAMESPACE namespace UG2d {
#define END_NAMESPACE }
#define USING_UG_NAMESPACES using namespace UG2d;
#define NS_PREFIX UG2d::
#endif

#else
/* normal C-compiler, no namespace-stuff */
# define START_NAMESPACE
# define END_NAMESPACE
# define NS_PREFIX
# define USING_UG_NAMESPACES
#endif

/* check if the required symbols exist */
#if !defined(NS_PREFIX) || !defined(START_NAMESPACE) || !defined(END_NAMESPACE)
# error missing symbol!
#endif

//#endif
