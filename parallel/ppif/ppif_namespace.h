// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef PPIF_NAMESPACE_H
#define PPIF_NAMESPACE_H

/****************************************************************************/
/*                                                                          */
/* Put ppif into its own namespace if it is compiled as c++                 */
/*                                                                          */
/****************************************************************************/

#ifdef __cplusplus
#define PPIF_NS_PREFIX PPIF::
#define USING_PPIF_NAMESPACE using namespace PPIF;
#else
#define PPIF_NS_PREFIX
#define USING_PPIF_NAMESPACE
#endif

#endif
