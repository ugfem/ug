// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef UG_PIXEL_H
#define UG_PIXEL_H

/*
   this is the common datastructure used by the bullet-plotter and the
   output devices

 */

typedef struct {
  unsigned char cindex;
  unsigned char intensity;
} PIXEL;

#endif
