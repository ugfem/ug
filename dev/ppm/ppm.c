// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ppm.c                                                                                                         */
/*																			*/
/* Purpose:   write ppm files for bullet plotter                                                        */
/*																			*/
/* Author:	  Michael Lampe                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   19.08.98 begin                                                                            */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "defaults.h"
#include "fileopen.h"
#include "ugdevices.h"
#include "initdev.h"
#include "general.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define COLORS 256

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

typedef struct {
  FILE  *file;
  INT header_length;
  INT width;
  INT height;
} PPM_WINDOW;

static OUTPUTDEVICE *ppm_OutputDevice;
static PPM_WINDOW   *ppm_Window;

/* color tables */
static short RedTab[COLORS];
static short GreenTab[COLORS];
static short BlueTab[COLORS];

/****************************************************************************/
/*																			*/
/*                        basic drawing functions                               */
/*																			*/
/****************************************************************************/

static void ppm_Move(SHORT_POINT point) {
  return;
}
static void ppm_Draw(SHORT_POINT point) {
  return;
}
static void ppm_Polyline(SHORT_POINT *points, INT n) {
  return;
}
static void ppm_InversePolyline(SHORT_POINT *points, INT n) {
  return;
}
static void ppm_Polygon(SHORT_POINT *points, INT n) {
  return;
}
static void ppm_ShadedPolygon(SHORT_POINT *points, INT n, DOUBLE intensity) {
  return;
}
static void ppm_InversePolygon(SHORT_POINT *points, INT n) {
  return;
}
static void ppm_ErasePolygon(SHORT_POINT *points, INT n) {
  return;
}
static void ppm_Marker(short n, short s, SHORT_POINT point) {
  return;
}
static void ppm_Polymark(short n, SHORT_POINT *points) {
  return;
}
static void ppm_InvMarker(short n, short s, SHORT_POINT point) {
  return;
}
static void ppm_InvPolymark(short n, SHORT_POINT *points) {
  return;
}
static void ppm_DrawText(const char *s, INT mode) {
  return;
}
static void ppm_CenteredText(SHORT_POINT point, const char *s, INT mode) {
  return;
}
static void ppm_ClearViewPort(void) {
  return;
}

/****************************************************************************/
/*																			*/
/*                             set functions                                    */
/*																			*/
/****************************************************************************/

static void ppm_SetLineWidth(short w) {
  return;
}
static void ppm_SetTextSize(short s) {
  return;
}
static void ppm_SetMarker(short s) {
  return;
}
static void ppm_SetMarkerSize(short s) {
  return;
}
static void ppm_SetColor(long index) {
  return;
}
static void ppm_SetPaletteEntry(long index, short r, short g, short b) {
  return;
}
static void ppm_GetPaletteEntry(long index, short *r, short *g, short *b) {
  return;
}
static void ppm_Flush(void) {
  return;
}

/* setup color map */
static void ppm_SetNewPalette(long x, long count, short *r, short *g, short *b)
{
  INT i;

  for (i = 0; i < count; i++) {
    RedTab[i]   = r[i];
    GreenTab[i] = g[i];
    BlueTab[i]  = b[i];
  }
  RedTab[1] = GreenTab[1] = BlueTab[1] = 0xD0;         /* set gray */
}

/****************************************************************************/
/*																			*/
/*                 copy the pixel buffer to the ppm file                    */
/*																			*/
/****************************************************************************/

static void ppm_PlotPixelBuffer(void *buffer, void *data, INT len,
                                int x, int y, int w, int h)
{
  unsigned char *p;
  int i, j;
  long offset;

  p = buffer;
  offset = ppm_Window->header_length + 3 * (y*ppm_Window->width + x);
  for (j = 0; j < h; j++) {
    fseek(ppm_Window->file, offset, SEEK_SET);
    for (i = 0; i < w; i++) {
      fputc(RedTab  [*p],   ppm_Window->file);
      fputc(GreenTab[*p],   ppm_Window->file);
      fputc(BlueTab [*p++], ppm_Window->file);
    }
    offset += 3 * ppm_Window->width;
  }
}

/****************************************************************************/
/*																			*/
/*                  operations for managing windows                             */
/*																			*/
/****************************************************************************/

static WINDOWID ppm_OpenWindow(const char *title, INT rename,
                               INT x, INT y, INT width, INT height,
                               INT *Global_LL, INT *Global_UR,
                               INT *Local_LL, INT *Local_UR, INT *error)
{
  PPM_WINDOW *ppm_window;
  FILE *file;
  char ppm_path[160], header[32];
  char white[3] = {255, 255, 255};
  INT i, n;

  /* allocate ppm window structure */
  *error = 0;
  ppm_window = (PPM_WINDOW *)malloc(sizeof(PPM_WINDOW));
  if (ppm_window == NULL) {
    *error = 1;
    return 0;
  }

  /* open ppm window */
  if (GetDefaultValue(DEFAULTSFILENAME, "ppmfilesdir", ppm_path) == 0)
    file = FileOpenUsingSearchPath_r(title, "w", ppm_path, rename);
  else
    file = fileopen(title, "w");
  if (file == NULL) {
    *error = 1;
    return 0;
  }
  ppm_window->file = file;

  /* fill in devices coordinate system */
  Global_LL[0] = 0; Global_LL[1] = height;
  Global_UR[0] = width; Global_UR[1] = height;
  Local_LL[0] = 0; Local_LL[1] = height;
  Local_UR[0] = width; Local_UR[1] = 0;

  ppm_window->width  = ++width;
  ppm_window->height = ++height;

  /* write the ppm file header */
  sprintf(header, "P6\n%d %d\n255\n", width, height);
  n = ppm_window->header_length = strlen(header);
  fwrite(header, 1, n, file);

  /* make a plain white canvas */
  for (i = 0; i < width*height; i++)
    fwrite(white, 3, 1, file);

  ppm_Window = ppm_window;

  return ((WINDOWID)ppm_window);
}

static INT ppm_CloseWindow(WINDOWID win)
{
  PPM_WINDOW *ppm_window;

  ppm_window = (PPM_WINDOW *)win;
  fclose(ppm_window->file);
  free(ppm_window);

  return 0;
}

static INT ppm_SetOutput(WINDOWID win)
{
  ppm_Window = (PPM_WINDOW *)win;
  return 0;
}

static INT ppm_UpdateOutput(WINDOWID win, INT tool)
{
  return 0;
}

/****************************************************************************/
/*																			*/
/*                            init ppm device                                   */
/*																			*/
/****************************************************************************/

static void ppm_InitInterface (OUTPUTDEVICE *thePort)
{
  /* init pointers to basic drawing functions */
  thePort->Move            = ppm_Move;
  thePort->Draw            = ppm_Draw;
  thePort->Polyline        = ppm_Polyline;
  thePort->Polygon         = ppm_Polygon;
  thePort->ShadedPolygon   = ppm_ShadedPolygon;
  thePort->InversePolygon  = ppm_InversePolygon;
  thePort->ErasePolygon    = ppm_ErasePolygon;
  thePort->Polymark        = ppm_Polymark;
  thePort->InvPolymark     = ppm_InvPolymark;
  thePort->DrawText        = ppm_DrawText;
  thePort->CenteredText    = ppm_CenteredText;
  thePort->ClearViewPort   = ppm_ClearViewPort;

  /* init pointers to set functions */
  thePort->SetLineWidth    = ppm_SetLineWidth;
  thePort->SetTextSize     = ppm_SetTextSize;
  thePort->SetMarker       = ppm_SetMarker;
  thePort->SetMarkerSize   = ppm_SetMarkerSize;
  thePort->SetColor        = ppm_SetColor;
  thePort->SetPaletteEntry = ppm_SetPaletteEntry;
  thePort->SetNewPalette   = ppm_SetNewPalette;

  /* init pointers to miscellaneous functions */
  thePort->GetPaletteEntry = ppm_GetPaletteEntry;
  thePort->Flush           = ppm_Flush;
  thePort->PlotPixelBuffer = ppm_PlotPixelBuffer;

  /* set palette */
  UgSetPalette(ppm_OutputDevice, COLOR_PALETTE);

  thePort->black = 255;
  thePort->gray = 1;
  thePort->white = 0;
  thePort->red = 254;
  thePort->green = 128;
  thePort->blue = 2;
  thePort->cyan = 65;
  thePort->orange = 220;
  thePort->yellow = 191;
  thePort->darkyellow = 205;
  thePort->magenta = 1;
  thePort->hasPalette = 1;
  thePort->range = 256;
  thePort->spectrumStart = 2;
  thePort->spectrumEnd = 254;
  thePort->signx = 1;
  thePort->signy = -1;
}

static OUTPUTDEVICE *ppm_InitDevice(void)
{
  /* create output device */
  if ((ppm_OutputDevice=CreateOutputDevice("ppm")) == NULL ) return(NULL);

  /* init output device   */
  ppm_OutputDevice->OpenOutput     = ppm_OpenWindow;
  ppm_OutputDevice->CloseOutput    = ppm_CloseWindow;
  ppm_OutputDevice->ActivateOutput = ppm_SetOutput;
  ppm_OutputDevice->UpdateOutput   = ppm_UpdateOutput;
  ppm_OutputDevice->v.locked       = 1;
  ppm_OutputDevice->PixelRatio     = 1.0;

  ppm_InitInterface(ppm_OutputDevice);

  UserWrite("output device 'ppm' created\n");

  return (ppm_OutputDevice);
}

INT InitPPMDevice(void)
{
  if ((ppm_InitDevice()) == NULL) return 1;

  return 0;
}
