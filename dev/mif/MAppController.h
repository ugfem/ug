// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  MAppController.h												*/
/*																			*/
/* Purpose:   Controller class for UG applications in OPENSTEP.				*/
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  Interdisziplin"ares Zentrum f"ur Wissenschaftliches Rechnen	*/
/*			  Universit"at Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69210 Heidelberg												*/
/*			  email: Volker.Reichenberger@IWR.Uni-Heidelberg.DE		        */
/*																			*/
/*	History:  June 4, 1999 begin (based on OPENSTEP code)					*/
/*																			*/
/****************************************************************************/

#import <AppKit/AppKit.h>

@interface MAppController : NSObject
{
  id inspector;
  id shell;

  NSString *_openPanelPath;
}
- (void)appendText:(id)sender;
- (void)closeGraphic:(id)sender;
- (void)executeCommand:(id)sender;
- (void)runScript:(id)sender;
- (void)newGraphic:(id)sender;
- (void)saveDocument:(id)sender;
- (void)saveDocumentAs:(id)sender;
- (void)showInfo:(id)sender;
- (void)showInspector:(id)sender;
- (void)showPreferences:(id)sender;
- (void)dealloc;
@end
