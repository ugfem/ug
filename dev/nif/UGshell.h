// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#import <appkit/appkit.h>

@interface UGshell : Object
{
  id shellWindow;
  id shellScrollView;
  id shellText;
  id infoPanel;
  id theMenu;
  id shellFont;
}

- setUp;
- shellWindow;
- shellText;
- appendToText:(const char *)val;
@end
