// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* GRAPE  5.1  header for sgi created */
#ifndef G_VERSION_ID
#define G_VERSION_ID "5.1"
#endif
#define G_HEADER_TIMESTAMP 823186875
#define G_HEADER_CODE 64291806
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#ifdef mem_alloc
#undef mem_alloc
#undef mem_free
#endif
#ifdef ardent
#include <sys/types.h>
#endif
#ifdef  __STDC__
#define ANSI_PARM(X) X
#else
#define ANSI_PARM(X) ()
#endif
#ifndef NIL
#define NIL 0L
#else
#undef NIL
#define NIL 0L
#endif
#ifndef ERR
#define ERR -1L
#else
#undef ERR
#define ERR -1L
#endif
typedef struct frame {
  struct frame *link;
  struct instance *object;
  struct method *method;
} FRAME;
typedef enum {
  msStatic   = 1,
  msOutdate  = 2,
  msUpdate   = 4
} METHSTAT;
typedef struct method {
  struct method *left;
  struct method *right;
  short balance;
  char *name;
  struct class *class;
  char *(*code)();
  int namelen;
  METHSTAT status;
  char *closure ;
} METHOD;
typedef struct class {
  struct class *itself;
  char *name;
  struct class *superclass;
  struct class *subclasses;
  struct class *next;
  struct instance *instances;
  struct method *methods;
  struct method *all_methods;
  struct method *last_method;
  char *last_called_name;
  int is_static_var;
  int isize;
} CLASS;
#ifdef INODE_TEST
#define INODE_ENTRY struct inode *inode;
#else
#define INODE_ENTRY
#endif
#ifdef G_DOUBLY_LINKED_INST
#define G_INST_NEXT next_inst_of_class, *prev_inst_of_class
#else
#define G_INST_NEXT next
#endif
#define INSTANCE_STRUCT \
  struct class *class; \
  char *name; \
  struct instance *G_INST_NEXT; \
  INODE_ENTRY \
  unsigned int refcount
typedef struct instance {
  INSTANCE_STRUCT;
} INSTANCE;
#ifdef __i386__
#define G_SORT_METHODS_BY_LEN
#endif
#ifdef G_SORT_METHODS_BY_LEN
#ifdef __i386__
#define G_METHCMP(aname, alen, bname, blen)\
  (((alen)-(blen)) ? ((alen)-(blen)) : memcmp((aname),(bname),(alen)))
#define G_METHCMPTMP(aname, alen, bname, blen, tmp)\
  ((tmp = ((alen) - (blen))) ? tmp : memcmp((aname), (bname),(alen)))
#else
#define G_METHCMP(aname, alen, bname, blen)\
  (((alen) - (blen)) ? ((alen) - (blen)) : strcmp ((aname), (bname)))
#define G_METHCMPTMP(aname, alen, bname, blen, tmp)\
  ((tmp = ((alen) - (blen))) ? tmp : strcmp ((aname), (bname)))
#endif
#else
#define G_METHCMP(aname, alen, bname, blen) (strcmp ((aname), (bname)))
#define G_METHCMPTMP(aname, alen, bname, blen, tmp)\
  (strcmp ((aname), (bname)))
#endif
#define G_METHSTRCMP(aname, bname)\
  G_METHCMP ((aname), strlen (aname), (bname), strlen (bname))
typedef double VEC2[2];
typedef double VEC3[3];
typedef double VEC4[4];
typedef float FVEC2[2];
typedef float FVEC3[3];
typedef float FVEC4[4];
typedef int INT2[2];
typedef int INT3[3];
typedef int INT4[4];
typedef double MATRIX33[3][3];
typedef double MATRIX44[4][4];
#define ASSIGN(p,q) ((p) = (q),(p)->refcount++,(p))
#define G_EMPTY_ARG
#ifndef __FILE__
#  define __FILE__ NULL
#endif
#ifndef __LINE__
#  define __LINE__ 0
#endif
#define ASSURE(condition,message,error_exit) \
  if ( ! (condition))  { \
    grape_error ( message, __FILE__, __LINE__, 0 ); \
    error_exit; \
  }
#ifdef DEBUG
#define MESSAGE(S) printf S
#else
#define MESSAGE(S)
#endif
extern void *mem_alloc ANSI_PARM((size_t));
extern void mem_free ANSI_PARM((void *, size_t));
extern void *mem_realloc ANSI_PARM((void *, size_t, size_t));
#define G_BLOCKSIZE (~((size_t)17))
extern int grape_error ANSI_PARM((char *, char *, int, int));
extern int _trace;
extern INSTANCE *start_method ANSI_PARM((int));
extern void end_method ANSI_PARM((void *));
extern char *(*grape())();
#define G_CLASS         1
#define G_INSTANCE      2
#define START_METHOD(mode) start_method((mode))
#define END_METHOD(val) { end_method((val)); return ((val)); }
#define GRAPE(obj,meth) (*grape((obj),(meth)))
typedef struct {void *class, *object;} G_CLASS_CAST_DESCRIPTOR;
extern G_CLASS_CAST_DESCRIPTOR g_class_cast_descriptor;
#define G_CAST_CLASS(obj,newclass)\
  (g_class_cast_descriptor.object = (obj),\
   g_class_cast_descriptor.class = (newclass),\
   &g_class_cast_descriptor)
typedef struct G_Node {
  struct G_Node *next;
  struct G_Node *prev;
  void *obj;
} G_NODE;
#define G_LIST_STRUCT \
  INSTANCE_STRUCT; \
  char *type; \
  G_NODE *head; \
  G_NODE *tail; \
  G_NODE *iter
typedef struct g_list {
  G_LIST_STRUCT;
} G_LIST;
int g_list_add_head(G_LIST *, void *);
int g_list_add_tail(G_LIST *, void *);
int g_list_add_pre(G_LIST *, void *);
int g_list_add_post(G_LIST *, void *);
int g_list_insert(G_LIST *, void *, void *);
void *g_list_rem_head(G_LIST *);
void *g_list_rem_tail(G_LIST *);
void *g_list_rem_curr(G_LIST *);
void *g_list_remove(G_LIST *, void *);
void *g_list_first(G_LIST *);
void *g_list_last(G_LIST *);
void *g_list_next(G_LIST *);
void *g_list_previous(G_LIST *);
void *g_list_current(G_LIST *);
int g_list_set_current(G_LIST *, void *);
#ifdef INODE_TEST
typedef enum {
  isUptodate = 1,
  isUpdating = 2
} INODESTAT;
typedef struct inode {
  INSTANCE *obj;
  G_LIST *ingr;
  G_LIST *dep;
  G_LIST *update;
  INODESTAT status;
} INODE;
typedef enum {
  psUptodate = 1,
  psUpdating = 2
} IPORTSTAT;
typedef struct iport {
  INODE *node;
  struct iport *conn;
  char *id;
  CLASS *type;
  IPORTSTAT status;
} IPORT;
#endif
typedef struct lightsource_dev {
  int type;
  int on_off;
  double color[3];
  double position[3];
  double direction[3];
  double open_angle;
} LIGHTSOURCE_DEV;
typedef struct clipplane_dev {
  int on_off;
  double normal[3];
  double distance;
} CLIPPLANE_DEV;
typedef struct suprop_dev {
  double emission[3];
  double ambient[3];
  double diffuse[3];
  double specular[3];
  double specular_exp;
  double transparency;
} SUPROP_DEV;
typedef struct {
  long x, y;
} RESOLUTION_TYPE;
typedef double RGB_TYPE[3];
typedef double NORMAL_TYPE[3];
typedef struct {
  long n;
  float parray[10][3];
  unsigned short iarray[10];
} GOURAUD_PATCH_TYPE;
typedef struct {
  int on_off;
  double zmin,zmax;
  double dark;
} DEPTHCUE_TYPE;
#define DEVICE_STRUCT \
  INSTANCE_STRUCT; \
  unsigned RW_flags
typedef struct device {
  DEVICE_STRUCT;
} DEVICE;
#define GRAPHICDEVICE_STRUCT \
  DEVICE_STRUCT; \
  int *info; \
  void (*update)(void); \
  void (*clear)(void); \
  void (*message)(long, char *); \
  int (*attribute)(int, int, void *); \
  void (*transform)(int, int, MATRIX44); \
  void (*move)(VEC3); \
  void (*draw)(VEC3); \
  void (*begin_patch)(void); \
  void (*patch_vertex)(VEC3); \
  void (*patch_vertex_index)(int); \
  void (*patch_normal)(VEC3); \
  void (*patch_color)(VEC3); \
  void (*end_patch)(void); \
  void (*move_direct)(VEC3); \
  void (*draw_direct)(VEC3); \
  void (*begin_patch_direct)(void); \
  void (*patch_vertex_direct)(VEC3); \
  void (*patch_normal_direct)(VEC3); \
  void (*patch_color_direct)(VEC3); \
  void (*end_patch_direct)(void); \
  void (*begin)(unsigned int); \
  void (*vertex)(VEC3); \
  void (*normal)(VEC3); \
  void (*color)(VEC3); \
  void (*color_alpha)(VEC3); \
  void (*end)(void); \
  void (*texture_init)(unsigned char *,unsigned int, int, int); \
  void (*texture_coord)(double, double); \
  int (*generate_lightsource)(int); \
  void (*lightsource)(int, int, LIGHTSOURCE_DEV *); \
  void (*delete_lightsource)(int); \
  int (*generate_clipplane)(void); \
  void (*clipplane)(int, int, CLIPPLANE_DEV *); \
  void (*delete_clipplane)(int); \
  void (*text)(VEC3, char *); \
  void (*set_full_screen)(int); \
  void (*video_record_function)(void); \
  int grid_patch; \
  int interactive; \
  int global_grid_patch
typedef struct graphicdevice {
  GRAPHICDEVICE_STRUCT;
} GRAPHICDEVICE;
#define G_GRID  0x000003b8
#define G_PATCH 0x000003d8
#define G_READONLY      1
#define G_WRITEONLY     2
#define G_READWRITE (G_READONLY | G_WRITEONLY)
#define G_MODE_QUERY    0
#define G_MODE_GET      1
#define G_MODE_SET      2
#define G_MODE_PRE      5
#define G_MODE_POST     6
#define G_MODE_SWAP     7
#define G_MATRIX_NONE    0
#define G_MATRIX_PROJECT 1
#define G_MATRIX_VIEW    2
#define G_MATRIX_MODEL   3
#define G_MATRIX_PVM     4
#define G_MATRIX_IVM     5
#define G_LIGHT_FREE    0
#define G_LIGHT_AMBIENT 1
#define G_LIGHT_DIRECT  2
#define G_LIGHT_POINT   3
#define G_LIGHT_SPOT    4
#define G_RASTER_ABILITY         1
#define G_COLOR_ABILITY          2
#define G_ZBUFF_ABILITY          3
#define G_RESOLUTION             4
#define G_ASPECTRATIO            5
#define G_BACKGROUND_COLOR       6
#define G_LINE_COLOR             7
#define G_LINE_TYPE              8
#define G_ANTIALIAS              9
#define G_DEPTHCUE              10
#define G_PATCH_SUPROP          13
#define G_TEXT_COLOR            14
#define G_TEXT_TYPE             15
#define G_DOUBLEBUFFER          16
#define G_LIGHT_MODEL           17
#define G_LIGHT_NUMBER          18
#define G_HAS_DOUBLEBUFFER      19
#define G_PATCH_DRAW_EDGE       20
#define G_ORIGIN                21
#define G_TEXTURE               22
#define G_POINTS                0x0000
#define G_LINES                 0x0001
#define G_LINE_LOOP             0x0002
#define G_LINE_STRIP            0x0003
#define G_TRIANGLES             0x0004
#define G_TRIANGLE_STRIP        0x0005
#define G_TRIANGLE_FAN          0x0006
#define G_QUADS                 0x0007
#define G_QUAD_STRIP            0x0008
#define G_POLYGON               0x0009
#define G_COLOR_INDEX           0x1900
#define G_RED                   0x1903
#define G_GREEN                 0x1904
#define G_BLUE                  0x1905
#define G_ALPHA                 0x1906
#define G_RGB                   0x1907
#define G_RGBA                  0x1908
#define G_LUMINANCE             0x1909
#define G_LUMINANCE_ALPHA       0x190A
void graphicdevice_define_stddev(GRAPHICDEVICE *, void (* )());
#define G_DEPTHSORT_DRAW_CUTTING_LINES 1
#define G_DEPTHSORT_CHECK_CUT 2
#define G_DEPTHSORT_RGB_PER_VERTEX 4
#define G_DEPTHSORT_OMIT_OBSCURED 8
#define G_DEPTHSORT_RAW_CLIPPING 16
#define G_DEPTHSORT_NORMAL_ACTION 0
#define G_DEPTHSORT_SHADED_ZSORT  1
#define G_DEPTHSORT_ZSORT         2
typedef enum stateflags {
  sfShiftKey      =  0x00000001,
  sfLockKey       =  0x00000002,
  sfControlKey    =  0x00000004,
  sfMod1Key       =  0x00000008,
  sfMod2Key       =  0x00000010,
  sfMod3Key       =  0x00000020,
  sfMod4Key       =  0x00000040,
  sfMod5Key       =  0x00000080,
  sfLeftMouse     =  0x00000100,
  sfMiddleMouse   =  0x00000200,
  sfRightMouse    =  0X00000400,
  sfOpt1Mouse     =  0x00000800,
  sfOpt2Mouse     =  0x00001000,
  sfSomeMouse     =  0x00001F00
} STATEFLAGS;
typedef enum controlkeys {
  ckReturn,
  ckEscape,
  ckBackspace,
  ckTab,
  ckDelete,
  ckHome,
  ckLeftArrow,
  ckUpArrow,
  ckRightArrow,
  ckDownArrow,
  ckPageUp,
  ckPageDown,
  ckEnd,
  ckBegin
} CONTROLKEYS;
typedef enum evtype {
  evNothing       = 0x00000000,
  evMousePressed  = 0x00000001,
  evMouseReleased = 0x00000002,
  evMouseHold     = 0x00000004,
  evMouseDragged  = 0x00000008,
  evMouse         = 0x0000000F,
  evNormalKey     = 0x00000010,
  evCtrlKey       = 0x00000020,
  evKey           = 0x000000F0,
  evRedraw        = 0x00000100,
  evResize        = 0x00000200,
  evMessage       = 0x0000FF00
} EVTYPE;
typedef enum whotype {
  evForCtl        = 0x00000001,
  evForGrph       = 0x00000002
} WHOTYPE;
typedef struct event {
  EVTYPE what;
  unsigned int which;
  unsigned int state;
  double x, y;
  double deltax, deltay;
  WHOTYPE who;
  struct layer *layer;
} EVENT;
typedef struct {
  double x1, y1, x2, y2;
} GRECT;
typedef struct grect_list {
  GRECT *rect;
  struct grect_list *next;
  int already_used;
} GRECT_LIST;
typedef enum g_clip_type {
  g_ctNone        = 0x00000000,
  g_ctEnableRects = 0x00000001,
  g_ctDisableRects
} G_CLIP_TYPE;
typedef struct g_clip_descriptor {
  G_CLIP_TYPE which;
  int max_no_of_rects;
} G_CLIP_DESCRIPTOR;
#define CONTROLDEVICE_STRUCT \
  DEVICE_STRUCT; \
  int *info; \
  void (*update)(void); \
  int (*attribute)(int, int, int *); \
  void (*clear)(void); \
  void (*set_color)(int); \
  void (*set_rgb_color)(VEC3); \
  void (*set_text_color)(int); \
  void (*get_color)(int *); \
  void (*change_color)(int, VEC3); \
  void (*moveto)(double, double); \
  void (*drawto)(double, double); \
  void (*relative_draw)(double, double); \
  void (*draw_string)(double, double, char *); \
  void (*get_mouse)(double *, double *); \
  void (*fix_mouse)(double, double, double, double, double, double); \
  void (*free_mouse)(void); \
  void (*set_mouse)(double, double); \
  void (*rectangle)(double, double, double, double); \
  void (*circle)(double, double, double); \
  void (*reinit)(void); \
  void (*get_control_size)(double *, double *); \
  void (*wait)(double); \
  void (*clear_event_queue)(void); \
  void (*get_spaceball_data)(double *, double *, double *, double *, \
                             double *, double *, double *, int *); \
  int (*handle_events)(EVENT *, int); \
  int (*menu_button)(void); \
  int (*read_string)(double, double, char *, double); \
  int (*middle_button)(void); \
  int (*get_free_color_index)(void); \
  G_CLIP_DESCRIPTOR * \
  (*kind_of_clipping)(void); \
  void (*set_clip_rects)(GRECT_LIST *); \
  double (*get_unit_per_pixel)(void); \
  double (*stringwidth)(char *); \
  double (*stringheight)(char *)
typedef struct controldevice {
  CONTROLDEVICE_STRUCT;
} CONTROLDEVICE;
#define GRAPHICPS_STRUCT \
  GRAPHICDEVICE_STRUCT
typedef struct graphicps {
  GRAPHICPS_STRUCT;
} GRAPHICPS;
#define SUPROP_STRUCT \
  INSTANCE_STRUCT; \
  SUPROP_DEV suprop_dev; \
  int grid_patch
typedef struct suprop {
  SUPROP_STRUCT;
} SUPROP;
#define SCENE_STRUCT \
  INSTANCE_STRUCT; \
  struct scene *next_scene; \
  INSTANCE *object; \
  MATRIX44 object_trans; \
  int matrix_flag; \
  char *method_name; \
  SUPROP *suprop
typedef struct scene {
  SCENE_STRUCT;
} SCENE;
bool_t g_xdr_suprop_dev(XDR *, SUPROP_DEV *);
#define CHAIN_STRUCT \
  SCENE_STRUCT; \
  int flag
typedef struct chain {
  CHAIN_STRUCT;
} CHAIN;
#ifndef __HANDLE_H
#define __HANDLE_H
#define CTL_PATCH_X_MIN_SIZE     28
#define CTL_PATCH_Y_MIN_SIZE     18
#define CTL_PATCH_Y_STD_SIZE     27
#define MANAGER_REDRAW            0
#define MANAGER_DONT_REDRAW       1
#define MANAGER_DISPLIST_REDRAW   MANAGER_REDRAW
typedef enum {
  mfMouseRelative,
  mfNorthWest,
  mfNorthEast,
  mfEastNorth,
  mfEastSouth,
  mfSouthEast,
  mfSouthWest,
  mfWestSouth,
  mfWestNorth
} MENU_FILL_FLAG;
#define MENU_FILL_TOP             mfNorthWest
#define MENU_FILL_BOTTOM          mfSouthEast
#define MENU_MOUSE_RELATIVE       mfMouseRelative
#define MENU_FILL_LEFT            mfWestSouth
#define MENU_FILL_RIGHT           mfNorthEast
#define MENU_AUTO_POSITION        -100.0
#define MENU_MAX_SIZE             -100.0
#define LINE_FLAT                 0
#define LINE_ELEVATED             1
#endif
#define BLACK_COLOR               0
#define RED_COLOR                 1
#define GREEN_COLOR               2
#define YELLOW_COLOR              3
#define BLUE_COLOR                4
#define MAGENTA_COLOR             5
#define CYAN_COLOR                6
#define WHITE_COLOR               7
#define FOREGROUND_COLOR          8
#define BACKGROUND_COLOR          9
#define PRESSED_UPPER_COLOR      10
#define PRESSED_LOWER_COLOR      11
#define TEXT_COLOR               12
#define LINE_COLOR               13
#define LIGHT_ON_COLOR           14
#define NUMBER_OF_USED_COLORS    15
#define NUMBER_OF_COLORS         64
#define RULER_FLOAT               dfFLOAT
#define RULER_INTEGER             dfINT
#define RULER_DOUBLE              dfDOUBLE
#define PLANE_DOUBLE              dfDOUBLE
#define PLANE_INTEGER             dfINT
#define KNOB_DOUBLE               dfDOUBLE
#define KNOB_INTEGER              dfINT
#define SLIDER_DOUBLE             dfDOUBLE
#define SLIDER_INTEGER            dfINT
#define PRESSED                   0
#define UNPRESSED                 1
#define G_USER_SUPPLIED   30
#define G_SYSTEM_SUPPLIED  0
#define MARKED            23
#define NOT_MARKED        24
#define MXBUF                     256
#define MAX_STRLEN                256
#define MAX_STRING_LENGTH         MAX_STRLEN
#define NOT_ACTIVE 0
#define ACTIVE     1
#define INACTIVE   2
#define ON  1
#define OFF 0
#define NO  0
#define YES 1
#ifndef TRUE
#define FALSE  0
#define TRUE   1
#endif
#define G_CROSS_EMPTY   ' '
#define G_CROSS_PUSHED  'p'  /* instance is pushed,
                              * pop it to go there */
#define G_CROSS_DOWN    'd'  /* go down: current instance
                              * has to move left */
#define G_CROSS_NEXT    'n'  /* like next_scene: curr inst
                              * has to move up */
#define G_CROSS_CHAINED 'c'  /* next element in a doubly
                              * linked list, like timesteps.
                              * Since this instance will
                              * show curr inst in its follow,
                              * don't push current instance */
#define G_CHECKFIELD_CUSTOM 0
#define G_CHECKFIELD_RADIO  1
#define G_CHECKFIELD_CBOXES 2
#define G_CHECKFIELD_CHECK_AND 1
#define G_CHECKFIELD_CHECK_XOR 2
#define G_CHECKFIELD_UNCHECK_AND 4
#define G_CHECKFIELD_UNCHECK_XOR 8
typedef enum {
  mbfNone = 0,
  mbfYes = 1,
  mbfNo  = 2,
  mbfCancel = 4,
  mbfOk = 8,
  mbfAbort = 16,
  mbfRecursive = 128
} MSGBOX_FLAGS;
#define G_MSGBOX_BUTTON_YES    mbfYes
#define G_MSGBOX_BUTTON_NO     mbfNo
#define G_MSGBOX_BUTTON_CANCEL mbfCancel
#define G_MSGBOX_BUTTON_OK     mbfOk
#define G_MSGBOX_BUTTON_ABORT  mbfAbort
#define G_MSGBOX_RECURSIVE     mbfRecursive
#define ALERT(condition,message,error_exit) \
  if ( ! (condition))  { \
    g_errorbox (message, __FILE__, __LINE__, # condition); \
    error_exit; \
  }
typedef enum {
  rfOff = OFF,
  rfOn = ON,
  rfDisplist
} REDRAW_FLAG;
typedef enum {
  afCenter = 0,
  afTop = 1,
  afBottom = 2,
  afLeft = 4,
  afRight = 8
} ALIGN_FLAGS;
typedef enum border_flags {
  bfNoBorder = 0x00000000,
  bfBorder   = 0x00000001,
  bfTitle    = 0x00000002
} BORDER_FLAGS;
typedef unsigned long G_HANDLE_FLAGS;
#define G_IF_REL_X ((G_HANDLE_FLAGS)1)
#define G_IF_REL_Y ((G_HANDLE_FLAGS)2)
#define ITEM_STRUCT \
  INSTANCE_STRUCT; \
  struct item *next_entry; \
  struct group *menu_ptr; \
  double sizex,sizey; \
  double posx,posy; \
  double constr_posx,constr_posy; \
  double constr_sizex, \
         constr_sizey; \
  int refresh; \
  G_HANDLE_FLAGS item_flags; \
  int fill; \
  VEC3 color; \
  int color_index; \
  char *add_inter_action; \
  char *remove_inter_action; \
  char *help_text
typedef struct item {
  ITEM_STRUCT;
} ITEM;
#define RECTANGLE_STRUCT \
  ITEM_STRUCT
typedef struct rectangle {
  RECTANGLE_STRUCT;
} RECTANGLE;
#define TEXTMESSAGE_STRUCT \
  ITEM_STRUCT; \
  char *label; \
  int alignment_x, alignment_y; \
  BORDER_FLAGS border_flags
typedef struct textmessage {
  TEXTMESSAGE_STRUCT;
} TEXTMESSAGE;
#define STATICTEXT_STRUCT \
  TEXTMESSAGE_STRUCT
typedef struct statictext {
  STATICTEXT_STRUCT;
} STATICTEXT;
#define LINE_STRUCT \
  ITEM_STRUCT; \
  int type
typedef struct line {
  LINE_STRUCT;
} LINE;
#define INTERACTIVE_STRUCT \
  ITEM_STRUCT; \
  struct side_effect *effects; \
  char *action; \
  char *method; \
  INSTANCE *inst; \
  char *label; \
  struct layer_list *layers; \
  int active; \
  REDRAW_FLAG redraw_flag
typedef struct interactive {
  INTERACTIVE_STRUCT;
} INTERACTIVE;
typedef struct layer_list {
  struct layer *layer;
  struct layer_list *next;
} LAYER_LIST;
typedef struct side_effect {
  struct side_effect *next;
  INTERACTIVE *entry;
} SIDE_EFFECT;
#define BUTTON_STRUCT \
  INTERACTIVE_STRUCT; \
  int on_off; \
  char *draw_method
typedef struct button {
  BUTTON_STRUCT;
} BUTTON;
typedef union value_pointer {
  int *Int;
  float *Float;
  double *Double;
  char *String;
} VALUE_POINTER;
typedef enum declaration_format {
  dfInt,
  dfFloat,
  dfDouble,
  dfString,
  dfGroup,
  dfRegexp
} DECLARATION_FORMAT;
#define dfINT dfInt
#define dfFLOAT dfFloat
#define dfDOUBLE dfDouble
#define dfSTRING dfString
#define dfREGULAR_EXPRESSION dfRegexp
#define INPUT_STRUCT \
  INTERACTIVE_STRUCT; \
  char *evaluation_method; \
  char *inv_evaluation_method; \
  DECLARATION_FORMAT prompt_type; \
  VALUE_POINTER varp; \
  char *prompt_entry; \
  int max_size; \
  int cursor_pos; \
  int offset; \
  int alignment; \
  char *print_format; \
  char *tab_action; \
  int is_keyactive; \
  int new_input
typedef struct input {
  INPUT_STRUCT;
} INPUT;
#define PROMPT_STRUCT \
  BUTTON_STRUCT; \
  char *prompt_display; \
  char *prompt_entry; \
  double reg_posx; \
  double reg_posy; \
  double reg_sizex; \
  double reg_sizey
typedef struct prompt {
  PROMPT_STRUCT;
} PROMPT;
#define XLIST_BUTTON_STRUCT \
  BUTTON_STRUCT
typedef struct xlist_button {
  XLIST_BUTTON_STRUCT;
} XLIST_BUTTON;
typedef struct string_list {
  struct string_list *next;
  char *str;
  int key;
} STRING_LIST;
#define PURELIST_STRUCT \
  INTERACTIVE_STRUCT; \
  int no_of_entries; \
  STRING_LIST *first_entry; \
  int first_row; \
  int last_chosen; \
  double max_width, \
         max_height
typedef struct purelist {
  PURELIST_STRUCT;
} PURELIST;
#define BIGBUTTON_STRUCT \
  BUTTON_STRUCT
typedef struct bigbutton {
  BIGBUTTON_STRUCT;
} BIGBUTTON;
#define STEP_BUTTON_STRUCT \
  BUTTON_STRUCT
typedef struct step_button {
  STEP_BUTTON_STRUCT;
} STEP_BUTTON;
typedef struct {
  int value;
  char *label;
} CYCLE_LABEL;
#define CYCLE_BUTTON_STRUCT \
  BUTTON_STRUCT; \
  int current_label; \
  CYCLE_LABEL *label_list
typedef struct cycle_button {
  CYCLE_BUTTON_STRUCT;
} CYCLE_BUTTON;
typedef struct instance_list {
  INSTANCE *inst;
  struct instance_list *next;
  char type;
} INSTANCE_LIST;
typedef struct method_list {
  void (*funct)();
  int no_of_item;
  struct method_list *next;
  struct item **item;
} METHOD_LIST;
typedef struct inst_meth_list {
  struct inst_meth_list *next;
  INSTANCE *inst;
  char *meth;
  int called;
} INST_METH_LIST;
#define MANAGER_STRUCT \
  INSTANCE_STRUCT; \
  struct group *menu; \
  struct graphicdevice *grdev; \
  struct controldevice *ctldev; \
  METHOD_LIST *meth_list; \
  INSTANCE_LIST *inst_list; \
  INSTANCE *current_instance; \
  char *current_display_method; \
  INST_METH_LIST *cross_dep; \
  int cross_dep_changed; \
  int redraw; \
  int draw_on_off; \
  int item_removed; \
  int new_handle_flag; \
  int button_display_flag; \
  INTERACTIVE *active; \
  INTERACTIVE *keyactive; \
  ITEM *deleted; \
  G_LIST *hotkeys; \
  int ask_on_exit
typedef struct manager {
  MANAGER_STRUCT;
} MANAGER;
#define GROUP_STRUCT \
  INTERACTIVE_STRUCT; \
  MANAGER *mgr; \
  ITEM *entries; \
  ITEM *dirty_entry; \
  struct side_effect \
  *replace_to; \
  struct side_effect \
  *replacer_of; \
  struct replaced_list \
  *replaced_items; \
  GRECT_LIST *free_list; \
  BORDER_FLAGS border
typedef struct group {
  GROUP_STRUCT;
} GROUP;
typedef struct replaced_list {
  struct replaced_list *next;
  INTERACTIVE *entry;
  GROUP *group;
} REPLACED_LIST;
#define LAYER_STRUCT \
  GROUP_STRUCT; \
  struct layer *under; \
  struct layer *over; \
  INTERACTIVE *creator; \
  int new_position; \
  G_HANDLE_FLAGS layer_flags; \
  double last_size_x, last_size_y
typedef struct layer {
  LAYER_STRUCT;
} LAYER;
#define RADIO_STRUCT \
  GROUP_STRUCT; \
  BUTTON *pressed_button
typedef struct radio {
  RADIO_STRUCT;
} RADIO;
#define SELECTOR_STRUCT \
  RADIO_STRUCT; \
  GROUP *selected_group
typedef struct selector {
  SELECTOR_STRUCT;
} SELECTOR;
typedef enum orientations {
  orHorizontal, orVertical
} ORIENTATIONS;
typedef enum bartypes {
  btRuler, btSlider, btBoundedRuler,
  btFunctionRuler, btScrollBar
}
BARTYPES;
typedef enum {
  spUndef   = (int) 0x00000000,
  spLoaded  = (int) 0x00000001,
  spChanged = (int) 0x00000002
} SPLINE_STATUS;
#define SPLINE1D_STRUCT \
  INSTANCE_STRUCT; \
  int degree, max_degree; \
  int num_of_intervals, max_num_of_intervals; \
  double *u, *f, *b, *tang; \
  SPLINE_STATUS status
typedef struct spline1d {
  SPLINE1D_STRUCT;
} SPLINE1D;
#define BEZIER1D_STRUCT \
  INSTANCE_STRUCT; \
  int degree, max_degree; \
  double *u, *b
typedef struct bezier1d {
  BEZIER1D_STRUCT;
} BEZIER1D;
#define SPLINE_EDITOR_STRUCT \
  GROUP_STRUCT; \
  struct puresplineeditor *pure; \
  struct button *help_button, *fix_x_button; \
  struct layer *help_layer
typedef struct spline_editor {
  SPLINE_EDITOR_STRUCT;
} SPLINE_EDITOR;
#define PURESPLINEEDITOR_STRUCT \
  INTERACTIVE_STRUCT; \
  double positionx, positiony; \
  double graph_win_xmin, graph_win_xmax, graph_win_ymin, graph_win_ymax; \
  double umin, umax, fmin, fmax; \
  double mouse_radius; \
  int draw_discr; \
  double unit; \
  double size; \
  double tang_length; \
  SPLINE1D *spline1d; \
  int active_knot, active_part; \
  int active_action; \
  int recompute_window; \
  int draw_knots; \
  double old_zoom_factor, zoom_factor; \
  int x_coord_fixed, y_coord_fixed
typedef struct puresplineeditor {
  PURESPLINEEDITOR_STRUCT;
} PURESPLINEEDITOR;
typedef unsigned long CHECKFIELD_VAR;
#define CHECKBOX_STRUCT \
  BUTTON_STRUCT; \
  int *the_flag; \
  struct checkfield *checkfield; \
  CHECKFIELD_VAR checkfield_check_and_mask, checkfield_check_xor_mask; \
  CHECKFIELD_VAR checkfield_uncheck_and_mask, checkfield_uncheck_xor_mask; \
  CHECKFIELD_VAR checkfield_mask
typedef struct checkbox {
  CHECKBOX_STRUCT;
} CHECKBOX;
#define CHECKFIELD_STRUCT \
  GROUP_STRUCT; \
  CHECKFIELD_VAR *varp; \
  CHECKFIELD_VAR used_numbers; \
  int mask_mode; \
  CHECKFIELD_VAR checkfield_check_and_mask, checkfield_check_xor_mask; \
  CHECKFIELD_VAR checkfield_uncheck_and_mask, checkfield_uncheck_xor_mask
typedef struct checkfield {
  CHECKFIELD_STRUCT;
} CHECKFIELD;
#define PURERULER_STRUCT \
  INTERACTIVE_STRUCT; \
  double position; \
  VEC3 bcolor; \
  int bcolor_index
typedef struct pureruler {
  PURERULER_STRUCT;
} PURERULER;
#define PURESLIDER_STRUCT \
  PURERULER_STRUCT
typedef struct pureslider {
  PURESLIDER_STRUCT;
} PURESLIDER;
#define PURESCROLLBAR_STRUCT \
  PURESLIDER_STRUCT; \
  double bar_size
typedef struct purescrollbar {
  PURESCROLLBAR_STRUCT;
} PURESCROLLBAR;
#define BAR1D_STRUCT \
  GROUP_STRUCT; \
  int has_spline_editor; \
  int has_arrows; \
  int has_prompt; \
  int has_layer; \
  ORIENTATIONS orientation; \
  BARTYPES type_of_ruler; \
  PURERULER *ruler; \
  double (*eval)(double); \
  double (*inv_eval)(double); \
  char *value_changed; \
  DECLARATION_FORMAT var_type; \
  VALUE_POINTER var_adr; \
  double offset; \
  double step_size; \
  double scale; \
  double min_value, \
         max_value; \
  struct layer *config_layer; \
  G_HANDLE_FLAGS bar1d_flags
typedef struct bar1d {
  BAR1D_STRUCT;
} BAR1D;
#define RULER_STRUCT \
  BAR1D_STRUCT
typedef struct ruler {
  RULER_STRUCT;
} RULER;
#define SLIDER_STRUCT \
  RULER_STRUCT
typedef struct slider {
  SLIDER_STRUCT;
} SLIDER;
#define SCROLLBAR_STRUCT \
  BAR1D_STRUCT; \
  double visible, gesamt
typedef struct scrollbar {
  SCROLLBAR_STRUCT;
} SCROLLBAR;
#define FUNCTION_BAR1D_STRUCT \
  SLIDER_STRUCT; \
  SPLINE1D *spline; \
  LAYER *edit_layer; \
  BUTTON *edit_button; \
  double *function_arg; \
  int effect_on_off; \
  int last_on_off; \
  double last_val
typedef struct function_bar1d {
  FUNCTION_BAR1D_STRUCT;
} FUNCTION_BAR1D;
#define FUNCTION_RULER_STRUCT \
  FUNCTION_BAR1D_STRUCT
typedef struct function_ruler {
  FUNCTION_RULER_STRUCT;
} FUNCTION_RULER;
#define FUNCTION_SLIDER_STRUCT \
  FUNCTION_BAR1D_STRUCT
typedef struct function_slider {
  FUNCTION_SLIDER_STRUCT;
} FUNCTION_SLIDER;
#define BOUNDED_RULER_STRUCT \
  RULER_STRUCT
typedef struct bounded_ruler {
  BOUNDED_RULER_STRUCT;
} BOUNDED_RULER;
typedef enum sort_method {
  smUnsorted,
  smAlphabetic,
  smFileList,
  smMethodList
} SORT_METHOD;
#define XLIST_STRUCT \
  GROUP_STRUCT; \
  PURELIST *list; \
  SCROLLBAR *scrollbar; \
  int no_visible; \
  char *chosen; \
  char *edt_txt; \
  int edt_key, chosen_key; \
  SORT_METHOD sorting
typedef struct xlist {
  XLIST_STRUCT;
} XLIST;
typedef enum open_mode {
  omRead,
  omWrite,
  omWriteRead
} OPEN_MODE;
#define FILELIST_STRUCT \
  XLIST_STRUCT; \
  char *home_dir; \
  char **actuell_dir; \
  char *own_dir; \
  char *wildcard; \
  OPEN_MODE open_mode
typedef struct filelist {
  FILELIST_STRUCT;
} FILELIST;
#define METHLIST_STRUCT \
  XLIST_STRUCT; \
  CLASS *current_class; \
  CLASS *initial_class; \
  char *meth_end; \
  int no_of_adds; \
  char **adds
typedef struct methlist {
  METHLIST_STRUCT;
} METHLIST;
typedef enum {
  mltDisplay = 1,
  mltSend    = 2,
  mltAny     = 4
} METHLAYER_TYPE;
#define METHLAYER_STRUCT \
  LAYER_STRUCT; \
  METHLAYER_TYPE type
typedef struct methlayer {
  METHLAYER_STRUCT;
} METHLAYER;
#define PLANE_STRUCT \
  INTERACTIVE_STRUCT; \
  char *evaluation_method; \
  char *inv_evaluation_method; \
  void *varx, *vary; \
  double tmp_varx, tmp_vary; \
  double offsetx, offsety; \
  double positionx, positiony; \
  DECLARATION_FORMAT type
typedef struct plane {
  PLANE_STRUCT;
} PLANE;
#define COLOR_SEL_STRUCT \
  GROUP_STRUCT; \
  char *evaluation_method; \
  char *inv_evaluation_method; \
  double *rgb_value; \
  VEC3 tmp_var, last_rgb_value; \
  BAR1D *component[3]; \
  BAR1D *hsv_component[3]; \
  BAR1D *rgb_component[3]; \
  int color_model
typedef struct color_sel {
  COLOR_SEL_STRUCT;
} COLOR_SEL;
#define SPHERE_STRUCT \
  INTERACTIVE_STRUCT; \
  char *evaluation_method; \
  char *inv_evaluation_method; \
  double *varx, *vary; \
  double tmp_varx, tmp_vary; \
  double offsetx, offsety; \
  double positionx, positiony; \
  double centerx, centery, centerz
typedef struct sphere {
  SPHERE_STRUCT;
} SPHERE;
#define GD_SPHERE_STRUCT \
  SPHERE_STRUCT
typedef struct gd_sphere {
  GD_SPHERE_STRUCT;
} GD_SPHERE;
typedef struct Mat_Edt_vals {
  double trans[3];
  double scale_dist;
  double rot[3];
  double spx, spy;
  VEC3 point;
  double plxy[2];
} MAT_EDT_VALS;
typedef enum {
  mefScale = 0x00000001,
  mefDist  = 0x00000002,
  mefJustRot = 0x00000004
} MAT_EDIT_FLAGS;
#define MATRIX_EDIT_STRUCT \
  GROUP_STRUCT; \
  VEC4 *matrix; \
  MAT_EDT_VALS new_vals; \
  MAT_EDT_VALS old_vals; \
  MAT_EDIT_FLAGS tflags; \
  GROUP *trans_ruls; \
  RULER *scale_dist_rul; \
  GROUP *rot_ruls; \
  GROUP *rots; \
  SELECTOR *trans_sel; \
  BUTTON *reset_butt; \
  SPHERE *sphere
typedef struct matrix_edit {
  MATRIX_EDIT_STRUCT;
} MATRIX_EDIT;
#define CAMERA_STRUCT \
  INSTANCE_STRUCT; \
  VEC3 eyepoint; \
  VEC3 lookpoint; \
  VEC3 up; \
  double focal; \
  double far,near
typedef struct camera {
  CAMERA_STRUCT;
} CAMERA;
extern double g_maintime;
extern int g_batch_mode;
int g_has_superclass(INSTANCE *, CLASS *);
void g_hide_layer (ITEM *);
int g_msgbox (char *, char *, MSGBOX_FLAGS);
#define msgbox g_msgbox
void g_infobox (char *, char *);
#define infobox g_infobox
void g_errorbox (char *, char *, int, char *);
double g_eval_ident_fn (double x);
double g_eval_angle_fn (double x);
double g_eval_exp_fn (double x);
double g_eval_log_fn (double x);
#define TEXTDEVICE_STRUCT \
  DEVICE_STRUCT; \
  FILE *fp; \
  int term_flags; \
  void (*query)(); \
  int (*seek )ANSI_PARM(( long int, int )); \
  int (*tell )ANSI_PARM(( long int * )); \
  int (*write_char   )ANSI_PARM(( char )); \
  int (*write_bytes  )ANSI_PARM(( size_t, void * )); \
  int (*write_string )ANSI_PARM(( char * )); \
  int (*write_short  )ANSI_PARM(( short int )); \
  int (*write_int    )ANSI_PARM(( int )); \
  int (*write_long   )ANSI_PARM(( long int )); \
  int (*write_float  )ANSI_PARM(( double )); \
  int (*write_double )ANSI_PARM(( double )); \
  int (*read_char   )ANSI_PARM(( char * )); \
  int (*read_bytes  )ANSI_PARM(( size_t, void * )); \
  int (*read_string )ANSI_PARM(( char * )); \
  int (*read_short  )ANSI_PARM(( short int * )); \
  int (*read_int    )ANSI_PARM(( int * )); \
  int (*read_long   )ANSI_PARM(( long int * )); \
  int (*read_float  )ANSI_PARM(( float * )); \
  int (*read_double )ANSI_PARM(( double * ))
typedef struct textdevice {
  TEXTDEVICE_STRUCT;
} TEXTDEVICE;
#define INPUTFILE_STRUCT \
  TEXTDEVICE_STRUCT
typedef struct inputfile {
  INPUTFILE_STRUCT;
} INPUTFILE;
#define OUTPUTFILE_STRUCT \
  TEXTDEVICE_STRUCT
typedef struct outputfile {
  OUTPUTFILE_STRUCT;
} OUTPUTFILE;
#define INOUTFILE_STRUCT \
  TEXTDEVICE_STRUCT
typedef struct inoutfile {
  INOUTFILE_STRUCT;
} INOUTFILE;
#define RESOURCEDEV_STRUCT \
  DEVICE_STRUCT; \
  char *resource; \
  G_LIST *filenames; \
  G_LIST *rcpath
typedef struct resourcedev {
  RESOURCEDEV_STRUCT;
} RESOURCEDEV;
typedef struct g_resource_descr {
  char *ident;
  void *value;
  DECLARATION_FORMAT type;
  int index;
  int count;
  int sep, last_sep;
} G_RESOURCE_DESCR;
void g_remove_resourcelist (G_LIST *);
void g_resourcelist_print (G_LIST *);
int g_resourcelist_find (G_LIST *, G_RESOURCE_DESCR *, int *);
G_LIST *g_rclist_add (G_LIST *, char *, void *,
                      DECLARATION_FORMAT);
void g_resourcedev_not_found (INSTANCE *, G_RESOURCE_DESCR *);
void g_resourcedev_parse (INSTANCE *, char *);
typedef enum {
  I_End = -99999,
  I_Name,
  I_Size,
  I_SizeX,
  I_SizeY,
  I_RSize,
  I_RSizeX,
  I_RSizeY,
  I_Pos,
  I_PosX,
  I_PosY,
  I_Color,
  I_ColorRGB,
  I_ColorIndex,
  I_FillMode,
  I_AddAction,
  I_DrawMethod,
  I_State,
  I_Align,
  I_Text,
  I_Action,
  I_Instance,
  I_Self,
  I_Method,
  I_Label,
  I_HelpText,
  I_Redraw,
  I_BarColor,
  I_BarColorRGB,
  I_BarColorIndex,
  I_EvalFn,
  I_MinMax,
  I_Min,
  I_Max,
  I_Offset,
  I_Scale,
  I_StepSize,
  I_ValueChanged,
  I_PrintFormat,
  I_Var,
  I_TabAction,
  I_Eval,
  I_Type,
  I_Item,
  I_Border,
  I_ButtonList,
  I_RulerList,
  I_SliderList,
  I_FixSize,
  I_Effect
} ITEM_TAG;
ITEM *new_item(CLASS *, ...);
typedef enum {
  P_End = -88888,
  P_Project,
  P_ProjectList,
  P_Class,
  P_ClassList,
  P_ClassCN,
  P_Method,
  P_MethodList,
  P_MethodCN,
  P_Item,
  P_ItemOpt,
  P_ButtonList,
  P_RulerList,
  P_SliderList,
  P_InitMethod,
  P_ExitMethod,
  P_AddMethod,
  P_Version,
  P_Date,
  P_Author
} PROJECT_TAG;
typedef struct {
  char *name;
  char *setup;
  void *(*func)();
  int hidden;
} PRTYPE;
typedef struct {
  char *name;
} PROJECT_DESCR;
typedef struct {
  CLASS **superclass;
  CLASS **class;
  char *name;
  int size;
} CLASS_DESCR;
typedef struct {
  CLASS **class;
  char *name;
  void *(*func)();
} METHOD_DESCR;
typedef struct {
  char *methname;
  INSTANCE **obj;
  char *txt;
  double l, h, x, y;
  int optmenu;
} BUTTON_DESCR;
typedef struct {
  void *var;
  char *txt;
  double x, y;
  char *action;
  int type;
  int optmenu;
} RULER_DESCR;
typedef struct {
  void *var;
  char *txt;
  double x, y;
  char *action;
  int type;
  int optmenu;
} SLIDER_DESCR;
#define PROJECT_STRUCT \
  INSTANCE_STRUCT; \
  int status; \
  char *setup; \
  int used; \
  int added; \
  G_LIST *projects; \
  G_LIST *classes; \
  G_LIST *methods; \
  G_LIST *items; \
  char *init; \
  char *exit; \
  char *add; \
  INSTANCE *object; \
  char *path; \
  int version; \
  char *date; \
  char *author; \
  void *load_descr
typedef struct project {
  PROJECT_STRUCT;
} PROJECT;
void g_project_include(PRTYPE *);
PROJECT *g_project_use(char *, int);
PROJECT *g_project_unuse(char *);
PROJECT *g_project_add(char *);
#define TRIANG1D_STRUCT \
  INSTANCE_STRUCT; \
  int number_of_points; \
  int max_number_of_points; \
  double *x, *y, *z
typedef struct triang1d {
  TRIANG1D_STRUCT;
} TRIANG1D;
typedef struct plane_parm_type {
  double x,y,z,d;
} PLANE_PARM;
typedef struct trace_dat {
  struct trace_dat *next;
  double coord[3];
  int element;
  double bary[4];
} TRACE_DAT;
typedef struct bnd_function {
  double (*g)();
  void *par;
} BND_FUNCTION;
#define TRIANG2D_STRUCT \
  INSTANCE_STRUCT; \
  int number_of_points; \
  int max_number_of_points; \
  double *x, *y, *z; \
  int number_of_elements; \
  int max_number_of_elements; \
  INT3 *vertex; \
  INT3 *neighbour
typedef struct triang2d {
  TRIANG2D_STRUCT;
} TRIANG2D;
#define FE2D_STRUCT \
  TRIANG2D_STRUCT; \
  int dimension_of_value; \
  int polynomial_order; \
  void (*f)(struct fe2d*, int, VEC3, double*); \
  int size_of_data; \
  char *data
typedef struct fe2d {
  FE2D_STRUCT;
} FE2D;
#define COLORBAR_STRUCT \
  INSTANCE_STRUCT; \
  INSTANCE (*funct)(); \
  INSTANCE *object; \
  char name_of_function[100]; \
  VEC3 xyz; \
  double height; \
  double width; \
  int display_on_off; \
  double eps; \
  double offset; \
  double color; \
  double min,max; \
  int number_of_values; \
  double values[100]; \
  double input_value; \
  char double_display[MAX_STRLEN]; \
  XLIST *xlist; \
  GROUP *group; \
  LAYER  *layer_options; \
  LAYER  *layer_values
typedef struct colorbar {
  COLORBAR_STRUCT;
} COLORBAR;
typedef void (f2_func)(FE2D*, int, double*, double*);
VEC2 *g_fe2d_grad ANSI_PARM ((FE2D*,int,double*,f2_func*,int,VEC2*));
int make_hsv_color(double , VEC3 );
#define TRIANG3D_STRUCT \
  INSTANCE_STRUCT; \
  int number_of_points; \
  int max_number_of_points; \
  double *x, *y, *z; \
  int number_of_elements; \
  int max_number_of_elements; \
  INT4 *vertex; \
  INT4 *neighbour
typedef struct triang3d {
  TRIANG3D_STRUCT;
} TRIANG3D;
#define FE3D_STRUCT \
  TRIANG3D_STRUCT; \
  int dimension_of_value; \
  int polynomial_order; \
  void (*f)(struct fe3d *, int, VEC4, double*); \
  int size_of_data; \
  char *data
typedef struct fe3d {
  FE3D_STRUCT;
} FE3D;
typedef struct clip3d_par {
  TRIANG3D *t;
  int it;
  int n;
  VEC3 v[5];
  double bary[5][4];
  double (*f)();
  char *var;
  char *scal;
  int flag ;
} CLIP3D_PAR;
typedef void (f3_func)(FE3D*, int, double*, double*);
VEC3 *g_fe3d_gradient ANSI_PARM ((FE3D*,int,double*,int,VEC3*));
VEC3 *g_fe3d_grad ANSI_PARM ((FE3D*,int,double*,f3_func*,int,VEC3*));
#define TIMESCENE_STRUCT \
  SCENE_STRUCT; \
  double time,object_time; \
  INSTANCE *dynamic; \
  char *dynamic_method_name; \
  int sync
typedef struct timescene {
  TIMESCENE_STRUCT;
} TIMESCENE;
#define TIMESTEP_STRUCT \
  INSTANCE_STRUCT; \
  double time; \
  struct timestep *pre_step,*post_step; \
  INSTANCE *pre_object,*post_object
typedef struct timestep {
  TIMESTEP_STRUCT;
} TIMESTEP;
typedef struct {
  double x, y;
} COMPLEX;
typedef COMPLEX CVEC3[3];
void g_dmatrix44_mult(MATRIX44, MATRIX44, MATRIX44);
void g_dmatrix44_mult_dvec3(MATRIX44, VEC3);
void g_dmatrix44_mult_dvec4(MATRIX44, VEC4);
void g_dmatrix44_mult_dvec4_dvec4(MATRIX44, VEC4, VEC4);
void g_dmatrix44_set_identity(MATRIX44);
void g_dmatrix44_transp(MATRIX44, MATRIX44);
int g_vec3_get_normal ANSI_PARM((VEC3, VEC3, VEC3));
int g_vec3_get_normal_to_plane ANSI_PARM((VEC3, VEC3, VEC3, VEC3));
int g_vec3_proj_onto_line ANSI_PARM((VEC3 ,VEC3));
int g_vec3_proj_onto_plane ANSI_PARM((VEC3 , VEC3));
void g_vec_change_value_int3 ANSI_PARM((INT3 *, int, int, int));
double g_vec3_sqrabs ANSI_PARM((VEC3));
double g_vec3_sqrdist ANSI_PARM((VEC3, VEC3));
void g_int3_assign ANSI_PARM((INT3, INT3));
void g_int3_list_assign ANSI_PARM((INT3 *, INT3 *, int));
void g_vec3_list_assign ANSI_PARM((VEC3 *, VEC3 *, int));
void g_vec3_set ANSI_PARM((VEC3, double, double, double));
void g_vec3_invert ANSI_PARM((VEC3));
void g_vec3_set_zero ANSI_PARM((VEC3));
void g_vec3_assign ANSI_PARM((VEC3, VEC3));
void g_vec3_sub ANSI_PARM((VEC3, VEC3, VEC3));
void g_vec3_add ANSI_PARM((VEC3, VEC3, VEC3));
void g_vec3_normalize ANSI_PARM((VEC3));
void g_vec3_set_length ANSI_PARM((VEC3 , double));
void g_vec3_crossprod ANSI_PARM((VEC3, VEC3, VEC3));
void g_vec3_mult ANSI_PARM((VEC3, double));
double g_vec3_abs ANSI_PARM((VEC3));
double g_vec3_skp ANSI_PARM((VEC3, VEC3));
double g_vec3_dist ANSI_PARM((VEC3, VEC3));
int g_vec3_cmp ANSI_PARM((VEC3 , VEC3, double));
double g_vec3_angle ANSI_PARM((VEC3, VEC3, VEC3));
double g_vec3_ctg ANSI_PARM((VEC3, VEC3));
void g_swap_int ANSI_PARM((int *, int *));
void g_swap_double ANSI_PARM((double *, double *));
void g_swap_vec3 ANSI_PARM((VEC3 , VEC3));
double g_matrix33_det ANSI_PARM((MATRIX33));
void g_matrix44_set_identity ANSI_PARM((MATRIX44));
void g_matrix44_mult ANSI_PARM((MATRIX44, MATRIX44, MATRIX44));
void g_matrix44_transp ANSI_PARM((MATRIX44, MATRIX44));
void g_matrix44_assign ANSI_PARM((MATRIX44, MATRIX44));
void g_matrix44_mult_vec3 ANSI_PARM((MATRIX44, VEC3));
double g_matrix44_det ANSI_PARM((MATRIX44));
int g_matrix44_inv ANSI_PARM((MATRIX44, MATRIX44));
int g_matrix44_invert_gauss ANSI_PARM((MATRIX44, MATRIX44));
int g_vec_cmp_double ANSI_PARM((double * , double *, int , double));
void g_vec_change_int3 ANSI_PARM((INT3 *, int, int *, int *, int));
void g_vec_change_int ANSI_PARM((int *, int, int *, int *, int));
void g_vec_change_value_int ANSI_PARM((int *, int, int, int));
VEC3 *g_vec_append_vec3 ANSI_PARM((VEC3 **, int *, int *, VEC3 *, int));
INT3 *g_vec_append_int3 ANSI_PARM((INT3 **, int *, int *, INT3 *, int));
double *g_vec_append_double ANSI_PARM((double **, int *, int *, double *, int));
int *g_vec_append_int ANSI_PARM((int **, int *, int *, int *, int));
INSTANCE **g_vec_append_instance
ANSI_PARM((INSTANCE ***, int *, int *, INSTANCE **, int));
int g_get_vecpos_int ANSI_PARM((int, int *, int));
double g_vec_max_double ANSI_PARM((double *, int));
double g_vec_min_double ANSI_PARM((double *, int));
int g_vec_max_int ANSI_PARM((int *, int));
int g_vec_min_int ANSI_PARM((int *, int));
void g_get_vecpos ANSI_PARM ((double *, int, double, int *, double *));
#define LIST_OF_INST_STRUCT \
  INSTANCE_STRUCT; \
  int num_of_objects; \
  int max_num_of_objects; \
  int *flags; \
  INSTANCE **objects
typedef struct list_of_inst {
  LIST_OF_INST_STRUCT;
} LIST_OF_INST;
#define LIST_OF_DOUBLE_STRUCT \
  INSTANCE_STRUCT; \
  int *flag; \
  int max_num_of_double; \
  int num_of_double; \
  double *vec
typedef struct list_of_double {
  LIST_OF_DOUBLE_STRUCT;
} LIST_OF_DOUBLE;
#define LIST_OF_INT_STRUCT \
  INSTANCE_STRUCT; \
  int *flag; \
  int max_num_of_int; \
  int num_of_int; \
  int *vec
typedef struct list_of_int {
  LIST_OF_INT_STRUCT;
} LIST_OF_INT;
#define LIST_OF_VEC3_STRUCT \
  INSTANCE_STRUCT; \
  int *flag; \
  int max_num_of_vec3; \
  int num_of_vec3; \
  VEC3 *vec
typedef struct list_of_vec3 {
  LIST_OF_VEC3_STRUCT;
} LIST_OF_VEC3;
LIST_OF_INST **list_of_inst_p_alloc ANSI_PARM((size_t));
LIST_OF_INST **list_of_inst_p_free  ANSI_PARM((LIST_OF_INST **,\
                                               size_t));
LIST_OF_INST **list_of_inst_p_realloc ANSI_PARM((LIST_OF_INST **,\
                                                 size_t, size_t));
double g_vec2_dist ANSI_PARM((double *p, double *q));
double g_vec2_sqrdist ANSI_PARM((double *p, double *q));
double g_vec2_sqrabs ANSI_PARM((double *v));
void g_vec2_assign ANSI_PARM((VEC2 d, VEC2 s));
void g_vec2_dec ANSI_PARM((VEC2 d, VEC2 s));
void g_vec2_sub ANSI_PARM((VEC2 d, VEC2 a, VEC2 b));
double g_vec2_skp ANSI_PARM((VEC2 a, VEC2 b));
double g_vec2_abs ANSI_PARM((VEC2 a));
double g_vec2_sub_skp ANSI_PARM((VEC2 a, VEC2 b, VEC2 n));
int g_vec2_normalize ANSI_PARM((VEC2 n));
int g_vec2_get_normal ANSI_PARM((VEC2 normal, VEC2 a, VEC2 b));
double g_dvec3_abs ANSI_PARM((VEC3 p));
double g_dvec3_circle ANSI_PARM((VEC3 c, VEC3 p[3]));
double g_dvec3_dist ANSI_PARM((VEC3 p, VEC3 q));
double g_dvec3_skp ANSI_PARM((VEC3 p, VEC3 q));
double g_dvec3_sqrabs ANSI_PARM((VEC3 p));
double g_dvec3_sqrdist ANSI_PARM((VEC3 p, VEC3 q));
int g_d_line_rotate ANSI_PARM((VEC3 pkt[2], double fi, MATRIX44 dmatrix));
int g_dvec3_get_normal ANSI_PARM((VEC3 normal, VEC3 v, VEC3 w));
int g_dvec3_get_normal_to_plane ANSI_PARM((VEC3 normal, VEC3 p, VEC3 q, VEC3 r));
int g_dvec3_normalize ANSI_PARM((VEC3 p));
int g_dvec3_proj_onto_line ANSI_PARM((VEC3 v, VEC3 line));
int g_dvec3_proj_onto_plane ANSI_PARM((VEC3 v, VEC3 normal));
int g_dvec3_set_length ANSI_PARM((VEC3 p, double length));
void g_dvec3_add ANSI_PARM((VEC3 r, VEC3 p, VEC3 q));
void g_dvec3_assign ANSI_PARM((VEC3 r, VEC3 p));
void g_dvec3_crossprod ANSI_PARM((VEC3 r, VEC3 p, VEC3 q));
void g_dvec3_sub ANSI_PARM((VEC3 r, VEC3 p, VEC3 q));
void g_matrix44_mult_dvec3 ANSI_PARM((MATRIX44 m, VEC3 p));
double *g_dvec3_mult ANSI_PARM((VEC3 v, double r));
double *g_dvec3_cv_vec3 ANSI_PARM((VEC3 v, VEC3 w));
double *g_vec3_cv_dvec3 ANSI_PARM((VEC3 v, VEC3 w));
int g_gcd ANSI_PARM((int , int));
double g_vec3_circle ANSI_PARM((VEC3 center, VEC3 p[3]));
#ifndef G_MAX
#define G_MAX(A,B)        ( (A) > (B) ? (A) : (B) )
#endif
#ifndef G_MIN
#define G_MIN(A,B)        ( (A) < (B) ? (A) : (B) )
#endif
#ifndef G_ABS
#define G_ABS(A)        ( 0 < (A) ? (A) : (-(A)) )
#endif
#ifndef G_SGN
#define G_SGN(A)       ( (A) < 0 ? (-1) : 1 )
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#define M_PI_2 1.570796326794896619
#endif
void *g_smart_mem_realloc ANSI_PARM((void *, size_t, size_t, size_t));
double *double_alloc ANSI_PARM((size_t));
double *double_free  ANSI_PARM((double *, size_t));
double *double_realloc ANSI_PARM((double *, size_t, size_t));
double *g_smart_double_realloc ANSI_PARM((double *, size_t, \
                                          size_t, size_t));
double **double_p_alloc ANSI_PARM((size_t));
double **double_p_free  ANSI_PARM((double **, size_t));
double **double_p_realloc ANSI_PARM((double **, size_t, size_t));
double *double_alloc ANSI_PARM((size_t));
double *double_free  ANSI_PARM((double *, size_t));
double *double_realloc ANSI_PARM((double *, size_t, size_t));
double **double_p_alloc ANSI_PARM((size_t));
double **double_p_free  ANSI_PARM((double **, size_t));
double **double_p_realloc ANSI_PARM((double **, size_t, size_t));
char *char_alloc ANSI_PARM((size_t));
char *char_free  ANSI_PARM((char *, size_t));
char *char_realloc ANSI_PARM((char *, size_t, size_t));
char **char_p_alloc ANSI_PARM((size_t));
char **char_p_free  ANSI_PARM((char **, size_t));
char **char_p_realloc ANSI_PARM((char **, size_t, size_t));
int *int_alloc ANSI_PARM((size_t));
int *int_free  ANSI_PARM((int *, size_t));
int *int_realloc ANSI_PARM((int *, size_t, size_t));
int **int_p_alloc ANSI_PARM((size_t));
int **int_p_free  ANSI_PARM((int **, size_t));
int **int_p_realloc ANSI_PARM((int **, size_t, size_t));
INT2 *int2_alloc ANSI_PARM((size_t));
INT2 *int2_free  ANSI_PARM((INT2 *, size_t));
INT2 *int2_realloc ANSI_PARM((INT2 *, size_t, size_t));
INT2 **int2_p_alloc ANSI_PARM((size_t));
INT2 **int2_p_free  ANSI_PARM((INT2 **, size_t));
INT2 **int2_p_realloc ANSI_PARM((INT2 **, size_t, size_t));
INT3 *int3_alloc ANSI_PARM((size_t));
INT3 *int3_free  ANSI_PARM((INT3 *, size_t));
INT3 *int3_realloc ANSI_PARM((INT3 *, size_t, size_t));
INT3 *g_smart_int3_realloc ANSI_PARM((INT3 *, size_t, size_t, size_t));
INT3 **int3_p_alloc ANSI_PARM((size_t));
INT3 **int3_p_free  ANSI_PARM((INT3 **, size_t));
INT3 **int3_p_realloc ANSI_PARM((INT3 **, size_t, size_t));
INT4 *int4_alloc ANSI_PARM((size_t));
INT4 *int4_free  ANSI_PARM((INT4 *, size_t));
INT4 *int4_realloc ANSI_PARM((INT4 *, size_t, size_t));
INT4 *g_smart_int4_realloc ANSI_PARM((INT4 *, size_t, size_t, size_t));
INT4 **int4_p_alloc ANSI_PARM((size_t));
INT4 **int4_p_free  ANSI_PARM((INT4 **, size_t));
INT4 **int4_p_realloc ANSI_PARM((INT4 **, size_t, size_t));
VEC2 *vec2_alloc ANSI_PARM((size_t));
VEC2 *vec2_free  ANSI_PARM((VEC2 *, size_t));
VEC2 *vec2_realloc ANSI_PARM((VEC2 *, size_t, size_t));
VEC2 **vec2_p_alloc ANSI_PARM((size_t));
VEC2 **vec2_p_free  ANSI_PARM((VEC2 **, size_t));
VEC2 **vec2_p_realloc ANSI_PARM((VEC2 **, size_t, size_t));
VEC3 *vec3_alloc ANSI_PARM((size_t));
VEC3 *vec3_free  ANSI_PARM((VEC3 *, size_t));
VEC3 *vec3_realloc ANSI_PARM((VEC3 *, size_t, size_t));
VEC3 **vec3_p_alloc ANSI_PARM((size_t));
VEC3 **vec3_p_free  ANSI_PARM((VEC3 **, size_t));
VEC3 **vec3_p_realloc ANSI_PARM((VEC3 **, size_t, size_t));
VEC3 *dvec3_alloc ANSI_PARM((size_t));
VEC3 *dvec3_free  ANSI_PARM((VEC3 *, size_t));
VEC3 *dvec3_realloc ANSI_PARM((VEC3 *, size_t, size_t));
VEC3 **dvec3_p_alloc ANSI_PARM((size_t));
VEC3 **dvec3_p_free  ANSI_PARM((VEC3 **, size_t));
VEC3 **dvec3_p_realloc ANSI_PARM((VEC3 **, size_t, size_t));
VEC4 *vec4_alloc ANSI_PARM((size_t));
VEC4 *vec4_free  ANSI_PARM((VEC4 *, size_t));
VEC4 *vec4_realloc ANSI_PARM((VEC4 *, size_t, size_t));
VEC4 **vec4_p_alloc ANSI_PARM((size_t));
VEC4 **vec4_p_free  ANSI_PARM((VEC4 **, size_t));
VEC4 **vec4_p_realloc ANSI_PARM((VEC4 **, size_t, size_t));
VEC4 *dvec4_alloc ANSI_PARM((size_t));
VEC4 *dvec4_free  ANSI_PARM((VEC4 *, size_t));
VEC4 *dvec4_realloc ANSI_PARM((VEC4 *, size_t, size_t));
VEC4 **dvec4_p_alloc ANSI_PARM((size_t));
VEC4 **dvec4_p_free  ANSI_PARM((VEC4 **, size_t));
VEC4 **dvec4_p_realloc ANSI_PARM((VEC4 **, size_t, size_t));
double (**dfnp_alloc ANSI_PARM((size_t))) ();
double (**dfnp_free  ANSI_PARM((double (**)(), size_t))) ();
double (**dfnp_realloc ANSI_PARM((double (**)(), size_t, size_t))) ();
INSTANCE **instance_p_alloc ANSI_PARM((size_t));
INSTANCE **instance_p_free  ANSI_PARM((INSTANCE **, size_t));
INSTANCE **instance_p_realloc ANSI_PARM((INSTANCE **, size_t, size_t));
CHAIN **chain_p_alloc ANSI_PARM((size_t));
CHAIN **chain_p_free  ANSI_PARM((CHAIN **, size_t));
CHAIN **chain_p_realloc ANSI_PARM((CHAIN **, size_t, size_t));
SPLINE1D **spline1d_p_alloc ANSI_PARM((size_t));
SPLINE1D **spline1d_p_free  ANSI_PARM((SPLINE1D **, size_t));
SPLINE1D **spline1d_p_realloc ANSI_PARM((SPLINE1D **, size_t, size_t));
TRIANG1D **triang1d_p_alloc ANSI_PARM((size_t));
TRIANG1D **triang1d_p_free  ANSI_PARM((TRIANG1D **, size_t));
TRIANG1D **triang1d_p_realloc ANSI_PARM((TRIANG1D **, size_t, size_t));
TRIANG2D **triang2d_p_alloc ANSI_PARM((size_t));
TRIANG2D **triang2d_p_free  ANSI_PARM((TRIANG2D **, size_t));
TRIANG2D **triang2d_p_realloc ANSI_PARM((TRIANG2D **, size_t, size_t));
TRIANG1D ***triang1d_pp_alloc ANSI_PARM((size_t));
TRIANG1D ***triang1d_pp_free  ANSI_PARM((TRIANG1D ***, size_t));
TRIANG1D ***triang1d_pp_realloc ANSI_PARM((TRIANG1D ***, size_t, size_t));
TRIANG2D ***triang2d_pp_alloc ANSI_PARM((size_t));
TRIANG2D ***triang2d_pp_free  ANSI_PARM((TRIANG2D ***, size_t));
TRIANG2D ***triang2d_pp_realloc ANSI_PARM((TRIANG2D ***, size_t, size_t));
TRIANG1D ****triang1d_ppp_alloc ANSI_PARM((size_t));
TRIANG1D ****triang1d_ppp_free  ANSI_PARM((TRIANG1D ****, size_t));
TRIANG1D ****triang1d_ppp_realloc ANSI_PARM((TRIANG1D ****, size_t, size_t));
TRIANG2D ****triang2d_ppp_alloc ANSI_PARM((size_t));
TRIANG2D ****triang2d_ppp_free  ANSI_PARM((TRIANG2D ****, size_t));
TRIANG2D ****triang2d_ppp_realloc ANSI_PARM((TRIANG2D ****, size_t, size_t));
FUNCTION_RULER **function_ruler_p_alloc ANSI_PARM((size_t));
FUNCTION_RULER **function_ruler_p_free  ANSI_PARM((FUNCTION_RULER **, size_t));
FUNCTION_RULER **function_ruler_p_realloc ANSI_PARM((FUNCTION_RULER **, size_t, size_t));
BOUNDED_RULER **bounded_ruler_p_alloc ANSI_PARM((size_t));
BOUNDED_RULER **bounded_ruler_p_free ANSI_PARM((BOUNDED_RULER **, size_t));
BOUNDED_RULER **bounded_ruler_p_realloc ANSI_PARM((BOUNDED_RULER **, size_t, size_t));
RULER **ruler_p_alloc ANSI_PARM((size_t));
RULER **ruler_p_free  ANSI_PARM((RULER **, size_t));
RULER **ruler_p_realloc ANSI_PARM((RULER **, size_t, size_t));
void g_show_bnd_box(INSTANCE *);
double g_dvec4_abs ANSI_PARM((VEC4));
double g_dvec4_sqrabs ANSI_PARM((VEC4));
double g_dvec4_skp ANSI_PARM((VEC4, VEC4));
double g_dvec4_dist ANSI_PARM((VEC4, VEC4));
double g_dvec4_sqrdist ANSI_PARM((VEC4, VEC4));
void g_dvec4_sub ANSI_PARM((VEC4, VEC4, VEC4));
void g_dvec4_add ANSI_PARM((VEC4, VEC4, VEC4));
void g_dvec4_mult ANSI_PARM((VEC4, double));
int g_dvec4_normalize ANSI_PARM((VEC4));
int g_dvec4_set_length ANSI_PARM((VEC4 , double));
void g_dvec4_assign ANSI_PARM((VEC4, VEC4));
void g_dvec4_list_assign ANSI_PARM ((VEC4 *, VEC4 *, int));
void g_dvec4_swap ANSI_PARM ((VEC4, VEC4));
double g_dvec4_angle ANSI_PARM((VEC4, VEC4, VEC4, VEC4));
double g_dvec4_ctg ANSI_PARM((VEC4, VEC4, VEC4, VEC4));
double g_dvec4_abs_sign ANSI_PARM((VEC4, VEC4));
double g_dvec4_sqrabs_sign ANSI_PARM((VEC4, VEC4));
double g_dvec4_skp_sign ANSI_PARM((VEC4, VEC4, VEC4));
double g_dvec4_dist_sign ANSI_PARM((VEC4, VEC4, VEC4));
double g_dvec4_sqrdist_sign ANSI_PARM((VEC4, VEC4, VEC4));
int g_dvec4_normalize_sign ANSI_PARM((VEC4, VEC4));
int g_dvec4_set_length_sign ANSI_PARM((VEC4 , double, VEC4));
void g_matrix44_mult_dvec4 ANSI_PARM((MATRIX44, VEC4));
void g_matrix44_mult_dvec4_dvec4 ANSI_PARM((MATRIX44, VEC4, VEC4));
void g_make_header ANSI_PARM ((INSTANCE *, char *, TEXTDEVICE *, long *));
void g_make_bottom ANSI_PARM ((TEXTDEVICE *, long));
int g_position ANSI_PARM ((char *, TEXTDEVICE *));
int g_vec3_get_arbitr_normal    ANSI_PARM ((VEC3 , VEC3));
int g_vec3_assign_triang_node   ANSI_PARM ((VEC3 , TRIANG2D *, int));
int g_get_triangle_with_edge    ANSI_PARM ((TRIANG2D *,int , int ));
int g_get_opposite_vertex       ANSI_PARM ((TRIANG2D *, int , int ));
int g_get_opp_vertex            ANSI_PARM ((INT3 , int , int ));
int g_get_opp_vertex_ind        ANSI_PARM ((INT3 , int , int ));
int g_triang2d_point_equal      ANSI_PARM ((TRIANG2D *, int , int ,double));
int g_triang2d_get_cls_pnt      ANSI_PARM ((TRIANG2D *, VEC3));
double g_triang2d_get_size      ANSI_PARM ((TRIANG2D *));
double g_triang2d_point_sqrdist ANSI_PARM ((TRIANG2D *, int , int ));
double g_triang2d_get_edge_length  ANSI_PARM ((TRIANG2D *, int , int ));
double g_triang2d_get_vertex_angle ANSI_PARM ((TRIANG2D *, int , int ));
int g_line_plane_intersect      ANSI_PARM ((TRIANG2D *, int, int, int, int, int, VEC3));
TRIANG2D *get_icosaeder         ANSI_PARM ((double *));
double g_triang1d_get_length      ANSI_PARM ((TRIANG1D *));
int g_triang1d_get_vec      ANSI_PARM ((TRIANG1D *, int, VEC3));
int g_line_reflect ANSI_PARM ((VEC3 [], MATRIX44));
int g_line_rotate ANSI_PARM ((VEC3 [], double, MATRIX44));
int g_plane_reflect ANSI_PARM ((VEC3 [], MATRIX44));
int g_point_reflect ANSI_PARM ((VEC3 , MATRIX44));
int g_translate ANSI_PARM ((VEC3 , MATRIX44));
int g_circle ANSI_PARM ((VEC3, VEC3, VEC3, VEC3, double, VEC3 *));
int g_helix ANSI_PARM ((VEC3, VEC3, VEC3, VEC3, double, VEC3 *));
void g_draw_circle(VEC3, double, VEC3, int, GRAPHICDEVICE *);
void g_draw_circle_orient(VEC3, double, VEC3, int, GRAPHICDEVICE *);
int g_draw_cyl(VEC3 *, VEC3 *, VEC3 *, VEC3 *, int, int, int, VEC3, GRAPHICDEVICE *);
int g_draw_arrow(VEC3, VEC3, double, int, int, int, VEC3, GRAPHICDEVICE *);
int g_solve4 ANSI_PARM ((double[][4], double[], double[]));
int g_solve3 ANSI_PARM ((VEC3 [], VEC3, VEC3));
int g_solve2 ANSI_PARM ((VEC2 [], VEC2, VEC2));
char *g_strdup ANSI_PARM ((const char *));
char *g_strfree ANSI_PARM ((char *));
char *g_strchange ANSI_PARM ((char *, const char *));
char *g_strcat ANSI_PARM ((char *, const char *));
char *g_strcut ANSI_PARM ((char *, size_t));
int g_strempty ANSI_PARM ((char *));
void g_check_bounds (double *, const double, const double);
#define NONE (-1)
#define ARRAY_STRUCT \
  INSTANCE_STRUCT; \
  int num; \
  int max; \
  char *data
typedef struct array {
  ARRAY_STRUCT;
} ARRAY;
bool_t g_xdr_instance(XDR *, INSTANCE **);
bool_t g_xdr_superclass(XDR *, INSTANCE *);
int g_xdr_add_type(char *, int, xdrproc_t);
#define g_add_array_type g_xdr_add_type
bool_t g_xdr_array(XDR *, void **, int, int, char *);
bool_t g_xdr_string(XDR *, char **);
bool_t g_xdr_version(XDR *, int *);
int g_xdr_get_version(void);
int g_xdr_get_revision(void);
int g_xdr_start(XDR *);
XDR *g_xdr_open_file(char *, OPEN_MODE);
int g_xdr_end(void);
int g_xdr_close_file(XDR *);
bool_t g_xdr_int2(XDR *, INT2);
bool_t g_xdr_int3(XDR *, INT3);
bool_t g_xdr_int4(XDR *, INT4);
bool_t g_xdr_vec2(XDR *, VEC3);
bool_t g_xdr_vec3(XDR *, VEC3);
bool_t g_xdr_vec4(XDR *, VEC4);
bool_t g_xdr_fvec3(XDR *, FVEC3);
bool_t g_xdr_matrix44(XDR *, MATRIX44);
#define REFINE 1
#define COARSE 2
#define G_REFINE 0x1
#define G_COARSE 0x2
#define G_MOVE   0x4
#define G_DELETE 0x8
#define G_NOT_REFINE 0x10
#define G_NOT_COARSE 0x20
#define G_MOVED   0x40
#define G_DELETED 0x80
typedef struct refine2d {
  unsigned char refine;
  unsigned char edge;
  unsigned char friend;
  unsigned char level;
  unsigned char vtx_level[3];
} REFINE2D;
#define ADAPT2D_STRUCT \
  FE2D_STRUCT; \
  REFINE2D *refine; \
  int number_of_bnd_functions; \
  int max_number_of_bnd_functions; \
  BND_FUNCTION *bnd_function; \
  void (*interpol)(), (*restrict)(); \
  void (*element_compress)(), (*point_compress)()
typedef struct adapt2d {
  ADAPT2D_STRUCT;
} ADAPT2D;
bool_t g_xdr_refine2d(XDR *, REFINE2D *);
typedef struct refine3d {
  char refine;
  char level[4];
  char edge[4];
} REFINE3D;
#define ADAPT3D_STRUCT \
  FE3D_STRUCT; \
  REFINE3D *refine; \
  int number_of_boundaries, max_number_of_boundaries; \
  BND_FUNCTION *bnd_function; \
  void (*interpol)(),(*restrict)(); \
  void (*element_compress)(),(*point_compress)()
typedef struct adapt3d {
  ADAPT3D_STRUCT;
} ADAPT3D;
bool_t g_xdr_refine3d(XDR *, REFINE3D *);
#define TRIANG0D_STRUCT \
  INSTANCE_STRUCT; \
  double x,y,z
typedef struct triang0d {
  TRIANG0D_STRUCT;
} TRIANG0D;
#define TRACE_STRUCT \
  TRIANG1D_STRUCT; \
  double *t
typedef struct trace {
  TRACE_STRUCT;
} TRACE;
#define STREAKLINE_STRUCT \
  TRACE_STRUCT; \
  int number_of_steps; \
  int max_number_of_steps
typedef struct streakline {
  STREAKLINE_STRUCT;
} STREAKLINE;
#define CLOUD_STRUCT \
  INSTANCE_STRUCT; \
  int number_of_points; \
  VEC3 *x
typedef struct cloud {
  CLOUD_STRUCT;
} CLOUD;
typedef struct mpoint {
  double *x,*y,*z;
  int it,old_it;
  double bary[4];
} MPOINT;
typedef struct moving_info {
  int maxn,n,i;
  INSTANCE **pre,**post;
  int *new_discr;
  double *time;
  int dim,space;
  int nsteps;
  double rulstime,ruletime;
  double runge_delta;
  double rulstreaktime;
  double stime,etime;
  double ctime,ntime;
  double step_delta;
  int maxnop,nop;
  MPOINT *p;
  int (*global_search)(),(*local_search)();
  double smu,cmu;
  void (*f)();
  double ref_eps ;
  double vertical_texture ;
  double horizontal_texture ;
} MOVING_INFO;
typedef struct optipar {
  double (*F1)(),(*DF1)(),(*F2)(),(*DF2)();
  double (*concentration)();
  int itmax;
  double eps;
} OPTIPAR;
typedef struct element2d_description ELEMENT2D_DESCRIPTION;
typedef struct element2d ELEMENT2D;
typedef struct f_data2d F_DATA2D;
typedef struct f_el_info2d F_EL_INFO2D;
struct element2d_description
{
  int number_of_vertices;
  int dimension_of_coord;
  float     **coord;
  int parametric_degree;
  int (*world_to_coord)(ELEMENT2D *,float *, float *);
  void (*coord_to_world)(ELEMENT2D *,float *, float *);
  int (*check_inside)(ELEMENT2D *,float *);
  ELEMENT2D *(*neighbour)(ELEMENT2D *, int, int, float *, float *);
  int (*boundary)(ELEMENT2D *, int);
};
struct element2d
{
  struct mesh2d         *mesh;
  float                 **vertex;
  int                   *vindex;
  int eindex;
  ELEMENT2D_DESCRIPTION *descr;
  int size_of_user_data;
  void                  *user_data;
};
#define MESH2D_STRUCT \
  INSTANCE_STRUCT; \
  ELEMENT2D     *(*first_element)(struct mesh2d *); \
  ELEMENT2D     *(*next_element)(ELEMENT2D *); \
  ELEMENT2D     *(*copy_element)(ELEMENT2D *); \
  void (*free_element)(ELEMENT2D *); \
  int dimension_of_world; \
  int max_dimension_of_coord; \
  int max_eindex; \
  int max_vindex; \
  int max_number_of_vertices; \
  F_DATA2D      *f_data; \
  int size_of_user_data; \
  void          *user_data
typedef struct mesh2d {
  MESH2D_STRUCT;
} MESH2D;
struct f_data2d
{
  char     *name;
  int dimension_of_value;
  int continuous_data;
  void (*f)(ELEMENT2D *, int,
            float[], float[]);
  void (*f_el_info)(ELEMENT2D *,
                    F_EL_INFO2D *);
  int size_of_user_data;
  void     *user_data;
  F_DATA2D *last, *next;
};
struct f_el_info2d
{
  int polynomial_degree;
};
#ifndef INSIDE
#define INSIDE -1
#define EXACT_NEIGHBOUR -1
#define FIRST_NEIGHBOUR 0
#define NEXT_NEIGHBOUR  1
#endif
int g_dbl_normal_to_plane(VEC3 , FVEC3 , FVEC3 , FVEC3 );
#define g_fvec3_to_vec3(f, d) \
  (d)[0]=(double)(f)[0]; (d)[1]=(double)(f)[1]; (d)[2]=(double)(f)[2];
typedef struct element3d_description ELEMENT3D_DESCRIPTION;
typedef struct helement3d_description HELEMENT3D_DESCRIPTION;
typedef struct element3d ELEMENT3D;
typedef struct helement3d HELEMENT3D;
typedef struct f_data3d F_DATA3D;
typedef struct f_hdata3d F_HDATA3D;
typedef struct f_el_info3d F_EL_INFO3D;
typedef struct f_hel_info3d F_HEL_INFO3D;
typedef struct vinherit VINHERIT;
struct element3d_description
{
  int number_of_vertices;
  int number_of_polygons;
  int       *polygon_length;
  int       **polygon_vertex;
  int       **polygon_neighbour;
  int dimension_of_coord;
  double     **coord;
  int parametric_degree;
  int (*world_to_coord)(ELEMENT3D *,float *, double *);
  void (*coord_to_world)(ELEMENT3D *,double *, float *);
  int (*check_inside)(ELEMENT3D *,double *);
  ELEMENT3D *(*neighbour)(ELEMENT3D *, int, int, double *, float *);
  int (*boundary)(ELEMENT3D *,  int);
};
struct element3d
{
  struct mesh3d          *mesh;
  float                  **vertex;
  int                    *vindex;
  int eindex;
  ELEMENT3D_DESCRIPTION  *descr;
  int size_of_user_data;
  void                   *user_data;
};
struct f_data3d
{
  char     *name;
  int dimension_of_value;
  int continuous_data;
  void (*f)(ELEMENT3D *, int,
            double[], double[]);
  void (*f_el_info)(ELEMENT3D *,
                    F_EL_INFO3D *);
  int size_of_user_data;
  void     *user_data;
  F_DATA3D *last, *next;
};
struct f_el_info3d
{
  int polynomial_degree;
};
#define MESH3D_DATA F_DATA3D
#define MESH3D_STRUCT \
  INSTANCE_STRUCT; \
  ELEMENT3D     *(*first_element)(); \
  ELEMENT3D     *(*next_element)(); \
  ELEMENT3D     *(*copy_element)(); \
  void (*free_element)(); \
  int max_dimension_of_coord; \
  int max_eindex; \
  int max_vindex; \
  int max_number_of_vertices; \
  MESH3D_DATA   *f_data; \
  int size_of_user_data; \
  void          *user_data
typedef struct mesh3d {
  MESH3D_STRUCT;
} MESH3D;
#undef MESH3D_DATA
struct vinherit {
  int np;
  int *pindex;
  double *pweight;
};
struct f_hdata3d
{
  char     *name;
  int dimension_of_value;
  int continuous_data;
  void (*f)(HELEMENT3D *, int,
            double[], double[]);
  void (*f_el_info)(HELEMENT3D *,
                    F_HEL_INFO3D *);
  int size_of_user_data;
  void     *user_data;
  F_DATA3D *last, *next;
  void (*get_bounds)(HELEMENT3D *,
                     double* ,double*);
  double (*get_estimate)(HELEMENT3D*,
                         double*);
  double est_bound;
};
struct f_hel_info3d
{
  int polynomial_degree;
};
struct helement3d_description
{
  int number_of_vertices;
  int number_of_polygons;
  int       *polygon_length;
  int       **polygon_vertex;
  int       **polygon_neighbour;
  int dimension_of_coord;
  double     **coord;
  int parametric_degree;
  int (*world_to_coord)(HELEMENT3D *,float *, double *);
  void (*coord_to_world)(HELEMENT3D *,double *, float *);
  int (*check_inside)(HELEMENT3D *,double *);
  HELEMENT3D *(*neighbour)(HELEMENT3D *, int, int, double *, float *);
  int (*boundary)(HELEMENT3D *,  int);
};
struct helement3d
{
  struct hmesh3d         *mesh;
  float                  **vertex;
  int                    *vindex;
  int eindex;
  HELEMENT3D_DESCRIPTION *descr;
  int size_of_user_data;
  void                   *user_data;
  HELEMENT3D             *parent;
  VINHERIT               *vinh;
};
#define MESH3D_DATA F_HDATA3D
#define HMESH3D_STRUCT \
  MESH3D_STRUCT; \
  HELEMENT3D    *(*first_child)(HELEMENT3D *); \
  HELEMENT3D    *(*next_child)(HELEMENT3D *); \
  HELEMENT3D    *(*first_macro)(); \
  HELEMENT3D    *(*next_macro)(HELEMENT3D *); \
  int max_level ; \
  int level_of_interest
typedef struct hmesh3d {
  HMESH3D_STRUCT;
} HMESH3D;
#undef MESH3D_DATA
typedef struct clipm3d_par {
  ELEMENT3D *e;
  int n;
  VEC3 *v;
  double *coord;
  double *fval;
  double (*f)();
  char *var;
  char *scal;
  int flag ;
  int dimension_of_value;
} CLIPM3D_PAR;
#ifndef INSIDE
#define INSIDE -1
#define EXACT_NEIGHBOUR -1
#define FIRST_NEIGHBOUR 0
#define NEXT_NEIGHBOUR  1
#endif
extern int g_probe_fineness;
#define GEOM2D_ACTIVE     0x1
#define GEOM2D_DISPLAY    0x2
#define GEOM2D_USERNORM   0x4
#define GEOM2D_DISP_BND   0x8
#define GEOM2D_DISP_HOR   0x10
#define GEOM2D_DISP_VER   0x20
#define GEOM2D_STRUCT \
  ADAPT2D_STRUCT; \
  int flag; \
  VEC3 *normal; \
  VEC3 *color_vertex; \
  VEC3 *color_patch; \
  LIST_OF_INST *bnd_curve; \
  CHAIN *hor, *ver
typedef struct geom2d {
  GEOM2D_STRUCT;
} GEOM2D;
#define CONS_STRUCT \
  INSTANCE_STRUCT; \
  TRIANG2D *geom; \
  int number_of_points; \
  int max_number_of_points; \
  int *vertex; \
  int flag
typedef struct cons {
  CONS_STRUCT;
} CONS;
#define FIX_CIRCLE    00004
#define FIX_CURVE     00010
#define FIX_STRAIGHT  00020
#define FIX_PLANE     00040
#define FIX_THREAD    00100
#define EDGE_STRAIGHT 00200
#define EDGE_CIRCLE   00400
#define EDGE_HELIX    01000
#define EDGE_SPLINE   02000
#define BND_CURVE_STRUCT \
  TRIANG1D_STRUCT; \
  int flag; \
  TRIANG2D *geom; \
  int *vertex; \
  int *neighbour; \
  struct { \
    void (*f)(struct bnd_curve *); \
    int size; \
    double *data; \
  } restrict
typedef struct bnd_curve {
  BND_CURVE_STRUCT;
} BND_CURVE;
#define LIST_OF_STRING_STRUCT \
  INSTANCE_STRUCT; \
  char *string; \
  int string_number; \
  struct list_of_string *next_string
typedef struct list_of_string {
  LIST_OF_STRING_STRUCT;
} LIST_OF_STRING;
#define LIST_OF_DESCRIPTOR_STRUCT \
  INSTANCE_STRUCT; \
  char *descriptor; \
  int descriptor_number; \
  int list_length; \
  LIST_OF_STRING *list_start, *list
typedef struct list_of_descriptor {
  LIST_OF_DESCRIPTOR_STRUCT;
} LIST_OF_DESCRIPTOR;
#define TIME_PARM_STRUCT \
  INSTANCE_STRUCT; \
  BOUNDED_RULER **int_rul; \
  FUNCTION_RULER **dbl_rul; \
  int num_int_parm; \
  int num_double_parm; \
  int *int_parm; \
  double *double_parm
typedef struct time_parm {
  TIME_PARM_STRUCT;
} TIME_PARM;
int g_time_parm_get_var_list(TIME_PARM *, char *);
void g_time_parm_get_var_values(TIME_PARM *, double, double *);
typedef enum {
  tskGeometry  = 0x00000001,
  tskGauss     = 0x00000002,
  tskConjugate = 0x00000004,
  tskDomain    = 0x00000008,
  tskNewton    = 0x00000010,
  tskExplicit  = 0x00000020,
  tskAny       = 0x000000ff
} TASK_TYPE;
#define TIME_OBJECT_STRUCT \
  INSTANCE_STRUCT; \
  TASK_TYPE task, task_mask; \
  int frame, num_of_frames; \
  int *status; \
  double *time, *time_step; \
  TIME_PARM *time_parm; \
  CHAIN *geometry, *inter; \
  SCENE *tree
typedef struct time_object {
  TIME_OBJECT_STRUCT;
} TIME_OBJECT;
TASK_TYPE g_read_task(FILE *, int);
char *g_task_name(TASK_TYPE);
#define REFLECTION_CONTROL_STRUCT \
  INSTANCE_STRUCT; \
  int type; \
  int obj; \
  int last_obj_index; \
  int object_index, bdy_index, bdy2_index
typedef struct reflection_control {
  REFLECTION_CONTROL_STRUCT;
} REFLECTION_CONTROL;
#define SURFACE_STRUCT \
  TIME_OBJECT_STRUCT; \
  CHAIN *model, *domain, *gauss, *conjugate, *conj_inter; \
  double *asso; \
  REFLECTION_CONTROL *reflect
typedef struct surface {
  SURFACE_STRUCT;
} SURFACE;
typedef union {
  struct {
    char c;
    int i;
  } op;
  double d;
} t_fun;
char  *fn_errstr ANSI_PARM((void));
char  *fn_errpos ANSI_PARM((void));
t_fun *fn_func   ANSI_PARM((char *, char *));
double fn_eval   ANSI_PARM((t_fun *, double *));
char  *fn_string ANSI_PARM((t_fun *, char *));
t_fun **t_fun_alloc   ANSI_PARM((int));
t_fun **t_fun_free    ANSI_PARM((t_fun **, int));
int get_t_fun         ANSI_PARM((t_fun **, char *, FILE *));
double tfun_min       ANSI_PARM((t_fun *, double, double));
double tfun_max       ANSI_PARM((t_fun *, double, double));
SPLINE1D *tfun2spline ANSI_PARM((t_fun *, double, double, int, char *, SPLINE1D *));
int funerr            ANSI_PARM((char *));
typedef struct domain_parm_rul {
  FUNCTION_RULER **u_min;
  FUNCTION_RULER **v_min;
  FUNCTION_RULER **u_max;
  FUNCTION_RULER **v_max;
  FUNCTION_RULER **u_dstrbn;
  FUNCTION_RULER **v_dstrbn;
} DOMAIN_PARM_RUL;
#define DOMAIN_PARM_STRUCT \
  INSTANCE_STRUCT; \
  DOMAIN_PARM_RUL rul; \
  int num_parts; \
  int max_parts; \
  int *flag; \
  int *num_u_lines; \
  int *num_v_lines; \
  int *num_u_sub; \
  int *num_v_sub; \
  SPLINE1D **u_dstrbn; \
  SPLINE1D **v_dstrbn; \
  t_fun **u_min; \
  t_fun **v_min; \
  t_fun **u_max; \
  t_fun **v_max; \
  double **u_min_d; \
  double **v_min_d; \
  double **u_max_d; \
  double **v_max_d
typedef struct domain_parm {
  DOMAIN_PARM_STRUCT;
} DOMAIN_PARM;
#define EXPLICIT_STRUCT \
  SURFACE_STRUCT; \
  DOMAIN_PARM *domain_parm; \
  COMPLEX (*compute_arg)(); \
  struct { \
    t_fun *x, *y, *z; \
  } coord_f
typedef struct explicit {
  EXPLICIT_STRUCT;
} EXPLICIT;
typedef struct cnode {
  VEC3 xyz;
  VEC3 normal;
  COMPLEX z, w, gauss;
  CVEC3 cinti, cintj;
  int np;
  int i, j;
  int ipoints;
} CNODE;
#define AMANDUS_STRUCT \
  EXPLICIT_STRUCT; \
  CHAIN *newton; \
  struct { \
    int type; \
    COMPLEX (*x)(COMPLEX), (*y)(COMPLEX), (*z)(COMPLEX); \
    COMPLEX (*g)(), (*dh)(), (*dg)(); \
    COMPLEX (*poly)(COMPLEX, COMPLEX); \
    COMPLEX (*dwpoly)(COMPLEX, COMPLEX); \
    COMPLEX (*dzpoly)(COMPLEX, COMPLEX); \
    COMPLEX (*newton_init)(COMPLEX); \
    void (*correct)(); \
  } weier_f; \
  CNODE *cnode, *icnode
typedef struct amandus {
  AMANDUS_STRUCT;
} AMANDUS;
#define KARL            37
#define HERMANN         38
#define ORBIT_STRUCT \
  TRIANG1D_STRUCT; \
  double time; \
  double length; \
  int ind; \
  double ind_pos
typedef struct orbit {
  ORBIT_STRUCT;
} ORBIT;
#define UFO_STRUCT \
  INSTANCE_STRUCT; \
  VEC3 mid
typedef struct ufo {
  UFO_STRUCT;
} UFO;
#define G_DSO_HANDLER_INIT g_scan_libs ()
#define JVC_STRUCT \
  INSTANCE_STRUCT; \
  int number
typedef struct jvc {
  JVC_STRUCT;
} JVC;
#define LVR_STRUCT \
  INSTANCE_STRUCT; \
  int number
typedef struct lvr {
  LVR_STRUCT;
} LVR;
#define GRAPHICGT_STRUCT \
  GRAPHICDEVICE_STRUCT
typedef struct graphicgt {
  GRAPHICGT_STRUCT;
} GRAPHICGT;
#define G_GRAPHICGT_EXISTS
extern char *g_output_envar;
void g_default_graphicdevice( CLASS **, char **);
#define GRAPHICX11_STRUCT \
  GRAPHICDEVICE_STRUCT
typedef struct graphicx11 {
  GRAPHICX11_STRUCT;
} GRAPHICX11;
void g_default_controldevice( CLASS **, char ** );
#define CONTROLX11_STRUCT \
  CONTROLDEVICE_STRUCT
typedef struct controlx11 {
  CONTROLX11_STRUCT;
} CONTROLX11;
#define ZERO_MACHINE 1.e-10
#define g_sqrt(x)   sqrt(x)
#define g_fabs(x)   fabs(x)
#define g_pow(x,y)  pow(x,y)
#define g_cos(x)    cos(x)
#define g_sin(x)    sin(x)
#define g_tan(x)    tan(x)
#define g_acos(x)   acos(x)
#define g_asin(x)   asin(x)
#define g_atan(x)   atan(x)
#define g_atan2(x)  atan2(x)
#define g_floor(x)  floor(x)
#define g_fmod(x,y) fmod(x,y)
#define g_ceil(x)   ceil(x)
#define g_exp(x)    exp(x)
#define g_log(x)    log(x)
#define g_log10(x)  log10(x)
#define g_sqrtf(x)   sqrtf(x)
#define g_fabsf(x)   fabsf(x)
#define g_powf(x,y)  powf(x,y)
#define g_cosf(x)    cosf(x)
#define g_sinf(x)    sinf(x)
#define g_tanf(x)    tanf(x)
#define g_acosf(x)   acosf(x)
#define g_asinf(x)   asinf(x)
#define g_atanf(x)   atanf(x)
#define g_atan2f(x)  atan2f(x)
#define g_floorf(x)  floorf(x)
#define g_fmodf(x,y) fmodf(x,y)
#define g_ceilf(x)   ceilf(x)
#define g_expf(x)    expf(x)
#define g_logf(x)    logf(x)
#define g_log10f(x)  log10f(x)
extern CLASS *Adapt2d;
extern CLASS *Adapt3d;
extern CLASS *Amandus;
extern CLASS *Array;
extern CLASS *Bar1d;
extern CLASS *Bezier1d;
extern CLASS *Bigbutton;
extern CLASS *Bnd_Curve;
extern CLASS *Bounded_Ruler;
extern CLASS *Button;
extern CLASS *Camera;
extern CLASS *Chain;
extern CLASS *CheckField;
extern CLASS *Checkbox;
extern CLASS *Cloud;
extern CLASS *Color_Sel;
extern CLASS *Colorbar;
extern CLASS *Cons;
extern CLASS *ControlDevice;
extern CLASS *ControlX11;
extern CLASS *Cycle_Button;
extern CLASS *Device;
extern CLASS *Domain_Parm;
extern CLASS *Explicit;
extern CLASS *Fe2d;
extern CLASS *Fe3d;
extern CLASS *FileList;
extern CLASS *Function_Bar1d;
extern CLASS *Function_Ruler;
extern CLASS *Function_Slider;
extern CLASS *G_List;
extern CLASS *Gd_Sphere;
extern CLASS *Geom2d;
extern CLASS *GraphicDevice;
extern CLASS *GraphicGt;
extern CLASS *GraphicPs;
extern CLASS *GraphicX11;
extern CLASS *Group;
extern CLASS *HMesh3d;
extern CLASS *InOutFile;
extern CLASS *Input;
extern CLASS *InputFile;
extern CLASS *Interactive;
extern CLASS *Item;
extern CLASS *Jvc;
extern CLASS *Layer;
extern CLASS *Line;
extern CLASS *List_Of_Descriptor;
extern CLASS *List_Of_Double;
extern CLASS *List_Of_Inst;
extern CLASS *List_Of_Int;
extern CLASS *List_Of_String;
extern CLASS *List_Of_Vec3;
extern CLASS *Lvr;
extern CLASS *Manager;
extern CLASS *Matrix_Edit;
extern CLASS *Mesh2d;
extern CLASS *Mesh3d;
extern CLASS *MethLayer;
extern CLASS *MethList;
extern CLASS *Orbit;
extern CLASS *OutputFile;
extern CLASS *Plane;
extern CLASS *Project;
extern CLASS *Prompt;
extern CLASS *PureList;
extern CLASS *PureRuler;
extern CLASS *PureScrollbar;
extern CLASS *PureSlider;
extern CLASS *PureSplineEditor;
extern CLASS *Radio;
extern CLASS *Rectangle;
extern CLASS *Reflection_Control;
extern CLASS *ResourceDev;
extern CLASS *Root;
extern CLASS *Ruler;
extern CLASS *Scene;
extern CLASS *Scrollbar;
extern CLASS *Selector;
extern CLASS *Slider;
extern CLASS *Sphere;
extern CLASS *Spline1d;
extern CLASS *Spline_Editor;
extern CLASS *StaticText;
extern CLASS *Step_Button;
extern CLASS *StreakLine;
extern CLASS *Suprop;
extern CLASS *Surface;
extern CLASS *TextDevice;
extern CLASS *Textmessage;
extern CLASS *TimeScene;
extern CLASS *TimeStep;
extern CLASS *Time_Object;
extern CLASS *Time_Parm;
extern CLASS *Trace;
extern CLASS *Triang0d;
extern CLASS *Triang1d;
extern CLASS *Triang2d;
extern CLASS *Triang3d;
extern CLASS *Ufo;
extern CLASS *XList;
extern CLASS *XList_Button;
