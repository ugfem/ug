// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#define FLEX_SCANNER
#define YY_FLEX_MAJOR_VERSION 2
#define YY_FLEX_MINOR_VERSION 5

#include <stdio.h>


/* cfront 1.2 defines "c_plusplus" instead of "__cplusplus" */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif


#ifdef __cplusplus

#include <stdlib.h>
#include <unistd.h>

/* Use prototypes in function declarations. */
#define YY_USE_PROTOS

/* The "const" storage-class-modifier is valid. */
#define YY_USE_CONST

#else   /* ! __cplusplus */

#if __STDC__

#define YY_USE_PROTOS
#define YY_USE_CONST

#endif  /* __STDC__ */
#endif  /* ! __cplusplus */

#ifdef __TURBOC__
 #pragma warn -rch
 #pragma warn -use
#include <io.h>
#include <stdlib.h>
#define YY_USE_CONST
#define YY_USE_PROTOS
#endif

#ifdef YY_USE_CONST
#define ngconst const
#else
#define ngconst
#endif


#ifdef YY_USE_PROTOS
#define YY_PROTO(proto) proto
#else
#define YY_PROTO(proto) ()
#endif

/* Returned upon end-of-file. */
#define YY_NULL 0

/* Promotes a possibly negative, possibly signed char to an unsigned
 * integer for use as an array index.  If the signed char is negative,
 * we want to instead treat it as an 8-bit unsigned char, hence the
 * double cast.
 */
#define YY_SC_TO_UI(c) ((unsigned int) (unsigned char) c)

/* Enter a start condition.  This macro really ought to take a parameter,
 * but we do it the disgusting crufty way forced on us by the ()-less
 * definition of BEGIN.
 */
#define BEGIN ng_start = 1 + 2 *

/* Translate the current start state into a value that can be later handed
 * to BEGIN to return to the state.  The YYSTATE alias is for lex
 * compatibility.
 */
#define YY_START ((ng_start - 1) / 2)
#define YYSTATE YY_START

/* Action number for EOF rule of a given start state. */
#define YY_STATE_EOF(state) (YY_END_OF_BUFFER + state + 1)

/* Special action meaning "start processing a new file". */
#define YY_NEW_FILE ngrestart( ngin )

#define YY_END_OF_BUFFER_CHAR 0

/* Size of default input buffer. */
#define YY_BUF_SIZE 16384

typedef struct ng_buffer_state *YY_BUFFER_STATE;

extern int ngleng;
extern FILE *ngin, *ngout;

#define EOB_ACT_CONTINUE_SCAN 0
#define EOB_ACT_END_OF_FILE 1
#define EOB_ACT_LAST_MATCH 2

/* The funky do-while in the following #define is used to turn the definition
 * int a single C statement (which needs a semi-colon terminator).  This
 * avoids problems with code like:
 *
 *      if ( condition_holds )
 *		ngless( 5 );
 *	else
 *		do_something_else();
 *
 * Prior to using the do-while the compiler would get upset at the
 * "else" because it interpreted the "if" statement as being all
 * done when it reached the ';' after the ngless() call.
 */

/* Return all but the first 'n' matched characters back to the input stream. */

#define ngless(n) \
  do \
  { \
    /* Undo effects of setting up ngtext. */ \
    *ng_cp = ng_hold_char; \
    YY_RESTORE_YY_MORE_OFFSET \
      ng_c_buf_p = ng_cp = ng_bp + n - YY_MORE_ADJ; \
    YY_DO_BEFORE_ACTION;             /* set up ngtext again */ \
  } \
  while ( 0 )

#define unput(c) ngunput( c, ngtext_ptr )

/* The following is because we cannot portably get our hands on size_t
 * (without autoconf's help, which isn't available because we want
 * flex-generated scanners to compile on their own).
 */
typedef unsigned int ng_size_t;


struct ng_buffer_state
{
  FILE *ng_input_file;

  char *ng_ch_buf;                      /* input buffer */
  char *ng_buf_pos;                     /* current position in input buffer */

  /* Size of input buffer in bytes, not including room for EOB
   * characters.
   */
  ng_size_t ng_buf_size;

  /* Number of characters read into ng_ch_buf, not including EOB
   * characters.
   */
  int ng_n_chars;

  /* Whether we "own" the buffer - i.e., we know we created it,
   * and can realloc() it to grow it, and should free() it to
   * delete it.
   */
  int ng_is_our_buffer;

  /* Whether this is an "interactive" input source; if so, and
   * if we're using stdio for input, then we want to use getc()
   * instead of fread(), to make sure we stop fetching input after
   * each newline.
   */
  int ng_is_interactive;

  /* Whether we're considered to be at the beginning of a line.
   * If so, '^' rules will be active on the next match, otherwise
   * not.
   */
  int ng_at_bol;

  /* Whether to try to fill the input buffer when we reach the
   * end of it.
   */
  int ng_fill_buffer;

  int ng_buffer_status;
#define YY_BUFFER_NEW 0
#define YY_BUFFER_NORMAL 1
  /* When an EOF's been seen but there's still some text to process
   * then we mark the buffer as YY_EOF_PENDING, to indicate that we
   * shouldn't try reading from the input source any more.  We might
   * still have a bunch of tokens to match, though, because of
   * possible backing-up.
   *
   * When we actually see the EOF, we change the status to "new"
   * (via ngrestart()), so that the user can continue scanning by
   * just pointing ngin at a new input file.
   */
#define YY_BUFFER_EOF_PENDING 2
};

static YY_BUFFER_STATE ng_current_buffer = 0;

/* We provide macros for accessing buffer states in case in the
 * future we want to put the buffer states in a more general
 * "scanner state".
 */
#define YY_CURRENT_BUFFER ng_current_buffer


/* ng_hold_char holds the character lost when ngtext is formed. */
static char ng_hold_char;

static int ng_n_chars;          /* number of characters read into ng_ch_buf */


int ngleng;

/* Points to current character in buffer. */
static char *ng_c_buf_p = (char *) 0;
static int ng_init = 1;         /* whether we need to initialize */
static int ng_start = 0;        /* start state number */

/* Flag which is used to allow ngwrap()'s to do buffer switches
 * instead of setting up a fresh ngin.  A bit of a hack ...
 */
static int ng_did_buffer_switch_on_eof;

void ngrestart YY_PROTO(( FILE *input_file ));

void ng_switch_to_buffer YY_PROTO(( YY_BUFFER_STATE new_buffer ));
void ng_load_buffer_state YY_PROTO(( void ));
YY_BUFFER_STATE ng_create_buffer YY_PROTO(( FILE *file, int size ));
void ng_delete_buffer YY_PROTO(( YY_BUFFER_STATE b ));
void ng_init_buffer YY_PROTO(( YY_BUFFER_STATE b, FILE *file ));
void ng_flush_buffer YY_PROTO(( YY_BUFFER_STATE b ));
#define YY_FLUSH_BUFFER ng_flush_buffer( ng_current_buffer )

YY_BUFFER_STATE ng_scan_buffer YY_PROTO(( char *base, ng_size_t size ));
YY_BUFFER_STATE ng_scan_string YY_PROTO(( ngconst char *ng_str ));
YY_BUFFER_STATE ng_scan_bytes YY_PROTO(( ngconst char *bytes, int len ));

static void *ng_flex_alloc YY_PROTO(( ng_size_t ));
static void *ng_flex_realloc YY_PROTO(( void *, ng_size_t ));
static void ng_flex_free YY_PROTO(( void * ));

#define ng_new_buffer ng_create_buffer

#define ng_set_interactive(is_interactive) \
  { \
    if ( ! ng_current_buffer ) \
      ng_current_buffer = ng_create_buffer( ngin, YY_BUF_SIZE ); \
    ng_current_buffer->ng_is_interactive = is_interactive; \
  }

#define ng_set_bol(at_bol) \
  { \
    if ( ! ng_current_buffer ) \
      ng_current_buffer = ng_create_buffer( ngin, YY_BUF_SIZE ); \
    ng_current_buffer->ng_at_bol = at_bol; \
  }

#define YY_AT_BOL() (ng_current_buffer->ng_at_bol)

typedef unsigned char YY_CHAR;
FILE *ngin = (FILE *) 0, *ngout = (FILE *) 0;
typedef int ng_state_type;
extern char *ngtext;
#define ngtext_ptr ngtext

static ng_state_type ng_get_previous_state YY_PROTO(( void ));
static ng_state_type ng_try_NUL_trans YY_PROTO(( ng_state_type current_state ));
static int ng_get_next_buffer YY_PROTO(( void ));
static void ng_fatal_error YY_PROTO(( ngconst char msg[] ));

/* Done after the current pattern has been matched and before the
 * corresponding action - sets up ngtext.
 */
#define YY_DO_BEFORE_ACTION \
  ngtext_ptr = ng_bp; \
  ngleng = (int) (ng_cp - ng_bp); \
  ng_hold_char = *ng_cp; \
  *ng_cp = '\0'; \
  ng_c_buf_p = ng_cp;

#define YY_NUM_RULES 13
#define YY_END_OF_BUFFER 14
static ngconst short int ng_accept[27] =
{   0,
    0,    0,   14,   12,    1,    2,    3,   12,   12,    4,
    11,    7,    9,    6,    8,   10,    1,    3,    0,    5,
    5,    4,    0,    0,    5,    0} ;

static ngconst int ng_ec[256] =
{   0,
    1,    1,    1,    1,    1,    1,    1,    1,    2,    3,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    2,    1,    1,    4,    1,    1,    1,    1,    1,
    1,    1,    5,    1,    6,    7,    1,    8,    8,    8,
    8,    8,    8,    8,    8,    8,    8,    1,    9,    1,
    1,    1,    1,    1,    1,   10,    1,    1,   11,    1,
    1,    1,   12,    1,    1,   13,    1,    1,    1,    1,
    1,    1,   14,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,

    15,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,

    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1} ;

static ngconst int ng_meta[16] =
{   0,
    1,    1,    2,    1,    3,    3,    1,    4,    1,    1,
    5,    1,    1,    1,    5} ;

static ngconst short int ng_base[30] =
{   0,
    0,    0,   40,   41,   37,   41,    0,    9,   30,   11,
    41,   41,   41,   41,   41,   41,   35,    0,   28,   13,
    27,   15,   26,   25,   17,   41,   23,   25,   28} ;

static ngconst short int ng_def[30] =
{   0,
    26,    1,   26,   26,   26,   26,   27,   26,   26,   26,
    26,   26,   26,   26,   26,   26,   26,   27,   26,   26,
    28,   26,   29,   26,   26,    0,   26,   26,   26} ;

static ngconst short int ng_nxt[57] =
{   0,
    4,    5,    6,    7,    4,    8,    9,   10,   11,   12,
    13,   14,   15,   16,    4,   19,   20,   19,   22,   19,
    20,   19,   22,   18,   25,   18,   18,   18,   23,   23,
    24,   24,   25,   25,   21,   21,   17,   21,   17,   26,
    3,   26,   26,   26,   26,   26,   26,   26,   26,   26,
    26,   26,   26,   26,   26,   26} ;

static ngconst short int ng_chk[57] =
{   0,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    8,    8,   10,   10,   20,
    20,   22,   22,   27,   25,   27,   27,   27,   28,   28,
    29,   29,   24,   23,   21,   19,   17,    9,    5,    3,
    26,   26,   26,   26,   26,   26,   26,   26,   26,   26,
    26,   26,   26,   26,   26,   26} ;

static ng_state_type ng_last_accepting_state;
static char *ng_last_accepting_cpos;

/* The intent behind this definition is that it'll catch
 * any uses of REJECT which flex missed.
 */
#define REJECT reject_used_but_not_detected
#define ngmore() ngmore_used_but_not_detected
#define YY_MORE_ADJ 0
#define YY_RESTORE_YY_MORE_OFFSET
char *ngtext;
#line 1 "ngin.l"
#define INITIAL 0
#line 2 "ngin.l"
/****************************************************************************/
/*                                                                          */
/* File:      grid.l                                                        */
/*                                                                          */
/* Purpose:   lexer for gridfiles                                           */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen                             */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: ug@ica3.uni-stuttgart.de                            */
/*                                                                          */
/* History:   19.2.98 begin,                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include "ng2d.h"
#include "ngin-yacc.h"

static int noline=1;

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);


#line 410 "lex.ng.c"

/* Macros after this point can all be overridden by user definitions in
 * section 1.
 */

#ifndef YY_SKIP_YYWRAP
#ifdef __cplusplus
extern "C" int ngwrap YY_PROTO(( void ));
#else
extern int ngwrap YY_PROTO(( void ));
#endif
#endif

#ifndef YY_NO_UNPUT
static void ngunput YY_PROTO(( int c, char *buf_ptr ));
#endif

#ifndef ngtext_ptr
static void ng_flex_strncpy YY_PROTO(( char *, ngconst char *, int ));
#endif

#ifdef YY_NEED_STRLEN
static int ng_flex_strlen YY_PROTO(( ngconst char * ));
#endif

#ifndef YY_NO_INPUT
#ifdef __cplusplus
static int nginput YY_PROTO(( void ));
#else
static int input YY_PROTO(( void ));
#endif
#endif

#if YY_STACK_USED
static int ng_start_stack_ptr = 0;
static int ng_start_stack_depth = 0;
static int *ng_start_stack = 0;
#ifndef YY_NO_PUSH_STATE
static void ng_push_state YY_PROTO(( int new_state ));
#endif
#ifndef YY_NO_POP_STATE
static void ng_pop_state YY_PROTO(( void ));
#endif
#ifndef YY_NO_TOP_STATE
static int ng_top_state YY_PROTO(( void ));
#endif

#else
#define YY_NO_PUSH_STATE 1
#define YY_NO_POP_STATE 1
#define YY_NO_TOP_STATE 1
#endif

#ifdef YY_MALLOC_DECL
YY_MALLOC_DECL
#else
#if __STDC__
#ifndef __cplusplus
#include <stdlib.h>
#endif
#else
/* Just try to get by without declaring the routines.  This will fail
 * miserably on non-ANSI systems for which sizeof(size_t) != sizeof(int)
 * or sizeof(void*) != sizeof(int).
 */
#endif
#endif

/* Amount of stuff to slurp up with each read. */
#ifndef YY_READ_BUF_SIZE
#define YY_READ_BUF_SIZE 8192
#endif

/* Copy whatever the last rule matched to the standard output. */

#ifndef ECHO
/* This used to be an fputs(), but since the string might contain NUL's,
 * we now use fwrite().
 */
#define ECHO (void) fwrite( ngtext, ngleng, 1, ngout )
#endif

/* Gets input and stuffs it into "buf".  number of characters read, or YY_NULL,
 * is returned in "result".
 */
#ifndef YY_INPUT
#define YY_INPUT(buf,result,max_size) \
  if ( ng_current_buffer->ng_is_interactive ) \
  { \
    int c = '*', n; \
    for ( n = 0; n < max_size && \
          (c = getc( ngin )) != EOF && c != '\n'; ++n ) \
      buf[n] = (char) c; \
    if ( c == '\n' ) \
      buf[n++] = (char) c; \
    if ( c == EOF && ferror( ngin ) ) \
      YY_FATAL_ERROR( "input in flex scanner failed" ); \
    result = n; \
  } \
  else if ( ((result = fread( buf, 1, max_size, ngin )) == 0) \
            && ferror( ngin ) ) \
    YY_FATAL_ERROR( "input in flex scanner failed" );
#endif

/* No semi-colon after return; correct usage is to write "ngterminate();" -
 * we don't want an extra ';' after the "return" because that will cause
 * some compilers to complain about unreachable statements.
 */
#ifndef ngterminate
#define ngterminate() return YY_NULL
#endif

/* Number of entries by which start-condition stack grows. */
#ifndef YY_START_STACK_INCR
#define YY_START_STACK_INCR 25
#endif

/* Report a fatal error. */
#ifndef YY_FATAL_ERROR
#define YY_FATAL_ERROR(msg) ng_fatal_error( msg )
#endif

/* Default declaration of generated scanner - a define so the user can
 * easily add parameters.
 */
#ifndef YY_DECL
#define YY_DECL int nglex YY_PROTO(( void ))
#endif

/* Code executed at the beginning of each rule, after ngtext and ngleng
 * have been set up.
 */
#ifndef YY_USER_ACTION
#define YY_USER_ACTION
#endif

/* Code executed at the end of each rule. */
#ifndef YY_BREAK
#define YY_BREAK break;
#endif

#define YY_RULE_SETUP \
  YY_USER_ACTION

YY_DECL
{
  register ng_state_type ng_current_state;
  register char *ng_cp, *ng_bp;
  register int ng_act;

#line 44 "ngin.l"


#line 564 "lex.ng.c"

  if ( ng_init )
  {
    ng_init = 0;

#ifdef YY_USER_INIT
    YY_USER_INIT;
#endif

    if ( ! ng_start )
      ng_start = 1;                     /* first start state */

    if ( ! ngin )
      ngin = stdin;

    if ( ! ngout )
      ngout = stdout;

    if ( ! ng_current_buffer )
      ng_current_buffer =
        ng_create_buffer( ngin, YY_BUF_SIZE );

    ng_load_buffer_state();
  }

  while ( 1 )                   /* loops until end-of-file is reached */
  {
    ng_cp = ng_c_buf_p;

    /* Support of ngtext. */
    *ng_cp = ng_hold_char;

    /* ng_bp points to the position in ng_ch_buf of the start of
     * the current run.
     */
    ng_bp = ng_cp;

    ng_current_state = ng_start;
ng_match:
    do
    {
      register YY_CHAR ng_c = ng_ec[YY_SC_TO_UI(*ng_cp)];
      if ( ng_accept[ng_current_state] )
      {
        ng_last_accepting_state = ng_current_state;
        ng_last_accepting_cpos = ng_cp;
      }
      while ( ng_chk[ng_base[ng_current_state] + ng_c] != ng_current_state )
      {
        ng_current_state = (int) ng_def[ng_current_state];
        if ( ng_current_state >= 27 )
          ng_c = ng_meta[(unsigned int) ng_c];
      }
      ng_current_state = ng_nxt[ng_base[ng_current_state] + (unsigned int) ng_c];
      ++ng_cp;
    }
    while ( ng_base[ng_current_state] != 41 );

ng_find_action:
    ng_act = ng_accept[ng_current_state];
    if ( ng_act == 0 )
    {                     /* have to back up */
      ng_cp = ng_last_accepting_cpos;
      ng_current_state = ng_last_accepting_state;
      ng_act = ng_accept[ng_current_state];
    }

    YY_DO_BEFORE_ACTION;


do_action:      /* This label is used only to access EOF actions. */


    switch ( ng_act )
    {     /* beginning of action switch */
    case 0 :                    /* must back up */
      /* undo the effects of YY_DO_BEFORE_ACTION */
      *ng_cp = ng_hold_char;
      ng_cp = ng_last_accepting_cpos;
      ng_current_state = ng_last_accepting_state;
      goto ng_find_action;

    case 1 :
      YY_RULE_SETUP
#line 46 "ngin.l"
      ;
      YY_BREAK
    case 2 :
      YY_RULE_SETUP
#line 47 "ngin.l"
      {noline++;}
      YY_BREAK
    case 3 :
      YY_RULE_SETUP
#line 48 "ngin.l"
      ;
      YY_BREAK
    case 4 :
      YY_RULE_SETUP
#line 49 "ngin.l"
      {nglval.ival=atol((const char *)ngtext); return (INT_VALUE);}
      YY_BREAK
    case 5 :
      YY_RULE_SETUP
#line 50 "ngin.l"
      {nglval.dval=strtod((const char *)ngtext,NULL); return (DOUBLE_VALUE);}
      YY_BREAK
    case 6 :
      YY_RULE_SETUP
#line 51 "ngin.l"
      {return (INODE);}
      YY_BREAK
    case 7 :
      YY_RULE_SETUP
#line 52 "ngin.l"
      {return (BNODE);}
      YY_BREAK
    case 8 :
      YY_RULE_SETUP
#line 53 "ngin.l"
      {return (LINE);}
      YY_BREAK
    case 9 :
      YY_RULE_SETUP
#line 54 "ngin.l"
      {return (ELEM);}
      YY_BREAK
    case 10 :
      YY_RULE_SETUP
#line 55 "ngin.l"
      {return (SIDE);}
      YY_BREAK
    case 11 :
      YY_RULE_SETUP
#line 56 "ngin.l"
      {return (TEND);}
      YY_BREAK
    case 12 :
      YY_RULE_SETUP
#line 57 "ngin.l"
      {ngerror(NULL);}
      YY_BREAK
    case 13 :
      YY_RULE_SETUP
#line 59 "ngin.l"
      ECHO;
      YY_BREAK
#line 712 "lex.ng.c"
    case YY_STATE_EOF(INITIAL) :
      ngterminate();

    case YY_END_OF_BUFFER :
    {
      /* Amount of text matched not including the EOB char. */
      int ng_amount_of_matched_text = (int) (ng_cp - ngtext_ptr) - 1;

      /* Undo the effects of YY_DO_BEFORE_ACTION. */
      *ng_cp = ng_hold_char;
      YY_RESTORE_YY_MORE_OFFSET

      if ( ng_current_buffer->ng_buffer_status == YY_BUFFER_NEW )
      {
        /* We're scanning a new file or input source.  It's
         * possible that this happened because the user
         * just pointed ngin at a new source and called
         * nglex().  If so, then we have to assure
         * consistency between ng_current_buffer and our
         * globals.  Here is the right place to do so, because
         * this is the first action (other than possibly a
         * back-up) that will match for the new input source.
         */
        ng_n_chars = ng_current_buffer->ng_n_chars;
        ng_current_buffer->ng_input_file = ngin;
        ng_current_buffer->ng_buffer_status = YY_BUFFER_NORMAL;
      }

      /* Note that here we test for ng_c_buf_p "<=" to the position
       * of the first EOB in the buffer, since ng_c_buf_p will
       * already have been incremented past the NUL character
       * (since all states make transitions on EOB to the
       * end-of-buffer state).  Contrast this with the test
       * in input().
       */
      if ( ng_c_buf_p <= &ng_current_buffer->ng_ch_buf[ng_n_chars] )
      {                   /* This was really a NUL. */
        ng_state_type ng_next_state;

        ng_c_buf_p = ngtext_ptr + ng_amount_of_matched_text;

        ng_current_state = ng_get_previous_state();

        /* Okay, we're now positioned to make the NUL
         * transition.  We couldn't have
         * ng_get_previous_state() go ahead and do it
         * for us because it doesn't know how to deal
         * with the possibility of jamming (and we don't
         * want to build jamming into it because then it
         * will run more slowly).
         */

        ng_next_state = ng_try_NUL_trans( ng_current_state );

        ng_bp = ngtext_ptr + YY_MORE_ADJ;

        if ( ng_next_state )
        {
          /* Consume the NUL. */
          ng_cp = ++ng_c_buf_p;
          ng_current_state = ng_next_state;
          goto ng_match;
        }

        else
        {
          ng_cp = ng_c_buf_p;
          goto ng_find_action;
        }
      }

      else switch ( ng_get_next_buffer() )
        {
        case EOB_ACT_END_OF_FILE :
        {
          ng_did_buffer_switch_on_eof = 0;

          if ( ngwrap() )
          {
            /* Note: because we've taken care in
             * ng_get_next_buffer() to have set up
             * ngtext, we can now set up
             * ng_c_buf_p so that if some total
             * hoser (like flex itself) wants to
             * call the scanner after we return the
             * YY_NULL, it'll still work - another
             * YY_NULL will get returned.
             */
            ng_c_buf_p = ngtext_ptr + YY_MORE_ADJ;

            ng_act = YY_STATE_EOF(YY_START);
            goto do_action;
          }

          else
          {
            if ( ! ng_did_buffer_switch_on_eof )
              YY_NEW_FILE;
          }
          break;
        }

        case EOB_ACT_CONTINUE_SCAN :
          ng_c_buf_p =
            ngtext_ptr + ng_amount_of_matched_text;

          ng_current_state = ng_get_previous_state();

          ng_cp = ng_c_buf_p;
          ng_bp = ngtext_ptr + YY_MORE_ADJ;
          goto ng_match;

        case EOB_ACT_LAST_MATCH :
          ng_c_buf_p =
            &ng_current_buffer->ng_ch_buf[ng_n_chars];

          ng_current_state = ng_get_previous_state();

          ng_cp = ng_c_buf_p;
          ng_bp = ngtext_ptr + YY_MORE_ADJ;
          goto ng_find_action;
        }
      break;
    }

    default :
      YY_FATAL_ERROR(
        "fatal flex scanner internal error--no action found" );
    }     /* end of action switch */
  }               /* end of scanning one token */
}         /* end of nglex */


/* ng_get_next_buffer - try to read in a new buffer
 *
 * Returns a code representing an action:
 *	EOB_ACT_LAST_MATCH -
 *	EOB_ACT_CONTINUE_SCAN - continue scanning from current position
 *	EOB_ACT_END_OF_FILE - end of file
 */

static int ng_get_next_buffer()
{
  register char *dest = ng_current_buffer->ng_ch_buf;
  register char *source = ngtext_ptr;
  register int number_to_move, i;
  int ret_val;

  if ( ng_c_buf_p > &ng_current_buffer->ng_ch_buf[ng_n_chars + 1] )
    YY_FATAL_ERROR(
      "fatal flex scanner internal error--end of buffer missed" );

  if ( ng_current_buffer->ng_fill_buffer == 0 )
  {               /* Don't try to fill the buffer, so this is an EOF. */
    if ( ng_c_buf_p - ngtext_ptr - YY_MORE_ADJ == 1 )
    {
      /* We matched a single character, the EOB, so
       * treat this as a final EOF.
       */
      return EOB_ACT_END_OF_FILE;
    }

    else
    {
      /* We matched some text prior to the EOB, first
       * process it.
       */
      return EOB_ACT_LAST_MATCH;
    }
  }

  /* Try to read more data. */

  /* First move last chars to start of buffer. */
  number_to_move = (int) (ng_c_buf_p - ngtext_ptr) - 1;

  for ( i = 0; i < number_to_move; ++i )
    *(dest++) = *(source++);

  if ( ng_current_buffer->ng_buffer_status == YY_BUFFER_EOF_PENDING )
    /* don't do the read, it's not guaranteed to return an EOF,
     * just force an EOF
     */
    ng_current_buffer->ng_n_chars = ng_n_chars = 0;

  else
  {
    int num_to_read =
      ng_current_buffer->ng_buf_size - number_to_move - 1;

    while ( num_to_read <= 0 )
    {                     /* Not enough room in the buffer - grow it. */
#ifdef YY_USES_REJECT
      YY_FATAL_ERROR(
        "input buffer overflow, can't enlarge buffer because scanner uses REJECT" );
#else

      /* just a shorter name for the current buffer */
      YY_BUFFER_STATE b = ng_current_buffer;

      int ng_c_buf_p_offset =
        (int) (ng_c_buf_p - b->ng_ch_buf);

      if ( b->ng_is_our_buffer )
      {
        int new_size = b->ng_buf_size * 2;

        if ( new_size <= 0 )
          b->ng_buf_size += b->ng_buf_size / 8;
        else
          b->ng_buf_size *= 2;

        b->ng_ch_buf = (char *)
                       /* Include room in for 2 EOB chars. */
                       ng_flex_realloc( (void *) b->ng_ch_buf,
                                        b->ng_buf_size + 2 );
      }
      else
        /* Can't grow it, we don't own it. */
        b->ng_ch_buf = 0;

      if ( ! b->ng_ch_buf )
        YY_FATAL_ERROR(
          "fatal error - scanner input buffer overflow" );

      ng_c_buf_p = &b->ng_ch_buf[ng_c_buf_p_offset];

      num_to_read = ng_current_buffer->ng_buf_size -
                    number_to_move - 1;
#endif
    }

    if ( num_to_read > YY_READ_BUF_SIZE )
      num_to_read = YY_READ_BUF_SIZE;

    /* Read in more data. */
    YY_INPUT( (&ng_current_buffer->ng_ch_buf[number_to_move]),
              ng_n_chars, num_to_read );

    ng_current_buffer->ng_n_chars = ng_n_chars;
  }

  if ( ng_n_chars == 0 )
  {
    if ( number_to_move == YY_MORE_ADJ )
    {
      ret_val = EOB_ACT_END_OF_FILE;
      ngrestart( ngin );
    }

    else
    {
      ret_val = EOB_ACT_LAST_MATCH;
      ng_current_buffer->ng_buffer_status =
        YY_BUFFER_EOF_PENDING;
    }
  }

  else
    ret_val = EOB_ACT_CONTINUE_SCAN;

  ng_n_chars += number_to_move;
  ng_current_buffer->ng_ch_buf[ng_n_chars] = YY_END_OF_BUFFER_CHAR;
  ng_current_buffer->ng_ch_buf[ng_n_chars + 1] = YY_END_OF_BUFFER_CHAR;

  ngtext_ptr = &ng_current_buffer->ng_ch_buf[0];

  return ret_val;
}


/* ng_get_previous_state - get the state just before the EOB char was reached */

static ng_state_type ng_get_previous_state()
{
  register ng_state_type ng_current_state;
  register char *ng_cp;

  ng_current_state = ng_start;

  for ( ng_cp = ngtext_ptr + YY_MORE_ADJ; ng_cp < ng_c_buf_p; ++ng_cp )
  {
    register YY_CHAR ng_c = (*ng_cp ? ng_ec[YY_SC_TO_UI(*ng_cp)] : 1);
    if ( ng_accept[ng_current_state] )
    {
      ng_last_accepting_state = ng_current_state;
      ng_last_accepting_cpos = ng_cp;
    }
    while ( ng_chk[ng_base[ng_current_state] + ng_c] != ng_current_state )
    {
      ng_current_state = (int) ng_def[ng_current_state];
      if ( ng_current_state >= 27 )
        ng_c = ng_meta[(unsigned int) ng_c];
    }
    ng_current_state = ng_nxt[ng_base[ng_current_state] + (unsigned int) ng_c];
  }

  return ng_current_state;
}


/* ng_try_NUL_trans - try to make a transition on the NUL character
 *
 * synopsis
 *	next_state = ng_try_NUL_trans( current_state );
 */

#ifdef YY_USE_PROTOS
static ng_state_type ng_try_NUL_trans( ng_state_type ng_current_state )
#else
static ng_state_type ng_try_NUL_trans( ng_current_state )
ng_state_type ng_current_state;
#endif
{
  register int ng_is_jam;
  register char *ng_cp = ng_c_buf_p;

  register YY_CHAR ng_c = 1;
  if ( ng_accept[ng_current_state] )
  {
    ng_last_accepting_state = ng_current_state;
    ng_last_accepting_cpos = ng_cp;
  }
  while ( ng_chk[ng_base[ng_current_state] + ng_c] != ng_current_state )
  {
    ng_current_state = (int) ng_def[ng_current_state];
    if ( ng_current_state >= 27 )
      ng_c = ng_meta[(unsigned int) ng_c];
  }
  ng_current_state = ng_nxt[ng_base[ng_current_state] + (unsigned int) ng_c];
  ng_is_jam = (ng_current_state == 26);

  return ng_is_jam ? 0 : ng_current_state;
}


#ifndef YY_NO_UNPUT
#ifdef YY_USE_PROTOS
static void ngunput( int c, register char *ng_bp )
#else
static void ngunput( c, ng_bp )
int c;
register char *ng_bp;
#endif
{
  register char *ng_cp = ng_c_buf_p;

  /* undo effects of setting up ngtext */
  *ng_cp = ng_hold_char;

  if ( ng_cp < ng_current_buffer->ng_ch_buf + 2 )
  {               /* need to shift things up to make room */
                  /* +2 for EOB chars. */
    register int number_to_move = ng_n_chars + 2;
    register char *dest = &ng_current_buffer->ng_ch_buf[
      ng_current_buffer->ng_buf_size + 2];
    register char *source =
      &ng_current_buffer->ng_ch_buf[number_to_move];

    while ( source > ng_current_buffer->ng_ch_buf )
      *--dest = *--source;

    ng_cp += (int) (dest - source);
    ng_bp += (int) (dest - source);
    ng_current_buffer->ng_n_chars =
      ng_n_chars = ng_current_buffer->ng_buf_size;

    if ( ng_cp < ng_current_buffer->ng_ch_buf + 2 )
      YY_FATAL_ERROR( "flex scanner push-back overflow" );
  }

  *--ng_cp = (char) c;


  ngtext_ptr = ng_bp;
  ng_hold_char = *ng_cp;
  ng_c_buf_p = ng_cp;
}
#endif  /* ifndef YY_NO_UNPUT */


#ifdef __cplusplus
static int nginput()
#else
static int input()
#endif
{
  int c;

  *ng_c_buf_p = ng_hold_char;

  if ( *ng_c_buf_p == YY_END_OF_BUFFER_CHAR )
  {
    /* ng_c_buf_p now points to the character we want to return.
     * If this occurs *before* the EOB characters, then it's a
     * valid NUL; if not, then we've hit the end of the buffer.
     */
    if ( ng_c_buf_p < &ng_current_buffer->ng_ch_buf[ng_n_chars] )
      /* This was really a NUL. */
      *ng_c_buf_p = '\0';

    else
    {                     /* need more input */
      int offset = ng_c_buf_p - ngtext_ptr;
      ++ng_c_buf_p;

      switch ( ng_get_next_buffer() )
      {
      case EOB_ACT_LAST_MATCH :
        /* This happens because ng_g_n_b()
         * sees that we've accumulated a
         * token and flags that we need to
         * try matching the token before
         * proceeding.  But for input(),
         * there's no matching to consider.
         * So convert the EOB_ACT_LAST_MATCH
         * to EOB_ACT_END_OF_FILE.
         */

        /* Reset buffer status. */
        ngrestart( ngin );

      /* fall through */

      case EOB_ACT_END_OF_FILE :
      {
        if ( ngwrap() )
          return EOF;

        if ( ! ng_did_buffer_switch_on_eof )
          YY_NEW_FILE;
#ifdef __cplusplus
        return nginput();
#else
        return input();
#endif
      }

      case EOB_ACT_CONTINUE_SCAN :
        ng_c_buf_p = ngtext_ptr + offset;
        break;
      }
    }
  }

  c = *(unsigned char *) ng_c_buf_p;            /* cast for 8-bit char's */
  *ng_c_buf_p = '\0';           /* preserve ngtext */
  ng_hold_char = *++ng_c_buf_p;


  return c;
}


#ifdef YY_USE_PROTOS
void ngrestart( FILE *input_file )
#else
void ngrestart( input_file )
FILE *input_file;
#endif
{
  if ( ! ng_current_buffer )
    ng_current_buffer = ng_create_buffer( ngin, YY_BUF_SIZE );

  ng_init_buffer( ng_current_buffer, input_file );
  ng_load_buffer_state();
}


#ifdef YY_USE_PROTOS
void ng_switch_to_buffer( YY_BUFFER_STATE new_buffer )
#else
void ng_switch_to_buffer( new_buffer )
YY_BUFFER_STATE new_buffer;
#endif
{
  if ( ng_current_buffer == new_buffer )
    return;

  if ( ng_current_buffer )
  {
    /* Flush out information for old buffer. */
    *ng_c_buf_p = ng_hold_char;
    ng_current_buffer->ng_buf_pos = ng_c_buf_p;
    ng_current_buffer->ng_n_chars = ng_n_chars;
  }

  ng_current_buffer = new_buffer;
  ng_load_buffer_state();

  /* We don't actually know whether we did this switch during
   * EOF (ngwrap()) processing, but the only time this flag
   * is looked at is after ngwrap() is called, so it's safe
   * to go ahead and always set it.
   */
  ng_did_buffer_switch_on_eof = 1;
}


#ifdef YY_USE_PROTOS
void ng_load_buffer_state( void )
#else
void ng_load_buffer_state()
#endif
{
  ng_n_chars = ng_current_buffer->ng_n_chars;
  ngtext_ptr = ng_c_buf_p = ng_current_buffer->ng_buf_pos;
  ngin = ng_current_buffer->ng_input_file;
  ng_hold_char = *ng_c_buf_p;
}


#ifdef YY_USE_PROTOS
YY_BUFFER_STATE ng_create_buffer( FILE *file, int size )
#else
YY_BUFFER_STATE ng_create_buffer( file, size )
FILE *file;
int size;
#endif
{
  YY_BUFFER_STATE b;

  b = (YY_BUFFER_STATE) ng_flex_alloc( sizeof( struct ng_buffer_state ) );
  if ( ! b )
    YY_FATAL_ERROR( "out of dynamic memory in ng_create_buffer()" );

  b->ng_buf_size = size;

  /* ng_ch_buf has to be 2 characters longer than the size given because
   * we need to put in 2 end-of-buffer characters.
   */
  b->ng_ch_buf = (char *) ng_flex_alloc( b->ng_buf_size + 2 );
  if ( ! b->ng_ch_buf )
    YY_FATAL_ERROR( "out of dynamic memory in ng_create_buffer()" );

  b->ng_is_our_buffer = 1;

  ng_init_buffer( b, file );

  return b;
}


#ifdef YY_USE_PROTOS
void ng_delete_buffer( YY_BUFFER_STATE b )
#else
void ng_delete_buffer( b )
YY_BUFFER_STATE b;
#endif
{
  if ( ! b )
    return;

  if ( b == ng_current_buffer )
    ng_current_buffer = (YY_BUFFER_STATE) 0;

  if ( b->ng_is_our_buffer )
    ng_flex_free( (void *) b->ng_ch_buf );

  ng_flex_free( (void *) b );
}


#ifndef YY_ALWAYS_INTERACTIVE
#ifndef YY_NEVER_INTERACTIVE
extern int isatty YY_PROTO(( int ));
#endif
#endif

#ifdef YY_USE_PROTOS
void ng_init_buffer( YY_BUFFER_STATE b, FILE *file )
#else
void ng_init_buffer( b, file )
YY_BUFFER_STATE b;
FILE *file;
#endif


{
  ng_flush_buffer( b );

  b->ng_input_file = file;
  b->ng_fill_buffer = 1;

#if YY_ALWAYS_INTERACTIVE
  b->ng_is_interactive = 1;
#else
#if YY_NEVER_INTERACTIVE
  b->ng_is_interactive = 0;
#else
  b->ng_is_interactive = file ? (isatty( fileno(file) ) > 0) : 0;
#endif
#endif
}


#ifdef YY_USE_PROTOS
void ng_flush_buffer( YY_BUFFER_STATE b )
#else
void ng_flush_buffer( b )
YY_BUFFER_STATE b;
#endif

{
  if ( ! b )
    return;

  b->ng_n_chars = 0;

  /* We always need two end-of-buffer characters.  The first causes
   * a transition to the end-of-buffer state.  The second causes
   * a jam in that state.
   */
  b->ng_ch_buf[0] = YY_END_OF_BUFFER_CHAR;
  b->ng_ch_buf[1] = YY_END_OF_BUFFER_CHAR;

  b->ng_buf_pos = &b->ng_ch_buf[0];

  b->ng_at_bol = 1;
  b->ng_buffer_status = YY_BUFFER_NEW;

  if ( b == ng_current_buffer )
    ng_load_buffer_state();
}


#ifndef YY_NO_SCAN_BUFFER
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE ng_scan_buffer( char *base, ng_size_t size )
#else
YY_BUFFER_STATE ng_scan_buffer( base, size )
char *base;
ng_size_t size;
#endif
{
  YY_BUFFER_STATE b;

  if ( size < 2 ||
       base[size-2] != YY_END_OF_BUFFER_CHAR ||
       base[size-1] != YY_END_OF_BUFFER_CHAR )
    /* They forgot to leave room for the EOB's. */
    return 0;

  b = (YY_BUFFER_STATE) ng_flex_alloc( sizeof( struct ng_buffer_state ) );
  if ( ! b )
    YY_FATAL_ERROR( "out of dynamic memory in ng_scan_buffer()" );

  b->ng_buf_size = size - 2;            /* "- 2" to take care of EOB's */
  b->ng_buf_pos = b->ng_ch_buf = base;
  b->ng_is_our_buffer = 0;
  b->ng_input_file = 0;
  b->ng_n_chars = b->ng_buf_size;
  b->ng_is_interactive = 0;
  b->ng_at_bol = 1;
  b->ng_fill_buffer = 0;
  b->ng_buffer_status = YY_BUFFER_NEW;

  ng_switch_to_buffer( b );

  return b;
}
#endif


#ifndef YY_NO_SCAN_STRING
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE ng_scan_string( ngconst char *ng_str )
#else
YY_BUFFER_STATE ng_scan_string( ng_str )
ngconst char *ng_str;
#endif
{
  int len;
  for ( len = 0; ng_str[len]; ++len )
    ;

  return ng_scan_bytes( ng_str, len );
}
#endif


#ifndef YY_NO_SCAN_BYTES
#ifdef YY_USE_PROTOS
YY_BUFFER_STATE ng_scan_bytes( ngconst char *bytes, int len )
#else
YY_BUFFER_STATE ng_scan_bytes( bytes, len )
ngconst char *bytes;
int len;
#endif
{
  YY_BUFFER_STATE b;
  char *buf;
  ng_size_t n;
  int i;

  /* Get memory for full buffer, including space for trailing EOB's. */
  n = len + 2;
  buf = (char *) ng_flex_alloc( n );
  if ( ! buf )
    YY_FATAL_ERROR( "out of dynamic memory in ng_scan_bytes()" );

  for ( i = 0; i < len; ++i )
    buf[i] = bytes[i];

  buf[len] = buf[len+1] = YY_END_OF_BUFFER_CHAR;

  b = ng_scan_buffer( buf, n );
  if ( ! b )
    YY_FATAL_ERROR( "bad buffer in ng_scan_bytes()" );

  /* It's okay to grow etc. this buffer, and we should throw it
   * away when we're done.
   */
  b->ng_is_our_buffer = 1;

  return b;
}
#endif


#ifndef YY_NO_PUSH_STATE
#ifdef YY_USE_PROTOS
static void ng_push_state( int new_state )
#else
static void ng_push_state( new_state )
int new_state;
#endif
{
  if ( ng_start_stack_ptr >= ng_start_stack_depth )
  {
    ng_size_t new_size;

    ng_start_stack_depth += YY_START_STACK_INCR;
    new_size = ng_start_stack_depth * sizeof( int );

    if ( ! ng_start_stack )
      ng_start_stack = (int *) ng_flex_alloc( new_size );

    else
      ng_start_stack = (int *) ng_flex_realloc(
        (void *) ng_start_stack, new_size );

    if ( ! ng_start_stack )
      YY_FATAL_ERROR(
        "out of memory expanding start-condition stack" );
  }

  ng_start_stack[ng_start_stack_ptr++] = YY_START;

  BEGIN(new_state);
}
#endif


#ifndef YY_NO_POP_STATE
static void ng_pop_state()
{
  if ( --ng_start_stack_ptr < 0 )
    YY_FATAL_ERROR( "start-condition stack underflow" );

  BEGIN(ng_start_stack[ng_start_stack_ptr]);
}
#endif


#ifndef YY_NO_TOP_STATE
static int ng_top_state()
{
  return ng_start_stack[ng_start_stack_ptr - 1];
}
#endif

#ifndef YY_EXIT_FAILURE
#define YY_EXIT_FAILURE 2
#endif

#ifdef YY_USE_PROTOS
static void ng_fatal_error( ngconst char msg[] )
#else
static void ng_fatal_error( msg )
char msg[];
#endif
{
  (void) fprintf( stderr, "%s\n", msg );
  exit( YY_EXIT_FAILURE );
}



/* Redefine ngless() so it works in section 3 code. */

#undef ngless
#define ngless(n) \
  do \
  { \
    /* Undo effects of setting up ngtext. */ \
    ngtext[ngleng] = ng_hold_char; \
    ng_c_buf_p = ngtext + n; \
    ng_hold_char = *ng_c_buf_p; \
    *ng_c_buf_p = '\0'; \
    ngleng = n; \
  } \
  while ( 0 )


/* Internal utility routines. */

#ifndef ngtext_ptr
#ifdef YY_USE_PROTOS
static void ng_flex_strncpy( char *s1, ngconst char *s2, int n )
#else
static void ng_flex_strncpy( s1, s2, n )
char *s1;
ngconst char *s2;
int n;
#endif
{
  register int i;
  for ( i = 0; i < n; ++i )
    s1[i] = s2[i];
}
#endif

#ifdef YY_NEED_STRLEN
#ifdef YY_USE_PROTOS
static int ng_flex_strlen( ngconst char *s )
#else
static int ng_flex_strlen( s )
ngconst char *s;
#endif
{
  register int n;
  for ( n = 0; s[n]; ++n )
    ;

  return n;
}
#endif


#ifdef YY_USE_PROTOS
static void *ng_flex_alloc( ng_size_t size )
#else
static void *ng_flex_alloc( size )
ng_size_t size;
#endif
{
  return (void *) malloc( size );
}

#ifdef YY_USE_PROTOS
static void *ng_flex_realloc( void *ptr, ng_size_t size )
#else
static void *ng_flex_realloc( ptr, size )
void *ptr;
ng_size_t size;
#endif
{
  /* The cast to (char *) in the following accommodates both
   * implementations that use char* generic pointers, and those
   * that use void* generic pointers.  It works with the latter
   * because both ANSI C and C++ allow castless assignment from
   * any pointer type to void*, and deal with argument conversions
   * as though doing an assignment.
   */
  return (void *) realloc( (char *) ptr, size );
}

#ifdef YY_USE_PROTOS
static void ng_flex_free( void *ptr )
#else
static void ng_flex_free( ptr )
void *ptr;
#endif
{
  free( ptr );
}

#if YY_MAIN
int main()
{
  nglex();
  return 0;
}
#endif
#line 59 "ngin.l"


NP_Error (int *line, char *text)
{
  *line=noline;
  strcpy(text,ngtext);
  return;
}
