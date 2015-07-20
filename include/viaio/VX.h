/*
 * $Id: VX.h 726 2004-03-08 13:12:45Z lohmann $
 *
 * Definitions associated with the VX library.
 */

#ifndef V_VX_h
#define V_VX_h 1

/*
 *  Copyright 1993, 1994 University of British Columbia
 *
 *  Permission to use, copy, modify, distribute, and sell this software and its
 *  documentation for any purpose is hereby granted without fee, provided that
 *  the above copyright notice appears in all copies and that both that
 *  copyright notice and this permission notice appear in supporting
 *  documentation. UBC makes no representations about the suitability of this
 *  software for any purpose. It is provided "as is" without express or
 *  implied warranty.
 *
 *  Author: Daniel Ko, UBC Laboratory for Computational Intelligence
 */

/* From Vista library: */
#include "viaio/Vlib.h"
#include "viaio/VEdges.h"
#include "viaio/VImage.h"

/* From standard C library: */
#include <stdio.h>

/* From Motif's VaSimple.h, and necessary under Motif V1.1.2 and X11R5
   to get Xm/Xm.h to compile under K&R C: */
#if !defined(__STDC__) && !defined(_NO_PROTO)
#define _NO_PROTO
#endif
#if defined(__STDC__) && defined(_NO_PROTO)
#undef _NO_PROTO
#endif

/* From Motif: */
#include <Xm/Xm.h>

/* For portability: */
#include <X11/Xfuncproto.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Symbolic constants: */
#define VXNmessageAreaNLines "messageAreaNLines"
#define VXCMessageAreaNLines "MessageAreaNLines"


/*
 *  Declarations of data structures.
 */

/*
 * VXInputType
 */
typedef enum {
    VXIkeyPress	      = 0,	/* when a key is pressed	   */
    VXIbuttonPress    = 1,	/* when a mouse button is pressed  */
    VXIbuttonRelease  = 2,	/* when a mouse button is released */
    VXIpointerMotion  = 3,	/* when the pointer is moved	   */
    VXnInputTypes     = 4
} VXInputType;

/*
 * VXModifierMask
 */
typedef enum {
    VXMshift   = ShiftMask,	/* the Shift key    */
    VXMctrl    = ControlMask,	/* the Ctrl key	    */
    VXMbutton1 = Button1Mask,	/* pointer button 1 */
    VXMbutton2 = Button2Mask,	/* pointer button 2 */
    VXMbutton3 = Button3Mask,	/* pointer button 3 */
    VXMbutton4 = Button4Mask,	/* pointer button 4 */
    VXMbutton5 = Button5Mask,	/* pointer button 5 */
    VXMmod1    = Mod1Mask,	/* modifier 1	    */
    VXMmod2    = Mod2Mask,	/* modifier 2	    */
    VXMmod3    = Mod3Mask,	/* modifier 3	    */
    VXMmod4    = Mod4Mask,	/* modifier 4	    */
    VXMmod5    = Mod5Mask	/* modifier 5	    */
} VXModifierMask;

/*
 * VXInputData
 */
typedef struct {

    VXInputType input_type;	/*
				 * - type of input that triggered the
				 *   callback
				 */

    int value;			/*
				 * - the ASCII code of the key being pressed
				 * - or the ID of the button being pressed
				 */

    VXModifierMask modifiers;	/*
				 * - the state of the buttons and
				 *   modifier keys when the input
				 *   event occurs
				 * - bit 1 means pressed and
				 *   bit 0 means un-pressed
				 */

    int row;			/*
				 * - row-coordinate of the location
				 *   where the input event occurs
				 */

    int column;			/*
				 * - column-coordinate of the location
				 *   where the input event occurs
				 */

} VXInputDataRec, *VXInputData;

/*
 * VXInputCallback
 */
typedef void (*VXInputCallback)(
#if NeedFunctionPrototypes
    VXInputData		/* input_data */,
    VPointer		/* client_data */
#endif
);

/*
 * VXMenuCallback
 */
typedef void (*VXMenuCallback)(
#if NeedFunctionPrototypes
    VPointer		/* client_data */
#endif
);

/*
 * VXAnswer
 */
typedef enum {
    VXAyes,			/* answer is yes    */
    VXAno,			/* answer is no	    */
    VXAcancel			/* answer is cancel */
} VXAnswer;


/*
 *  VXWarning is obsolete -- use VWarning instead.
 */

#define VXWarning VWarning


/*
 *  Declarations of library routines.
 */

/* From VXDialog.c: */

extern void VXPopupMessageBox (
#if NeedFunctionPrototypes
    VStringConst	/* title */,
    VStringConst	/* message */
#endif
);

extern VString VXPopupInputBox (
#if NeedFunctionPrototypes
    VStringConst	/* title */,
    VStringConst	/* prompt */,
    VStringConst	/* text */
#endif
);

extern VXAnswer VXPopupYesNoBox (
#if NeedFunctionPrototypes
    VStringConst	/* title */,
    VStringConst	/* question */
#endif
);

extern void VXPopupTextBox (
#if NeedFunctionPrototypes
    int			/* nrows */,
    int			/* ncolumns */,
    VStringConst	/* title */,
    VStringConst	/* text */
#endif
);

extern VString VXPopupFileBox (
#if NeedFunctionPrototypes
    VStringConst	/* title */
#endif
);

/* From VXImage.c: */

extern VBoolean VXSetImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    double		/* zoom */,
    int			/* row_center */,
    int			/* column_center */
#endif
);

/* From VXInit.c: */

extern void VXInit (
#if NeedFunctionPrototypes
    VStringConst	/* class */,
    VStringConst *	/* default_res */,
    int *		/* argc */,
    char **		/* argv */
#endif
);

extern void VXAppMainLoop (
#if NeedFunctionPrototypes
   void
#endif
);

extern void VXReportValidOptions (
#if NeedFunctionPrototypes
   void
#endif
);

/* From VXInput.c: */

extern void VXAddInputCallback (
#if NeedFunctionPrototypes
    VXInputType		/* input_type */,
    VXInputCallback	/* callback */,
    VPointer		/* client_data */
#endif
);

/* From VXLine.c: */

extern VBoolean VXSetLineColor (
#if NeedFunctionPrototypes
    VStringConst	/* color_name */
#endif
);

extern void VXSetLineWidth (
#if NeedFunctionPrototypes
    double		/* width */
#endif
);

extern void VXClearLines (
#if NeedFunctionPrototypes
    void
#endif
);

extern VBoolean VXDrawLine (
#if NeedFunctionPrototypes
    double		/* r1 */,
    double		/* c1 */,
    double		/* r2 */,
    double		/* c2 */
#endif
);

extern VBoolean VXDrawEdges (
#if NeedFunctionPrototypes
   VEdges		/* edges */
#endif
);

/* From VXMenu.c: */

extern void VXAddMenu (
#if NeedVarargsPrototypes
    VStringConst	/* menu_name */,
    ...
#endif
);

/* From VXMisc.c: */

extern void VXDisplayMessage (
#if NeedVarargsPrototypes
    VBooleanPromoted	/* overwrite */,
    VStringConst	/* format */,
    ...
#endif
);

extern void VXShowMessageArea (
#if NeedFunctionPrototypes
    void
#endif
);

extern void VXHideMessageArea (
#if NeedFunctionPrototypes
    void
#endif
);

extern Widget VXGetImageViewWidget (
#if NeedFunctionPrototypes
    void
#endif
);

extern Widget VXGetApplicationShell (
#if NeedFunctionPrototypes
    void
#endif
);

extern VBoolean VXIsColorDisplay (
#if NeedFunctionPrototypes
    void
#endif
);

/* From VXOverlays.c: */

extern void VXStoreOverlays (
#if NeedFunctionPrototypes
    void
#endif
);

extern void VXRestoreOverlays (
#if NeedFunctionPrototypes
    void
#endif
);

/* From VXText.c: */

extern VBoolean VXSetTextFont (
#if NeedFunctionPrototypes
    VStringConst	/* fontname */
#endif
);

extern VBoolean VXSetTextColor (
#if NeedFunctionPrototypes
    VStringConst	/* color_name */
#endif
);

extern void VXClearTexts (
#if NeedFunctionPrototypes
    void
#endif
);

extern VBoolean VXDrawText (
#if NeedFunctionPrototypes
    VStringConst	/* str */,
    double		/* r */,
    double		/* c */
#endif
);

#ifdef __cplusplus
}
#endif

#endif	/* V_VX_h */
