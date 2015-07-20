/*
 *  $Id: VXPrivate.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  This file contains private declarations for VX routines.
 */

#ifndef V_VXPrivate_h
#define V_VXPrivate_h 1

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

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/colormap.h"
#include "viaio/mu.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* From X11R5 Xt and Motif: */
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <Xm/Xm.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Macro(s): */
#define XtVCMW  XtVaCreateManagedWidget
#define XtVCW	XtVaCreateWidget
#define XtVGV	XtVaGetValues
#define XtVSV	XtVaSetValues


/*
 *  Data structure declarations.
 */

/*
 * Type: VRec
 *
 * Record for storing image related data.
 */
typedef struct {
    VImage image;		/* Image currently managed by VX. */
    VBand band;			/* Band of image currently managed by VX. */
    float row_scale;		/* # of screen pixel per image pixel in
				   row-dimension. */
    float zoom_level;		/* Zoom level in VXSetImage */
} VRec;


/*
 * Type: XRec
 *
 * Record for storing widget related data.
 */
typedef struct {
    XtAppContext appContext;	/* Application context. */
    VColormap vcolormap;	/* colormap, depth, visual, and colors */

    Widget topLevel;		/* Application shell. */
    Widget encloseAll;
    Widget mainWindow;		/* Main window. */
    Widget menuBar;		/* Menu bar. */
    Widget imageViewFrame;	/* Frame of imageView. */
    Widget imageView;		/* VImageView widget. */
    Widget msgAreaFrame;	/* Frame of msgArea (+ scrolled window). */
    Widget msgArea;		/* Text widget used as message area. */

    int init_width;		/* Init width of imageView. */
    int init_height;		/* Init height of imageView. */
    int cur_width;		/* Current width of imageView. */
    int cur_height;		/* Current height of imageView. */

    int msg_area_nlines;	/* Number of lines visible in message area. */

    Window busyWindow;		/* For displaying busy cursor. */
} XRec;


/*
 * Type: ORec
 *
 * Record for storing overlay related data.
 */
typedef struct {
    VBoolean pixmap_consistent;	/* pixmap for storing overlays is valid */
    GC gc;			/* graphic context ID */
    Pixmap pixmap;		/* pixmap for storing overlays */
    VBoolean pixmap_allocated;	/* TRUE iff pixmap is allocated */
} ORec;


/*
 * Type: AppRec
 *
 * Record for storing VX application data.
 */
typedef struct {
    VBoolean initialized;	/* TRUE: VXInit() has been called */
    VBoolean in_main_loop;	/* TRUE: VXAppMainLoop() is called */
    VRec v;			/* Vista image related data */
    XRec x;			/* X widget related data */
    ORec o;			/* overlay related data */
} AppRec;


/*
 * Global variables.
 */

extern AppRec VX_App;		/* application data */


/*
 *  Declarations of routines.
 */

/* From the menu module: */
extern VBoolean VX_InitMenu (void);

/* From the input module: */
extern VBoolean VX_InitInput (void);

/* From the image module: */
extern void VX_Zoomed (Widget, XtPointer, XtPointer);

/* From the line module: */
extern VBoolean VX_InitLine (void);
extern void VX_GetLineGC (void);
extern void VX_RedrawLines (void);
extern void VX_StoreLines (void);
extern void VX_RestoreLines(void);

/* From the text module: */
extern VBoolean VX_InitText (void);
extern void VX_GetTextGC (void);
extern void VX_RedrawTexts (void);
extern void VX_StoreTexts (void);
extern void VX_RestoreTexts(void);

/* From the overlays module: */
extern void VX_RedrawOverlays (Widget, XtPointer, XtPointer);

/* From the dialog module: */
extern VBoolean VX_InitDialog (void);
extern void VX_Warning (VStringConst);

#ifdef __cplusplus
}
#endif

#endif /* V_VXPrivate_h */
