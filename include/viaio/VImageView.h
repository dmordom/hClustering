/*
 *  $Id: VImageView.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  This file defines the public interface to the VImageView widget.
 */

#ifndef V_VImageView_h
#define V_VImageView_h 1

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
 *  Authors: Arthur Pope, Daniel Ko, Dan Razzell,
 *	     UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/VImage.h"

/* From Xt: */
#include <X11/Intrinsic.h>

/* For portability: */
#include <X11/Xfuncproto.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Resources:

 Name		     Class		RepType		Default Value
 ----		     -----		-------		-------------
 (from Object)
 destroyCallback     Callback		Pointer		NULL

 (from RectObj)
 ancestorSensitive   Sensitive		Boolean		True
 borderWidth	     BorderWidth	Dimension	1
 height		     Height		Dimension	0
 sensitive	     Sensitive		Boolean		True
 width		     Width		Dimension	0
 x		     Position		Position	0
 y		     Position		Position	0

 (from Core)
 accelerators	     Accelerators	AcceleratorTabl	NULL
 background	     Background		Pixel		White
 backgroundPixmap    Pixmap		Pixmap		XtUnspecifiedPixmap
 borderColor	     BorderColor	Pixel		XtDefaultForeground
 colormap	     Colormap		Colormap	CopyFromParent
 depth		     Depth		Int		CopyFromParent
 mappedWhenManaged   MappedWhenManaged	Boolean		True
 screen		     Screen		Screen		XtCopyScreen
 translations	     Translations	TranslationTabl

 (from VImageView)
 absolute	     Absolute		Boolean		TRUE
 band		     Band		Int		0
 columnCenter        ColumnCenter       Int             0
 cursor		     Cursor		Cursor		None
 exposeCallback	     Callback		Pointer	        NULL
 image		     Image		Pointer		NULL
 inputCallback	     Callback		Pointer		NULL
 moveZoomCenterCallback Callback        Pointer         NULL
 proportion	     Proportion		Boolean		True
 resize		     Resize		Boolean		False
 rowCenter           RowCenter          Int             0
 usePixmap	     UsePixmap		Boolean		TRUE
 vColormap	     VColormap		VColormap	NULL
 zoomInCallback      Callback           Pointer         NULL
 zoomLevel           ZoomLevel          Int             100
 zoomOutCallback     Callback           Pointer         NULL
*/

/* Resource names specific to VImageView: */
#define VxNabsolute "absolute"
#define VxNband "band"
#define VxNcolumnCenter "columnCenter"
#define VxNcursor "cursor"
#define VxNexposeCallback "exposeCallback"
#define VxNforeground "foreground"
#define VxNimage "image"
#define VxNinputCallback "inputCallback"
#define VxNmoveZoomCenterCallback "moveZoomCenterCallback"
#define VxNproportion "proportion"
#define VxNresize "resize"
#define VxNrowCenter "rowCenter"
#define VxNusePixmap "usePixmap"
#define VxNvColormap "vColormap"
#define VxNzoomInCallback "zoomInCallback"
#define VxNzoomOutCallback "zoomOutCallback"
#define VxNzoomLevel "zoomLevel"

/* Class names specific to VImageView: */
#define VxCAbsolute "Absolute"
#define VxCBand "Band"
#define VxCCallback "Callback"
#define VxCColumnCenter "ColumnCenter"
#define VxCCursor "Cursor"
#define VxCImage "Image"
#define VxCProportion "Proportion"
#define VxCResize "Resize"
#define VxCRowCenter "RowCenter"
#define VxCUsePixmap "UsePixmap"
#define VxCVColormap "VColormap"
#define VxCZoomLevel "ZoomLevel"

/* Specific VImageViewWidget class and instance datatypes: */
typedef struct V_ImageViewClassRec *VImageViewWidgetClass;
typedef struct V_ImageViewRec *VImageViewWidget;

/* The class constant: */
extern WidgetClass vImageViewWidgetClass;

/* Convenience routines: */
extern VBoolean VImageViewWindowToImage (
#if NeedFunctionPrototypes
    Widget 		/* w */,
    int			/* x */,
    int			/* y */,
    double *		/* rowp */,
    double *		/* columnp */
#endif
);

extern VBoolean VImageViewClipToImage (
#if NeedFunctionPrototypes
    Widget 		/* w */,
    int			/* x */,
    int			/* y */,
    double *		/* rowp */,
    double *		/* columnp */
#endif
);

extern VBoolean VImageViewImageToWindow (
#if NeedFunctionPrototypes
    Widget 		/* w */,
    double 		/* row */,
    double 		/* column */,
    int	*		/* xp */,
    int	*		/* yp */
#endif
);

extern VBoolean VImageViewPixelSize (
#if NeedFunctionPrototypes
    Widget		/* w */,
    double *		/* width */,
    double *		/* height */
#endif
);

extern void VImageViewRedraw (
#if NeedFunctionPrototypes
    Widget              /* w */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_VImageView_h */
