/*
 *  $Id: colormap.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  This file contains definitions for Vista library support of X Windows
 *  colormap allocation.
 */

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
 *  Authors: Dan Razzell, Art Pope
 *	     UBC Laboratory for Computational Intelligence
 */

#ifndef V_colormap_h
#define V_colormap_h 1

/* From X Windows libraries: */
#include <X11/Xlib.h>
#include <X11/Intrinsic.h>

/* For portability: */
#include <X11/Xfuncproto.h>

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Colormap type definition.
 */

typedef struct V_ColormapRec {
    Screen *screen;
    Atom property;
    XVisualInfo vinfo;

    /* A standard colormap describes the color mapping: */
    XStandardColormap stdcmap;

    /* If the colors allocated are not contiguous, the pixel values computed
       using stdcmap are used to index the following table, from which the
       actual pixel values are obtained: */
    unsigned long *indcmap;
    VBoolean *indcmap_alloced;		/* which need freeing */

    /* All grayscale rendering uses this ramp of pixel values: */
    int ngrays;
    unsigned long *invgmap;
    VBoolean *invgmap_alloced;		/* which need freeing */
} *VColormap;


/*
 *  Colormap accessors.
 */

#define VColormapColormap(vc)	((vc)->stdcmap.colormap)
#define VColormapDepth(vc)	((vc)->vinfo.depth)
#define VColormapProperty(vc)	((vc)->property)
#define VColormapVisual(vc)	((vc)->vinfo.visual)


/*
 *  Declarations of library routines.
 */

extern VColormap VCreateColormap (
#if NeedFunctionPrototypes
    Screen *		/* screen */,
    Atom		/* property */,
    long		/* vinfo_mask */,
    XVisualInfo *	/* vinfo_template */
#endif
);

extern void VDestroyColormap (
#if NeedFunctionPrototypes
    VColormap		/* vcolormap */
#endif
);

extern void VColormapRGBPixel (
#if NeedFunctionPrototypes
    VColormap		/* vcolormap */,
    XColor *		/* color */
#endif
);

extern void VColormapGrayPixel (
#if NeedFunctionPrototypes
    VColormap		/* vcolormap */,
    XColor *		/* color */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_colormap_h */
