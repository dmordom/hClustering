/*
 *  $Id: VGraph.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  Definitions associated with Graph: their representation in files and
 *  in memory, and operations that can be performed with them.
 */

#ifndef V_VGraph_h
#define V_VGraph_h 1

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
 *  Author: David Lowe, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/file.h"

/* For portability: */
#include <X11/Xfuncproto.h>

#ifdef __cplusplus
extern "C" {
#endif


/*  Structures for representing Graph in memory.
 *  -------------------------------------------
 *
 *  A Graph is stored as a table of node records. Nodes are connected by an
 *  adjancency list. The actual data of a node is stored in an opaque field
 *  of the node; the size of this field is given in the Graph structure.
 *  Thus, the representation of the node data is unknown to this structure.
 */

typedef struct V_GraphRec {
    int nnodes;			/* number of nodes */
    int nfields;		/* size of fields in a node�s private area */
    VRepnKind node_repn;	/* data representation in a node */
    VAttrList attributes;	/* list of other attributes */
    struct VNodestruct **table;	/* node table of Graph */
    int size;			/* number of places in table */
    int lastUsed;		/* last entry used in table */
    int iter;			/* iteration counter in sequential access */
    int useWeights;		/* TRUE iff weights are used */
} VGraphRec, *VGraph;

typedef struct VNodebaseStruct {
    unsigned int hops: 31;	/* numbor of hops in this node */
    unsigned int visited: 1;	/* true iff seen before */
    VFloat weight;		/* weight of this node */
    struct VAdjstruct *head;
} VNodeBaseRec, *VNodeBase;

typedef struct VNodestruct {
    VNodeBaseRec base;
    char data[1];		/* private data area of node starts here */
} VNodeRec, *VNode;

typedef struct VAdjstruct {
    unsigned int id;		/* node reference */
    VFloat weight;		/* weight of this node */
    struct VAdjstruct *next;	/* list of adjacent nodes */
} VAdjRec, *VAdjacency;

/*
 *  Attributes used to represent Graph
 *  ----------------------------------
 */

/* Attribute type names: */
#define VGraphAttr		"Graph"
#define VNGraphNodesAttr	"nnodes"
#define VNNodeFieldsAttr	"nfields"
#define VNNodeWeightsAttr	"useWeights"


/*
 *  Macros for accessing Graph attributes in memory.
 *  -----------------------------------------------
 */

#define VGraphNNodes(graph)		(graph->nnodes)

#define VGraphNFields(graph)		(graph->nfields)

#define VGraphNSize(graph)		(graph->size)

#define VGraphAttrList(graph)		(graph->attributes)

#define VGraphGetNode(graph, nid)	(graph->table[nid-1])

#define VGraphNodeIsFree(graph, nid)	(graph->table[nid-1] == 0)

#define VNodeRepn(graph)		(graph->node_repn)

#define VNodeSize(graph) \
	(sizeof(VNodeBaseRec) + (graph->nfields * VRepnPrecision(graph->node_repn)) / 8)

#define VNodeTestVisit(node)	(((VNodeBase)node)->visited == TRUE)

#define VNodeSetVisit(node)	(((VNodeBase)node)->visited = TRUE)

#define VNodeClearVisit(node)	(((VNodeBase)node)->visited = FALSE)

/*
 *  Declarations of library routines.
 */

/* From Graph.c: */

extern VGraph VCreateGraph (
#if NeedFunctionPrototypes
    int			/* nnodes */,
    int			/* nfields */,
    VRepnKind		/* node_repn */,
    int 		/* save weights */
#endif
);

extern VGraph VCopyGraph (
#if NeedFunctionPrototypes
    VGraph		/* Graph */
#endif
);

extern void VDestroyGraph (
#if NeedFunctionPrototypes
    VGraph		/* Graph */
#endif
);

extern int VReadGraphs (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList *		/* attributes */,
    VGraph **		/* graphs */
#endif
);

extern VBoolean VWriteGraphs (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList		/* attributes */,
    int			/* ngraphs */,
    VGraph *		/* graphs */
#endif
);

extern int VGraphLookupNode (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    VNode		/*  node */
#endif
);

extern int VGraphAddNode (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    VNode		/*  node */
#endif
);

extern int VGraphAddNodeAt (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    VNode		/*  node */,
    int			/*  position */
#endif
);

extern int VGraphLinkNodes (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  a */,
    int			/*  b */
#endif
);

extern int VGraphUnlinkNodes (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  a */,
    int			/*  b */
#endif
);

extern VPointer VGraphFirstNode (
#if NeedFunctionPrototypes
    VGraph		/*  graph */
#endif
);

extern VPointer VGraphNextNode (
#if NeedFunctionPrototypes
    VGraph		/*  graph */
#endif
);

extern void VGraphClearVisit (
#if NeedFunctionPrototypes
    VGraph		/*  graph */
#endif
);

extern int VGraphResizeFields (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  nfields */
#endif
);

extern int VGraphNCycles (
#if NeedFunctionPrototypes
    VGraph		/*  graph */
#endif
);

extern void VGraphToggleNodesFrom (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  i */
#endif
);

extern void VDestroyNode (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  i */
#endif
);

extern void VGraphDestroyNodesFrom (
#if NeedFunctionPrototypes
    VGraph		/*  graph */,
    int			/*  i */
#endif
);

extern void VGraphClearHops (
#if NeedFunctionPrototypes
    VGraph		/*  graph */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_VGraph_h */
