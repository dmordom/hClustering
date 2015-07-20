//---------------------------------------------------------------------------
//
// Project: OpenWalnut ( http://www.openwalnut.org )
//
// Copyright 2009 OpenWalnut Community, BSV-Leipzig and CNCF-CBS
// For more information see http://www.openwalnut.org/copying
//
// This file is part of OpenWalnut.
//
// OpenWalnut is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// OpenWalnut is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenWalnut. If not, see <http://www.gnu.org/licenses/>.
//
// This file is also part of the
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno
//
//---------------------------------------------------------------------------

#ifndef WHTREEPROCESSER_H
#define WHTREEPROCESSER_H

// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <cstdlib>

// hClustering
#include "WHtree.h"

typedef enum
{
    HTPR_SIZERATIO, HTPR_JOINSIZE, HTPR_JOINLEVEL
} HTPROC_MODE;

typedef enum
{
    HTPR_C_CONSTANT, HTPR_C_LINEAR, HTPR_C_SQ
} HTPROC_COLLAPSE;

/**
 * class implements methods for dendrogram processing, such as pruning, node collapsing and decimation
 * it operates over the class WHtree
 *
 */
class WHtreeProcesser
{
public:
    //! Constructor
    explicit WHtreeProcesser( WHtree* const tree );

    //! Destructor
    ~WHtreeProcesser();

    // === PUBLIC MEMBER FUNCTIONS ===


    //! Prunes the tree according to the options specified
    //! \param condition value of the condition to choose when pruning (size of joining cluster, size ratio or joining level)
    //! \param safeSize maximum size that a cluster can have in order to be pruned, bigger clusters wont be prune regardless of the condition
    //! \param pruneType type of pruning method (max size, maxlevel, sizeratio)
    //! \return pair with the number of nodes and leaves eliminated
    std::pair< size_t, size_t > pruneTree( float condition, size_t safeSize, const HTPROC_MODE pruneType );

    //! Prunes the tree from randomly chosen leaves
    //! \param numberPruned vnumber of leaves to be pruned
    //! \param seed seed for the random number generator to choose the leaves to be pruned
    //! \return pair with the number of nodes and leaves eliminated
    std::pair< size_t, size_t > pruneRandom( const size_t numberPruned, unsigned int seed = 0 );

    //! Collapses the nodes of the tree that have branch lenghts equal or lower to the one indicated
    //! \param flatGap branch lenght limit to be collapsed
    //! \param distLevelLimit distance level under which no collapsing will take place
    //! \return number of nodes eliminated
    size_t collapseTree( const dist_t flatGap, const dist_t distLevelLimit, const bool keepBaseNodes );

    //! Collapses the nodes of the tree that have branch lenghts smoller than the distance level multiplied by the coefficient indicated
    //! \param coefficient coefficient for collapsing depending on distance level
    //! \return number of nodes eliminated
    size_t collapseTreeLinear( const dist_t coefficient, const bool keepBaseNodes );

    //! Collapses the nodes of the tree that have branch lenghts smoller than the squared distance level multiplied by the coefficient indicated
    //! \param coefficient coefficient for collapsing depending on distance level
    //! \return number of nodes eliminated
    size_t collapseTreeSquare( const dist_t coefficient, const bool keepBaseNodes );

    //! Collapses the nodes of the subtree that have branch lenghts equal or lower to the one indicated
    //! \param flatGap branch lenght limit to be collapsed
    //! \param distLevelLimit distance level under which no collapsing will take place
    //! \param root root node of the subtree
    //! \return number of nodes eliminated
    size_t collapseBranch( const dist_t flatGap, const dist_t distLevelLimit, size_t root, const bool keepBaseNodes );

    //! performs the tree smoothin
    //! \param safeSize minimum size of the clusters to be pruned in the first stage
    //! \param smoothSize size of the clusters looked for
    //! \param smoothGap maximum gap to be smoothed
    //! \return pair with the number of nodes and leaves eliminated
//    std::pair< size_t, size_t > smoothTree( const size_t safeSize, const size_t smoothSize, const dist_t smoothGap );

    //! Collapses all the nodes of the subtrees indicated
    //! \param selection vector with the root nodes of the subtree
    //! \return number of nodes eliminated
    size_t flattenSelection( const std::list< size_t > selection, const bool keepBaseNodes );
    size_t flattenSelection( const std::vector< size_t > &selection, const bool keepBaseNodes );
    size_t flattenSelection( const std::vector< nodeID_t > &selection, const bool keepBaseNodes );

    //! Prunes all of the subtrees indicated
    //! \param selection vector with the nodes to be pruned
    //! \return pair with the number of nodes and leaves eliminated
    std::pair< size_t, size_t > pruneSelection( const std::vector< size_t > &selection );
    std::pair< size_t, size_t > pruneSelection( const std::vector< nodeID_t > &selection );

    //! flags for pruning all of the subtrees indicated
    //! \param selection vector with the nodes to be pruned
    void flagSelection( const std::vector< size_t > &selection );
    void flagSelection( const std::vector< nodeID_t > &selection );

    void flagLeaves( const std::vector< size_t > &selection );


    //! Reduces the grid size by the ratio specified
    //! \param coarseRatio coarsing ratio
    void coarseTree( const unsigned int coarseRatio );

    //! converts the base nodes into leaves, pruning out all excess leaves and information
    //! \return new number of leaves of the tree
    size_t baseNodes2Leaves();




    void forceMonotonicityUp();
    void forceMonotonicityDown();
    void forceMonotonicity();

    size_t debinarize( const bool keepBaseNodes );


private:
    // === PRIVATE MEMBER DATA ===


    //! tree object
    WHtree &m_tree;

    // === PRIVATE MEMBER FUNCTIONS ===


    //! Collapses all the nodes of the subtree
    //! \param root root node of the subtree
    //! \return number of nodes eliminated
    size_t flattenBranch( size_t root, const bool keepBaseNodes );

    //! eliminates branchings in the tree, raising them to the parent level
    //! \param thisNodeID id of the sub-branch node where to start the flatteing
    //! \param coefficient coefficient for collapsing depending on distance level
    //! \param collapseMode collapsing mode
    void collapseNode( const size_t thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode = HTPR_C_CONSTANT );

    //! eliminates non monotonic steps in the tree, raising them to the parent level
    //! \param thisNodeID id of the sub-branch node where to start the flatteing
    //! \param coefficient coefficient for collapsing depending on distance level
    void collapseNode( const nodeID_t &thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode = HTPR_C_CONSTANT );
};

#endif  // WHTREEPROCESSER_H
