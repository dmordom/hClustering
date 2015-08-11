//---------------------------------------------------------------------------
//
// Project: OpenWalnut ( http://www.openwalnut.org )
//
// Copyright 2009 OpenWalnut Community, BSV@Uni-Leipzig and CNCF@MPI-CBS
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
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
// This file is also part of OpenWalnut ( http://www.openwalnut.org ).
//
// For more reference on the underlying algorithm and research they have been used for refer to:
// - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
//   A hierarchical method for whole-brain connectivity-based parcellation.
//   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528
// - Moreno-Dominguez, D. (2014).
//   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
//   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
//   ISBN 978-3-941504-45-5
//
// hClustering is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
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

// type of tree pruning
typedef enum
{
    HTPR_SIZERATIO, // based on size ratio of joining nodes
    HTPR_JOINSIZE, // based on absolute size values of joinign nodes
    HTPR_JOINLEVEL // based on joining distance level
} HTPROC_MODE;

// type of of node collapsing
typedef enum
{
    HTPR_C_CONSTANT, // collapse all node pairs with distances lower than a constant distance difference
    HTPR_C_LINEAR, // normalize node distance difference linearly depending on their distance level
    HTPR_C_SQ // normalize node distance difference squarely depending on their distance level
} HTPROC_COLLAPSE;

/**
 * this class implements methods for dendrogram processing, such as pruning, node collapsing and decimation
 * it operates over the class WHtree.
 * As a result of its operations the tree structure is changed and previous structure cannot be recovered unless a class copy is previously made.
 */
class WHtreeProcesser
{
public:
    /**
     * Constructor
     * \param tree a pointer to the tree to be processed
     */
    explicit WHtreeProcesser( WHtree* const tree );

    //! Destructor
    ~WHtreeProcesser();

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Prunes the tree according to the options specified
     * \param condition value of the condition to choose when pruning (size of joining cluster, size ratio or joining level)
     * \param safeSize maximum size that a cluster can have in order to be pruned, bigger clusters wont be prune regardless of the condition
     * \param pruneType type of pruning method (max size, maxlevel, sizeratio)
     * \return pair with the number of nodes and leaves eliminated
     */
    std::pair< size_t, size_t > pruneTree( float condition, size_t safeSize, const HTPROC_MODE pruneType );

    /**
     * Prunes the tree from randomly chosen leaves
     * \param numberPruned number of leaves to be pruned
     * \param seed seed for the random number generator to choose the leaves to be pruned
     * \return pair with the number of nodes and leaves eliminated
     */
    std::pair< size_t, size_t > pruneRandom( const size_t numberPruned, unsigned int seed = 0 );

    /**
     * Collapses the nodes of the tree that have branch lenghts equal or lower to the one indicated
     * \param flatGap branch lenght limit to be collapsed
     * \param distLevelLimit distance level under which no collapsing will take place
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t collapseTree( const dist_t flatGap, const dist_t distLevelLimit, const bool keepBaseNodes );

    /**
     * Collapses the nodes of the tree that have branch lenghts smoller than the distance level multiplied by the coefficient indicated
     * \param coefficient coefficient for collapsing depending on distance level
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t collapseTreeLinear( const dist_t coefficient, const bool keepBaseNodes );

    /**
     * Collapses the nodes of the tree that have branch lenghts smoller than the squared distance level multiplied by the coefficient indicated
     * \param coefficient coefficient for collapsing depending on distance level
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t collapseTreeSquare( const dist_t coefficient, const bool keepBaseNodes );

    /**
     * Collapses the nodes of the subtree that have branch lenghts equal or lower to the one indicated
     * \param flatGap branch lenght limit to be collapsed
     * \param distLevelLimit distance level under which no collapsing will take place
     * \param root root node of the subtree
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t collapseBranch( const dist_t flatGap, const dist_t distLevelLimit, size_t root, const bool keepBaseNodes );

    /**
     * Collapses all the nodes of the subtrees indicated (using only non-leaf node ID)
     * \param selection list with the root nodes of the subtree
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t flattenSelection( const std::list< size_t > selection, const bool keepBaseNodes );

    /**
     * Collapses all the nodes of the subtrees indicated (using only non-leaf node ID)
     * \param selection vector with the root nodes of the subtree
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t flattenSelection( const std::vector< size_t > &selection, const bool keepBaseNodes );

    /**
     * Collapses all the nodes of the subtrees indicated (using full boolean-integer node ID)
     * \param selection vector with the root nodes of the subtree
     * \param keepBaseNodes flag to avoid collapsing base nodes (those that have only single-leaf children)
     * \return number of nodes eliminated
     */
    size_t flattenSelection( const std::vector< nodeID_t > &selection, const bool keepBaseNodes );

    /**
     * Prunes all of the subtrees indicated (using only non-leaf node ID)
     * \param selection vector with the nodes to be pruned
     * \return pair with the number of nodes and leaves eliminated
     */
    std::pair< size_t, size_t > pruneSelection( const std::vector< size_t > &selection );

    /**
     * Prunes all of the subtrees indicated (using full boolean-integer node ID)
     * \param selection vector with the nodes to be pruned
     * \return pair with the number of nodes and leaves eliminated
     */
    std::pair< size_t, size_t > pruneSelection( const std::vector< nodeID_t > &selection );

    /**
     * set flags for pruning all of the subtrees indicated (using only non-leaf node ID)
     * \param selection vector with the nodes to be pruned
     */
    void flagSelection( const std::vector< size_t > &selection );

    /**
     * set flags for pruning all of the subtrees indicated (using full boolean-integer node ID)
     * \param selection vector with the nodes to be pruned
     */
    void flagSelection( const std::vector< nodeID_t > &selection );

    /**
     * set flags for pruning all of the leaves indicated
     * \param selection vector with the leaves to be pruned
     */
    void flagLeaves( const std::vector< size_t > &selection );

    /**
     * Reduces the tree grid size by the ratio specified (as if it had been obtained from an image of lower resolution)
     * \param coarseRatio coarsing ratio
     */
    void coarseTree( const unsigned int coarseRatio );

    /**
     * converts the base nodes into leaves, pruning out all excess leaves and information
     * \return new number of leaves of the tree
     */
    size_t baseNodes2Leaves();

    /**
     * calls on private WHtree member to eliminate non monotonic steps in the tree while keeping as much information as possible.
     * it searches bottom-up for non-monotonic steps. When detected a corrected value for the parent node is computed
     * this value consists of a average between the non-monotonic children levels weighted by their sizes in number of leaves and the parent level weighted by the sizes of the remaining monotonic children
     * in orther to avoid infinite loops an error tolerance is included. default value is 1*E-5
     * \param errorMult multiplier to increase the base error tolerance (default value is 1 makign the tolerance 1*E-5, maximum value is 100 making the tolerance 1*E-3)
     */
    void forceMonotonicity( double errorMult = 1 );

    /**
     * calls on private WHtree member to eliminate non monotonic steps in the tree bottom-up, lowering the level of children involved in the non-monotonic steps to the parent level
     */
    void forceMonotonicityUp();

    /**
     * calls on private WHtree member to eliminate non monotonic steps in the tree top-down, raising the level of parents involved in the non-monotonic steps to the children level
     */
    void forceMonotonicityDown();

    /**
     * calls on private WHtree member to merge nodes joining binary at the same level into non-binary structures
     * \param keepBaseNodes if set base nodes ( nodes with only single leaves as children) will not be merged
     * \return number of nodes collapsed
     */
    size_t debinarize( const bool keepBaseNodes );


    /**
     * Writes a file with the base nodes (meta-leaves) obtained at the end of the homgeneous merging initial stage
     * \param baseNodes a vector with the base node IDs
     * \param filename the path to the file where the information will be written to
     */
    void writeBases( const std::vector< size_t >& baseNodes, const std::string& filename ) const;

    //    /**
    //     * performs the tree smoothin
    //     * \param safeSize minimum size of the clusters to be pruned in the first stage
    //     * \param smoothSize size of the clusters looked for
    //     * \param smoothGap maximum gap to be smoothed
    //     * \return pair with the number of nodes and leaves eliminated
    //     */
    //    std::pair< size_t, size_t > smoothTree( const size_t safeSize, const size_t smoothSize, const dist_t smoothGap );

private:
    // === PRIVATE MEMBER DATA ===

    WHtree &m_tree; //!< tree object

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Collapses all the nodes of the subtree
     * \param root root node of the subtree
     * \param keepBaseNodes if set base nodes ( nodes with only single leaves as children) will not be merged
     * \return number of nodes eliminated
     */
    size_t flattenBranch( size_t root, const bool keepBaseNodes );

    /**
     * eliminates branchings in the tree, raising them to the parent level
     * \param thisNodeID id of the sub-branch node where to start the flatteing
     * \param coefficient coefficient for collapsing depending on distance level
     * \param collapseMode collapsing mode
     */
    void collapseNode( const size_t thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode = HTPR_C_CONSTANT );

    /**
     * eliminates non monotonic steps in the tree, raising them to the parent level
     * \param thisNodeID id of the sub-branch node where to start the flatteing
     * \param coefficient coefficient for collapsing depending on distance level
     * \param collapseMode collapsing mode
     */
    void collapseNode( const nodeID_t &thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode = HTPR_C_CONSTANT );
};

#endif  // WHTREEPROCESSER_H
