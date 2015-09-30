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
// - Moreno-Dominguez, D., Anwander, A., & Knoesche, T. R. (2014).
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


#ifndef WHTREE_H
#define WHTREE_H

// std library
#include <utility>
#include <vector>
#include <list>
#include <string>


#include <boost/filesystem.hpp>

// hClustering
#include "WHcoord.h"
#include "WHnode.h"
#include "WFileParser.h"


/**
 * this class implements a hierarchical tree and provides functions for creation, partitioning, selection and navigation
 */
class WHtree
{
public:
    /**
     * Constructor
     */
    explicit WHtree();

    /**
     * Constructor
     * \param filename string containing the name of the file with the tree data
     */
    explicit WHtree( std::string filename );

    /**
     * Constructor
     * \param treeName name of the tree
     * \param datasetGridInit type of grid coordinate space
     * \param datasetSizeInit size of the brain dataset
     * \param numStreamlinesInit number of streamlines per voxel on the tractograms that were processed to build this tree
     * \param logFactorInit logarithmic normalization factor that was applied on the tracts when building the tree
     * \param leavesInit vector of leaves
     * \param nodesInit vector of nodes
     * \param trackidInit vector of seed track ids
     * \param coordInit vector of seed voxel coordinates
     * \param discardInit vector of discarded seed voxel coordinates
     * \param cpccInit cpcc value of the tree
     */
    explicit WHtree( std::string treeName, HC_GRID datasetGridInit, WHcoord datasetSizeInit, size_t numStreamlinesInit, float logFactorInit,
                     std::vector<WHnode> leavesInit, std::vector<WHnode> nodesInit, std::vector<size_t> trackidInit, std::vector<WHcoord> coordInit,
                     std::list<WHcoord> discardInit, float cpccInit = 0 );

    /**
     * Constructor
     * \param object tree object
     */
    explicit WHtree( const WHtree &object );

    /**
     * Destructor
     */
    ~WHtree();



    // === PUBLIC MEMBER FUNCTIONS ===



    /**
     * Returns the tree name
     * \return value containing the tree name
     */
    std::string getName() const;

    /**
     * Returns the loading status of the tree
     * \return load indicator
     */
    bool isLoaded() const;

    /**
     * Returns the number of leaves of the tree
     * \return number of leaves
     */
    size_t getNumLeaves() const;

    /**
     * Returns the number of nodes of the tree
     * \return number of leaves
     */
    size_t getNumNodes() const;

    /**
     * Returns the leaves of the tree
     * \return leaves vector
     */
    std::vector<WHnode> getLeaves() const;

    /**
     * Returns the  nodes of the tree
     * \return nodes vector
     */
    std::vector<WHnode> getNodes() const;

    /**
     * Returns the number of leaves of the tree
     * \return number of leaves
     */
    size_t getNumDiscarded() const;

    /**
     * Returns the dataset size
     * \return coordinate containing the dataset dimensions
     */
    WHcoord getDataSize() const;

    /**
     * Returns the dataset grid type
     * \return enum containing the grid type
     */
    HC_GRID getDataGrid() const;

    /**
     * Returns the cpcc value of the tree
     * \return value containing the cpcc
     */
    float getCpcc() const;

    /**
     * Returns the roi coordinate vector
     * \return vector containing the roi coordinates
     */
    std::vector<WHcoord> getRoi() const;

    /**
     * Returns the tractogram IDs vector
     * \return vector containing the tractogram ids per leaf
     */
    std::vector<size_t> getTrackids() const;

    /**
     * Returns the discarded list
     * \return list containing the discarded coordinates
     */
    std::list<WHcoord> getDiscarded() const;

    /**
     * Returns the cutting values for the selected partitions
     * \return vector containing the values for the selected partitions
     */
    std::vector<float> getSelectedValues() const;

    /**
     * Returns the cutting values for the selected partition
     * \param index index of the desired partition
     * \return cutting values for the selected partition
     */
    float getSelectedValues( size_t index ) const;

    /**
     * Returns the set of the selected partitions
     * \return vector containing the selected partitions
     */
    std::vector< std::vector< size_t > > getSelectedPartitions() const;

    /**
     * Returns the indicated selected partition from the saved set
     * \param index index of the desired partition
     * \return vector containing the selected partition
     */
    std::vector< size_t > getSelectedPartitions( size_t index ) const;

    /**
     * Returns the set of the selected partition colors
     * \return vector containing the selected partitions colors
     */
    std::vector< std::vector< WHcoord > > getSelectedColors() const;

    /**
     * Returns the selected partition colors from the saved set
     * \param index index of the desired partition
     * \return vector containing the selected partition colors
     */
    std::vector< WHcoord > getSelectedColors( size_t index ) const;

    /**
     * writes out a small report status of the tree
     * \param longMsg if set a longer two-line message will be output
     * \return string with the status message
     */
    std::string getReport( bool longMsg = true ) const;

    /**
     * checks tree integrity
     * \return tree integrity indicator
     */
    bool check() const;

    /**
     * returns a const pointer to the node indicated by the ID
     * \param thisNode full id of the required node
     * \return const reference to desired node/leaf object
     */
    const WHnode& getNode( const nodeID_t &thisNode ) const;

    /**
     * returns a const pointer to the node indicated by the ID
     * \param thisNode id of the required node
     * \return const reference to desired node object
     */
    const WHnode& getNode( const size_t thisNode ) const;

    /**
     * returns a const pointer to the leaf indicated by the ID
     * \param thisLeaf id of the required leaf
     * \return const reference to desired leaf object
     */
    const WHnode& getLeaf( const size_t thisLeaf ) const;

    /**
     * returns a const reference to the root node
     * \return a const reference to the root node

     */
    const WHnode& getRoot() const;

    /**
     * returns the leaf ID of a coordinate in the tree
     * \param thisCoord coordinate of the required leaf
     * \return ID of the leaf
     */
    size_t getLeafID( const WHcoord &thisCoord ) const;

    /**
     * returns the corresponding track ID to a leaf ID on the tree
     * \param leafID input leaf ID
     * \return ID of seed tractogram corresponding to that leaf
     */
    size_t getTrackID( const size_t &leafID ) const;

    /**
     * Returns a vector with all the leaf IDs contained in that cluster
     * \param nodeID id of the selected node
     * \return vector with all the leaf IDs contained in that cluster
     */
    std::vector<size_t> getLeaves4node( const size_t nodeID ) const;

    /**
     * Returns a vector with all the leaf IDs contained in that cluster
     * \param nodeID full id of the selected node
     * \return vector with all the leaf IDs contained in that cluster
     */
    std::vector<size_t> getLeaves4node( const nodeID_t &nodeID ) const;

    /**
     * Returns a vector with all the node IDs contained in that branch
     * \param nodeID id of the selected node
     * \return vector with all the node IDs contained in that brnach
     */
    std::vector<size_t> getBranchNodes( const size_t nodeID ) const;

    /**
     * the coordinates of a particular leaf
     * \param leafID id of the selected leaf
     * \return coordinate object with the coordinates of that leaf
     */
    WHcoord getCoordinate4leaf( const size_t leafID ) const;

    /**
     * Returns a vector with all the seed voxel coordinates contained in that cluster
     * \param nodeID id of the selected node
     * \return vector with all the seed voxel coordinates contained in that cluster
     */
    std::vector<WHcoord> getCoordinates4node( const size_t nodeID ) const;

    /**
     * Returns a vector with all the seed voxel coordinates contained in that cluster
     * \param nodeID full id of the selected node
     * \return vector with all the seed voxel coordinates contained in that cluster
     */
    std::vector<WHcoord> getCoordinates4node( const nodeID_t &nodeID ) const;

    /**
     * Returns the mean coordinate of all the seed voxel coordinates contained in that cluster
     * \param nodeID id of the selected node
     * \return coordinate object with the mean coordinates of that node
     */
    WHcoord getMeanCoordinate4node( const size_t nodeID ) const;

    /**
     * Returns the mean coordinate of all the seed voxel coordinates contained in that cluster
     * \param nodeID full id of the selected node
     * \return coordinate object with the mean coordinates of that node
     */
    WHcoord getMeanCoordinate4node( const nodeID_t &nodeID ) const;

    /**
     * retrieves the the first common ancestor node of two given nodes
     * \param nodeID1 first node id
     * \param nodeID2 second node id
     * \return ID of ancestor
     */
    size_t getCommonAncestor( const size_t nodeID1, const size_t nodeID2 ) const;

    /**
     * retrieves the the first common ancestor node of two given nodes
     * \param nodeID1 first node/leaf full id
     * \param nodeID2 second node/leaf full id
     * \return ID of ancestor
     */
    nodeID_t getCommonAncestor( const nodeID_t &nodeID1, const nodeID_t &nodeID2 ) const;

    /**
     * node route to root node
     * \param nodeID full id of the selected node
     * \return nodeId vector of nodes in the route
     */
    std::vector<nodeID_t> getRoute2Root( const nodeID_t &nodeID ) const;

    /**
     * checks the joining order of a given triplet
     * \param nodeIDa first node/leaf of the triplet
     * \param nodeIDb second node/leaf of the triplet
     * \param nodeIDc second node/leaf of the triplet
     * \return joining order: 0="non resolved" (a,b,c join at the same time), 1="ab before c", 2="ac before b", 3="bc before a"
     */
    unsigned int getTripletOrder( const nodeID_t &nodeIDa, const nodeID_t &nodeIDb, const nodeID_t &nodeIDc ) const;

    /**
     * returns a vector with all the nodes that have direct link to leaves fro the indicated subtree
     * \param root id of the selected subtree root node
     * \return nodeId vector of base nodes
     */
    std::vector<size_t> getBaseNodes( const size_t root ) const;

    /**
     * returns a vector with all the nodes that have direct link to leaves from the indicated subtree
     * \param root full id of the selected subtree root node
     * \return nodeId vector of base nodes
     */
    std::vector<nodeID_t> getBaseNodes( const nodeID_t &root ) const;

    /**
     * returns a vector with all the nodes that have direct link to leaves from the main root
     * \return nodeId vector of base nodes
     */
    std::vector<size_t> getRootBaseNodes() const;

    /**
     * tests if all base nodes have only leaf children
     * \return true if all base nodes have only leaf children
     */
    bool testRootBaseNodes() const;

    /**
     * computes the cophenetic distance between 2 nodes
     * \param nodeID1 first node to compute distance from
     * \param nodeID2 second node to compute distance from
     * \return cophenetic distance
     */
    dist_t getDistance( const size_t nodeID1, const size_t nodeID2 ) const;

    /**
     * computes the cophenetic distance between 2 leaf/nodes
     * \param nodeID1 first node/leaf to compute distance from
     * \param nodeID2 second node/leaf to compute distance from
     * \return cophenetic distance
     */
    dist_t getDistance( const nodeID_t &nodeID1, const nodeID_t &nodeID2 ) const;

    /**
     * computes the cophenetic distance between 2 leaves
     * \param coord1 first coordinate to compute distance from
     * \param coord2 second ncoordinate to compute distance from
     * \return cophenetic distance
     */
    dist_t getDistance( const WHcoord &coord1, const WHcoord &coord2 ) const;

    /**
     * computes the cophenetic distance between 2 leaves
     * \param leafID1 first leaf to compute distance from
     * \param leafID2 second leaf to compute distance from
     * \return cophenetic distance
     */
    dist_t getLeafDistance( const size_t leafID1, const size_t leafID2 ) const;



    /**
     * Reorders the nodes in the vector according to their size
     * \param nodeVector vector of nodes to be reordered
     */
    void sortBySize( std::vector<size_t>* const nodeVector ) const;

    /**
     * Reorders the nodes in the vector according to their size
     * \param nodeVector vector of nodes to be reordered
     */
    void sortBySize( std::vector<nodeID_t>* const nodeVector ) const;

    /**
     * Reorders the nodes in the list according to their size
     * \param nodeList list of nodes to be reordered
     */
    void sortBySize( std::list<size_t>* const nodeList ) const;

    /**
     * Reorders the nodes in the list according to their size
     * \param nodeList vector of nodes to be reordered
     */
    void sortBySize( std::list<nodeID_t>* const nodeList ) const;


    /**
     * Reorders the nodes in the vector according to their hierarchical level
     * \param nodeVector vector of nodes to be reordered
     */
    void sortByHLevel( std::vector<size_t>* const nodeVector ) const;

    /**
     * Reorders the nodes in the list according to their hierarchical level
     * \param nodeList list of nodes to be reordered
     */
    void sortByHLevel( std::list<size_t>* const nodeList ) const;


    /**
     * Reorders the nodes in the vector according to their hierarchical level
     * \param nodeVector vector of nodes to be reordered
     */
    void sortByHLevel( std::vector<nodeID_t>* const nodeVector ) const;

    /**
     * Reorders the nodes in the list according to their hierarchical level
     * \param nodeList vector of nodes to be reordered
     */
    void sortByHLevel( std::list<nodeID_t>* const nodeList ) const;


    /**
     * stores the contained leaves for each node
     */
    void loadContainedLeaves();

    /**
     * frees the memory of m_containedLeaves
     */
    void clearContainedLeaves();

    /**
     * converts all the coordinates to the grid format indicated
     * \param newGrid coordinate space in which to convert
     * \return flag indicating if the change was performed
     */
    bool convert2grid( const HC_GRID newGrid );

    /**
     * Reads the the tree data from a file and creates the structure, containing it in the leaves, nodes and coordinates vectors
     * \param filename string containing std::vector<hNode> m_leavesthe name of the file with the tree data
     * \return success indicator
     */
    bool readTree( const std::string &filename );

    /**
     * writes the current tree data
     * \param filename string containing the name of the file to be written
     * \param isNifti a flag indicating whether the tree should be written in nifti coordinates
     * \return success indicator
     */
    bool writeTree( const std::string &filename, const bool isNifti = true ) const;

    /**
     * writes all redundant tree data for debugging purposes
     * \param filename string containing the name of the file to be written
     * \return success indicator
     */
    bool writeTreeDebug( const std::string &filename ) const;

    /**
     * writes dta for working with the old walnut module
     * \param filename string containing the name of the file to be written
     * \return success indicator
     */
    bool writeTreeOldWalnut( const std::string &filename ) const;

    /**
     * writes only join data without coordinates or reading labels
     * \param filename string containing the name of the file to be written
     * \return success indicator
     */
    bool writeTreeSimple( const std::string &filename ) const;

    /**
     * insert a set of partitions into the tree data
     * \param selectedPartitions a set of partitions
     * \param selectedValues the quality values for those partitions
     * \param selectedColors the colors for each cluster of each partitions
     */
    void insertPartitions( const std::vector<std::vector<size_t> > &selectedPartitions,  const std::vector< float > &selectedValues,
                           const std::vector<std::vector<WHcoord> > &selectedColors = std::vector<std::vector<WHcoord> >() );

    /**
     * insert a set of partitions colors into the tree data in case a set of saved partitions is already present
     * \param selectedColors the set of new colors for the included partitions
     */
    void insertPartColors( const std::vector<std::vector<WHcoord> > &selectedColors );

    /**
     * delete all data related to saved partitions
     */
    void clearPartitions();

    /**
     * delete only color data of the saved partitions
     */
    void clearPartColors();

    /**
     * gets all possble sets of sub-partitions in a tree given an initial partition and a depth level (where leaves/nodes are identified with boolean-integer pair)
     \param thisPartition the starting partition
     \param depthLevel the number of levels down to look for sub-partitions
     \param excludeLeaves if set (1 ), the the algorithm will not further subdivide a base-node (where all its children are single leaves)
     \retval partitionSet the vector of returned sub-partitions
     \return a set of integers identifying each of the clusters from the resulting sub-partitions with its parent cluster from the starting partition
     */
    std::vector< std::vector< unsigned int > > getBranching( const std::vector < nodeID_t > &thisPartition,
                                                             size_t depthLevel,
                                                             std::vector< std::vector < nodeID_t > > *partitionSet,
                                                             const bool excludeLeaves );

    /**
     * gets all possble sets of sub-partitions in a tree given an initial partition and a depth level (where only non-leaf nodes are considered and are identified by a single integer)
     \param thisPartition the starting partition
     \param depthLevel the number of levels down to look for sub-partitions
     \retval partitionSet the vector of returned sub-partitions
     \return a set of integers identifying each of the clusters from the resulting sub-partitions with its parent cluster from the starting partition
     */
    std::vector< std::vector< unsigned int > > getBranching( const std::vector < size_t > &thisPartition,
                                                             size_t depthLevel,
                                                             std::vector< std::vector < size_t > > *partitionSet );



    friend class treeManager;
    friend class WHtreeProcesser;
    friend class WHtreePartition;

    friend class CnbTreeBuilder;
    friend class graphTreeBuilder;
    friend class randCnbTreeBuilder;
    friend class treeComparer;
    friend class pruneTree;
    friend class image2treeBuilder;


protected:
private:
    // === PRIVATE MEMBER DATA ===


    //! Stores all the leaves of the tree (represents a single voxel, for several purposes a leaf also counts as a node
    bool m_loadStatus;

    //! Stores the size of the dataset this tree was built from
    WHcoord m_datasetSize;

    //! Stores the the type of coordinate grid of the dataset
    HC_GRID m_datasetGrid;

    //! number of streamlines generated form each seed voxel in the tractograms used to build the tree
    size_t m_numStreamlines;

    //! logarithmic normalization factor used in the tractograms when building the tree
    float m_logFactor;

    //! Stores the cpcc value of the tree
    float m_cpcc;

    //! Stores the name of the tree file
    std::string m_treeName;

    //! Stores all the leaves of the tree (represents a single voxel, for several purposes a leaf also counts as a node
    std::vector<WHnode> m_leaves;

    //! Stores all the nodes (non-leaves) of the tree
    std::vector<WHnode> m_nodes;

    //! Stores the coordinates of the leaf voxels
    std::vector<WHcoord> m_coordinates;

    //! Stores the ids of the seed tracts correesponding to each leaf
    std::vector<size_t> m_trackids;

    //! Stores the coordinates of the voxels that were discarded during the process of building the tree
    std::list<WHcoord> m_discarded;

    //! Stores the contained leafs of each node
    std::vector<std::vector<size_t> > m_containedLeaves;

    //! Stores a set of selected partitions
    std::vector< std::vector<size_t> > m_selectedPartitions;

    //! Stores the quality values of the selected partitions
    std::vector< float > m_selectedValues;

    //! Stores the colors of the selected partitions
    std::vector< std::vector< WHcoord > > m_selectedColors;

    // === PRIVATE MEMBER FUNCTIONS ===


    /**
     * returns a pointer to the node indicated by the ID
     * \param thisNode full id of the required node
     * \return pointer to desired node/leaf object
     */
    WHnode* fetchNode( const nodeID_t &thisNode );

    /**
     * returns a pointer to the node indicated by the ID
     * \param thisNode id of the required node
     * \return pointer to desired node/leaf object
     */
    WHnode* fetchNode( const size_t thisNode );

    /**
     * returns a pointer to the leaf indicated by the ID
     * \param thisLeaf id of the required leaf
     * \return pointer to desired leaf object
     */
    WHnode* fetchLeaf( const size_t thisLeaf );

    /**
     * returns a pointer to the root node
     * \return a pointer to root node object
     */
    WHnode* fetchRoot();

    /**
     * cleans leaves/nodes set to be pruned
     * \param outLookup pointer to a vector where to store the lookup table of changes done to nodes
     * \return pair with the number of leaves and nodes pruned
     */
    std::pair<size_t, size_t> cleanup( std::vector<size_t> *outLookup = 0 );

    /**
     * merges nodes joining binary at the same level into non-binary structures
     * \param keepBaseNodes if set basenodes will not be marged
     * \return number of nodes collapsed
     */
    size_t debinarize( const bool keepBaseNodes = false );

    /**
     * eliminates non monotonic steps in the tree while keeping as much information as possible.
     * it searches bottom-up for non-monotonic steps. When detected a corrected value for the parent node is computed
     * this value consists of a average between the non-monotonic children levels weighted by their sizes in number of leaves and the parent level weighted by the sizes of the remaining monotonic children
     * in orther to avoid infinite loops an error tolerance is included. default value is 1*E-5
     * \param errorMult multiplier to increase the base error tolerance (default value is 1 makign the tolerance 1*E-5, maximum value is 100 making the tolerance 1*E-3)
     */
    void forceMonotonicity( double errorMult = 1 );

    /**
     * eliminates non monotonic steps in the tree bottom-up, lowering the level of children involved in the non-monotonic steps to the parent level
     */
    void forceMonotonicityUp();

    /**
     * eliminates non monotonic steps in the tree top-down, raising the level of parents involved in the non-monotonic steps to the children level
     */
    void forceMonotonicityDown();

    //! implements a compare operator for nodes inside a tree structure, depending on the Size
    struct compSize
    {
        const WHtree* const m_tree; //!< stores pointer to tree we work on

        /**
         * structure constructor, stores the tree pointer into the data member
         * \param tree the tree pointer
         */
        explicit compSize( const WHtree* const tree ): m_tree( tree )
        {
        }

        /**
         * the comparison operator result when using nodeID_t types
         * \param lhs the left hand side object
         * \param rhs the right hand side object
         * \return the result of the comparison operation
         */
        bool operator()( const nodeID_t &lhs, const nodeID_t &rhs ) const
        {
            return ( ( ( m_tree->getNode( lhs ) ).getSize() ) < ( ( m_tree->getNode( rhs ) ).getSize() ) );
        }

        /**
         * the comparison operator result when using size_t types
         * \param lhs the left hand side object
         * \param rhs the right hand side object
         * \return the result of the comparison operation
         */
        bool operator()( const size_t lhs, const size_t rhs ) const
        {
            return ( ( ( m_tree->getNode( lhs ) ).getSize() ) < ( ( m_tree->getNode( rhs ) ).getSize() ) );
        }
    };

    //! implements a compare operator for nodes, depending on the hierarchical level
    struct compHLevel
    {
        const WHtree* const m_tree; //!< stores pointer to tree we work on

        /**
         * structure constructor, stores the tree pointer into the data member
         * \param tree the tree pointer
         */
        explicit compHLevel( const WHtree* const tree ): m_tree( tree )
        {
        }

        /**
         * the comparison operator result when using nodeID_t types
         * \param lhs the left hand side object
         * \param rhs the right hand side object
         * \return the result of the comparison operation
         */
        bool operator()( const nodeID_t &lhs, const nodeID_t &rhs ) const
        {
            if( ( ( m_tree->getNode( lhs ) ).getHLevel() ) == ( ( m_tree->getNode( rhs ) ).getHLevel() ) )
            {
                return( lhs < rhs );
            }
            else
            {
                return( ( ( m_tree->getNode( lhs ) ).getHLevel() ) < ( ( m_tree->getNode( rhs ) ).getHLevel() ) );
            }
        }

        /**
         * the comparison operator result when using size_t types
         * \param lhs the left hand side object
         * \param rhs the right hand side object
         * \return the result of the comparison operation
         */
        bool operator()( const size_t lhs, const size_t rhs ) const
        {
            if( ( ( m_tree->getNode( lhs ) ).getHLevel() ) == ( ( m_tree->getNode( rhs ) ).getHLevel() ) )
            {
                return( lhs < rhs );
            }
            else
            {
                return( ( ( m_tree->getNode( lhs ) ).getHLevel() ) < ( ( m_tree->getNode( rhs ) ).getHLevel() ) );
            }
        }
    };
};

// === IN-LINE PUBLIC MEMBERS ===

inline std::string WHtree::getName() const
{
    return m_treeName;
}

inline bool WHtree::isLoaded() const
{
    return m_loadStatus;
}

inline size_t WHtree::getNumLeaves() const
{
    return m_leaves.size();
}

inline size_t WHtree::getNumNodes() const
{
    return m_nodes.size();
}

inline std::vector<WHnode> WHtree::getLeaves() const
{
    return m_leaves;
}

inline std::vector<WHnode> WHtree::getNodes() const
{
    return m_nodes;
}
inline size_t WHtree::getNumDiscarded() const
{
    return m_discarded.size();
}

inline WHcoord WHtree::getDataSize() const
{
    return m_datasetSize;
}

inline HC_GRID WHtree::getDataGrid() const
{
    return m_datasetGrid;
}

inline float WHtree::getCpcc() const
{
    return m_cpcc;
}

inline std::vector<WHcoord> WHtree::getRoi() const
{
    return m_coordinates;
}

inline std::vector<size_t> WHtree::getTrackids() const
{
    return m_trackids;
}

inline std::list<WHcoord> WHtree::getDiscarded() const
{
    return m_discarded;
}

inline std::vector<float> WHtree::getSelectedValues() const
{
    return m_selectedValues;
}

inline float WHtree::getSelectedValues( size_t index ) const
{
    if( index >= m_selectedValues.size() )
    {
        return -1;
    }
    else
    {
        return m_selectedValues[ index ];
    }
}

inline std::vector< std::vector< size_t > > WHtree::getSelectedPartitions() const
{
    return m_selectedPartitions;
}

inline std::vector< size_t > WHtree::getSelectedPartitions( size_t index ) const
{
    if( index >= m_selectedPartitions.size() )
    {
        return std::vector< size_t >();
    }
    else
    {
        return m_selectedPartitions[ index ];
    }
}

inline std::vector< std::vector< WHcoord > > WHtree::getSelectedColors() const
{
    return m_selectedColors;
}

inline std::vector< WHcoord > WHtree::getSelectedColors( size_t index ) const
{
    if( index >= m_selectedColors.size() )
    {
        return std::vector< WHcoord >();
    }
    else
    {
        return m_selectedColors[ index ];
    }
}

#endif  // WHTREE_H
