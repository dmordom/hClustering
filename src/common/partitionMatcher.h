//---------------------------------------------------------------------------
//
// Project: hCLustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
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
// hCLustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------


#ifndef TREECOMPARER_H
#define TREECOMPARER_H


// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <sstream>


// hClustering
#include "WHcoord.h"
#include "WHtree.h"


#define TREE1 true
#define TREE2 false

/**
 * This class implements methods in order to read partitions from a reference tree and find the best matching corresponding partitions in a target tree
 * The reference and found target partitions will be color-matched as best as possible.
 * There is also the possibility of only color-matching predefined partitions of the target tree to predefined partitions of the reference tree.
 */
class partitionMatcher
{
public:
    /**
     * Constructior
     * \param refTree a pointer to the reference tree
     * \param targetTree a pointer to the target tree
     * \param matchFilename the filepath to the trees base-matching file
     * \param verbose the verbose output flag
     */
    partitionMatcher( WHtree* const refTree, WHtree* const targetTree, std::string matchFilename, bool verbose );

    //! Destructor.
    ~partitionMatcher();

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Analyzes the number and size of the base nodes of each tree and returns a report message with the results
     * \return the report message
     */
    std::string reportBaseNodes() const;

    /**
     * Find matching partitions from the reference tree in the target tree
     * \param lambda the weighting factor for the signature-matrix matching, if less than 0 bidirectional seed voxel overlap matchign will be used
     * \param levleDepth the hierarchical search depth for best matchign partitions, if not set or set to 0 the search depth will be adaptively decreased as the number of clusters in the search increases in order to optimize computing time
     */
    void findMatchingPartitions( const float lambda, const size_t levleDepth = 0);

    /**
     * Match colors across the partitions present in both trees. Colors of reference tree will be matched onto target tree.
     * \param exclusive if set, the unmatched clusters of the target tree will be set to white color, otherwise, they will be given random coloring.
     * \return reference colors change flag. If true, some clusters of the reference tree were recolored to improve color matching (in cases of multiple-to-1 match)
     */
    bool matchColors( bool exclusive = false );

private:
    // === PRIVATE DATA MEMBERS ===

    bool                  m_verbose;             //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    WHtree&               m_refTree;             //!< A reference to the reference tree
    WHtree&               m_targetTree;          //!< A reference to the target tree
    std::ofstream*        m_logfile;             //!< A pointer to the output log file stream
    std::vector< size_t > m_refBaseNodes;        //!< The native order reference tree base nodes
    std::vector< size_t > m_targetBaseNodes;     //!< The native order target tree base nodes
    std::vector< size_t > m_refMatchedBases;     //!< The matched order reference tree base nodes
    std::vector< size_t > m_targetMatchedBases;  //!< The matched order target tree base nodes
    std::vector< size_t > m_fullCorrespondence;  //!< The correspondence base nodes lookup table across trees
    std::vector< std::vector< size_t > > m_refMatchedBases4Node;    //!< A list of corresponding target tree nodes for each reference tree node
    std::vector< std::vector< size_t > > m_targetMatchedBases4Node; //!< A list of corresponding reference tree nodes for each target tree node

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Loads and tests the integrity of the base nodes for both trees
     */
    void testBaseNodes();

    /**
     * Loads the base-nodes correspondence lookup table from file
     * \param matchFilename the path to the correspondence table file
     */
    void loadCorrespondence( const std::string &matchFilename );

    /**
     * Retrieves the relative base-node ID of a base node (its position in the base-node vector) from the absolute node ID on the tree structure
     * \param absoluteID the ID of the node in the tree structure
     * \param baseNodes the base node vector
     * \return the relative ID of the node in the base node vector
     */
    size_t findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t >& baseNodes ) const;

    /**
     * Returns the signature matrix of a partition in respect to one of the trees.
     * The signature matrix defines a value for each pair of base-nodes: 1 the base nodes are found in the same cluster, 0 if otherwise
     * (the more correlated two signature matrix of two partitions across trees are, the best matching those partitions are)
     * \param partition the vector of node IDs that define the partition
     * \param forRefTree a boolean value defining to which tree whe are referring, if true=refTree, if false=targetTree
     * \return the computed signature matrix
     */
    std::vector< std::vector< bool > > getSignatureMatrix( const std::vector< size_t>& partition, const bool forRefTree ) const;

    /**
     * Evaluate the matching degree of two partitions across trees through their respective signature matrices
     * \param lambda the size-difference correction factor. defines the weighting on how much thedifference in number of clusters between partitions will affect the matching value.
     *        if lambda is set to 0, cluster number difference will not affect at all. If set to 1, it will have as much weight as signature correlation.
     * \param refSize size of the reference tree partition (in number of clusters)
     * \param refSignature signature matrix of the reference tree partitio
     * \param targetSize size of the target tree partition (in number of clusters)
     * \param targetSignature signature matrix of the target tree partition
     * \return the parition match value
     */
    float evalSignaturePartMatch( const float lambda, size_t refSize, std::vector< std::vector< bool > >& refSignature, size_t targetSize,  std::vector< std::vector< bool > >& targetSignature ) const;

    /**
     * Returns the number of matched-base-nodes overlapping across two clusters
     * \param refCluster the reference cluster
     * \param targetCluster the target cluster
     * \return the number of matchign base nodes
     */
    size_t getClusterOverlap( const size_t refCluster, const size_t targetCluster ) const;

    /**
     * Returns a matrix with the number of overlapping base-nodes of each cluster in partition1 to each cluster in partition2 (unidirectional)
     * \param partition1 the reference partition
     * \param partition2 the target partition
     * \return the overlap matrix
     */
    std::pair< std::vector< std::vector< size_t > >, std::vector< std::vector< size_t > > > getClusterOverlapMatrix( const std::vector< size_t > &partition1, const std::vector< size_t > &partition2 ) const;

    /**
     * Computes a unidirectionl best matching lookup table based on the overlap matrix
     * \param matchMatrix theoverlap match matrix
     * \return the unidirectional matching lookup table
     */
    std::pair< std::vector< size_t >, std::vector< size_t > > getClusterMatchTable( const std::vector< std::vector< size_t > > &matchMatrix ) const;

    /**
     * Computes an overlap-based bidirectional best matching lookup table of two partitions across trees
     * \param partition1 the reference tree partition
     * \param partition2 the target tree partition
     * \retval tablePoint1 a pointer to where the first unidirectional lookup table will be saved
     * \retval tablePoint2 a pointer to where the second unidirectional lookup table will be saved
     * \return the quality value of the matching
     */
    double evalOverlapPartMatch( const std::vector< size_t > &partition1, const std::vector< size_t > &partition2,
                       std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint1,
                       std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint2) const;

    /**
     * Assigns a pre-defined static depth-search value  (for partition matching) based on the size of the current target partition
     * Bigger partitions (in number of clusters) are assigned lower search depth in order to optimize computation time
     * \param partSize the current target partition size
     * \return the assigned search depth level
     */
    size_t assignDepth( const size_t partSize );

    /**
     * Shifts the current color to a similar hue/intensity according to a shift index
     * used to color differently but similiarly multiple clusters from a partition that match a single cluster in the matched partition.
     * \param color the original color to be shifted
     * \param shiftIndex the shift index (different indexes will produce different shifted colors)
     * \return the resulting shifted color
     */
    WHcoord shiftColor( const WHcoord &color, const size_t shiftIndex ) const;

    /**
     * Given a partition of the hierarchical tree, recursively scans and returns all partitions derived from parforming branchings on the nodes of the original partition
     * taht yield new partitions smaller or equal (in number of clusters) to a max Partition sized accepted as parameter
     * \param maxPartSize the maximum partition size, partitions bigger than this size will not be returned
     * \param thisPart the partition to be scanned and sub-partitioned
     * \retval partVectorPtr a pointer to the vector of partitions where the obtained sub-partitions will be saved to
     */
    void searchPartition( size_t maxPartSize, const std::vector< size_t > &thisPart, std::vector< std::vector< size_t > >* partVectorPtr);
};

#endif  // TREECOMPARER_H
