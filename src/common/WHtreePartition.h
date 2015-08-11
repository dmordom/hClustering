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


#ifndef WHTREEPARTITION_H
#define WHTREEPARTITION_H

// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <cstdlib>
#include <limits>


// hClustering
#include "WHtree.h"

// type of classic partitioning
typedef enum
{
    HTP_HOZ, // horizontal cut (distance level)
    HTP_SIZE, // biggest cluster size limit
    HTP_HLEVEL // hierarchical level limit
}
HT_PARTMODE;

// type of optimized partitioning
typedef enum
{
    HTP2_OPT, // optimized partition spread-separation index
    HTP2_CSD, // Minimum cluster size difference
    HTP2_MIAD, // Maximum allowed intra cluster distance
    HTP2_WIAD, // Maximum allowed weighted intra cluster distance
    HTP2_MIRD, // Minimum required inter cluster distance
    HTP2_WIRD // Minimum required weighted inter cluster distance
}
HT_PARTMODE2;

// type of condition to use for partitioning
typedef enum
{
    HTC_VALUE, // value (i.e, distance level)
    HTC_CNUM // number of clusters
}
HT_CONDITION;


/**
 * this class operates over the class WHtree to extract partitions and analyze partition quality
 * it may change some non-essential members of the tree class such as saved selected partitions
 * but the tree structure is never modified.
 */
class WHtreePartition
{
public:
    /**
     * Constructor
     * \param tree a pointer to the tree to be partitioned
     */
    explicit WHtreePartition( WHtree* const tree );

    //! Destructor
    ~WHtreePartition();

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Obtain the partition quality values for the optimized (inter vs weighted intra cluster distance) partitions at each granulairty of the tree
     * \param levelDepth number of levels down where to search for the best branching decision
     * \retval partitionValuesPointer pointer to the vector where to save the partition quality values
     * \retval patitionVectorPointer pointer to the vector where to save the partition actual partitions
     * \param verbose flag to activate debug output
     */
    void scanOptimalPartitions( const size_t levelDepth,
                                std::vector< float >* partitionValuesPointer,
                                std::vector< std::vector< size_t> >* patitionVectorPointer,
                                const bool verbose = false );

    /**
     * Obtain the partition quality values for the classic partitions at each granulairty of the tree
     * \retval partitionValuesPointer pointer to the vector where to save the partition quality values
     * \retval patitionVectorPointer pointer to the vector where to save the partition actual partitions
     * \param verbose flag to activate debug output
     */
    void scanHozPartitions( std::vector< float >* partitionValuesPointer,
                            std::vector< std::vector< size_t> >* patitionVectorPointer,
                            const bool verbose = false );

    /**
     * Filter a set of partitions across granularity levels in orther to keep only local maxima within a certain granularity "radius"
     * \param filterRadius radius (in granularity steps) within to keep the local maxima
     * \retval partValPoint pointer to the partition value vector to be decimated
     * \retval partVectorPoint  pointer to the partition vector to be decimated
     * \return the index of the partition with the absolute maximum of quality vaule
     */
    size_t filterMaxPartitions( const unsigned int filterRadius,
                                std::vector< float > *partValPoint,
                                std::vector< std::vector< size_t> > *partVectorPoint );

    /**
     * write a set of partitions and their quality values to a text file
     * \param partFileName file name where to write  the partition selection
     * \param partitionValues partition value vector to be written
     * \param patitionVector partition vector to be written
     */
    void writePartitionSet( std::string partFileName,
                            const std::vector< float > &partitionValues,
                            const std::vector< std::vector< size_t> > &patitionVector );


    /**
     * obtain the partition corresponding to the second lowest level of granularity of a tree (using full boolean-integer node ID)
     * that is, those those clusters correspondign to nodes having only one branching in between them and single leaves
     * \param partition vector where the return partition is saved
     */
    void level2granularity( std::vector<nodeID_t>* partition ) const;

    /**
     * obtain the partition corresponding to the second lowest level of granularity of a tree (using only non-leaf node ID)
     * that is, those those clusters correspondign to nodes having only one branching in between them and single leaves
     * \param partition vector where the return partition is saved
     */
    void level2granularity( std::vector<size_t>* partition ) const;

    /**
     * Gets a classic partition for the tree
     * (those where clusters are further subdivided until a condition is met, based on number of clusters, level or biggest cluster size)
     * \param compValue comparison value to select partition
     * \retval partition vector where to save the partition
     * \param mode tipe of data on which to base the partition (distance, size h. level)
     * \param condition condition on which to build the partition (condition value or number of clusters)
     * \param excludeLeaves if set base nodes will not be further partitioned
     * \param root root of the subtree where to find the aprtition
     * \return value at wich partition was found
     */
    float partitionClassic( float compValue, std::vector<nodeID_t>* const partition, const HT_PARTMODE mode,
                            const HT_CONDITION condition, const bool excludeLeaves, const size_t root ) const;

    /**
     * Gets an optimized partition for the tree
     * (those where a partition characteristic must be cosnidered for several partitions and levels at a time to decide which subset of the tree to continue navigating,
     * until a condition is met)
     * \param compValue comparison value to select partition
     * \retval partition vector where to save the partition
     * \param mode tipe of data on which to base the partition ( spread-separation index, minimum size difference,
     * maximum allowed intra cluster distance, maximum allowed weighted intra cluster distance, minimum required inter cluster distance,
     * minimum required weighted inter cluster distance)
     * \param condition condition on which to build the partition (condition value or number of clusters)
     * \param excludeLeaves if set base nodes will not be further partitioned
     * \param root root of the subtree where to find the aprtition
     * \param levelDepth number of levels down where to search for the best branching decision
     * \return value at wich partition was found
     */
    float partitionOptimized( float compValue, std::vector<nodeID_t>* const partition, const HT_PARTMODE2 mode,
                            const HT_CONDITION condition, const bool excludeLeaves, const size_t root, const size_t levelDepth ) const;

    /**
     * Finds clusters with sharp boundaries (long branches) the tree
     * \param compValue comparison value to select partition
     * \retval partition vector where to save the partition
     * \param excludeLeaves if set base nodes will not be further partitioned
     * \param root root of the subtree where to find the aprtition
     * \param normalized flag to indicate whether the absolute length of the brancghwes should be considered of if they should be normalized by distance level in the tree
     * \return value at wich partition was found
     */
    float partitionSharp( const float compValue,
                          std::vector<nodeID_t>* const partition,
                          const bool excludeLeaves,
                          const size_t root,
                          const bool normalized  ) const;

    /**
     * Finds regions without sharp inner boundaries the tree
     * \param compValue comparison value to select partition
     * \retval partition vector where to save the partition
     * \param excludeLeaves if set base nodes will not be further partitioned
     * \param root root of the subtree where to find the aprtition
     * \return value at wich partition was found
     */
    float partitionSmooth( const float compValue, std::vector<nodeID_t>* const partition, const bool excludeLeaves, const size_t root ) const;


private:
    // === PRIVATE MEMBER DATA ===

    //! tree object
    WHtree &m_tree;

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Evaluate partition cluster size difference (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return cluster size difference value of the partition
     */
    std::pair< float, float > evalPartClustSizeDiff( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate partition cluster size difference (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return cluster size difference value of the partition
     */
    std::pair< float, float > evalPartClustSizeDiff( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate partition quality: spread-separation index  (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return quality value of the partition
     */
    float evalPartOptimal( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate partition quality: spread-separation index (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return quality value of the partition
     */
    float evalPartOptimal( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate spread-separation index for selected partition (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return spread-separation index value of the partition
     */
    float evalSSindex( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate spread-separation index for selected partition (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return spread-separation index value of the partition
     */
    float evalSSindex( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate intra cluster distance factor for selected partition (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return intra cluster distance factor value of the partition
     */
    float evalPartIntraDist( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate intra cluster distance factor for selected partition (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return intra cluster distance factor value of the partition
     */
    float evalPartIntraDist( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate weighted intra cluster distance factor for selected partition (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return weighted intra cluster distance factor value of the partition
     */
    float evalPartIntraDistWeighted( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate weighted intra cluster distance factor for selected partition (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return weighted intra cluster distance factor value of the partition
     */
    float evalPartIntraDistWeighted( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate approximated inter cluster distance factor for selected partition (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return approximated inter cluster distance factor value of the partition
     */
    float evalPartBranchDist( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate approximated inter cluster distance factor for selected partition (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return approximated inter cluster distance factor value of the partition
     */
    float evalPartBranchDist( const std::vector<nodeID_t> &partition ) const;

    /**
     * Evaluate approximated weighted inter cluster distance factor for selected partition (using only non-leaf node ID)
     * \param partition vector with the partition to be evaluated
     * \return approximated weighted inter cluster distance value of the partition
     */
    float evalPartBranchDistWeighted( const std::vector<size_t> &partition ) const;

    /**
     * Evaluate approximated weighted inter cluster distance factor for selected partition (using full boolean-integer node ID)
     * \param partition vector with the partition to be evaluated
     * \return approximated weighted inter cluster distance value of the partition
     */
    float evalPartBranchDistWeighted( const std::vector<nodeID_t> &partition ) const;

    // DEPRECATED MEMEBER FUNCTIONS

    /**
     * Obtain pairwise intercluster distance matrix for the selected partition, or derive new icd matrix from an old one following the subdivision of one of the clusters (DEPRECATED)
     * \param oldPartition vector with the partition to be evaluated or the previously evaluated partition
     * \param branchPos index of the cluster that was further subdivided if considerign a previously calculated matrix, 0 if no previous calculation exists (default)
     * \param branch vector with sub-clusters to replace original cluster after subdivision, empty if no previous calculation exists (default)
     * \param oldMatrix icd matrix from previous calculation, empty if no previous calculation exists (default)
     * \return pairwise inter cluster distance matrix for the selected (or derived) parititon
     */
    std::vector< std::vector< dist_t > > getICDmatrix( const std::vector<size_t> &oldPartition,
                                                       size_t branchPos = 0,
                                                       std::vector<size_t> branch = std::vector<size_t>(),
                                                       std::vector< std::vector< dist_t > > oldMatrix = std::vector< std::vector< dist_t > >() );

    /**
     * Evaluate full inter cluster distance factor for selected partition (DEPRECATED)
     * \param partition vector with the partition to be evaluated
     * \param icdMatrix pairwise inter-cluster distance matrix for the selected partition
     * \return full inter cluster distance factor value of the partition
     */
    float evalPartInterDist( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const;

    /**
     * Evaluate full weighted inter cluster distance factor for selected partition (DEPRECATED)
     * \param partition vector with the partition to be evaluated
     * \param icdMatrix pairwise inter-cluster distance matrix for the selected partition
     * \return full weighted inter cluster distance factor value of the partition
     */
    float evalPartInterDistWeighted( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const;
};

#endif  // WHTREEPARTITION_H
