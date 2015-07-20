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

typedef enum
{
    HTP_HOZ,
    HTP_SIZE,
    HTP_HLEVEL
}
HT_PARTMODE;

typedef enum
{
    HTP2_CSD,
    HTP2_MIAD,
    HTP2_WIAD,
    HTP2_MIRD,
    HTP2_WIRD,
    HTP2_OPT
}
HT_PARTMODE2;

typedef enum
{
    HTC_VALUE,
    HTC_CNUM
}
HT_CONDITION;


/**
 * it operates over the class WHtree
 *
 */
class WHtreePartition
{
public:
    //! Constructor
    explicit WHtreePartition( WHtree* const tree );

    //! Destructor
    ~WHtreePartition();

    // === PUBLIC MEMBER FUNCTIONS ===

    void scanOptimalPartitions( const size_t levelDepth, std::vector< float > &partitionValues,
                                std::vector< std::vector< size_t> > &patitionVector, const bool verbose = false );


    void scanHozPartitions(std::vector< float > &partitionValues, std::vector< std::vector< size_t> > &patitionVector, const bool verbose = false );

    size_t filterMaxPartitions(const unsigned int filterRadius, std::vector< float > *partValPoint, std::vector< std::vector< size_t> > *partVectorPoint );


    void writePartitionSet( std::string partFileName, const std::vector< float > &partitionValues, const std::vector< std::vector< size_t> > &patitionVector );

    void level2granularity( std::vector<nodeID_t>* const partition ) const;
    void level2granularity( std::vector<size_t>* const partition ) const;

    //! Gets a partition for the tree
    //! \param compValue comparison value to select partition
    //! \param partition vector where to save the partition
    //! \param mode tipe of data on whoch to base the partition (distance, size h. level)
    //! \param condition condition on which to build the partition (condition value or number of clusters)
    //! \param excludeLeaves if set base nodes will not be further partitioned
    //! \param root root of the subtree where to find the aprtition
    //! \return value at wich partition was found
    float partitionClassic( float compValue, std::vector<nodeID_t>* const partition, const HT_PARTMODE mode,
                            const HT_CONDITION condition, const bool excludeLeaves, const size_t root ) const;

    //! Gets an optimized partition for the tree
    //! \param compValue comparison value to select partition
    //! \param partition vector where to save the partition
    //! \param mode tipe of data on whoch to base the partition (distance, size h. level)
    //! \param condition condition on which to build the partition (condition value or number of clusters)
    //! \param excludeLeaves if set base nodes will not be further partitioned
    //! \param root root of the subtree where to find the aprtition
    //! \return value at wich partition was found
    float partitionOptimized( float compValue, std::vector<nodeID_t>* const partition, const HT_PARTMODE2 mode,
                            const HT_CONDITION condition, const bool excludeLeaves, const size_t root, const size_t levelDepth ) const;

    //! Finds clusters with sharp boundaries (long branches) the tree
    //! \param compValue comparison value to select partition
    //! \param partition vector where to save the partition
    //! \param excludeLeaves if set base nodes will not be further partitioned
    //! \param root root of the subtree where to find the aprtition
    //! \return value at wich partition was found
    float partitionSharp( const float compValue, std::vector<nodeID_t>* const partition, const bool excludeLeaves, const size_t root, const bool normalized  ) const;

    //! Finds clusters without inner boundaries the tree
    //! \param compValue comparison value to select partition
    //! \param partition vector where to save the partition
    //! \param excludeLeaves if set base nodes will not be further partitioned
    //! \param root root of the subtree where to find the aprtition
    //! \return value at wich partition was found
    float partitionSmooth( const float compValue, std::vector<nodeID_t>* const partition, const bool excludeLeaves, const size_t root ) const;


private:
    // === PRIVATE MEMBER DATA ===


    //! tree object
    WHtree &m_tree;

    // === PRIVATE MEMBER FUNCTIONS ===


    //! Gets an optimized partition for the tree
    //! \param partition vector with the partition to be evaluated
    //! \return quality value of the partition
    std::pair< float, float > evalPartClustSizeDiff( const std::vector<size_t> &partition ) const;
    std::pair< float, float > evalPartClustSizeDiff( const std::vector<nodeID_t> &partition ) const;


    float evalPartOptimal( const std::vector<size_t> &partition ) const;
    float evalPartOptimal( const std::vector<nodeID_t> &partition ) const;

    float evalSSindex( const std::vector<size_t> &partition ) const;
    float evalSSindex( const std::vector<nodeID_t> &partition ) const;


    //! Gets an optimized partition for the tree
    //! \param partition vector with the partition to be evaluated
    //! \return quality value of the partition
    float evalPartIntraDist( const std::vector<size_t> &partition ) const;
    float evalPartIntraDist( const std::vector<nodeID_t> &partition ) const;


    //! Gets an optimized partition for the tree
    //! \param partition vector with the partition to be evaluated
    //! \return quality value of the partition

    float evalPartIntraDistWeighted( const std::vector<size_t> &partition ) const;
    float evalPartIntraDistWeighted( const std::vector<nodeID_t> &partition ) const;


    float evalPartBranchDist( const std::vector<size_t> &partition) const;
    float evalPartBranchDist( const std::vector<nodeID_t> &partition) const;


    float evalPartBranchDistWeighted( const std::vector<size_t> &partition) const;
    float evalPartBranchDistWeighted( const std::vector<nodeID_t> &partition) const;




    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    std::vector< std::vector< dist_t > > getICDmatrix( const std::vector<size_t> &oldPartition, size_t branchPos = 0, std::vector<size_t> branch = std::vector<size_t>(),
                                                      std::vector< std::vector< dist_t > > oldMatrix = std::vector< std::vector< dist_t > >() );

    //! Gets an optimized partition for the tree
    //! \param partition vector with the partition to be evaluated
    //! \return quality value of the partition
    float evalPartInterDist( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const;

    float evalPartInterDistWeighted( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const;



};

#endif  // WHTREEPARTITION_H
