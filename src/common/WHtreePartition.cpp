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

#include <vector>
#include <list>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <limits>


#include "WStringUtils.h"

#include "WHtreePartition.h"



WHtreePartition::WHtreePartition( WHtree* const tree ) : m_tree( *tree )
{
}

WHtreePartition::~WHtreePartition()
{
    //Cleanup
}

// === PUBLIC MEMBER FUNCTIONS ===

void WHtreePartition::scanOptimalPartitions( const size_t levelDepth, std::vector< float >* partitionValuesPointer,
                                             std::vector< std::vector< size_t> >* partitionVectorPointer, const bool verbose )
{
    std::vector< float >& partitionValues = *partitionValuesPointer;
    std::vector< std::vector< size_t> >& partitionVector = *partitionVectorPointer;

    partitionValues.clear();
    partitionVector.clear();

    size_t match250( 0 ), mismatch250( 0 ), match500( 0 ), mismatch500( 0 ), match1000( 0 ), mismatch1000( 0 );
    size_t match2000( 0 ), mismatch2000( 0 ), match5000( 0 ), mismatch5000( 0 );


    std::vector< size_t > currentPartition, lastPartition;
    float currentValue;

    // do first step
    {
        const WHnode treeRoot( m_tree.getRoot() );
        std::vector< nodeID_t > firstKids( treeRoot.getChildren() );
        for( size_t i = 0; i < firstKids.size(); ++i )
        {
            // get kids of root as first partition set
            if( firstKids[i].first )
            {
                currentPartition.push_back( firstKids[i].second );
            }
        }

        currentValue = evalPartOptimal( currentPartition );

        partitionValues.push_back( currentValue );
        partitionVector.push_back( currentPartition );
    } // end first step

    lastPartition = currentPartition;
    if( verbose )
    {
        std::cout << "Step: " << 0 << ". Current partition size: " << currentPartition.size() << ". Current value: " << currentValue << std::flush;
    }

    size_t stepNr( 0 );

    while( true )
    {
        ++stepNr;
        std::cout << "\rStep: " << stepNr << ". Current partition size: ";
        std::cout << currentPartition.size() << ". Current value: ";
        std::cout << currentValue << "       " << std::flush;

        std::vector<float> currentLevels;
        currentLevels.reserve( currentPartition.size() );
        size_t highestLevelindex( 0 );
        for( size_t i = 0; i < currentPartition.size(); ++i )
        {
            currentLevels.push_back( m_tree.getNode( currentPartition[i] ).getDistLevel() );
            if( currentLevels[i] > currentLevels[highestLevelindex] )
            {
                highestLevelindex = i;
            }
        }

        lastPartition = currentPartition;

        std::vector< std::vector< size_t > > derivedPartitionSet;
        std::vector< std::vector< unsigned int > > derivedIndexes( m_tree.getBranching( currentPartition, levelDepth, &derivedPartitionSet ) );

        std::vector< double > derivedPartitionValues( derivedPartitionSet.size(), -1 );
        bool stopLoop( true );

        // get the partitions and values for all possible branchings
        #pragma omp parallel for schedule( guided )
        for( size_t i = 0; i < derivedPartitionSet.size(); ++i )
        {
            // if there are possible branchings continue looping
            #pragma omp critical
            stopLoop = false;

            //get the value for this partition
            derivedPartitionValues[i] = evalPartOptimal( derivedPartitionSet[i] );
        } // endFor

        //if there are no more possible branchings stop the loop
        if( stopLoop )
        {
            break;
        }


        size_t bestPartitionIndex( 0 );
        float bestValue( -1 );
        for( size_t j = 0; j < derivedPartitionValues.size(); ++j )
        {
            //if value is better than previous ones, keep partition and matrix
            if( derivedPartitionValues[j] > bestValue )
            {
                bestValue = derivedPartitionValues[j];
                bestPartitionIndex = j;
            }
        }

        unsigned int firstBranchIndex( derivedIndexes[bestPartitionIndex].front() );

        // find the first branching partition corresponding to the best partition found
        if( derivedIndexes[bestPartitionIndex].size() == 1 )
        {
            // it was a first branch partition
            currentPartition = derivedPartitionSet[bestPartitionIndex];
        }
        else
        {
            size_t nextBestIndex( 0 );
            bool notFound( true );

            for( size_t j = 0; j < derivedIndexes.size(); ++j )
            {
                if( ( derivedIndexes[j].size() == 1 ) && ( derivedIndexes[j].front() ==  firstBranchIndex ) )
                {
                    nextBestIndex = j;
                    notFound = false;
                }
            }

            if( notFound )
            {
                throw std::runtime_error( "ERROR @ partitionMatcher::findMatchingPartitions(): best match index not found" );
            }

            currentPartition = derivedPartitionSet[nextBestIndex];
            currentValue = derivedPartitionValues[nextBestIndex];
        }

        //introduce best partition on the vector and update current vector and matrix
        partitionValues.push_back( currentValue );
        partitionVector.push_back( currentPartition );


        std::vector< nodeID_t > partHozFullID;
        partitionClassic( currentPartition.size(), &partHozFullID, HTP_HOZ, HTC_CNUM, true, m_tree.getRoot().getID() );
        std::vector< size_t > partHoz;
        partHoz.reserve( partHozFullID.size() );
        for( size_t i = 0; i < partHozFullID.size(); ++i )
        {
            if( partHozFullID[i].first )
            {
                partHoz.push_back( partHozFullID[i].second );
            }
        }

        std::vector< size_t > currentCopy( currentPartition );
        std::sort( partHoz.begin(), partHoz.end() );
        std::sort( currentCopy.begin(), currentCopy.end() );


        //bool matchcondition( firstBranchIndex == highestLevelindex );
        bool matchcondition( partHoz == currentCopy );


        if( currentPartition.size() <= 250 )
        {
            if( matchcondition )
            {
                ++match250;
            }
            else
            {
                ++mismatch250;
            }
        }
        if( currentPartition.size() <= 500 )
        {
            if( matchcondition )
            {
                ++match500;
            }
            else
            {
                ++mismatch500;
            }
        }
        if( currentPartition.size() <= 1000 )
        {
            if( matchcondition )
            {
                ++match1000;
            }
            else
            {
                ++mismatch1000;
            }
        }
        if( currentPartition.size() <= 2000 )
        {
            if( matchcondition )
            {
                ++match2000;
            }
            else
            {
                ++mismatch2000;
            }
        }
        if( currentPartition.size() <= 5000 )
        {
            if( matchcondition )
            {
                ++match5000;
            }
            else
            {
                ++mismatch5000;
            }
        }
    } // end infinite loop

    if( verbose )
    {
        std::cout << std::endl;
    }

    /*
    std::cout << "Are partition decisions the same as in horizontal partitioning???" << std::endl;
    std::cout << "250 % match: " << ( match250*100.0 )/( match250+mismatch250 );
    std::cout << ". Total matches: " << match250 << ". fails: " << mismatch250 << std::endl;
    std::cout << "500 % match: " << ( match500*100.0 )/( match500+mismatch500 );
    std::cout << ". Total matches: " << match500 << ". fails: " << mismatch500 << std::endl;
    std::cout << "1000 % match: " << ( match1000*100.0 )/( match1000+mismatch1000 );
    std::cout << ". Total matches: " << match1000 << ". fails: " << mismatch1000 << std::endl;
    std::cout << "2000 % match: " << ( match2000*100.0 )/( match2000+mismatch2000 );
    std::cout << ". Total matches: " << match2000 << ". fails: " << mismatch2000 << std::endl;
    std::cout << "5000 % match: " << ( match5000*100.0 )/( match5000+mismatch5000 );
    std::cout << ". Total matches: " << match5000 << ". fails: " << mismatch5000 << std::endl;
    */



    return;
} // end scanOptimalPartitions() -------------------------------------------------------------------------------------


void WHtreePartition::scanHozPartitions( std::vector< float >* partitionValuesPointer,
                                         std::vector< std::vector< size_t> >* partitionVectorPointer,
                                         const bool verbose )
{
    std::vector< float >& partitionValues = *partitionValuesPointer;
    std::vector< std::vector< size_t> >& partitionVector = *partitionVectorPointer;

    partitionValues.clear();
    partitionVector.clear();

    std::vector< size_t > currentPartition, lastPartition;
    float currentValue;

    // do first step
    {
        const WHnode treeRoot( m_tree.getRoot() );
        std::vector< nodeID_t > firstKids( treeRoot.getChildren() );
        for( size_t i = 0; i < firstKids.size(); ++i )
        {
            // get kids of root as first partition set
            if( firstKids[i].first )
            {
                currentPartition.push_back( firstKids[i].second );
            }
        }

        currentValue = evalPartOptimal( currentPartition );

        partitionValues.push_back( currentValue );
        partitionVector.push_back( currentPartition );
    } // end first step

    lastPartition = currentPartition;
    if( verbose )
    {
        std::cout << "Step: " << 0 << ". Current partition size: " << currentPartition.size() << ". Current value: " << currentValue << std::flush;
    }

    size_t stepNr( 0 );
    size_t partSizeCounter( currentPartition.size() );
    size_t lastSize( currentPartition.size() );

    while( partSizeCounter <= 5000 )
    {
        ++stepNr;
        ++partSizeCounter;
        std::cout << "\rStep: " << stepNr << ". Current partition size: ";
        std::cout << currentPartition.size() << ". Current value: ";
        std::cout << currentValue << "       " << std::flush;

        std::vector<nodeID_t> partHozFullID;
        partitionClassic( partSizeCounter, &partHozFullID, HTP_HOZ, HTC_CNUM, true, m_tree.getRoot().getID() );
        currentPartition.clear();
        currentPartition.reserve( partHozFullID.size() );
        for( size_t i = 0; i < partHozFullID.size(); ++i )
        {
            if( partHozFullID[i].first )
            {
                currentPartition.push_back( partHozFullID[i].second );
            }
        }

        currentValue = evalPartOptimal( currentPartition );
        partitionValues.push_back( currentValue );
        partitionVector.push_back( currentPartition );
        partSizeCounter = currentPartition.size();

        if( lastSize == partSizeCounter )
        {
            break;
        }

        lastSize = partSizeCounter;
    } // end infinite loop

    if( verbose )
    {
        std::cout << std::endl;
    }


    return;
} // end scanHozPartitions() -------------------------------------------------------------------------------------



size_t WHtreePartition::filterMaxPartitions( const unsigned int filterRadius,
                                             std::vector< float > *partValPoint,
                                             std::vector< std::vector< size_t> > *partVectorPoint )
{
    std::vector< float > &partValues( *partValPoint );
    std::vector< std::vector< size_t> > &partVector( *partVectorPoint );

    std::vector< float > filteredValues;
    std::vector< std::vector< size_t> > filteredPartitions;

    size_t partCount( 0 ), bestIndex( 0 );
    float bestValue( -1 );

    for( size_t i = 0; i < partValues.size(); ++i )
    {
            // get in a vector the allements within the filter raidzs, centre element is always first entry
        std::vector< float > currentValues;
            currentValues.reserve( 1+ ( 2*filterRadius ) );
            size_t jDownLimit( ( i > filterRadius ) ? ( i - filterRadius ) : 0 );
            size_t jUpLimit( std::min( i + filterRadius + 1, partValues.size() ) );


            size_t centreIndex( 0 );
            for( size_t j = jDownLimit; j < jUpLimit; ++j )
            {
                if( j == i )
                {
                    centreIndex = currentValues.size();
                }
                currentValues.push_back( partValues[j] );
            }

            // centre element is the biggest
            if( std::max_element( currentValues.begin(), currentValues.end() ) == ( currentValues.begin() + centreIndex ) )
            {
                filteredValues.push_back( partValues[i] );
                filteredPartitions.push_back( partVector[i] );

                if( partValues[i] > bestValue )
                {
                    bestValue = partValues[i];
                    bestIndex = partCount;
                }

                ++partCount;
            }
    }

    partValues.swap( filteredValues );
    partVector.swap( filteredPartitions );

    return bestIndex;
} // end filterMaxPartitions() -------------------------------------------------------------------------------------



void WHtreePartition::writePartitionSet( std::string partFileName,
                                         const std::vector< float > &partitionValues,
                                         const std::vector< std::vector< size_t> > &patitionVector )
{
    std::ofstream partFile( partFileName.c_str() );
    if( !partFile )
    {
        std::cerr << "ERROR: unable to open out file: \"" << partFileName << "\"" << std::endl;
        exit( -1 );
    }



    partFile << "#value size" << std::endl;

    for( size_t i = 0; i < partitionValues.size(); ++i )
    {
        partFile << partitionValues[i] << " " << patitionVector[i].size() << std::endl;
    }

    return;
} // end writePartitionSet() -------------------------------------------------------------------------------------


// === PRIVATE MEMBER FUNCTIONS ===


float WHtreePartition::evalPartOptimal( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartOptimal( partFullID );
} // end evalPartOptimal() -------------------------------------------------------------------------------------
float WHtreePartition::evalPartOptimal( const std::vector<nodeID_t> &partition ) const
{
//    float intraDvalue( evalPartIntraDistWeighted( partition ) );
//    float interDvalue( evalPartBranchDist ( partition ) );
//    return interDvalue / intraDvalue;

    return evalSSindex( partition );
} // end evalPartOptimal() -------------------------------------------------------------------------------------



float WHtreePartition::evalSSindex( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalSSindex( partFullID );
} // end evalSSindex() -------------------------------------------------------------------------------------
float WHtreePartition::evalSSindex( const std::vector<nodeID_t> &partition ) const
{
//    double ssSum( 0 );
//    double sizeSum( 0 );

      double spreadSum( 0 );
      double sepSum( 0 );
      double sizeSum( 0 );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        const WHnode &thisNode( m_tree.getNode( partition[i] ) );
        const WHnode &parentNode( m_tree.getNode( thisNode.getParent() ) );

//        ssSum += ( parentNode.getDistLevel() * thisNode.getSize() ) / thisNode.getDistLevel() ;
//        sizeSum += thisNode.getSize();

        spreadSum += thisNode.getSize() *thisNode.getDistLevel();
        sepSum += parentNode.getDistLevel();
        sizeSum += thisNode.getSize();
    }
//    return ssSum / sizeSum ;

    double part1( sizeSum/partition.size() );
    double part2( sepSum/spreadSum );

    return ( part1*part2 );
} // end evalSSindex() -------------------------------------------------------------------------------------


std::pair< float, float > WHtreePartition::evalPartClustSizeDiff( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartClustSizeDiff( partFullID );
} // end evalPartitClustSizeD() -------------------------------------------------------------------------------------
std::pair< float, float > WHtreePartition::evalPartClustSizeDiff( const std::vector<nodeID_t> &partition ) const
{
    double diffSqSum( 0 ), sizeSum( 0 );
    for( size_t i = 0; i < partition.size() - 1; ++i )
    {
        double size1( m_tree.getNode( partition[i] ).getSize() );
        for( size_t j = i+1; j < partition.size(); ++j )
        {
            double size2( m_tree.getNode( partition[j] ).getSize() );
            double sizeDif( size1 - size2 );
            diffSqSum += sizeDif*sizeDif;
        }
        sizeSum += size1;
    }
    sizeSum += m_tree.getNode( partition.back() ).getSize();
    double M( partition.size() * ( partition.size() -1 ) / 2.0 );
    return std::make_pair( sizeSum/partition.size(), diffSqSum/M);
} // end evalPartitClustSizeD() -------------------------------------------------------------------------------------

float WHtreePartition::evalPartIntraDist( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartIntraDist( partFullID );
} // end evalPartitWintraD() -------------------------------------------------------------------------------------

float WHtreePartition::evalPartIntraDist( const std::vector<nodeID_t> &partition ) const
{
     double iaDistSum( 0 );
     for( size_t i = 0; i < partition.size(); ++i )
     {
         const WHnode &thisNode( m_tree.getNode( partition[i] ) );
         iaDistSum += thisNode.getDistLevel();
     }
     return iaDistSum/partition.size();
} // end evalPartitWintraD() -------------------------------------------------------------------------------------


float WHtreePartition::evalPartIntraDistWeighted( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartIntraDistWeighted( partFullID );
} // end evalPartitIntraDweighted() -------------------------------------------------------------------------------------
float WHtreePartition::evalPartIntraDistWeighted( const std::vector<nodeID_t> &partition ) const
{
    double iaWDistSum( 0 );
    size_t sizeSum( 0 );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        const WHnode &thisNode( m_tree.getNode( partition[i] ) );
        iaWDistSum += thisNode.getDistLevel()*thisNode.getSize();
        sizeSum += thisNode.getSize();
    }
    return iaWDistSum/sizeSum;
} // end evalPartitIntraDweighted() -------------------------------------------------------------------------------------

float WHtreePartition::evalPartBranchDist( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartBranchDist( partFullID );
} // end evalPartitDnb() -------------------------------------------------------------------------------------

float WHtreePartition::evalPartBranchDist( const std::vector<nodeID_t> &partition ) const
{
    double iaDistSum( 0 );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        const WHnode &thisNode( m_tree.getNode( partition[i] ) );
        const WHnode &parentNode( m_tree.getNode( thisNode.getParent() ) );
         iaDistSum += parentNode.getDistLevel();
    }
    return iaDistSum/partition.size();
} // end evalPartitDnb() -------------------------------------------------------------------------------------


float WHtreePartition::evalPartBranchDistWeighted( const std::vector<size_t> &partition ) const
{
    std::vector<nodeID_t> partFullID;
    partFullID.reserve( partition.size() );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        partFullID.push_back( std::make_pair( true, partition[i] ) );
    }
    return evalPartBranchDistWeighted( partFullID );
} // end evalPartitDnb() -------------------------------------------------------------------------------------
float WHtreePartition::evalPartBranchDistWeighted( const std::vector<nodeID_t> &partition ) const
{
    double iaDistSum( 0 ), sizeSum( 0 );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        const WHnode &thisNode( m_tree.getNode( partition[i] ) );
        size_t thisSize( thisNode.getSize() );
        const WHnode &parentNode( m_tree.getNode( thisNode.getParent() ) );
         iaDistSum += parentNode.getDistLevel() * thisSize;
        sizeSum += thisSize;
    }
    return iaDistSum/sizeSum;
} // end evalPartitDnb() -------------------------------------------------------------------------------------





float WHtreePartition::partitionClassic( float compValue, std::vector<nodeID_t> * partition, const HT_PARTMODE mode, const HT_CONDITION condition,
                                         const bool excludeLeaves, const size_t root ) const
{
    partition->clear();
    if( root > m_tree.getRoot().getID() )
    {
        std::cerr << "ERROR @ WHtree::partitionClassic(): branch root ID is out of boundaries (ID: "
                  << root << ", # nodes: " << m_tree.getNumNodes() << " )." << std::endl;
        return 0;
    }

    std::list<nodeID_t> worklist;
    std::list<nodeID_t> storelist;
    worklist.push_back( std::make_pair( true, root ) );
    float currentValue( 0 );
    WHnode currentNode( m_tree.getNode( root ) );
    bool loopCondition( true );

    switch( mode )
    {
    case HTP_HOZ:
        currentValue = m_tree.getNode( root ).getDistLevel();
        if( condition == HTC_CNUM )
        {
            compValue = static_cast<int>( compValue );
        }
        break;
    case HTP_SIZE:
        currentValue = m_tree.getNode( root ).getSize();
        compValue = static_cast<int>( compValue );
        break;
    case HTP_HLEVEL:
        currentValue = m_tree.getNode( root ).getHLevel();
        compValue = static_cast<int>( compValue );
        break;
    default:
        std::cerr << "ERROR @ WHtree::partitionClassic(): mode not recognized " << std::endl;
        return 0;
        break;
    }

    switch( condition )
    {
    case HTC_VALUE:
        loopCondition = ( currentValue > compValue );
        break;
    case HTC_CNUM:
        loopCondition = ( worklist.size() + storelist.size() < compValue );
        break;
    default:
        std::cerr << "ERROR @ WHtree::partitionClassic(): condition not recognized " << std::endl;
        return 0;
        break;
    }

    while( loopCondition )
    {
        const WHnode current( m_tree.getNode( worklist.front() ) );
        worklist.pop_front();
        std::vector<nodeID_t> kids( current.getChildren() );
        for( size_t i = 0; i < kids.size(); ++i )
        {
            const WHnode thisKid( m_tree.getNode( kids[i] ) );
            if( thisKid.isLeaf() )
            {
                storelist.push_back( thisKid.getFullID() );
            }
            else if( thisKid.getHLevel() == 1 && excludeLeaves )
            {
                storelist.push_back( thisKid.getFullID() );
            }
            else
            {
                worklist.push_back( thisKid.getFullID() );
            }
        }

        if( worklist.empty() )
        {
            break;
        }

        switch( mode )
        {
        case HTP_HOZ:
            worklist.sort();
            worklist.reverse();
            currentNode = m_tree.getNode( worklist.front() );
            currentValue = currentNode.getDistLevel();
            break;
        case HTP_SIZE:
            m_tree.sortBySize( &worklist );
            worklist.reverse();
            currentNode = m_tree.getNode( worklist.front() );
            currentValue = currentNode.getSize();
            break;
        case HTP_HLEVEL:
            m_tree.sortByHLevel( &worklist );
            worklist.reverse();
            currentNode = m_tree.getNode( worklist.front() );
            currentValue = currentNode.getHLevel();
            break;
        default:
            break;
        }

        switch( condition )
        {
        case HTC_VALUE:
            loopCondition = ( currentValue > compValue );
            break;
        case HTC_CNUM:
            loopCondition = ( worklist.size() + storelist.size() < compValue );
            break;
        default:
            break;
        }
    }

    storelist.insert( storelist.end(), worklist.begin(), worklist.end() );
    worklist.clear();

    // assign output current value
    float output( 0 );
    switch( mode )
    {
    case HTP_HOZ:
        if( !currentNode.isRoot() )
        {
            output = ( currentNode.getDistLevel() + m_tree.getNode( currentNode.getID()+1 ).getDistLevel() ) * 0.5;
        }
        else
        {
            output = ( 1+currentNode.getDistLevel() )*0.5;
        }
        break;
    case HTP_SIZE:
        output = currentValue;
        break;
    case HTP_HLEVEL:
        output = currentValue;
        break;
    default:
        break;
    }

    storelist.sort();
    partition->reserve( storelist.size() );
    for( std::list<nodeID_t>::iterator it( storelist.begin() ); it != storelist.end(); ++it )
    {
        partition->push_back( *it );
    }
    return output;
} // end partitionClassic() -------------------------------------------------------------------------------------


float WHtreePartition::partitionOptimized( float compValue, std::vector<nodeID_t> * partition, const HT_PARTMODE2 mode, const HT_CONDITION condition,
                                           const bool excludeLeaves, const size_t root, const size_t levelDepth ) const
{
    std::cout << "DEBUG. level depth: " << levelDepth << std::endl;


    size_t OPTPARTLIMIT = 500;

    partition->clear();
    if( root > m_tree.getRoot().getID() )
    {
        std::cerr << "ERROR @ WHtree::partitionOptimized(): branch root ID is out of boundaries (ID: "
                  << root << ", # nodes: " << m_tree.getNumNodes() << " )." << std::endl;
        return 0;
    }

    std::vector< nodeID_t > bestPartition;
    float condValue( 0 );
    const WHnode currentNode( m_tree.getNode( root ) );
    bestPartition.push_back( currentNode.getFullID() );

    bool loopCondition( true );

    switch( mode )
    {
    case HTP2_CSD:
        condValue = m_tree.getNode( root ).getSize();
        compValue = static_cast<int>( compValue );
        break;
    case HTP2_MIAD:
    case HTP2_WIAD:
        condValue = 1;
        if( condition == HTC_CNUM )
        {
            compValue = static_cast<int>( compValue );
        }
        break;
    case HTP2_MIRD:
    case HTP2_WIRD:
    case HTP2_OPT:
        std::cout << "DEBUG: Spread-Separation" << std::endl;

        condValue = 0;
        if( condition == HTC_CNUM )
        {
            compValue = static_cast<int>( compValue );
        }
        else
        {
            partition->push_back( currentNode.getFullID() );
            return currentNode.getDistLevel();
        }
        break;
    default:
        return 0;
        break;
    }

    switch( condition )
    {
    case HTC_VALUE:
        loopCondition = ( condValue > compValue );
        break;
    case HTC_CNUM:
        std::cout << "DEBUG: by clust num" << std::endl;
        loopCondition = ( bestPartition.size()  < compValue );
        break;
    default:
        break;
    }

    while( loopCondition )
    {
        std::vector< std::vector< nodeID_t > > derivedPartitionSet;
        std::vector< std::vector< unsigned int > > derivedIndexes( m_tree.getBranching( bestPartition,
                                                                                        levelDepth,
                                                                                        &derivedPartitionSet,
                                                                                        excludeLeaves ) );
        std::vector< float > derivedPartitionValues;
        float bestValue( 0 );


        // initialize values
        switch( mode )
        {
        case HTP2_CSD:
        case HTP2_MIAD:
        case HTP2_WIAD:
            derivedPartitionValues.assign( derivedPartitionSet.size(), std::numeric_limits< double >::max() );
            bestValue =  std::numeric_limits< double >::max();
            break;
        case HTP2_MIRD:
        case HTP2_WIRD:
        case HTP2_OPT:
            derivedPartitionValues.assign( derivedPartitionSet.size(), -1 );
            bestValue = -1;
            break;
        default:
            break;
        }

        bool stopLoop( true );

        // get the partitions and values for all possible branchings
#pragma omp parallel for schedule( guided )
        for( size_t i = 0; i < derivedPartitionSet.size(); ++i )
        {
            // if there are possible branchings continue looping
#pragma omp critical
            stopLoop = false;

            derivedPartitionValues[i] = evalPartOptimal( derivedPartitionSet[i] );
            switch( mode )
            {
            case HTP2_CSD:
                derivedPartitionValues[i]=evalPartClustSizeDiff( derivedPartitionSet[i] ).second;
                break;
            case HTP2_MIAD:
                derivedPartitionValues[i]=evalPartIntraDist( derivedPartitionSet[i] );
                break;
            case HTP2_WIAD:
                derivedPartitionValues[i]=evalPartIntraDistWeighted( derivedPartitionSet[i] );
                break;
            case HTP2_MIRD:
                derivedPartitionValues[i] = ( evalPartBranchDist( derivedPartitionSet[i] ) );
                break;
            case HTP2_WIRD:
                derivedPartitionValues[i] = ( evalPartBranchDistWeighted( derivedPartitionSet[i] ) );
                break;
            case HTP2_OPT:
                derivedPartitionValues[i] = ( evalPartOptimal( derivedPartitionSet[i] ) );
                break;
            default:
                break;
            }
        } // endFor

        //if there are no more possible branchings stop the loop
        if( stopLoop )
        {
            break;
        }

        // search for the best partition
        size_t bestPartIndex( 0 );
        for( size_t i = 0; i < derivedPartitionValues.size(); ++i )
        {
            switch( mode )
            {
            case HTP2_CSD:
            case HTP2_MIAD:
            case HTP2_WIAD:
                if( derivedPartitionValues[i] < bestValue )
                {
                    bestValue = derivedPartitionValues[i];
                    bestPartIndex = i;
                }
                break;
            case HTP2_MIRD:
            case HTP2_WIRD:
            case HTP2_OPT:
                if( derivedPartitionValues[i] > bestValue )
                {
                    bestValue = derivedPartitionValues[i];
                    bestPartIndex = i;
                }
                break;
            default:
                break;
            }
        }

        // find the first branching partition corresponding to the best partition found
        if( derivedIndexes[bestPartIndex].size() == 1 )
        {
            // it was a first branch partition
            bestPartition = derivedPartitionSet[bestPartIndex];
        }
        else
        {
            unsigned int firstBranchIndex( derivedIndexes[bestPartIndex].front() );

            size_t nextBestIndex( 0 );
            bool notFound( true );

            for( size_t j = 0; j < derivedIndexes.size(); ++j )
            {
                if( ( derivedIndexes[j].size() == 1 ) && ( derivedIndexes[j].front() ==  firstBranchIndex ) )
                {
                    nextBestIndex = j;
                    notFound = false;
                }
            }

            if( notFound )
            {
                throw std::runtime_error( "ERROR @ partitionMatcher::findMatchingPartitions(): best match index not found" );
            }

            bestPartition = derivedPartitionSet[nextBestIndex];
            bestValue = derivedPartitionValues[nextBestIndex];
        }

        // update condition value
        switch( mode )
        {
        case HTP2_CSD:
            condValue = evalPartClustSizeDiff( bestPartition ).first;
            break;
        case HTP2_MIAD:
        case HTP2_WIAD:
        case HTP2_MIRD:
        case HTP2_WIRD:
        case HTP2_OPT:
            condValue = bestValue;
            break;
        default:
            break;
        }

        // test condition
        switch( condition )
        {
        case HTC_VALUE:
            if( bestPartition.size() > OPTPARTLIMIT )
            {
                loopCondition = false;
            }

            switch( mode )
            {
            case HTP2_CSD:
            case HTP2_MIAD:
            case HTP2_WIAD:
            case HTP2_MIRD:
            case HTP2_WIRD:
                loopCondition = ( condValue < compValue );
                break;
            case HTP2_OPT:
                loopCondition = ( condValue > compValue );
                break;
            default:
                break;
            }

            break;
        case HTC_CNUM:
            loopCondition = ( bestPartition.size() < compValue );
            break;
        default:
            break;
        }
    } // end while loop


    // assign output current value
    std::sort( bestPartition.begin(), bestPartition.end() );
    *partition = bestPartition;
    return condValue;
} // end partitionOptimized() -------------------------------------------------------------------------------------


float WHtreePartition::partitionSharp( const float compValue,
                                       std::vector<nodeID_t>* partition,
                                       const bool excludeLeaves,
                                       const size_t root,
                                       const bool normalized ) const
{
    partition->clear();
    float longestBranch( 0 );

    if( root > m_tree.getRoot().getID() )
    {
        std::cerr << "ERROR @ WHtree::partitionSharp(): branch root ID is out of boundaries (ID: "
                  <<  root << ", # nodes: " << m_tree.getNumNodes() << " )." << std::endl;
        return 0;
    }

    std::list<nodeID_t> worklist;
    std::list<nodeID_t> storelist;

    worklist.push_back( std::make_pair( true, root ) );
    while( !worklist.empty() )
    {
        const WHnode& current( m_tree.getNode( worklist.front() ) );
        worklist.pop_front();
        std::vector<nodeID_t> kids( current.getChildren() );
        for( size_t i = 0; i < kids.size(); ++i )
        {
            const WHnode& thisKid( m_tree.getNode( kids[i] ) );
            float branchValue( ( current.getDistLevel()-thisKid.getDistLevel() )  );
            if( normalized )
            {
                branchValue = branchValue / thisKid.getDistLevel();
            }
            if( branchValue >= compValue )
            {
                storelist.push_back( thisKid.getFullID() );
                if( branchValue > longestBranch )
                {
                    longestBranch = branchValue;
                }
            }
            if( thisKid.isLeaf() )
            {
                continue;
            }
            else if( thisKid.getHLevel() == 1 && excludeLeaves )
            {
                continue;
            }
            else
            {
                worklist.push_back( thisKid.getFullID() );
            }
        }
    }
    storelist.sort();

    partition->reserve( storelist.size() );
    for( std::list<nodeID_t>::iterator it( storelist.begin() ); it != storelist.end(); ++it )
    {
        partition->push_back( *it );
    }
    return longestBranch;
} // end partitionSharp() -------------------------------------------------------------------------------------


float WHtreePartition::partitionSmooth( const float compValue, std::vector<nodeID_t>* partition, const bool excludeLeaves, const size_t root ) const
{
    partition->clear();
    float longestGap( 0 );

    if( root > m_tree.getRoot().getID() )
    {
        std::cerr << "ERROR @ WHtree::partitionSmooth(): branch root ID is out of boundaries (ID: "
                  <<  root << ", # nodes: " << m_tree.getNumNodes() << " )." << std::endl;
        return 0;
    }

    std::vector<unsigned int> conditionMet( m_tree.getNumNodes(), 0 );

    // if base nodes are always included flag them
    if( excludeLeaves )
    {
        std::vector<size_t> bases( m_tree.getBaseNodes( root ) );
        for( size_t i = 0; i < bases.size(); ++i )
        {
            conditionMet[bases[i]] = 1;
        }
    }
    else
    {
        // mark all base nodes that fit condition
        for( size_t i = 0; i < m_tree.getNumLeaves(); ++i )
        {
            const WHnode& dad( m_tree.getNode( m_tree.getLeaf( i ).getParent() ) );
            if( dad.getDistLevel() <= compValue )
            {
                conditionMet[dad.getID()] = 1;
                if( dad.getDistLevel() > longestGap )
                {
                    longestGap = dad.getDistLevel();
                }
            }
        }
    }

    for( size_t i = 0; i < conditionMet.size(); ++i )
    {
        float level = m_tree.getNode( i ).getDistLevel();

        if( conditionMet[i] == 1 )
        {
            if( m_tree.getNode( i ).isRoot() )
            {
                partition->push_back( m_tree.getRoot().getFullID() );
                return longestGap;
            }
            const WHnode& dad( m_tree.getNode( m_tree.getNode( i ).getParent() ) );
            dist_t gap( dad.getDistLevel()-level );
            if( gap <= compValue )
            {
                if( gap > longestGap )
                {
                    longestGap = gap;
                }
                conditionMet[dad.getID()] = 1;
            }
        }
        else if( conditionMet[i] == 0 )
        {
            std::vector<nodeID_t> kids( m_tree.getNode( i ).getChildren() );
            for( size_t j = 0; j < kids.size(); ++j )
            {
                if( kids[j].first )
                {
                    if( conditionMet[kids[j].second] == 1 )
                    {
                        conditionMet[kids[j].second] = 2;
                    }
                }
            }
        }
        else
        {
            std::cerr << "ERROR @ WHtree::partitionSmooth()" << std::endl;
        }
    }

    for( size_t i = 0; i < conditionMet.size(); ++i )
    {
        if( conditionMet[i] == 2 )
        {
            partition->push_back( m_tree.getNode( i ).getFullID() );
        }
    }
    std::sort( partition->begin(), partition->end() );
    return longestGap;
} // end partitionSmooth() -------------------------------------------------------------------------------------


// DEPRECATED FUNCTIONS

std::vector< std::vector< dist_t > > WHtreePartition::getICDmatrix( const std::vector<size_t> &oldPartition,
                                                                    size_t branchPos, std::vector<size_t> branch,
                                                                    std::vector< std::vector< dist_t > > oldMatrix )
{
    std::vector< std::vector< dist_t > >newMatrix;
    newMatrix.reserve( oldPartition.size() + branch.size() );

    // do whole matrix new
    if( oldMatrix.empty() )
    {
        for( size_t i = 0; i < oldPartition.size(); ++i )
        {
            std::vector< dist_t > newLine;
            newLine.reserve( i );
            for( size_t j = 0; j < i; ++j )
            {
                newLine.push_back( m_tree.getNode( m_tree.getCommonAncestor( oldPartition[i], oldPartition[j] ) ).getDistLevel() );
            }
            newMatrix.push_back( newLine );
        }
    }
    else
    {
        for( size_t i = 0; i < oldPartition.size(); ++i )
        {
            // if i<branchPos the line is the same as in the old matrix
            if( i < branchPos )
            {
                newMatrix.push_back( oldMatrix[i] );
            }
            // if i == i the old line has to be replaced by as many lines as new nodesi in the branching
            else if( i == branchPos )
            {
                for( size_t j = 0; j < branch.size(); ++j )
                {
                    std::vector< float > newLine( oldMatrix[branchPos] );
                    std::vector< float > extra( j,  m_tree.getNode( branch[j] ).getDistLevel() );
                    newLine.insert( newLine.end(), extra.begin(), extra.end() );
                    newMatrix.push_back( newLine );
                }
            }
            else
            {
                std::vector< float > newLine( oldMatrix[i] );
                std::vector< float >::iterator insertPoint( newLine.begin() + branchPos );
                std::vector< float > extra( branch.size()-1, *insertPoint );
                newLine.insert( insertPoint, extra.begin(), extra.end() );
                newMatrix.push_back( newLine );
            }
        }
    }
    return newMatrix;
} // end getICDmatrix() -------------------------------------------------------------------------------------


float WHtreePartition::evalPartInterDist( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const
{
    if( partition.size() != icdMatrix.size() )
    {
        std::cerr << "ERROR: matrix size: " <<  icdMatrix.size() << ". PArtition size: " << partition.size() << std::endl;
        throw std::runtime_error( " ERROR @ WHtreePartition::evalPartitInterD(): partition and icd matrix dimensions dont match" );
    }

    double irDistSum( 0 );

    for( size_t i = 0; i < partition.size(); ++i )
    {
        if( icdMatrix[i].size() != i )
        {
            std::cerr << "ERROR: matrix row " << i << " size: " << icdMatrix[i].size() << std::endl;
            throw std::runtime_error( " ERROR @ WHtreePartition::evalPartitInterD(): matrix has wrong dimensions" );
        }
        for( size_t j = 0; j < i; ++j )
        {
            dist_t ieDist( icdMatrix[i][j] );
            irDistSum += ieDist;
        }
    }
    double M( partition.size() * ( partition.size() -1 ) / 2.0 );

    return irDistSum/M;
} // end evalPartitInter() -------------------------------------------------------------------------------------

float WHtreePartition::evalPartInterDistWeighted( const std::vector<size_t> &partition, const std::vector< std::vector < dist_t > > &icdMatrix ) const
{
    if( partition.size() != icdMatrix.size() )
    {
        std::cerr << "ERROR: matrix size: " <<  icdMatrix.size() << ". PArtition size: " << partition.size() << std::endl;
        throw std::runtime_error( " ERROR @ WHtreePartition::evalPartitInterD(): partition and icd matrix dimensions dont match" );
    }

    double irDistSum( 0 ), sizeSum( 0 );

    for( size_t i = 0; i < partition.size(); ++i )
    {
        const WHnode nodeI( m_tree.getNode( partition[i] ) );
        size_t sizeI( nodeI.getSize() );
        if( icdMatrix[i].size() != i )
        {
            std::cerr << "ERROR: matrix row " << i << " size: " << icdMatrix[i].size() << std::endl;
            throw std::runtime_error( " ERROR @ WHtreePartition::evalPartitInterD(): matrix has wrong dimensions" );
        }
        for( size_t j = 0; j < i; ++j )
        {
            const WHnode nodeJ( m_tree.getNode( partition[j] ) );
            size_t sizeJ( nodeJ.getSize() );
            dist_t ieDist( icdMatrix[i][j] );
            irDistSum += ieDist * ( sizeI+sizeJ );
            sizeSum += sizeI+sizeJ;
        }
    }

    return irDistSum/sizeSum;
} // end evalPartitInterDweighted() -------------------------------------------------------------------------------------


void WHtreePartition::level2granularity( std::vector<nodeID_t>* partitionPointer ) const
{
    std::vector<nodeID_t>& partition( *partitionPointer );
    partition.clear();
    std::vector<size_t> partitionst;
    level2granularity( &partitionst );
    partition.reserve( partitionst.size() );

    for( size_t i = 0; i < partitionst.size(); ++i )
    {
        partition.push_back( std::make_pair( true, partitionst[i] ) );
    }
    return;
} // end level2granularity() -------------------------------------------------------------------------------------

void WHtreePartition::level2granularity( std::vector<size_t>* partitionPointer ) const
{
    std::vector<size_t>& partition( *partitionPointer );
    partition.clear();

    std::vector<size_t> baseNodes( m_tree.getRootBaseNodes() );
    std::vector<bool> checkvector( m_tree.getNumNodes(), false );

    for( size_t i = 0; i < baseNodes.size(); ++i )
    {
        if( checkvector[ baseNodes[i] ] )
        {
            continue;
        }

        size_t parent( m_tree.getNode( baseNodes[i] ).getParent().second );
        if( m_tree.getNode( parent ).getHLevel() > 2 )
        {
            partition.push_back( baseNodes[i] );
            checkvector[ baseNodes[i] ] = true;
        }
        else
        {
            partition.push_back( parent );
            checkvector[ parent ] = true;

            std::vector<nodeID_t> kids( m_tree.getNode( parent ).getChildren() );
            for( size_t j = 0; j < kids.size(); ++j )
            {
                if( kids[j].first )
                {
                    checkvector[ kids[j].second ] = true;
                }
            }
        }
    }
    return;
} // end level2granularity() -------------------------------------------------------------------------------------
