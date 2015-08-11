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



// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <algorithm>

#include "WStringUtils.h"

#include "partitionMatcher.h"

#define WARNINGS false
#define DEBUG false



// PUBLIC member functions

partitionMatcher::partitionMatcher( WHtree* const refTree, WHtree* const targetTree, std::string matchFilename, bool verbose ) :
    m_refTree( *refTree ), m_targetTree( *targetTree ), m_verbose( verbose )
{

    if( m_verbose )
    {
        std::cout<< "Testing trees basenodes..."<<std::flush;
    }
    testBaseNodes();
    if( m_verbose )
    {
        std::cout<< ". OK" << std::endl << "Loading correspondence table..."<<std::flush;
    }
    loadCorrespondence(matchFilename);
    if( m_verbose )
    {
        std::cout<< ". OK"<<std::endl;
    }
    if( m_refTree.getSelectedValues().empty() )
    {
        throw std::runtime_error( "ERROR: Tree 1 has no saved partitions to be matched.");
    }
    if( !m_targetTree.getSelectedValues().empty() && m_targetTree.getSelectedValues().size() != m_refTree.getSelectedValues().size() )
    {
        std::cerr<< "WARNING: Partitions of tree 1 and tree 2 have different sizes. Tree 2 partitions were cleared."<<std::endl;
        m_targetTree.clearPartitions();
    }
} // end partitionMatcher::partitionMatcher() -------------------------------------------------------------------------------------


partitionMatcher::~partitionMatcher()
{
}

std::string partitionMatcher::reportBaseNodes() const
{
    std::stringstream message;
    {
        size_t bMax( 0 ), bMin( m_refTree.getNumLeaves() ), numBig( 0 ), numSmall( 0 );
        for( std::vector< size_t >::const_iterator iter( m_refBaseNodes.begin() ); iter != m_refBaseNodes.end(); ++iter )
        {
            size_t currentSize( m_refTree.getNode( *iter ).getSize() );
            if( currentSize > bMax )
            {
                bMax = currentSize;
            }
            bMin = std::min( bMin, currentSize );
            if( currentSize >= 100 )
            {
                ++numBig;
            }
            else if( currentSize <= 10 )
            {
                ++numSmall;
            }
        }
        message << "Tree1: " << m_refBaseNodes.size() << " base nodes. Biggest: " << bMax << ". Smallest: " << bMin << ". "
                        << numBig << " >= 100." << numSmall << " <= 10." << std::endl;
    }

    {
        size_t bMax( 0 ), bMin( m_targetTree.getNumLeaves() ), numBig( 0 ), numSmall( 0 );
        for( std::vector< size_t >::const_iterator iter( m_targetBaseNodes.begin() ); iter != m_targetBaseNodes.end(); ++iter )
        {
            size_t currentSize( m_targetTree.getNode( *iter ).getSize() );
            if( currentSize > bMax )
            {
                bMax = currentSize;
            }
            bMin = std::min( bMin, currentSize );
            if( currentSize >= 100 )
            {
                ++numBig;
            }
            else if( currentSize <= 10 )
            {
                ++numSmall;
            }
        }
        message << "Tree2: " << m_targetBaseNodes.size() << " base nodes. Biggest: " << bMax << ". Smallest: " << bMin << ". "
                        << numBig << " >= 100." << numSmall << " <= 10." << std::flush;
    }
    message << std::endl<<"Trees have " <<  m_refMatchedBases.size() << " Matched nodes "<< std::flush;

    std::string outMessage( message.str() );
    return outMessage;
} // end partitionMatcher::reportBaseNodes() -------------------------------------------------------------------------------------


bool partitionMatcher::matchColors(bool exclusive)
{

    std::cout<< "Matching colors"<<std::endl;

    std::vector< std::vector< WHcoord > > tree1partColors, tree2partColors;
    bool tree1changed( false );

    std::vector< WHcoord >testColors(m_refTree.getSelectedColors( 0 ));
    if(testColors.empty())
    {
        std::cout<< "Tree 1 has no saved colors for selected partitions. Ignoring color matching"<<std::endl;
        return false;
    }


    size_t numPartitions(0);
    {
        std::vector< std::vector< size_t > > tree1partitions( m_refTree.getSelectedPartitions() );
        std::vector< std::vector< size_t > > tree2partitions( m_targetTree.getSelectedPartitions() );
        if( tree1partitions.size() != tree2partitions.size() )
        {
            throw std::runtime_error( "ERROR: trees have different number of partitions.");
        }
        numPartitions = tree1partitions.size();
    }

    if( m_verbose )
    {
        std::cout<< "Matching parttion colors for " << numPartitions << " partitions." << std::endl;
    }

    for( size_t i = 0; i < numPartitions; ++i )
    {

        size_t reColored1(0), reColored2(0), noMatch1(0), noMatch2(0);
        size_t totalMatches(0);

        if( m_verbose )
        {
            std::cout<< std::endl<< "Partition " << i << ": Matching... " << std::endl;
        }

        std::vector< WHcoord > partColors1, partColors2;

        partColors1 = m_refTree.getSelectedColors( i );


        const std::vector< size_t > partition1( m_refTree.getSelectedPartitions( i ) );
        const std::vector< size_t > partition2( m_targetTree.getSelectedPartitions( i ) );
        size_t part1size( partition1.size() );
        size_t part2size( partition2.size() );

        std::pair< std::vector< size_t >, std::vector< size_t > > matchSet1;
        std::pair< std::vector< size_t >, std::vector< size_t > > matchSet2;
        std::vector< size_t > &matchTable1( matchSet1.first );
        std::vector< size_t > &matchTable2( matchSet2.first );
        std::vector< size_t > &matchValues1( matchSet1.second );
        std::vector< size_t > &matchValues2( matchSet2.second );


        float partMatchValue( evalOverlapPartMatch( partition1, partition2, &matchSet1, &matchSet2 ) );


        if( m_verbose )
        {
            std::cout<< "Matching done. Quality index: "<< partMatchValue  << std::endl;
            std::cout<< "Coloring... " << std::endl;
        }

        /// %%%%%% DEBUG %%%%%%
        if ( DEBUG )
        {
            std:: cout <<" == Table 1: =="<<std::endl;
            for( size_t j = 0; j < matchTable1.size(); ++j )
            {
                std::cout << "1:" << j << " (" << m_refTree.getNode(partition1[j]).getSize() << ";" <<  partColors1[j].m_x << "-" << partColors1[j].m_y << "-" << partColors1[j].m_z << ") -> 2:" ;

                if(matchTable1[j] >= part2size )
                {
                    std:: cout <<"*";
                }
                else
                {
                    std:: cout <<matchTable1[j];
                    std:: cout << " -> 1:";
                    if(matchTable2[matchTable1[j]] >= part1size )
                    {
                        std:: cout <<"*";
                    }
                    else
                    {
                        std:: cout <<matchTable2[matchTable1[j]];
                    }
                }

                std::cout  << std::endl;
            }

            std:: cout <<" == Table 2: =="<<std::endl;
            for( size_t j = 0; j < matchTable2.size(); ++j )
            {
                std::cout << "2:" << j <<";"<< m_targetTree.getNode(partition2[j]).getSize() << " -> 1:" ;

                if(matchTable2[j] >= part1size )
                {
                    std:: cout <<"*";
                }
                else
                {
                    std:: cout <<matchTable2[j];
                    std:: cout << " -> 2:";
                    if(matchTable1[matchTable2[j]] >= part2size )
                    {
                        std:: cout <<"*";
                    }
                    else
                    {
                        std:: cout <<matchTable1[matchTable2[j]];
                    }
                }

                std::cout  << std::endl;
            }
        }

        /// %%%%%% END DEBUG %%%%%%


        // to check if all clusters have been accounted for
        std::vector< bool > test1( part1size, false );
        std::vector< bool > test2( part2size, false );
        WHcoord blankColor( 333, 333, 333 );
        partColors2.resize( part2size, blankColor );

        for( size_t j = 0; j < matchTable2.size(); ++j )
        {
            size_t part1matched(matchTable2[j] );
            WHcoord matchedColor( partColors1[part1matched] );

            // it was already processed
            if( test2[j] == true )
            {
                continue;
            }
            // cluster j of part2 has no leaves in common with any cluster from part 1
            if( part1matched >= part1size )
            {
                test2[j] = true;
                ++noMatch2;
                if (exclusive)
                {
                    partColors2[j] = WHcoord(255,255,255);
                }
                else
                {
                    partColors2[j] = WHcoord();
                }
            }
            else
            {
                size_t matchCount( std::count( matchTable2.begin(), matchTable2.end(), part1matched ) );
                if( matchCount == 1 ) // only one cluster of partition 2 matches to this cluster of partition 1
                {
                    partColors2[j] = partColors1[ part1matched ];
                    ++totalMatches;
                    test2[j] = true;

                    // check if part1 cluster matches with it, otherwise  warn
                    if( WARNINGS && matchTable1[part1matched] != j )
                    {
                        std::cerr<< " WARNING: part2:" << j << " is the only one matched to part1:" << part1matched;
                        std::cerr<< " but this one is matched to part2: "<< matchTable1[part1matched] << std::endl;
                    }

                }
                else // several clusters of partition 2 match to the same cluster of partition 1
                {
                    std::vector< size_t > repClustersIndexes;

                    repClustersIndexes.reserve( matchCount );
;
                    float bestValue( 0 );
                    size_t bestIndex( 0 );

                    //recover matching values
                    for( size_t k = 0; k < matchTable2.size(); ++k )
                    {
                        if( matchTable2[k] == part1matched )
                        {
                            repClustersIndexes.push_back( k );
                            if( matchValues2[k] > bestValue )
                            {
                                bestValue = matchValues2[k];
                                bestIndex = k;
                            }
                        }
                    }

                    // assing matching color to biggest match, if it doesnt match with part 1 match give warning
                    partColors2[ bestIndex ] = matchedColor;
                    ++totalMatches;
                    if( WARNINGS && matchTable1[matchTable2[bestIndex]] != bestIndex )
                    {
                        std::cerr<< " WARNING: part2:" << bestIndex << " is the best match to part1:" << matchTable2[bestIndex];
                        std::cerr<< " but this one is matched to part2: "<< matchTable1[matchTable2[bestIndex]] << std::endl;
                    }
                    test2[bestIndex] = true;
                    size_t addIndex( 0 );

                    // for remaining clusters, cheack if any cluster of partition 1 was assigned to them
                    // if so give them that color, otherwise give them a similar color to the matched one
                    for( size_t k = 0; k < repClustersIndexes.size(); ++k )
                    {
                        if( repClustersIndexes[k] == bestIndex )
                        {
                            continue;
                        }
                        test2[repClustersIndexes[k]] = true;

                        std::vector<size_t>::const_iterator findIter(std::find(matchTable1.begin(), matchTable1.end(), repClustersIndexes[k] ) );
                        size_t candidateChosen(0);

                        // a cluster from partition 1 is assigned to it
                        if(findIter != matchTable1.end() )
                        {
                            size_t candidateCount( std::count(matchTable1.begin(), matchTable1.end(), repClustersIndexes[k]) );
                            if( candidateCount > 1)
                            {
                                std::vector< size_t > candidateIndexes;
                                candidateIndexes.reserve( candidateCount );
                                size_t bestCandidateIndex( 0 );
                                size_t bestCandidateValue( 0 );

                                //recover matching values
                                for( size_t l = 0; l < matchTable1.size(); ++l )
                                {
                                    if( matchTable1[l] == repClustersIndexes[k] )
                                    {
                                        candidateIndexes.push_back( l );
                                        if( matchValues1[l] > bestCandidateValue )
                                        {
                                            bestCandidateValue = matchValues1[l];
                                            bestCandidateIndex = l;
                                        }
                                    }
                                }
                                candidateChosen = bestCandidateIndex;
                            }
                            else
                            {
                                candidateChosen = findIter-matchTable1.begin();
                            }

                            std::cerr<< " WARNING: part2:" << repClustersIndexes[k] << " is the a second-best match to part1:" << part1matched;
                            std::cerr<< " bu part1: "<< candidateChosen << " is matched to it. updating matching table"  << std::endl;

                            matchTable2[repClustersIndexes[k]] = candidateChosen;
                            partColors2[repClustersIndexes[k]] = partColors1[candidateChosen];
                            ++totalMatches;


                        }
                        else
                        {
                            partColors2[ repClustersIndexes[k] ] = shiftColor( matchedColor, addIndex++ );
                        }
                    }
                    reColored2 += addIndex;
                }
            }
        }

        if( reColored2 > 0 )
        {
            std::cout<< reColored2 << " clusters of partition 2 were shifted-colored due to one-to-multiple matching." << std::endl;
        }
        if( noMatch2 > 0 )
        {
            std::cout<< noMatch2 << " clusters of partition 2 had no Match." << std::endl;
        }

        std::cout<<"Checking reverse table:" << std::endl;


        for( size_t j = 0; j < matchTable1.size(); ++j )
        {

            size_t part2matched(matchTable1[j]);
            WHcoord matchedColor( partColors2[ part2matched ] );



            // it was already processed
            if( test1[j] == true )
            {
                continue;
            }

            // cluster j of part1 has no leaves in common with any cluster from part 2
            if( part2matched >= part2size )
            {
                test1[j] = true;
                ++noMatch1;
                if (exclusive)
                {
                    partColors1[j] = WHcoord(255,255,255);
                    tree1changed = true;
                }
            }
            else
            {
                size_t matchCount( std::count( matchTable1.begin(), matchTable1.end(), part2matched ) );
                if( matchCount == 0 ) // should not happen
                {
                    throw std::runtime_error(" No matches for partition in part1");
                }
                else if( matchCount == 1 ) // only one cluster of partition 1 matches to this cluster of partition 2
                {
                    test1[j] = true;

                    // check if part2 cluster matches with it, otherwise  warn
                    if( WARNINGS && matchTable2[part2matched] != j )
                    {
                        std::cerr<< " WARNING: part1:" << j << " is the only one matched to part2:" << part2matched;
                        std::cerr<< " but this one is matched to part1: "<< matchTable2[part2matched] << std::endl;
                    }
                }
                else // several clusters of partition 1 match to the same cluster of partition 2
                {
                    std::vector< size_t > repClustersIndexes;

                    repClustersIndexes.reserve( matchCount );
                    size_t bestValue( 0 );
                    size_t bestIndex( 0 );


                    // check to which of the clusters was partition 2 cluster assigned


                    // recover matching values
                    for( size_t k = 0; k < matchTable1.size(); ++k )
                    {
                        if( matchTable1[k] == part2matched )
                        {
                            repClustersIndexes.push_back( k );

                            if( matchValues1[k] > bestValue )
                            {
                                bestValue = matchValues1[k];
                                bestIndex = k;
                            }
                        }
                    }


                    // If best match does not match with that cluster give a warning
                    test1[bestIndex] = true;
                    if( WARNINGS && bestIndex != matchTable2[ part2matched ] )
                    {
                        std::cerr<< " WARNING:  part1:" << bestIndex << " is biggest match to part2:";
                        std::cerr<< part2matched << " , but this was assigned to part1:" << matchTable2[ part2matched ]  << std::endl;
                    }

                    size_t addIndex( 0 );

                    // check other clusters
                    for( size_t k = 0; k < repClustersIndexes.size(); ++k )
                    {
                        if( repClustersIndexes[k] == bestIndex )
                        {
                            continue;
                        }

                        test1[repClustersIndexes[k]] = true;


                        std::vector<size_t>::const_iterator findIter(std::find( matchTable2.begin(), matchTable2.end(), repClustersIndexes[k] ));
                        // check if any cluster of partition 2 were assigned to these other clusters, infomr and ignore if so
                        if( WARNINGS && findIter != matchTable2.end() )
                        {
                            std::cerr<< " PSEUDO-WARNING: cluster part1:" << repClustersIndexes[k] << " is one of multiple ones assigned to part2:" << part2matched;
                            std::cerr<< " but part2: "<< findIter-matchTable2.begin() <<" is matched to it"<< std::endl;
                            continue;
                        }

                        // give similar colors to remaining clusters of partition 1
                        partColors1[ repClustersIndexes[k] ] = shiftColor( matchedColor, addIndex++ );
                        tree1changed = true;

                    }

                    reColored1 += addIndex;
                }
            }
        }

        std::cout<< totalMatches << " matched pairs." << std::endl;

        if( reColored1 > 0 )
        {
            std::cout<< reColored1 << " clusters of partition 1 were shifted-colored due to one-to-multiple matching." << std::endl;
        }
        if( noMatch1 > 0 )
        {
            std::cout<< noMatch1 << " clusters of partition 1 had no Match." << std::endl;
        }


        /// %%%%%% DEBUG %%%%%%

        if( DEBUG)
        {

            std:: cout <<" == Table 1: =="<<std::endl;
            for( size_t j = 0; j < matchTable1.size(); ++j )
            {
                std::cout << "1:" << j << " ("  << m_refTree.getNode(partition1[j]).getSize() << ";"<< partColors1[j].m_x << "-" << partColors1[j].m_y << "-" << partColors1[j].m_z << ") -> 2:" ;

                if(matchTable1[j] >= part2size )
                {
                    std:: cout <<"*";
                }
                else
                {
                    std:: cout <<matchTable1[j];
                    std:: cout << " -> 1:";
                    if(matchTable2[matchTable1[j]] >= part1size )
                    {
                        std:: cout <<"*";
                    }
                    else
                    {
                        std:: cout <<matchTable2[matchTable1[j]];
                    }
                }

                std::cout  << std::endl;
            }

            std:: cout <<" == Table 2: =="<<std::endl;
            for( size_t j = 0; j < matchTable2.size(); ++j )
            {
                std::cout << "2:" << j << " ("  << m_targetTree.getNode(partition2[j]).getSize() << ";"<< partColors2[j].m_x << "-" << partColors2[j].m_y << "-" << partColors2[j].m_z << ") -> 1:" ;

                if(matchTable2[j] >= part1size )
                {
                    std:: cout <<"*";
                }
                else
                {
                    std:: cout <<matchTable2[j];
                    std:: cout << " -> 2:";
                    if(matchTable1[matchTable2[j]] >= part2size )
                    {
                        std:: cout <<"*";
                    }
                    else
                    {
                        std:: cout <<matchTable1[matchTable2[j]];
                    }
                }

                std::cout  << std::endl;
            }
        }

        /// %%%%%% END DEBUG %%%%%%


        for( size_t j = 0; j < test1.size(); ++j )
        {
            if( !test1[j] )
            {
                std::cerr<< " WARNING: cluster " << j << " of partition 1 was not assigned." << std::endl;
            }
        }
        for( size_t j = 0; j < test2.size(); ++j )
        {
            if( !test2[j] )
            {
                std::cerr<< " WARNING: cluster " << j << " of partition 2 was not assigned." << std::endl;
            }
        }

        tree1partColors.push_back( partColors1 );
        tree2partColors.push_back( partColors2 );


    } // end partition for

    m_refTree.insertPartColors( tree1partColors );
    m_targetTree.insertPartColors( tree2partColors );

    return tree1changed;


} // end partitionMatcher::matchColors() -------------------------------------------------------------------------------------





void partitionMatcher::findMatchingPartitions( const float lambda, const size_t predefDepth )
{
    bool overlapMatching( lambda < 0);
    bool autoDepth( predefDepth == 0 );
    size_t levelDepth( 1 );

    if( m_verbose )
    {
        if( overlapMatching )
        {
            std::cout << "Bidireactional cluster overlap partition-matching." << std::endl;
        }
        else
        {
            std::cout << "Signature based partition-matching . Lambda: " << lambda << std::endl;
        }

        if( autoDepth )
        {
            std::cout << "level depth assigned automatically depending on size of partition" << std::endl;
        }
        else
        {
            std::cout << "fixed level depth for partition exploration: " << predefDepth << std::endl;
        }
    }

    if( !autoDepth )
    {
        levelDepth = predefDepth;
    }

    if( !m_targetTree.getSelectedValues().empty() )
    {
        std::cerr<< "WARNING @  partitionMatcher::findMatchingPartitions(): Tree 2 had partitions saved, they have been deleted."<<std::endl;
    }
    m_targetTree.clearPartitions();

    std::vector< std::vector< size_t > > tree1partitions( m_refTree.getSelectedPartitions() );
    std::vector< std::vector< size_t > > tree2partitions;
    tree2partitions.reserve( tree1partitions.size() );
    std::vector< float > tree2partValues;
    tree2partValues.reserve( tree1partitions.size() );



        /// %%%%%%%%% DEBUG %%%%%%%%%%%%
//    {
//        std::vector< std::vector< size_t > > smallpartitions;
//        std::vector< float > smallvalues;
//        std::vector< std::vector< WHcoord > > smallcolors;

//        for( size_t i = 0; i < tree1partitions.size(); ++i )
//        {

//            if ( tree1partitions[i].size() > 150 )
//            {
//                continue;
//            }
//            else
//            {
//                smallpartitions.push_back( tree1partitions[i] );
//                smallvalues.push_back( m_tree1.getSelectedValues( i ));
//                smallcolors.push_back( m_tree1.getSelectedColors( i ));
//            }
//        }
//        m_tree1.insertPartitions( smallpartitions , smallvalues, smallcolors );
//        tree1partitions.swap(smallpartitions);
//    }
        /// %%%%%%%%% END DEBUG %%%%%%%%%%%%


    std::cout << "Tree 1 has " << tree1partitions.size() << " saved partitions, finding matches in Tree 2... " << std::endl;

    for( size_t i = 0; i < tree1partitions.size(); ++i )
    {

        std::cout << "Getting best match for partition " << i << " in Tree 1 with " << tree1partitions[i].size() << " clusters.";
        if( autoDepth )
        {
            levelDepth = assignDepth( tree1partitions[i].size() );
            std::cout << " Search level depth: " << levelDepth;
        }
        std::cout << std::endl;


        std::vector< std::vector< bool > > tree1Signature;

        if( !overlapMatching )
        {
            tree1Signature = getSignatureMatrix( tree1partitions[i], TREE1 );
        }


        std::vector< size_t > lastPartition, keptPartition;
        double lastValue(0), keptValue;

        // do first step
        {
            const WHnode treeRoot( m_targetTree.getRoot() );
            lastPartition.push_back(treeRoot.getID());

            if( overlapMatching )
            {
                std::pair< std::vector< size_t >, std::vector< size_t > > matchSet1;
                std::pair< std::vector< size_t >, std::vector< size_t > > matchSet2;
                lastValue =  evalOverlapPartMatch( tree1partitions[i], lastPartition, &matchSet1, &matchSet2 );
            }
            else
            {
                std::vector< std::vector< bool > > tree2signature( getSignatureMatrix( lastPartition, TREE2 ) );
                lastValue = evalSignaturePartMatch( lambda, tree1partitions[i].size(), tree1Signature, lastPartition.size(),tree2signature );
            }
        }// end first step

        while ( true )
        {

            std::cout << "\rLast try: " << lastPartition.size() << " clusters. Match value: " << lastValue <<"            "<< std::flush;

            std::vector< std::vector< size_t > > derivedPartitionSet;
            std::vector< std::vector< unsigned int > > derivedIndexes( m_targetTree.getBranching(lastPartition, levelDepth, &derivedPartitionSet));

            std::vector< double > derivedPartitionValues( derivedPartitionSet.size(), -1 );
            bool stopLoop( true );

            // get the partitions and values for all possible branchings
            #pragma omp parallel for schedule(guided)
            for( size_t j = 0; j < derivedPartitionSet.size(); ++j )
            {
                // if there are possible branchings continue looping
                stopLoop = false;
                if( overlapMatching )
                {
                    std::pair< std::vector< size_t >, std::vector< size_t > > matchSet1;
                    std::pair< std::vector< size_t >, std::vector< size_t > > matchSet2;
                    derivedPartitionValues[j] =  evalOverlapPartMatch( tree1partitions[i], derivedPartitionSet[j], &matchSet1, &matchSet2 );
                }
                else
                {
                    std::vector< std::vector< bool > > derivedSignature( getSignatureMatrix( derivedPartitionSet[j], TREE2 ) );
                    derivedPartitionValues[j] =( evalSignaturePartMatch( lambda,  tree1partitions[i].size(), tree1Signature, derivedPartitionSet[j].size(), derivedSignature ) );
                }
            }// endFor

            //if there are no more possible branchings stop the loop
            if( stopLoop )
            {
                break;
            }

            size_t bestPartitionIndex;
            float bestValue( -1 );
            for( size_t j = 0; j < derivedPartitionValues.size(); ++j )
            {
                //if value is better than previous ones, keep partition and matrix
                if( derivedPartitionValues[j] > bestValue )
                {
                    bestValue = derivedPartitionValues[j];
                    bestPartitionIndex=j;
                }
            }

            // if no better partition was found and we went over the parttiion number, stop
            if( bestValue <= lastValue )
            {
                if( derivedPartitionSet[bestPartitionIndex].size() > ( tree1partitions[i].size() + ( tree1partitions[i].size() / 10 ) + 10 ) )
                {
                    break;
                }
            }

            // otherwise, continue

            // find the first branching partition corresponding to the best partition found
            if( derivedIndexes[bestPartitionIndex].size() == 1 )
            {
                // it was a first branch partition
                lastPartition = derivedPartitionSet[bestPartitionIndex];
            }
            else
            {
                unsigned int firstBranchIndex( derivedIndexes[bestPartitionIndex].front() );

                size_t nextBestIndex(0);
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

                lastPartition = derivedPartitionSet[nextBestIndex];
                bestValue = derivedPartitionValues[nextBestIndex];

            }

            if( bestValue > lastValue ) // a new best partition was found
            {
                keptValue = bestValue;
                keptPartition = lastPartition;
            }

            lastValue = bestValue;
        } // end infinite loop

        std::cout <<std::endl << "Best match in Tree 2 found to be one with: " << keptPartition.size() << " clusters and a partition distance of "<< keptValue << std::endl;

        std::sort(keptPartition.begin(),keptPartition.end());
        std::reverse(keptPartition.begin(),keptPartition.end());

        tree2partitions.push_back( keptPartition );
        tree2partValues.push_back( keptValue );

    } // end big for

    m_targetTree.insertPartitions( tree2partitions, tree2partValues );

} // end partitionMatcher::findMatchingPartitions() -------------------------------------------------------------------------------------



void partitionMatcher::testBaseNodes()
{
    m_refBaseNodes = m_refTree.getRootBaseNodes();
    m_targetBaseNodes = m_targetTree.getRootBaseNodes();

    if( m_refBaseNodes.empty() || m_refBaseNodes.empty() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::testBaseNodes(): base node vectors are empty" );
    }
    if( !m_refTree.testRootBaseNodes() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::testBaseNodes(): reference tree is not purely with meta-leaves" );
    }
    if( !m_targetTree.testRootBaseNodes() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::testBaseNodes(): target tree is not purely with meta-leaves" );
    }
    return;
} // end partitionMatcher::testBaseNodes() -------------------------------------------------------------------------------------



void partitionMatcher::loadCorrespondence( const std::string &matchFilename )
{
    // clear vectors
    {
        std::vector< size_t > emptyTable1;
        m_fullCorrespondence.swap( emptyTable1 );

        std::vector< size_t > emptyVector1;
        std::vector< size_t > emptyVector2;
        m_refMatchedBases.swap( emptyVector1 );
        m_targetMatchedBases.swap( emptyVector2 );

        std::vector< std::vector< size_t > > emptyList1;
        std::vector< std::vector< size_t > > emptyList2;
        m_refMatchedBases4Node.swap( emptyList1 );
        m_targetMatchedBases4Node.swap( emptyList2 );
    }

    WFileParser parser( matchFilename );
    if( !parser.readFile() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): Parser error" );
    }
    std::vector<std::string> lines = parser.getRawLines();
    if( lines.size() == 0 )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): File is empty" );
    }
    std::vector< std::vector< std::string> >matchStrings = parser.getLinesForTagSeparated( "correspondence" );
    if( matchStrings.size() == 0 )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): matching table was not found file" );
    }

    std::vector< std::pair < size_t, size_t > > fullIDtable;
    fullIDtable.reserve( matchStrings.size() );

    for( size_t i = 0; i < matchStrings.size(); ++i )
    {
        fullIDtable.push_back( std::make_pair( string_utils::fromString< size_t >( matchStrings[i][0] ), string_utils::fromString< size_t >( matchStrings[i][1] ) ) );
    }
    if( m_refBaseNodes.size() != fullIDtable.size()  )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): correspondance vector size does not match basenodes vector" );
    }
    m_fullCorrespondence.resize( fullIDtable.size(), 0 );

    for( size_t i = 0; i < fullIDtable.size(); ++i)
    {
        size_t rID1( findRelativeBasenodeID( fullIDtable[i].first, m_refBaseNodes ) );
        if( rID1 >= m_refBaseNodes.size())
        {
            throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): node from correspondence table was not found among tree 1 basenodes" );
        }
        size_t rID2( m_targetBaseNodes.size() + 1  );

        if( fullIDtable[i].second < m_targetTree.getNumNodes() )
        {
            rID2 = ( findRelativeBasenodeID( fullIDtable[i].second, m_targetBaseNodes ) );
            if( rID2 >= m_targetBaseNodes.size())
            {
                throw std::runtime_error( "ERROR @ partitionMatcher::loadMatchTable(): node from correspondence table was not found among tree 2 basenodes" );
            }
            m_refMatchedBases.push_back( fullIDtable[i].first );
            m_targetMatchedBases.push_back( fullIDtable[i].second );
        }
        m_fullCorrespondence[ rID1 ] = rID2;
    }


    // get a list of contained matched base nodes (relative IDs) for each node of each tree


    /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_refMatchedBases4Node.resize(m_refTree.getNumNodes());
    for( size_t i = 0; i < m_refMatchedBases.size(); ++i )
    {
        m_refMatchedBases4Node[m_refMatchedBases[i]].push_back(i);
    }
    for( size_t i = 0; i < m_refMatchedBases4Node.size(); ++i )
    {
        std::vector< nodeID_t > kids(m_refTree.getNode(i).getChildren());
        for( size_t j = 0; j < kids.size(); ++j )
        {
            if( kids[j].first )
            {
                m_refMatchedBases4Node[i].insert(m_refMatchedBases4Node[i].end(), m_refMatchedBases4Node[kids[j].second].begin(), m_refMatchedBases4Node[kids[j].second].end());
            }
        }
    }

    m_targetMatchedBases4Node.resize(m_targetTree.getNumNodes());
    for( size_t i = 0; i < m_targetMatchedBases.size(); ++i )
    {
        m_targetMatchedBases4Node[m_targetMatchedBases[i]].push_back(i);
    }
    for( size_t i = 0; i < m_targetMatchedBases4Node.size(); ++i )
    {
        std::vector< nodeID_t > kids(m_targetTree.getNode(i).getChildren());
        for( size_t j = 0; j < kids.size(); ++j )
        {
            if( kids[j].first )
            {
                m_targetMatchedBases4Node[i].insert(m_targetMatchedBases4Node[i].end(), m_targetMatchedBases4Node[kids[j].second].begin(), m_targetMatchedBases4Node[kids[j].second].end());
            }
        }
    }

    return;
} // end partitionMatcher::loadMatchTable() -------------------------------------------------------------------------------------


size_t partitionMatcher::findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t > &baseNodes ) const
{
    std::vector<size_t>::const_iterator findIter( std::find( baseNodes.begin(), baseNodes.end(), absoluteID ) );
    if( findIter == baseNodes.end() )
    {
        return baseNodes.size()+1;
    }
    else
    {
        return findIter - baseNodes.begin();
    }
} // end partitionMatcher::findRelativeBasenodeID() -------------------------------------------------------------------------------------

std::vector< std::vector< bool > > partitionMatcher::getSignatureMatrix( const std::vector< size_t> &partition, const bool forRefTree ) const
{
    std::vector< size_t >  membership( m_refMatchedBases.size(), 0 );
    const WHtree *treePointer;
    const std::vector< size_t > *baseNodePointer;

    if( forRefTree )
    {
        treePointer = &m_refTree;
        baseNodePointer = &m_refMatchedBases;
    }
    else
    {
        treePointer = &m_targetTree;
        baseNodePointer = &m_targetMatchedBases;
    }

    // for each partition cluster, find all their basenodes and update the basenode membership table
    const WHtree &tree( *treePointer );
    std::vector< size_t > baseNodes( *baseNodePointer );
    size_t doneBaseCount( 0 );
    for( size_t i = 0; i < partition.size(); ++i )
    {
        std::list <size_t> worklist;
        worklist.push_back( partition[i] );
        while( !worklist.empty() )
        {
            const WHnode currentNode( tree.getNode( worklist.front() ) );
            worklist.pop_front();
            if( currentNode.getHLevel() == 1 )
            {
                size_t relativeID;
                relativeID = findRelativeBasenodeID( currentNode.getID(), baseNodes );
                if( relativeID < baseNodes.size() )
                {
                    membership[relativeID]=i;
                    ++doneBaseCount;
                }
            }
            else
            {
                std::vector<nodeID_t> kids( currentNode.getChildren() );
                for( size_t j = 0; j < kids.size(); ++j )
                {
                    const WHnode thisKid( tree.getNode( kids[j] ) );
                    if( thisKid.isNode() )
                    {
                        worklist.push_back( thisKid.getID() );
                    }
                }
            }
        }
    }
    if( doneBaseCount != baseNodes.size() )
    {
        std::cerr << "doneBaseCount: "<< doneBaseCount <<". baseNodeSize: "<< baseNodes.size()<<std::endl;
        throw std::runtime_error( "ERROR @ partitionMatcher::getBaseNodeMembership(): not all bases where assigned membership" );
    }

    std::vector< std::vector< bool > > signature;
    signature.reserve(membership.size());
    for( size_t i = 0; i < membership.size(); ++i )
    {
        std::vector< bool > signatureLine;
        signatureLine.reserve(i);
        for( size_t j = 0; j < i; ++j )
        {
            signatureLine.push_back( membership[i]==membership[j]);
        }
        signature.push_back(signatureLine);
    }
    return signature;
} // end partitionMatcher::getBaseNodeMembership() -------------------------------------------------------------------------------------


float partitionMatcher::evalSignaturePartMatch(  const float lambda, size_t refSize, std::vector< std::vector< bool > >& refSignature,  size_t targetSize, std::vector< std::vector< bool > >& targetSignature ) const
{
    if( refSignature.size() != targetSignature.size() )
     {
        throw std::runtime_error( "ERROR @ partitionMatcher::evalPartitionMatch(): signature matrices do not have the same size" );
    }

    double sum1( 0 ), sum2( 0 ), sumProd( 0 );
    double M ((refSignature.size() * (refSignature.size() -1)) / 2.0);


    for( size_t i = 0; i < refSignature.size(); ++i )
    {
        for( size_t j = 0; j < i; ++j )
        {
            sum1 += refSignature[i][j];
            sum2 += targetSignature[i][j];
            sumProd += refSignature[i][j] * targetSignature[i][j];
        }
    }

    // do the final computations
    double mean1( sum1 / M );
    double mean2( sum2 / M );
    double numerator( ( sumProd / M ) - ( mean2 * mean1 ) );
    double denominator1( ( mean1 ) * ( 1 - mean1 ) );
    double denominator2( ( mean2 ) * ( 1 - mean2 ) );
    double matchDist( numerator / sqrt( denominator1 * denominator2 ) );

    double bigPart(std::max(refSize,targetSize));
    double smallPart(std::min(refSize,targetSize));

    double final = matchDist + (lambda*(smallPart/bigPart));

    return final;
} // end partitionMatcher::evalPartitionMatch() -------------------------------------------------------------------------------------

size_t partitionMatcher::getClusterOverlap( const size_t cluster1, const size_t cluster2 ) const
{
    size_t matchedLeaves( 0 );

    std::vector< size_t > bases1 (m_refMatchedBases4Node[cluster1]);
    std::vector< size_t > bases2 (m_targetMatchedBases4Node[cluster2]);


    for( size_t i = 0; i < bases1.size(); ++i )
    {
        if( std::find(bases2.begin(),bases2.end(),bases1[i]) != bases2.end() )
        {
            ++matchedLeaves;
            //matchedLeaves += tree.getNode( matchedBases[basesA[i]] ).getSize();
        }
    }
    return matchedLeaves;
} // end partitionMatcher::getClusterMatch() -------------------------------------------------------------------------------------




std::pair< std::vector< std::vector< size_t > >, std::vector< std::vector< size_t > > > partitionMatcher::getClusterOverlapMatrix(  const std::vector< size_t > &partition1, const std::vector< size_t > &partition2) const
{

    std::vector< std::vector< size_t > > clusterMatchMatrix1( partition1.size(), std::vector< size_t >( partition2.size(), 0 ) );
    std::vector< std::vector< size_t > > clusterMatchMatrix2( partition2.size(), std::vector< size_t >( partition1.size(), 0 ) );

    #pragma omp parallel for
    for( size_t i = 0; i < partition1.size(); ++i )
    {
        for( size_t j = 0; j < partition2.size(); ++j )
        {
            size_t thisMatch( getClusterOverlap( partition1[i], partition2[j] ) );
            clusterMatchMatrix1[i][j] =  thisMatch;
            clusterMatchMatrix2[j][i] =  thisMatch;

        }
    }
    return std::make_pair( clusterMatchMatrix1, clusterMatchMatrix2 );

} // end partitionMatcher::getClusterMatchMatrix() -------------------------------------------------------------------------------------

std::pair< std::vector< size_t >, std::vector< size_t > > partitionMatcher::getClusterMatchTable( const std::vector< std::vector< size_t > > &matchMatrix ) const
{
    std::vector< size_t > matchTable;
    std::vector< size_t > matchValues;
    matchTable.reserve( matchMatrix.size() );
    matchValues.reserve( matchMatrix.size() );


    for( size_t i = 0; i < matchMatrix.size(); ++i )
    {
        std::vector< size_t >::const_iterator maxIter( std::max_element( matchMatrix[i].begin(), matchMatrix[i].end() ) );
        float maxValue( *maxIter );
        if ( maxValue > 0 )
        {
            matchTable.push_back(  maxIter - matchMatrix[i].begin() );
        }
        else
        {
            matchTable.push_back( matchMatrix[i].size() + 1 );
        }
        matchValues.push_back(maxValue);
    }

    return std::make_pair( matchTable, matchValues );
} // end partitionMatcher::getMatchTable() -------------------------------------------------------------------------------------

double partitionMatcher::evalOverlapPartMatch( const std::vector< size_t > &partition1, const std::vector< size_t > &partition2,
                                    std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint1,
                                    std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint2 ) const
{

    std::pair< std::vector< std::vector< size_t > >, std::vector< std::vector< size_t > > > clusterMatchMatrices( getClusterOverlapMatrix( partition1, partition2 ) );
    std::vector< std::vector< size_t > > &clusterMatchMatrix1( clusterMatchMatrices.first );
    std::vector< std::vector< size_t > > &clusterMatchMatrix2( clusterMatchMatrices.second );

    if( clusterMatchMatrix1.empty() || clusterMatchMatrix1.empty() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::doPartMatch(): matrices are empty");
    }

    if( clusterMatchMatrix1.size() != clusterMatchMatrix2.front().size() || clusterMatchMatrix1.front().size() != clusterMatchMatrix2.size() )
    {
        throw std::runtime_error( "ERROR @ partitionMatcher::doPartMatch(): matrices have wrong dimensions");
    }

    /// DEBUG///
    if( DEBUG )
    {
        std::cout<< std::endl<<"Matrix 1->2 " << std::flush;
        for( size_t i = 0; i < clusterMatchMatrix1.size(); ++i )
        {
            std::cout<< std::endl;

            for( size_t j = 0; j < clusterMatchMatrix1[i].size(); ++j )
            {
                std::cout<< clusterMatchMatrix1[i][j]<<" ";
            }
            std::cout<< std::endl;
        }
        std::cout<<  std::endl<<"Matrix 2->2 " << std::flush;
        for( size_t i = 0; i < clusterMatchMatrix2.size(); ++i )
        {
            std::cout<< std::endl;

            for( size_t j = 0; j < clusterMatchMatrix2[i].size(); ++j )
            {
                std::cout<< clusterMatchMatrix2[i][j]<<" ";
            }
            std::cout<< std::endl;
        }
    }
    /// END DEBUG///

    std::pair< std::vector< size_t >, std::vector< size_t > > matchSet1( getClusterMatchTable( clusterMatchMatrix1) );
    std::pair< std::vector< size_t >, std::vector< size_t > > matchSet2( getClusterMatchTable( clusterMatchMatrix2) );
    std::vector< size_t > &matchTable1( matchSet1.first );
    std::vector< size_t > &matchTable2( matchSet2.first );
    std::vector< size_t > &matchValues1( matchSet1.second );
    std::vector< size_t > &matchValues2( matchSet2.second );

    std::vector< size_t > &finalMatchTable1( tablePoint1->first );
    std::vector< size_t > &finalMatchTable2( tablePoint2->first );
    finalMatchTable1 = ( matchTable1 );
    finalMatchTable2 = ( matchTable2 );

    std::vector< size_t > &finalMatchValues1( tablePoint1->second );
    std::vector< size_t > &finalMatchValues2( tablePoint2->second );
    finalMatchValues1 = ( matchValues1 );
    finalMatchValues2 = ( matchValues2 );

    std::vector< bool > checkTable1( matchTable1.size(), false );
    std::vector< bool > checkTable2( matchTable2.size(), false );

    size_t addedSize(0);
    size_t totalMatches(0);


    // first sweep (non matched clusters
    for( size_t i = 0; i < matchTable1.size(); ++i )
    {
        if( matchTable1[i] >= matchTable2.size() )
        {
            checkTable1[i] = true;
        }
    }
    for( size_t i = 0; i < matchTable2.size(); ++i )
    {
        if( matchTable2[i] >= matchTable1.size() )
        {
            checkTable2[i] = true;
        }
    }


    bool doContinue(true);
    bool isFirstRound(true);
    size_t reDone(0);


    while(doContinue)
    {
        doContinue = false;

        // first round
        for( size_t i = 0; i < matchTable1.size(); ++i )
        {
            if( checkTable1[i] ) // element was already matched (or had no match)
            {
                continue;
            }
            size_t match2( matchTable1[i] );
            size_t match1( matchTable2[match2] );

            if( match1 == i )   // tables coincide
            {
                if( checkTable2[match2] ) // element was already matched (or had no match)
                {
                    throw std::runtime_error( "ERROR @ partitionMatcher::doPartMatch(): this should not happen");
                }

                if( matchValues1[match1] != matchValues2[match2] )
                {
                    throw std::runtime_error( "ERROR @ partitionMatcher::getPartMatchIndex(): match values differ" );
                }
                addedSize += matchValues1[match1];
                ++totalMatches;
                finalMatchTable1[match1]=match2;
                finalMatchTable2[match2]=match1;
                finalMatchValues1[match1] = matchValues1[match1];
                finalMatchValues2[match2] = matchValues1[match1];
                checkTable1[match1] = true;
                checkTable2[match2] = true;
                if( !isFirstRound )
                {
                    ++reDone;
                }
            }
        }

        isFirstRound = false;

        // adjust round for table 1
        for( size_t i = 0; i < matchTable1.size(); ++i )
        {
            if( checkTable1[i] ) // element was already matched (or had no match)
            {
                continue;
            }
            else if( checkTable2[matchTable1[i]] ) // the element it is pointing at was matched
            {
                size_t newBestMatch(0);
                size_t newBestValue(0);
                // get the next best matching
                for( size_t j = 0; j<clusterMatchMatrix1[i].size(); ++j)
                {
                    if( checkTable2[j])
                    {
                        continue;
                    }
                    else
                    {
                        if( clusterMatchMatrix1[i][j] > newBestValue )
                        {
                            newBestValue = clusterMatchMatrix1[i][j];
                            newBestMatch = j;
                        }
                    }
                }

                if ( newBestValue > (finalMatchValues1[i]/2 ) ) // it found another possible match (must be at least half of the original best match
                {
                    matchValues1[i] = newBestValue;
                    matchTable1[i] = newBestMatch;
                    doContinue = true;
                }
            }
            else // it is a chain, gets messy, forget it
            {
                continue;
            }
        }

        // adjust round for table 2
        for( size_t i = 0; i < matchTable2.size(); ++i )
        {
            if( checkTable2[i] ) // element was already matched (or had no match)
            {
                continue;
            }
            else if( checkTable1[matchTable2[i]] ) // the element it is pointing at was matched
            {
                size_t newBestMatch(0);
                size_t newBestValue(0);
                // get the next best matching
                for( size_t j = 0; j<clusterMatchMatrix2[i].size(); ++j)
                {
                    if( checkTable1[j])
                    {
                        continue;
                    }
                    else
                    {
                        if( clusterMatchMatrix2[i][j] > newBestValue )
                        {
                            newBestValue = clusterMatchMatrix2[i][j];
                            newBestMatch = j;
                        }
                    }
                }

                if ( newBestValue > (finalMatchValues2[i]/2 ) ) // it found another possible match (must be at least half of the original best match
                {
                    matchValues2[i] = newBestValue;
                    matchTable2[i] = newBestMatch;
                    doContinue = true;
                }
            }
            else // it is a chain, gets messy, forget it
            {
                continue;
            }
        }
    }
    double addedSizeDouble( addedSize );
    double quality(addedSizeDouble / m_refMatchedBases.size());

//    std::cout << "Q: "<< quality <<". S: " << partition2.size() << ". Total matches: " << totalMatches;
//    if (  reDone > 0)
//    {
//        std::cout << " ("<< reDone << " fixed)";

//    }
//    std::cout << std::endl;





    return quality;

} // end partitionMatcher::getPartMatchIndex() -------------------------------------------------------------------------------------


size_t partitionMatcher::assignDepth( const size_t partSize )
{
    if( partSize < 40 )
    {
        return 5;
    }
    else if ( partSize < 90 )
    {
        return 4;
    }
    else if ( partSize < 200 )
    {
        return 3;
    }
    else if ( partSize < 350 )
    {
        return 2;
    }
    else
    {
        return 1;
    }
} // end partitionMatcher::assignDepth() -------------------------------------------------------------------------------------


WHcoord partitionMatcher::shiftColor( const WHcoord &color, const size_t shiftIndex ) const
{
    WHcoord outColor( color );
    size_t shiftedCoord( shiftIndex%3 );
    int amount( 30 * ( (shiftIndex / 3) +1  )  );
    coord_t value1(0), value2(0);

    switch( shiftedCoord )
    {
    case 0:
        value1 = outColor.m_x;
        value2 = outColor.m_y;

        break;
    case 1:
        value1 = outColor.m_y;
        value2 = outColor.m_z;

        break;
    case 2:
        value1 = outColor.m_x;
        value2 = outColor.m_z;
        break;
    default:
        break;
    }

    if( value1 >= 128 )
    {
        value1 -= amount;
    }
    else
    {
        value1 += amount;
    }

    if( value2 >= 128 )
    {
        value2 -= amount;
    }
    else
    {
        value2 += amount;
    }


    switch( shiftedCoord )
    {
    case 0:
        outColor.m_x = value1;
        outColor.m_y = value2;
        break;
    case 1:
        outColor.m_y = value1;
        outColor.m_z = value2;
        break;
    case 2:
        outColor.m_x = value1;
        outColor.m_z = value2;
        break;
    default:
        break;
    }

    return outColor;
} // end partitionMatcher::shiftColor() -------------------------------------------------------------------------------------


void partitionMatcher::searchPartition( const size_t maxPartSize, const std::vector< size_t > &thisPart, std::vector< std::vector< size_t > >* partVectorPtr)
{
    std::vector< std::vector< size_t > >& partVector( *partVectorPtr );
    if (thisPart.size() >= maxPartSize)
    {
        partVector.push_back( thisPart );
        return;
    }
    for( size_t j = 0; j < thisPart.size(); ++j )
    {
        std::vector< size_t > partBranched(thisPart);

        {

            const WHnode thisBranch( m_targetTree.getNode( thisPart[j] ) );
            // if this branch is a base node, dont consider it
            if( thisBranch.getHLevel() == 1 )
            {
                continue;
            }
            else
            {

                // get nodes of branching for this branch
                std::vector< size_t > branchNodes;
                {
                    std::vector< nodeID_t > kids( thisBranch.getChildren() );
                    for( size_t k = 0; k < kids.size(); ++k )
                    {
                        if( kids[k].first )
                        {
                            branchNodes.push_back( kids[k].second );
                        }
                    }
                }
                // get the corresponding partition considering this branching
                partBranched.erase( partBranched.begin() + j );
                partBranched.insert( partBranched.begin() + j, branchNodes.begin(), branchNodes.end() );
            }
        }
        searchPartition( maxPartSize, partBranched, &partVector );
    } // endFor
    return;
}
