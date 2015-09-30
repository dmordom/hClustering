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

#include "treeComparer.h"



treeComparer::treeComparer( WHtree* const tree1, WHtree* const tree2, bool verbose ) :
    m_tree1( *tree1 ),
    m_tree2( *tree2 ),
    m_maxPhysDist( 0 ),
    m_tractThreshold1( 0 ),
    m_tractThreshold2( 0 ),
    m_logfile( 0 ),
    m_realBaseNodes( false ),
    m_coordsFromFile( false ),
    m_meanTractsFromFile( false ),
    m_verbose( verbose )
{
    if( m_tree1.m_logFactor != m_tree2.m_logFactor )
    {
        std::cerr << "WARNING @ treeComparer::treeComparer: trees have different track log factors: tree 1: " << m_tree1.m_logFactor << ". tree 2: " << m_tree2.m_logFactor << std::endl;
        if( m_tree1.m_numStreamlines ==  m_tree2.m_numStreamlines )
        {
            if( m_tree1.m_numStreamlines != 0 )
            {
                m_tree1.m_logFactor = log10( m_tree1.m_numStreamlines );
                m_tree2.m_logFactor = log10( m_tree2.m_numStreamlines );
                std::cerr << "Trees have same tracking streamlines number. Recomputed track log factors to: tree 1: " << m_tree1.m_logFactor << ". tree 2: " << m_tree2.m_logFactor << std::endl;
            }
            else
            {
                m_tree1.m_logFactor = 0;
                m_tree2.m_logFactor = 0;
                std::cerr << "Both trees have tracking streamlines number = 0. Resetting both track log factors to 0" << std::endl;
            }
        }
        else if( m_tree1.m_logFactor == 0 )
        {
            if(  m_tree1.m_numStreamlines != 0 )
            {
                m_tree1.m_logFactor = log10( m_tree1.m_numStreamlines );
                std::cerr << "Reset tree 1 track log factor to: " << m_tree1.m_logFactor << std::endl;
            }
            else
            {
                m_tree2.m_logFactor = 0;
                std::cerr << "Reset tree 2 track log factor to: 0" << std::endl;
            }

        }
        else if( m_tree2.m_logFactor == 0)
        {
            if(  m_tree2.m_numStreamlines != 0 )
            {
                m_tree2.m_logFactor = log10( m_tree2.m_numStreamlines );
                std::cerr << "Reset tree 2 track log factor to: " << m_tree2.m_logFactor << std::endl;
            }
            else
            {
                m_tree1.m_logFactor = 0;
                std::cerr << "Reset tree 1 track log factor to: 0" << std::endl;
            }
        }
    }
    fetchBaseNodes( false );
    m_initialSizes = std::make_pair(m_tree1.getNumLeaves(), m_tree2.getNumLeaves());
}
treeComparer::treeComparer( WHtree* const tree1, WHtree* const tree2, treeComparer &comparer ):
  m_tree1(*tree1),
  m_tree2(*tree2),
  m_singleTractFolder1(comparer.m_singleTractFolder1),
  m_singleTractFolder2(comparer.m_singleTractFolder2),
  m_meanTractFolder1(comparer.m_meanTractFolder1),
  m_meanTractFolder2(comparer.m_meanTractFolder2),
  m_maxPhysDist(comparer.m_maxPhysDist),
  m_tractThreshold1(comparer.m_tractThreshold1),
  m_tractThreshold2(comparer.m_tractThreshold2),
  m_logfile(comparer.m_logfile),
  m_baseNodes1(comparer.m_baseNodes1),
  m_originalBaseNodes1(comparer.m_originalBaseNodes1),
  m_baseCoords1(comparer.m_baseCoords1),
  m_noiseLevels1(comparer.m_noiseLevels1),
  m_baseNodes2(comparer.m_baseNodes2),
  m_originalBaseNodes2(comparer.m_originalBaseNodes2),
  m_baseCoords2(comparer.m_baseCoords2),
  m_noiseLevels2(comparer.m_noiseLevels2),
  m_initialSizes(comparer.m_initialSizes),
  m_baseDistMatrix(comparer.m_baseDistMatrix),
  m_realBaseNodes(comparer.m_realBaseNodes),
  m_coordsFromFile(comparer.m_coordsFromFile),
  m_meanTractsFromFile(comparer.m_meanTractsFromFile),
  m_fullCorrespondence(comparer.m_fullCorrespondence),
  m_newCorrespondence(comparer.m_newCorrespondence),
  m_newCorrespReverse(comparer.m_newCorrespReverse),
  m_correspDistances(comparer.m_correspDistances),
  m_verbose( comparer.m_verbose )
{}
// PUBLIC member functions

void treeComparer::setRelativeThreshold( float thresholdRatio )
{
    if( thresholdRatio <= 0 || thresholdRatio >= 1 )
    {
        m_tractThreshold1 = 0;
        m_tractThreshold2 = 0;
        if( thresholdRatio != 0 )
        {
            std::cerr << "WARNING @ treeComparer::setRelativeThreshold(): threshold ratio provided (" << thresholdRatio << ") is out of bounds [0,1), using a value of 0.0 (no thresholding)" << std::endl;
        }
    }
    else
    {
        if ( m_tree1.m_logFactor == 0 )
        {
            m_tractThreshold1 = thresholdRatio; // if using natural units the normalized tract threshold is the same as the threshold ratio
        }
        else
        {
            m_tractThreshold1 = log10( m_tree1.m_numStreamlines * thresholdRatio ) / m_tree1.m_logFactor; // if using log units the threshold must be calculated
        }

        if ( m_tree2.m_logFactor == 0 )
        {
            m_tractThreshold2 = thresholdRatio; // if using natural units the normalized tract threshold is the same as the threshold ratio
        }
        else
        {
            m_tractThreshold2 = log10( m_tree2.m_numStreamlines * thresholdRatio ) / m_tree2.m_logFactor; // if using log units the threshold must be calculated
        }
    }
    if( m_verbose )
    {
        std::cout << "Relative threshold: " << thresholdRatio << ". Normalized threshold 1: " << m_tractThreshold1 << ". Normalized threshold 2: " << m_tractThreshold2 <<std::endl;
    }
    return;
}


std::pair< float, float >treeComparer::simpleTriplets( size_t sampleFreq ) const
{
    if( m_baseNodes1.size() != m_baseNodes2.size() )
    {
        throw std::runtime_error( "ERROR @ treeCompare::simpleTriplets(): base node vectors have different sizes" );
    }

    bool modeNodes( false );
    size_t loopLength( 0 );

    if( m_baseNodes1.empty() )
    {
        if( m_verbose )
        {
            std::cout << "Computing leaf-wise simple triplets comparison..." << std::endl;
        }
        if( m_tree1.getNumLeaves() != m_tree2.getNumLeaves() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::simpleTriplets(): trees have different sizes" );
        }

        if( m_tree1.m_coordinates != m_tree2.m_coordinates )
        {
            std::cerr << ( "WARNING @ treeCompare::simpleTriplets(): trees have different coordinates" ) << std::endl;
        }

        modeNodes = false;
        loopLength = m_tree1.getNumLeaves();
    }
    else
    {
        if( m_verbose )
        {
            std::cout << "Computing baseNode-wise simple triplets comparison..." << std::endl;
        }
        if( m_baseNodes1.size() != m_newCorrespondence.size() )
        {
            throw std::runtime_error(
                            "ERROR @ treeCompare::simpleTriplets(): correspondance vector size does not match basenodes vector" );
        }

        modeNodes = true;
        loopLength = m_newCorrespondence.size();
    }

    if( m_verbose && sampleFreq != 1 )
    {
        std::cout << "Subsampling frequency: " << sampleFreq << std::endl;
        loopLength = ( ( size_t )( loopLength / sampleFreq ) ) * sampleFreq;
    }

    std::vector< size_t > tripletSumVect( loopLength - 2, 0 );
    std::vector< size_t > tripletWeightedSumVect( loopLength - 2, 0 );
    std::vector< size_t > sizeSumVector( loopLength - 2, 0 );


    size_t doneCount( 0 );
    time_t loopStartTime( time( NULL ) ), lastTime( time( NULL ) );
    double totalTriplets( boost::math::binomial_coefficient< double >( ( loopLength / sampleFreq ), 3 ) );

    // loop through all possible triplets
#pragma omp parallel for schedule( dynamic, 100 )
    for( size_t i = 0; i < loopLength - 2; i += sampleFreq )
    {
        unsigned int result1( 0 ), result2( 0 );
        size_t sizeElement( 0 );

        for( size_t j = i + sampleFreq; j < loopLength; j += sampleFreq )
        {
            for( size_t k = j + sampleFreq; k < loopLength; k += sampleFreq )
            {
                if( modeNodes )
                {
                    result1 = m_tree1.getTripletOrder( std::make_pair( true, m_baseNodes1[i] ),
                                                       std::make_pair( true, m_baseNodes1[j] ),
                                                       std::make_pair( true, m_baseNodes1[k] ) );
                    result2 = m_tree2.getTripletOrder( std::make_pair( true, m_baseNodes2[m_newCorrespondence[i]] ),
                                                       std::make_pair( true, m_baseNodes2[m_newCorrespondence[j]] ),
                                                       std::make_pair( true, m_baseNodes2[m_newCorrespondence[k]] ) );

                    size_t size1 ( ( m_tree1.getNode( m_baseNodes1[i] ).getSize() ) +
                                   ( m_tree1.getNode( m_baseNodes1[j] ).getSize() ) +
                                   ( m_tree1.getNode( m_baseNodes1[k] ).getSize() ) );
                    size_t size2 ( ( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[i]] ).getSize() ) +
                                   ( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[j]] ).getSize() ) +
                                   ( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[k]] ).getSize() ) );
                    sizeElement = size1 + size2;
                }
                else
                {
                    result1 = m_tree1.getTripletOrder( m_tree1.getLeaf( i ).getFullID(),
                                                       m_tree1.getLeaf( j ).getFullID(),
                                                       m_tree1.getLeaf( k ).getFullID() );
                    result2 = m_tree2.getTripletOrder( m_tree2.getLeaf( i ).getFullID(),
                                                       m_tree2.getLeaf( j ).getFullID(),
                                                       m_tree2.getLeaf( k ).getFullID() );
                    sizeElement = 6;
                }

                sizeSumVector[i] += sizeElement;
                if( result1 == result2 )
                {
                    tripletWeightedSumVect[i] += sizeElement;
                    ++tripletSumVect[i];
                }

#pragma omp atomic
                ++doneCount;
            } // end for (k)

#pragma omp single nowait // only one thread executes output
            if( m_verbose )
            {
                time_t currentTime( time( NULL ) );
                if( currentTime - lastTime > 1 )
                {
                    lastTime = currentTime;
                    size_t localCount = doneCount;
                    float progress = ( localCount * 100. ) / ( totalTriplets );
                    size_t elapsedTime( difftime( currentTime, loopStartTime ) );
                    std::stringstream message;
                    message << "\r" << ( int )progress << " % completed. Expected remaining time: ";
                    if( progress > 0 )
                    {
                        size_t expected_remain( elapsedTime * ( ( 100. - progress ) / progress ) );
                        message << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' "
                                        << ( ( expected_remain % 3600 ) % 60 ) << "\". ";
                    }
                    message << "Elapsed time: ";
                    message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' "
                                    << ( ( elapsedTime % 3600 ) % 60 ) << "\". ";
                    std::cout << message.str() << std::flush;
                }
            } // end m_verbose
        } // end for (j)
    } // end parallel for (i)

    if( m_verbose )
    {
        std::cout << "\r100 % completed. Total triples: " << doneCount << "     " << std::endl;
    }

    if( totalTriplets != ( ( double )doneCount ) )
    {
        if( sampleFreq == 1 )
        {
            std::cerr << "Total triples by count: " << doneCount << std::endl;
            std::cerr << "Total triples by formula: " << ( size_t )totalTriplets << std::endl;
            throw std::runtime_error(
                            "ERROR @ treeCompare::simpleTriplets(): theoretical and calculated number of triplets dont match" );
        }
        else if( m_verbose )
        {
            std::cerr << "Total triples by count: " << doneCount << std::endl;
            std::cerr << "Total triples by formula: " << ( size_t )totalTriplets << std::endl;
            std::cerr << "WARNING @ treeCompare::simpleTriplets(): theoretical and calculated number of triplets dont match" << std::endl;
        }
    }

    double compareSum( 0 ), weightedSum( 0 ), sizeSum( 0 );
    for( size_t i = 0; i < tripletSumVect.size(); ++i )
    {
        compareSum += tripletSumVect[i];
        weightedSum += tripletWeightedSumVect[i];
        sizeSum += sizeSumVector[i];
    }

    double tripletsCoef( ( compareSum ) / ( double )doneCount );
    double weightedTripletsCoef( ( weightedSum ) / sizeSum );


    if( m_verbose )
    {
        std::cout << "unweighted STC: " << boost::lexical_cast< std::string >( tripletsCoef ) << std::endl;
        std::cout << "size-weighted STC: " << boost::lexical_cast< std::string >( weightedTripletsCoef ) << std::endl;
    }

    return std::make_pair( tripletsCoef, weightedTripletsCoef );
} // end treeComparer::simpleTriplets() -------------------------------------------------------------------------------------


std::pair< std::pair< float, float >, std::pair< float, float > > treeComparer::doTcpcc() const
{
    if( m_baseNodes1.size() != m_baseNodes2.size() )
    {
        throw std::runtime_error( "ERROR @ treeCompare::doCpct(): base node vectors have different sizes" );
    }
    if( !m_noiseLevels1.empty() )
    {
        if( m_noiseLevels1.size() != m_noiseLevels2.size() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::doCpct(): noise level vectors have different sizes" );
        }
        if( m_noiseLevels1.size() != m_baseNodes1.size() || m_noiseLevels2.size() != m_baseNodes2.size() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::doCpct(): noise level vectors have different sizes to base node vectors" );
        }
    }

    bool modeNodes( false );
    size_t loopLength( 0 );

    if( m_baseNodes1.empty() )
    {
        if( m_verbose )
        {
            std::cout << "Computing leaf-wise tree cophenetic correlation comparison..." << std::endl;
        }
        if( m_tree1.getNumLeaves() != m_tree2.getNumLeaves() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::doCpct(): trees have different sizes" );
        }
        if( m_tree1.m_coordinates != m_tree2.m_coordinates )
        {
            std::cerr << "WARNING @ treeManager::doCpct: Trees do not have the same seed voxels" << std::endl;
        }

        modeNodes = false;
        loopLength = m_tree1.getNumLeaves();
    }
    else
    {
        if( m_verbose )
        {
            std::cout << "Computing baseNode-wise cophenetic correlation comparison..." << std::endl;
        }
        if( m_baseNodes1.size() != m_newCorrespondence.size() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::doCpct(): correspondance vector size does not match basenodes vector" );
        }

        modeNodes = true;
        loopLength = m_newCorrespondence.size();
    }

    const size_t loopLengthMinusOne( loopLength - 1 );

    std::vector< double > sumT1vec( loopLengthMinusOne, 0 ), sumT2vec( loopLengthMinusOne, 0 ), sqT1vec( loopLengthMinusOne, 0 ),
                    sqT2vec( loopLengthMinusOne, 0 ), sumProdvec( loopLengthMinusOne, 0 );
    std::vector< double > wSumT1vec( loopLengthMinusOne, 0 ), wSumT2vec( loopLengthMinusOne, 0 ), wSqT1vec( loopLengthMinusOne, 0 ),
                    wSqT2vec( loopLengthMinusOne, 0 ), wSumProdvec( loopLengthMinusOne, 0 );
    std::vector< size_t > sumSize1Vec( loopLengthMinusOne, 0 ), sumSize2Vec( loopLengthMinusOne, 0 ), sumSizeProdVec( loopLengthMinusOne, 0 );
    std::vector< size_t > sumSqSize1Vec( loopLengthMinusOne, 0 ), sumSqSize2Vec( loopLengthMinusOne, 0 );

    size_t totalPairs( loopLength * loopLengthMinusOne / 2.0 );
    double N2( loopLength * loopLength );
    std::vector< size_t > usedPairsVect( loopLengthMinusOne, 0 );


    size_t doneCount( 0 );
    time_t loopStartTime( time( NULL ) ), lastTime( time( NULL ) );

#pragma omp parallel for schedule( dynamic, 100 )
    for( size_t i = 0; i < loopLengthMinusOne; ++i )
    {
        size_t localCount( 0 );


        for( size_t j = i + 1; j < loopLength; ++j )
        {
            double dist1( 0 ), dist2( 0 );
            size_t sizeComp1( 0 ), sizeComp2( 0 );

            if( modeNodes )
            {
                dist1 = m_tree1.getDistance( m_baseNodes1[i], m_baseNodes1[j] );
                dist2 = m_tree2.getDistance( m_baseNodes2[m_newCorrespondence[i]], m_baseNodes2[m_newCorrespondence[j]] );

                sizeComp1 = ( m_tree1.getNode( m_baseNodes1[i] ).getSize() ) +
                            ( m_tree1.getNode( m_baseNodes1[j] ).getSize() );
                sizeComp2 = ( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[i]] ).getSize() ) +
                            ( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[j]] ).getSize() );
            }
            else
            {
                dist1 = m_tree1.getLeafDistance( i, j );
                dist2 = m_tree2.getLeafDistance( i, j );
                sizeComp1 = 2;
                sizeComp2 = 2;
            }

//            float avrgDist( ( dist1 + dist2 ) / 2 );
//            float avrgMatch( ( m_correspDistances[i].first + m_correspDistances[j].first ) / 2 );
//            if( avrgDist <= avrgMatch )
//            {
//                continue;
//            }

            // test if the tree values being taken into account are higher than the distance between matched leaves, otherwise pair is not used
            if( !m_noiseLevels1.empty() && ( ( dist1 <= m_noiseLevels1[i] ) || ( dist1 <= m_noiseLevels1[j] ) ) )
            {
                continue;
            }
            if( !m_noiseLevels2.empty() && ( ( dist2 <= m_noiseLevels2[i] ) || ( dist2 <= m_noiseLevels2[j] ) ) )
            {
                continue;
            }



            sumT1vec[i] += dist1;
            sumT2vec[i] += dist2;
            sqT1vec[i] += ( dist1 * dist1 );
            sqT2vec[i] += ( dist2 * dist2 );
            sumProdvec[i] += ( dist1 * dist2 );

            double wDist1( dist1 * sizeComp1 );
            double wDist2( dist2 * sizeComp2 );
            wSumT1vec[i] += wDist1;
            wSumT2vec[i] += wDist2;
            wSqT1vec[i] += ( wDist1 * wDist1 );
            wSqT2vec[i] += ( wDist2 * wDist2 );
            wSumProdvec[i] += ( wDist1 * wDist2 );
            sumSize1Vec[i] += sizeComp1;
            sumSize2Vec[i] += sizeComp2;
            sumSqSize1Vec[i] += ( sizeComp1 * sizeComp1 );
            sumSqSize2Vec[i] += ( sizeComp2 * sizeComp2 );
            sumSizeProdVec[i] += ( sizeComp1 * sizeComp2 );
            usedPairsVect[i] += 1;
        } // end inner for (j)



#pragma omp critical
        {
            doneCount += ( loopLength - i - 1 );
            localCount = doneCount;
        }

#pragma omp single nowait // only one thread executes output
        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                lastTime = currentTime;

                float progress = ( localCount * 100. ) / ( totalPairs );
                size_t elapsedTime( difftime( currentTime, loopStartTime ) );
                std::stringstream message;
                message << "\r" << ( int )progress << " % completed. Expected remaining time: ";
                if( progress > 0 )
                {
                    size_t expected_remain( elapsedTime * ( ( 100. - progress ) / progress ) );
                    message << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                    % 3600 ) % 60 ) << "\". ";
                }
                message << "Elapsed time: ";
                message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' " << ( ( elapsedTime % 3600 ) % 60 )
                                << "\". ";
                std::cout << message.str() << std::flush;
            }
        } // end m_verbose
    } // end parallel for (i)



    if( m_verbose )
    {
        size_t elapsedTime( difftime( time( NULL ), loopStartTime ) );
        std::stringstream message;
        message << "\r 100 % completed. Elapsed time: ";
        message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' " << ( ( elapsedTime % 3600 ) % 60 ) << "\".  Doing vector sums...";
        std::cout << message.str() << std::endl;
    }

    double sumT1( 0 ), sumT2( 0 ), sqT1( 0 ), sqT2( 0 ), sumProd( 0 );
    double wSumT1( 0 ), wSumT2( 0 ), wSqT1( 0 ), wSqT2( 0 ), wSumProd( 0 );
    double sumSize1( 0 ), sumSize2( 0 ), sumSqSize1( 0 ), sumSqSize2( 0 ), sumSizeProd( 0 );
    size_t usedPairs( 0 );


#pragma omp parallel sections
    {
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sumT1 += sumT1vec[i];
                wSumT1 += wSumT1vec[i];
                usedPairs += usedPairsVect[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sumT2 += sumT2vec[i];
                wSumT2 += wSumT2vec[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sqT1 += sqT1vec[i];
                wSqT1 += wSqT1vec[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sqT2 += sqT2vec[i];
                wSqT2 += wSqT2vec[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sumProd += sumProdvec[i];
                wSumProd += wSumProdvec[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sumSize1 += sumSize1Vec[i];
                sumSize2 += sumSize2Vec[i];
                sumSizeProd += sumSizeProdVec[i];
            }
        }
#pragma omp section
        {
            for( size_t i = 0; i < loopLengthMinusOne; ++i )
            {
                sumSqSize1 += sumSqSize1Vec[i];
                sumSqSize2 += sumSqSize2Vec[i];
            }
        }
    }

    if( m_verbose )
        std::cout << "\rSums obtained, doing final caculations..." << std::flush;

    // do the final computations
    double meanT1( sumT1 / usedPairs );
    double meanT2( sumT2 / usedPairs );
    double numerator( ( sumProd / usedPairs ) - ( meanT2 * meanT1 ) );
    double denominator1( ( sqT1 / usedPairs ) - ( meanT1 * meanT1 ) );
    double denominator2( ( sqT2 / usedPairs ) - ( meanT2 * meanT2 ) );
    float sCPCC( numerator / sqrt( denominator1 * denominator2 ) );

    double wMeanT1( wSumT1 / sumSize1 );
    double wMeanT2( wSumT2 / sumSize2 );
    double wNumerator( ( wSumProd / sumSizeProd ) - ( wMeanT1 * wMeanT2 ) );
    double wDenominator1( ( wSqT1 / sumSqSize1 ) - ( wMeanT1 * wMeanT1 ) );
    double wDenominator2( ( wSqT2 / sumSqSize2 ) - ( wMeanT2 * wMeanT2 ) );
    float tCPCC( wNumerator / sqrt( wDenominator1 * wDenominator2 ) );

    float effectiveGran( N2 / ( N2 - (2*usedPairs) ) );


//    std::cout<< std::endl<< " SumT1: " << sumT1 << "  usedPairs: "<< usedPairs << std::endl;
//    std::cout<<  " SumT2: " << sumT2 << "  usedPairs: "<< usedPairs << std::endl;
//    std::cout<<  " wSumT1: " << wSumT1 << "  sumSize1: "<< sumSize1 << std::endl;
//    std::cout<<  " wSumT2: " << wSumT2 << "  sumSize2: "<< sumSize2 << std::endl;
//    std::cout<<  " sumProd: " << sumProd << "  wSumProd: "<< wSumProd << std::endl;
//    std::cout<<  " A: " << ( sumProd / usedPairs ) << "  wA: "<< ( wSumProd / sumSizeProd ) << std::endl;
//    std::cout<<  " meanT1: " << meanT1 << "  wMeanT1: "<< wMeanT1 << std::endl;
//    std::cout<< " meanT2: " << meanT2 << "  wMeanT2: "<< wMeanT2 << std::endl;
//    std::cout<<  " num: " << numerator << "  wNum: "<< wNumerator << std::endl;
//    std::cout<<  " den1: " << denominator1 << "  wDen1: "<< wDenominator1 << std::endl;
//    std::cout<<  " den2: " << denominator2 << "  wDen2: "<< wDenominator2 << std::endl;




    if( denominator1 == 0 || denominator2 == 0 )
    {
        std::cerr << "WARNING @ cpct(): one or two of the trees is completely flat, no structure... CPCT will be set to 0"
                        << std::endl;
        sCPCC = 0;
        tCPCC = 0;
    }

    if( m_verbose )
    {
        std::cout << std::endl;
        std::cout << "Weighted tCPCC: " << tCPCC << std::endl;
        std::cout << "Simple CPCC: " << sCPCC << std::endl;
        std::cout << "Used pairs (%): " << (usedPairs*100.0)/totalPairs << std::endl;
        std::cout << "Effective granularity: " << effectiveGran << std::endl;

    }

    return std::make_pair( std::make_pair( tCPCC, sCPCC ), std::make_pair( (usedPairs*1.0)/totalPairs, effectiveGran ) );
} // end "cpct()" -----------------------------------------------------------------


bool treeComparer::leafCorrespondence()
{
    m_baseNodes1.clear();
    m_baseNodes2.clear();
    m_baseCoords1.clear();
    m_baseCoords2.clear();

    if( m_tree1.getDataGrid() != m_tree2.getDataGrid() )
    {
        std::cout << "Trees are in different coordinate grids: ";
        if( m_tree1.convert2grid( HC_NIFTI ) )
        {
            std::cout << "Tree 1 was converted to nifti coordinates";
        }
        if( m_tree2.convert2grid( HC_NIFTI ) )
        {
            std::cout << "Tree 2 was converted to nifti coordinates";
        }
        std::cout << std::endl;
    }

    if( m_tree1.getDataGrid() != m_tree2.getDataGrid() )
    {
        throw std::runtime_error( "ERROR @ treeCompare::leafCorrespondence(): did not manage to convert trees to the smae grid" );
    }

    if( m_tree1.m_coordinates == m_tree2.m_coordinates )
    {
        return false;
    }

    for( size_t i = 0; i < m_tree1.m_coordinates.size(); ++i )
    {
        if( std::find( m_tree2.m_coordinates.begin(), m_tree2.m_coordinates.end(), m_tree1.m_coordinates[i] )
                        == m_tree2.m_coordinates.end() )
        {
            m_tree1.fetchLeaf( i )->setFlag( true );
        }
    }
    std::pair< size_t, size_t > pruned1( m_tree1.cleanup() );
    if( m_verbose )
    {
        std::cout << "Eliminated " << pruned1.first << " leaves and " << pruned1.second << " nodes of Tree 1" << std::endl;
    }

    for( size_t i = 0; i < m_tree2.m_coordinates.size(); ++i )
    {
        if( std::find( m_tree1.m_coordinates.begin(), m_tree1.m_coordinates.end(), m_tree2.m_coordinates[i] )
                        == m_tree1.m_coordinates.end() )
        {
            m_tree2.fetchLeaf( i )->setFlag( true );
        }
    }
    std::pair< size_t, size_t > pruned2( m_tree2.cleanup() );
    if( m_verbose )
    {
        std::cout << "Eliminated " << pruned2.first << " leaves and " << pruned2.second << " nodes of Tree 2" << std::endl;
    }

    if( m_tree1.m_coordinates != m_tree2.m_coordinates )
    {
        throw std::runtime_error( "ERROR @ treeCompare::leafCorrespondence(): failed to equalize the leaves" );
    }

    return true;
} // end "leafCorrespondence()" -----------------------------------------------------------------


void treeComparer::greedyCorrespondence( float dissimThreshold, bool redoCoords )
{

    if( dissimThreshold > 1 )
    {
        dissimThreshold = 1;
    }
    if( dissimThreshold < 0.1 )
    {
        dissimThreshold = 0.1;
    }
    std::vector< size_t > protoCorrespTable, correspTable;
    std::vector< std::pair< float, float > > correspDistances;
    std::vector< bool > isMatched1, isMatched2;

    fetchBaseNodes( false );
    if (m_baseNodes1.size() != m_baseCoords1.size() || m_baseNodes2.size() != m_baseCoords2.size())
    {
        if ( redoCoords )
        {
            if ( m_verbose )
            {
                std::cout<< "Getting cluster coordinate information..." <<std::endl;
            }
            fetchBaseNodes( true );
        }
        else
        {
            m_baseCoords1.clear();
            m_baseCoords1.resize( m_baseNodes1.size(), WHcoord() );
            m_baseCoords2.clear();
            m_baseCoords2.resize( m_baseNodes2.size(), WHcoord() );
        }
    }

    if ( m_baseDistMatrix.empty())
        throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): base node distance matrix is empty" );

    if( ( m_baseDistMatrix.size() !=  m_baseNodes1.size() ) || ( m_baseDistMatrix.front().size() != m_baseNodes2.size() ) )
        throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): base node distance matrix dimensions dont match base node vecotrs" );

    std::vector< std::vector< dist_t > > baseDistMatrix(m_baseDistMatrix);

    if( m_verbose )
    {
        std::cout << "Computing base-node distance table by simple greedy correspondence:" << std::endl;
    }

    const size_t nomatch( m_initialSizes.second );
    protoCorrespTable.resize( m_baseNodes1.size(), nomatch );
    isMatched1.resize( m_baseNodes1.size(), false );
    isMatched2.resize( m_baseNodes2.size(), false );

    std::vector< size_t > oldBaseNodes1( m_baseNodes1 ), oldBaseNodes2( m_baseNodes2 );


    std::list< size_t > leftNodes1( m_baseNodes1.begin(), m_baseNodes1.end() );
    std::list< size_t > leftNodes2( m_baseNodes2.begin(), m_baseNodes2.end() );



    // greedy matching
    while( !leftNodes1.empty() || !leftNodes2.empty() )
    {

        std::vector<size_t> matchVector( baseDistMatrix.size(), 0 );
        std::vector<dist_t> minDistVector( baseDistMatrix.size(), 2 );


        // search for lowest distance
        #pragma omp parallel for schedule(guided)
        for( size_t i = 0; i < baseDistMatrix.size(); ++i )
        {
            for( size_t j = 0; j < baseDistMatrix[i].size(); ++j )
            {
                if( baseDistMatrix[i][j] < minDistVector[i] )
                {
                    minDistVector[i] = baseDistMatrix[i][j];
                    matchVector[i] = j;
                }
            }
        }

        size_t match1( 0 ), match2( 0 );
        dist_t minDist( 2 );

        for( size_t i = 0; i < matchVector.size(); ++i )
        {
            if( minDistVector[i] < minDist )
            {

                minDist = minDistVector[i];
                match1 = i;
                match2 = matchVector[i];
            }
        }

        // if no more matches are found, stop
        if( minDist > dissimThreshold )
            break;

        // update correspondence table, reduce lists and clear row and column of distance table
        protoCorrespTable[match1] = match2;
        leftNodes1.remove( m_baseNodes1[match1] );
        leftNodes2.remove( m_baseNodes2[match2] );
        baseDistMatrix[match1].assign( baseDistMatrix[match1].size(), 2 );
        for( size_t i = 0; i < baseDistMatrix.size(); ++i )
        {
            baseDistMatrix[i][match2] = 2;
        }
        isMatched1[match1]=true;
        isMatched2[match2]=true;
    } // end inner while


    m_fullCorrespondence = protoCorrespTable;


    std::vector<size_t> nodeLookup1, nodeLookup2;
    bool deletion( false );

    // delete excess nodes or those without proper match
    if( !leftNodes1.empty() )
    {

        std::cout << "Removing " << leftNodes1.size() << " base nodes from tree1...";
        size_t sizeSum( 0 );
        for( std::list< size_t >::iterator iter( leftNodes1.begin() ); iter != ( leftNodes1.end() ); ++iter )
        {
            std::vector< size_t > leaves2prune( m_tree1.getLeaves4node( *iter ) );
            sizeSum += leaves2prune.size();

            for( size_t pruneIndex = 0; pruneIndex < leaves2prune.size(); ++pruneIndex )
            {
                m_tree1.fetchLeaf( leaves2prune[pruneIndex] )->setFlag( true );
            }
        }
        std::cout << "mean size: " << sizeSum/leftNodes1.size() << " leaves." << std::endl;
        m_tree1.cleanup( &nodeLookup1 );
        deletion = true;
    }
    else
    {
        nodeLookup1.clear();
        nodeLookup1.reserve( m_tree1.getNumNodes() );
        for( size_t i = 0; i < m_tree1.getNumNodes(); ++i )
        {
            nodeLookup1.push_back( i );
        }
    }


    if( !leftNodes2.empty() )
    {
        std::cout << "Removing " << leftNodes2.size() << " base nodes from tree2...";
        size_t sizeSum( 0 );
        for( std::list< size_t >::iterator iter( leftNodes2.begin() ); iter != ( leftNodes2.end() ); ++iter )
        {
            std::vector< size_t > leaves2prune( m_tree2.getLeaves4node( *iter ) );
            sizeSum += leaves2prune.size();

            for( size_t pruneIndex = 0; pruneIndex < leaves2prune.size(); ++pruneIndex )
            {
                m_tree2.fetchLeaf( leaves2prune[pruneIndex] )->setFlag( true );
            }
        }
        std::cout << "mean size: " << sizeSum/leftNodes2.size() << " leaves." << std::endl;
        m_tree2.cleanup( &nodeLookup2 );
        deletion = true;
    }
    else
    {
        nodeLookup2.clear();
        nodeLookup2.reserve( m_tree2.getNumNodes() );
        for( size_t i = 0; i < m_tree2.getNumNodes(); ++i )
        {
            nodeLookup2.push_back( i );
        }
    }

    if( deletion ) // one or both of the trees was pruned update ids in matching vector
    {
        fetchBaseNodes(false);

        std::cout << "Updating correspondence table..." << std::endl;
        if( nodeLookup1.size() < protoCorrespTable.size() || nodeLookup2.size() < protoCorrespTable.size() )
        {
            std::cerr << "Correspondence vector size: " << protoCorrespTable.size() << ". Lookup1: " << nodeLookup1.size()
                      << ". Lookup2: " << nodeLookup2.size() << std::endl;
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): lookups are smaller than correspondence" );
        }
        if( m_baseNodes1.size() != m_baseNodes2.size() )
        {
            std::cerr << "basenodes1 size: " << m_baseNodes1.size() << ". basenodes2 size: " << m_baseNodes2.size() << std::endl;
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): new basenodes dimensions dont match" );
        }

        correspTable.resize( m_baseNodes1.size(), 0 );
        correspDistances.resize( m_baseNodes1.size(), std::make_pair( 2, 99 ) );

        std::cout << "getting new IDs..." << std::endl;

        for( size_t i = 0; i < protoCorrespTable.size(); ++i )
        {
            size_t oldRelativeID1( i );
            size_t oldAbsoluteID1( oldBaseNodes1[oldRelativeID1] );
            size_t newAbsoluteID1( nodeLookup1[oldAbsoluteID1] );

            size_t oldRelativeID2( protoCorrespTable[oldRelativeID1] );
            if( oldRelativeID2 == nomatch)
            {
                continue; // the node the table referred to has been deleted
            }
            size_t oldAbsoluteID2( oldBaseNodes2[oldRelativeID2] );
            size_t newAbsoluteID2( nodeLookup2[oldAbsoluteID2] );

            if( newAbsoluteID1 >= m_tree1.getNumNodes() || newAbsoluteID2 >= m_tree2.getNumNodes() )
            {
                continue; // the node the table referred to has been deleted
            }

            size_t newRelativeID1( std::find( m_baseNodes1.begin(), m_baseNodes1.end(), newAbsoluteID1 ) - m_baseNodes1.begin() );
            size_t newRelativeID2( std::find( m_baseNodes2.begin(), m_baseNodes2.end(), newAbsoluteID2 ) - m_baseNodes2.begin() );

            if( newRelativeID1 >= m_baseNodes1.size() || newRelativeID2 >= m_baseNodes2.size() )
            {
                std::cerr << "new abs ID1: " << newAbsoluteID1 << ". new abs ID2: " << newAbsoluteID2 << std::endl;
                throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): new IDs dont match basenodes" );
            }

            correspTable[newRelativeID1] = newRelativeID2;

            float tractDist( m_baseDistMatrix[oldRelativeID1][oldRelativeID2] );
            float clusterEucDist(m_baseCoords1[oldRelativeID1].getPhysDist(m_baseCoords2[oldRelativeID2]) );

            correspDistances[newRelativeID1] = std::make_pair( tractDist , clusterEucDist);
        }

        if( std::count( correspTable.begin(), correspTable.end(), nomatch) != 0 )
        {
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): error in correspondence table" );
        }

        // eliminate deleted nodes from distance matrix
        {
            std::vector< std::vector< dist_t > > emptyMatrix;
            baseDistMatrix.swap(emptyMatrix);
        }
        baseDistMatrix.resize(correspTable.size(),std::vector<float>( correspTable.size(), 0 ) );

        size_t countRow( 0 );

        std::cout << "Cropping distance table..." << std::endl;


        for( size_t i = 0; i < isMatched1.size(); ++i )
        {
            if( isMatched1[i])
            {
                size_t countCol( 0 );
                for( size_t j = 0; j < isMatched2.size(); ++j )
                {
                    if( isMatched2[j])
                    {
                        baseDistMatrix[countRow][countCol] =  m_baseDistMatrix[i][j];
                        ++countCol;
                    }
                }
                ++countRow;
            }
        }
        m_baseDistMatrix.swap(baseDistMatrix);
    }
    else
    {
        correspTable = protoCorrespTable;

        for( size_t i = 0; i < protoCorrespTable.size(); ++i )
        {
            float tractDist( m_baseDistMatrix[i][protoCorrespTable[i]] );
            float clusterEucDist(m_baseCoords1[i].getPhysDist(m_baseCoords2[protoCorrespTable[i]]) );
            correspDistances.resize( m_baseNodes1.size(), std::make_pair( 2, 99 ) );
            correspDistances[i] = std::make_pair( tractDist , clusterEucDist);
        }
    }

    if( m_verbose )
    {
        std::cout << reportBaseNodes() << std::endl;
    }

    m_newCorrespondence = correspTable;

    m_newCorrespReverse.clear();
    m_newCorrespReverse.resize(m_newCorrespondence.size(),0);
    for( size_t i = 0; i < m_newCorrespondence.size(); ++i )
    {
        m_newCorrespReverse[m_newCorrespondence[i]]=i;
    }
    m_correspDistances =  correspDistances;

    return;
} // end treeComparer::simpleCorrespondence() -------------------------------------------------------------------------------------

std::vector<float> treeComparer::rateCorrespondence()
{
    std::vector<float> ratingResults;

    if( m_verbose )
    {
        std::cout << "Rating matching quality..." << std::endl;
    }

    if( m_baseNodes1.size() != m_baseNodes2.size() )
    {
        throw std::runtime_error( "ERROR @ treeCompare::rateCorrespondence(): base node vectors have different sizes" );
    }


    if( m_baseNodes1.empty() )
    {

        std::cerr << "WARNING @ treeCompare::rateCorrespondence(): leaf-wise matching performed, cannot be rated" << std::endl;
        ratingResults.assign(4, 0);
        return ratingResults;
    }


    if( m_baseNodes1.size() != m_newCorrespondence.size() )
    {
        throw std::runtime_error(
                    "ERROR @ treeCompare::rateCorrespondence(): correspondance vector size does not match basenodes vector" );
    }
    if( m_baseNodes1.size() != m_correspDistances.size() )
    {
        throw std::runtime_error(
                    "ERROR @ treeCompare::rateCorrespondence(): correspondance distance vector size does not match basenodes vector" );
    }

    size_t loopLength(m_newCorrespondence.size());


    size_t sizeSum1(0), sizeSum2(0), sizeSqSum1(0), sizeSqSum2(0), sizeProdSum(0);
    double distSum(0), distWeightSum(0), physDistSum(0), physDistWeightSum(0);
    float minDist( 1 ), maxDist( 0 );

    // loop
    for( size_t i = 0; i < loopLength; ++i )
    {
        size_t size1( m_tree1.getNode( m_baseNodes1[i] ).getSize() );
        size_t size2( m_tree2.getNode( m_baseNodes2[m_newCorrespondence[i]] ).getSize() );
        float distance(m_correspDistances[i].first);
        float physDist(m_correspDistances[i].second);


        if (distance != m_baseDistMatrix[i][m_newCorrespondence[i]])
        {
            std::cerr<< "ERROR @ treeComparer::rateCorrespondence(): distance in distance table does not correspond with distance in matrix" <<std::endl;
            std::cerr<< "Table: "<< distance << ". Matrix: " << m_baseDistMatrix[i][m_newCorrespondence[i]] <<std::endl;
        }

        if (minDist > distance )
        {
            minDist = distance;
        }
        if (maxDist < distance )
        {
            maxDist = distance;
        }

        sizeSum1 += size1;
        sizeSum2 += size2;
        sizeSqSum1 += ( size1 * size1 );
        sizeSqSum2 += ( size2 * size2 );
        sizeProdSum += ( size1 * size2 );

        distSum += distance;
        distWeightSum += ( distance * ( size1 + size2 ) );
        physDistSum += physDist;
        physDistWeightSum += ( physDist * ( size1 + size2 ) );

    } // end for
    double N( loopLength );
    double meanSize1( sizeSum1 / N );
    double meanSize2( sizeSum2 / N );

    double sizeCorrelNum( ( sizeProdSum /  N ) - ( meanSize1 * meanSize2 ) );
    double sizeCorrelDen1( ( sizeSqSum1 / N ) - ( meanSize1 * meanSize1 ) );
    double sizeCorrelDen2( ( sizeSqSum2 / N ) - ( meanSize2 * meanSize2 ) );

    double sizeCorrel( sizeCorrelNum / std::sqrt( sizeCorrelDen1 * sizeCorrelDen2 ) );
    ratingResults.push_back(sizeCorrel);

    double meanMatchDist( distSum / N );
    ratingResults.push_back(meanMatchDist);

    double weightedMatchDist( distWeightSum / ( sizeSum1 + sizeSum2 ) );
    ratingResults.push_back(weightedMatchDist);

    double amountMatched( ( 1.0 *( sizeSum1 + sizeSum2 ) ) / ( m_initialSizes.first + m_initialSizes.second ) );
    ratingResults.push_back(amountMatched);

    double meanPhysDist( physDistSum / N );
    ratingResults.push_back(meanPhysDist);

    double weightedPhysDist( physDistWeightSum / ( sizeSum1 + sizeSum2 ) );
    ratingResults.push_back(weightedPhysDist);

    if( m_verbose )
    {
        std::cout << "% of basenodes matched:\t " << 100.0*amountMatched  << std::endl;
        std::cout << "Size-Weighted Match Distance:\t " << weightedMatchDist << std::endl;
        std::cout << "Max dist: " << maxDist << ". min Dist: " << minDist << std::endl;
        std::cout << "Size-Weighted Euclidean Distance:\t " << weightedPhysDist << std::endl;

    }

    return ratingResults;

} // end treeComparer::rateCorrespondence() -------------------------------------------------------------------------------------


std::pair< float, float > treeComparer::applyNoiseBaseline( const float noiseAlpha )
{


    if( m_verbose )
    {
        std::cout << "Applying matching noise corrections to the trees. Alpha: "<< noiseAlpha << std::endl;
    }
    float numNodesTree1( m_tree1.getNumNodes() );
    float maxgran1( noiseBaseline( TREE1, noiseAlpha ) );
    float reductiontree1( (numNodesTree1-m_tree1.getNumNodes())/numNodesTree1 );

    float numNodesTree2( m_tree2.getNumNodes() );
    float maxgran2( noiseBaseline( TREE2, noiseAlpha ) );
    float reductiontree2( (numNodesTree2-m_tree2.getNumNodes())/numNodesTree2 );

    float averageNoisyLoss( ( reductiontree1 + reductiontree2 )  / 2.0 );
    if( m_verbose )
    {
        std::cout << "Tree1 lost " << 100.0 * reductiontree1 <<" % of its nodes" << std::endl;
        std::cout << "Tree2 lost " << 100.0 * reductiontree2 <<" % of its nodes" << std::endl;
        std::cout << "Overall % structure loss: "<< 100 * averageNoisyLoss << std::endl;
        std::cout << "Size of maxgran part 1: " << maxgran1 << std::endl;
        std::cout << "Size of maxgran part 2: " << maxgran2 << std::endl;
        std::cout << "Average maxgran size: "<< ( ( maxgran1 + maxgran2 )  / 2.0 ) << std::endl;
    }
    return std::make_pair(maxgran1,maxgran2);
} // end treeComparer::applyNoiseBaseline() -------------------------------------------------------------------------------------


bool treeComparer::fetchBaseNodes( bool doGetCoords )
{
    m_baseNodes1.clear();
    m_baseNodes2.clear();
    m_baseNodes1 = m_tree1.getRootBaseNodes();
    m_baseNodes2 = m_tree2.getRootBaseNodes();

    m_realBaseNodes = true;

    if( !m_tree1.testRootBaseNodes() )
    {
        m_realBaseNodes = false;
    }
    if( !m_tree2.testRootBaseNodes() )
    {
        m_realBaseNodes = false;
    }

    if ( m_originalBaseNodes1.empty() )
    {
        m_originalBaseNodes1 = m_baseNodes1;
    }
    if ( m_originalBaseNodes2.empty() )
    {
        m_originalBaseNodes2 = m_baseNodes2;
    }

    if ( doGetCoords )
    {
        getBaseCoords();
    }
    return m_realBaseNodes;
} // end treeComparer::fetchBaseNodes() -------------------------------------------------------------------------------------


void treeComparer::getBaseCoords()
{

    std::cout << "Obtaining base node coordinates..." << std::endl;

    m_baseCoords1.clear();
    m_baseCoords1.resize( m_baseNodes1.size() );


//    #pragma omp parallel for
    for( size_t i = 0; i < m_baseNodes1.size(); ++i )
    {
        if( m_coordsFromFile )
        {
            fileManagerFactory clusterFMF( m_meanTractFolder1 );
            fileManager& clusterMaskmanager( clusterFMF.getFM() );
            std::string clusterMaskFilename( clusterMaskmanager.getClusterMaskFilename( m_baseNodes1[i] ) );
            clusterMaskmanager.loadMaskImage( clusterMaskFilename );
            m_baseCoords1[i]=( clusterMaskmanager.meanCoordFromMask() );
        }
        else
        {
            m_baseCoords1[i]=( m_tree1.getMeanCoordinate4node( m_baseNodes1[i] ) );
        }
    }

    m_baseCoords2.clear();
    m_baseCoords2.resize( m_baseNodes2.size() );



//    #pragma omp parallel for
    for( size_t i = 0; i < m_baseNodes2.size(); ++i )
    {
        if( m_coordsFromFile )
        {
            fileManagerFactory clusterFMF( m_meanTractFolder2 );
            fileManager& clusterMaskmanager( clusterFMF.getFM() );
            std::string clusterMaskFilename( clusterMaskmanager.getClusterMaskFilename( m_baseNodes2[i] ) );
            clusterMaskmanager.loadMaskImage( clusterMaskFilename );
            m_baseCoords2[i]=( clusterMaskmanager.meanCoordFromMask() );
        }
        else
        {
            m_baseCoords2[i]=( m_tree2.getMeanCoordinate4node( m_baseNodes2[i] ) );
        }
    }


    return;
} // end treeComparer::getBaseNodeCoords() -------------------------------------------------------------------------------------


std::string treeComparer::reportBaseNodes() const
{
    std::stringstream message;

    {
        size_t bMax( 0 ), bMin( m_tree1.getNumLeaves() ), numBig( 0 ), numSmall( 0 );

        for( std::vector< size_t >::const_iterator iter( m_baseNodes1.begin() ); iter != m_baseNodes1.end(); ++iter )
        {
            size_t currentSize( m_tree1.getNode( *iter ).getSize() );
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

        message << "Tree1: " << m_baseNodes1.size() << " base nodes. Biggest: " << bMax << ". Smallest: " << bMin << ". "
                        << numBig << " >= 100." << numSmall << " <= 10." << std::endl;
    }

    {
        size_t bMax( 0 ), bMin( m_tree2.getNumLeaves() ), numBig( 0 ), numSmall( 0 );

        for( std::vector< size_t >::const_iterator iter( m_baseNodes2.begin() ); iter != m_baseNodes2.end(); ++iter )
        {
            size_t currentSize( m_tree2.getNode( *iter ).getSize() );
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

        message << "Tree2: " << m_baseNodes2.size() << " base nodes. Biggest: " << bMax << ". Smallest: " << bMin << ". "
                        << numBig << " >= 100." << numSmall << " <= 10." << std::flush;
    }
    std::string outMessage( message.str() );

    return outMessage;
} // end treeComparer::reportBaseNodes() -------------------------------------------------------------------------------------


void treeComparer::getBaseDistMatrix()
{
    if( m_meanTractFolder1.empty() || m_meanTractFolder2.empty() )
        throw std::runtime_error( "ERROR @ treeCompare::getBaseDistMatrix(): Location of mean tract folders is invalid" );

    if( m_baseNodes1.empty() || m_baseNodes2.empty() )
        throw std::runtime_error( "ERROR @ treeComparer::getBaseDistMatrix(): one (or both) of the base node vectors is empty" );


    m_baseDistMatrix.clear();

    std::vector< std::vector< dist_t > > baseDistMatrix;
    {
        std::vector< dist_t > fillvect( m_baseNodes2.size(), 1 );
        baseDistMatrix.resize( m_baseNodes1.size(), fillvect );
    }

    if( m_verbose )
    {
        std::cout << "Obtaining base node information..." << std::endl;
    }
    fetchBaseNodes( true );
    if( m_verbose )
    {
        std::cout << "base node information obtained" << std::endl;
    }


    if( m_meanTractsFromFile )
    {
        if( m_verbose )
        {
            std::cout << "Mean tracts will be read from files" << std::endl;
        }
    }
    else
    {
        if( m_singleTractFolder1.empty() || m_singleTractFolder2.empty() )
        {
            throw std::runtime_error( "ERROR @ treeCompare::getBaseDistMatrix(): Location of single tracts folders is invalid" );
        }

        if( m_verbose )
        {
            std::cout << "Calculating and writing base node mean tracts" << std::endl;
        }

        treeManager manager1( &m_tree1, m_verbose );
        manager1.setSingleTractFolder( m_singleTractFolder1 );
        manager1.setMeanTractFolder( m_meanTractFolder1 );
        manager1.writeMeanTracts( m_baseNodes1 );

        treeManager manager2( &m_tree2, m_verbose );
        manager2.setSingleTractFolder( m_singleTractFolder2 );
        manager2.setMeanTractFolder( m_meanTractFolder2 );
        manager2.writeMeanTracts( m_baseNodes2 );
    }

    size_t progCount( 0 );
    time_t lastTime( time( NULL ) ), startTime( time( NULL ) );

    if( m_verbose )
        std::cout << "Calculating distance matrix" << std::endl;

    // loop through base nodes of tree 1
    #pragma omp parallel for schedule(guided)
    for( size_t i = 0; i < m_baseNodes1.size(); ++i )
    {
        // file managers
        fileManagerFactory nodeFileMF1( m_meanTractFolder1 );
        fileManager& nodeFM1( nodeFileMF1.getFM() );
        nodeFM1.readAsLog();
        nodeFM1.readAsUnThres();
        fileManagerFactory nodeFileMF2( m_meanTractFolder2 );
        fileManager& nodeFM2( nodeFileMF2.getFM() );
        nodeFM2.readAsLog();
        nodeFM2.readAsUnThres();

        // get mean position for base node
        WHcoord baseCoord1( m_baseCoords1[i] );

        // base node tract
        compactTract baseTract1;
        nodeFM1.readNodeTract( m_baseNodes1[i], &baseTract1 );
        baseTract1.threshold( m_tractThreshold1 );
        baseTract1.computeNorm();

        for( size_t j = 0; j < m_baseNodes2.size(); ++j )
        {
            float pDist( baseCoord1.getPhysDist( m_baseCoords2[j] ) );
            if( pDist > m_maxPhysDist && m_maxPhysDist > 0 )
            {
                #pragma omp atomic
                ++progCount;
                continue; // if areas are too far away from each other, leave as maximum distance (1)
            }

            compactTract baseTract2;
            nodeFM2.readNodeTract( m_baseNodes2[j], &baseTract2 );
            baseTract2.threshold( m_tractThreshold2 );
            baseTract2.computeNorm();
            baseDistMatrix[i][j] = baseTract1.tractDistance( baseTract2 );
        }

        #pragma omp atomic
        ++progCount;
        #pragma omp single nowait // only one thread executes output
        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                lastTime = currentTime;
                size_t currentCount( progCount );
                float progress = ( currentCount ) * 100. / ( m_baseNodes1.size() * m_baseNodes2.size() );
                size_t elapsedTime( difftime( currentTime, startTime ) );
                std::stringstream message;
                message << "\r" << ( int )progress << " % completed. Expected remaining time: ";
                if( progress > 0 )
                {
                    int expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                    message << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                    % 3600 ) % 60 ) << "\". ";
                }
                message << "Elapsed time: ";
                message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' " << ( ( elapsedTime % 3600 ) % 60 )
                                << "\". ";
                std::cout << message.str() << std::flush;
            }
        } // end m_verbose
    } // end parallel loop
    std::cout << "\r100 % Completed (" << m_baseNodes1.size() << "x" << m_baseNodes2.size() << " distance matrix)" << std::endl;

    m_baseDistMatrix = baseDistMatrix;
    return;
} // end treeComparer::getBaseDistMatrix() -------------------------------------------------------------------------------------

void treeComparer::writeBaseDistMatrix( std::string matrixFilename )
{
    fileManagerFactory matrixFileMF;
    fileManager& matrixFM( matrixFileMF.getFM() );
    matrixFM.writeInFloat();
    matrixFM.storeZipped();
    matrixFM.writeMatrix(matrixFilename, VTFloat32 ,m_baseDistMatrix );
    return;
} // end treeComparer::writeBaseDistMatrix() -------------------------------------------------------------------------------------

void treeComparer::readBaseDistMatrix( std::string matrixFilename )
{
    fileManagerFactory matrixFileMF;
    fileManager& matrixFM( matrixFileMF.getFM() );
    matrixFM.readMatrix( matrixFilename, &m_baseDistMatrix );
    return;
} // end treeComparer::readBaseDistMatrix() -------------------------------------------------------------------------------------

void treeComparer::randomCorrespondence()
{

    if (m_baseNodes1.empty() || m_baseCoords1.empty() || m_baseNodes2.empty() || m_baseCoords2.empty() )
    {
        std::cout<<"Fetching nodes and coordinates"<<std::endl;
        fetchBaseNodes( true );
    }

    std::vector< size_t > oldBaseNodes1( m_baseNodes1 ), oldBaseNodes2( m_baseNodes2 );

    const size_t nomatch( m_initialSizes.second );
    std::vector< size_t > protoCorrespondence, correspondence;
    protoCorrespondence.resize( m_baseNodes1.size(), nomatch );

    std::vector< size_t > values( m_baseNodes2.size(), 0 );
    for( size_t i = 0; i < values.size(); ++i )
    {
        values[i] = i;
    }

    // This is the underlying integer random number generator
    boost::mt19937 igen;
    igen.seed( time( NULL ) );

    // uniformly distributed random number generator between 0 and 1
    boost::variate_generator< boost::mt19937, boost::uniform_01< > > u01Prng( igen, boost::uniform_01< >() );

    std::list< size_t > leftNodes1;
    std::list< size_t > leftNodes2;

    // do it for every base node
    for (size_t i = 0; i < m_baseNodes1.size(); ++i )
    {
        std::vector< size_t > candidateIndexes;
        candidateIndexes.reserve(m_baseNodes2.size());

        // obtain a list of indexes to basenodes in tree 2 in the vicinity of this basenode of tree 1
        for(size_t j = 0; j < values.size(); ++j )
        {
            float pDist( m_baseCoords1[i].getPhysDist( m_baseCoords2[values[j]] ) );
            if( pDist <= m_maxPhysDist || m_maxPhysDist == 0 )
            {
                candidateIndexes.push_back(j);
            }
        }

        // if there are no nodes left on the vicinity, this node will not be matched
        if( candidateIndexes.empty())
        {
            leftNodes1.push_back(m_baseNodes1[i]);
            continue;
        }

        // randomly choose a candidate index   ( if there is more than one )
        if( candidateIndexes.size() == 1 )
        {
            protoCorrespondence[i] = ( values[candidateIndexes.front()] );
            values.erase( values.begin() + candidateIndexes.front() );
        }
        else
        {
            size_t chosenIndex( std::floor( u01Prng() * candidateIndexes.size() ) );
            if( chosenIndex == candidateIndexes.size() )
                --chosenIndex;

            protoCorrespondence[i] = ( values[candidateIndexes[chosenIndex]] );
            values.erase( values.begin() + candidateIndexes[chosenIndex] );
        }
    }

    m_fullCorrespondence = protoCorrespondence;


    //nodes of tree 2 without a match
    for (size_t i = 0; i < values.size(); ++i )
    {
        leftNodes2.push_back(m_baseNodes2[values[i]]);
    }


    ///////////////////////////////////////////////////


    std::vector<size_t> nodeLookup1, nodeLookup2;
    bool deletion( false );

    // delete excess nodes or those without proper match
    if( !leftNodes1.empty() )
    {

        std::cout << "Removing " << leftNodes1.size() << " base nodes from tree1...";
        size_t sizeSum( 0 );
        for( std::list< size_t >::iterator iter( leftNodes1.begin() ); iter != ( leftNodes1.end() ); ++iter )
        {
            std::vector< size_t > leaves2prune( m_tree1.getLeaves4node( *iter ) );
            sizeSum += leaves2prune.size();

            for( size_t pruneIndex = 0; pruneIndex < leaves2prune.size(); ++pruneIndex )
            {
                m_tree1.fetchLeaf( leaves2prune[pruneIndex] )->setFlag( true );
            }
        }
        std::cout << "mean size: " << sizeSum/leftNodes1.size() << " leaves." << std::endl;
        m_tree1.cleanup( &nodeLookup1 );
        deletion = true;
    }
    else
    {
        nodeLookup1.clear();
        nodeLookup1.reserve( m_tree1.getNumNodes() );
        for( size_t i = 0; i < m_tree1.getNumNodes(); ++i )
        {
            nodeLookup1.push_back( i );
        }
    }


    if( !leftNodes2.empty() )
    {
        std::cout << "Removing " << leftNodes2.size() << " base nodes from tree2...";
        size_t sizeSum( 0 );
        for( std::list< size_t >::iterator iter( leftNodes2.begin() ); iter != ( leftNodes2.end() ); ++iter )
        {
            std::vector< size_t > leaves2prune( m_tree2.getLeaves4node( *iter ) );
            sizeSum += leaves2prune.size();

            for( size_t pruneIndex = 0; pruneIndex < leaves2prune.size(); ++pruneIndex )
            {
                m_tree2.fetchLeaf( leaves2prune[pruneIndex] )->setFlag( true );
            }
        }
        std::cout << "mean size: " << sizeSum/leftNodes2.size() << " leaves." << std::endl;
        m_tree2.cleanup( &nodeLookup2 );
        deletion = true;
    }
    else
    {
        nodeLookup2.clear();
        nodeLookup2.reserve( m_tree2.getNumNodes() );
        for( size_t i = 0; i < m_tree2.getNumNodes(); ++i )
        {
            nodeLookup2.push_back( i );
        }
    }

    if( deletion ) // one or both of the trees was pruned update ids in matching vector
    {
        fetchBaseNodes(false);

        std::cout << "Updating correspondence table..." << std::endl;
        if( nodeLookup1.size() < protoCorrespondence.size() || nodeLookup2.size() < protoCorrespondence.size() )
        {
            std::cerr << "Correspondence vector size: " << protoCorrespondence.size() << ". Lookup1: " << nodeLookup1.size()
                      << ". Lookup2: " << nodeLookup2.size() << std::endl;
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): lookups are smaller than correspondence" );
        }
        if( m_baseNodes1.size() != m_baseNodes2.size() )
        {
            std::cerr << "basenodes1 size: " << m_baseNodes1.size() << ". basenodes2 size: " << m_baseNodes2.size() << std::endl;
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): new basenodes dimensions dont match" );
        }

        correspondence.resize( m_baseNodes1.size(), 0 );

        std::cout << "getting new IDs..." << std::endl;

        for( size_t i = 0; i < protoCorrespondence.size(); ++i )
        {
            size_t oldRelativeID1( i );
            size_t oldAbsoluteID1( oldBaseNodes1[oldRelativeID1] );
            size_t newAbsoluteID1( nodeLookup1[oldAbsoluteID1] );

            size_t oldRelativeID2( protoCorrespondence[oldRelativeID1] );
            if( oldRelativeID2 == nomatch)
            {
                continue; // the node the table referred to has been deleted
            }
            size_t oldAbsoluteID2( oldBaseNodes2[oldRelativeID2] );
            size_t newAbsoluteID2( nodeLookup2[oldAbsoluteID2] );

            if( newAbsoluteID1 >= m_tree1.getNumNodes() || newAbsoluteID2 >= m_tree2.getNumNodes() )
            {
                continue; // the node the table referred to has been deleted
            }

            size_t newRelativeID1( std::find( m_baseNodes1.begin(), m_baseNodes1.end(), newAbsoluteID1 ) - m_baseNodes1.begin() );
            size_t newRelativeID2( std::find( m_baseNodes2.begin(), m_baseNodes2.end(), newAbsoluteID2 ) - m_baseNodes2.begin() );

            if( newRelativeID1 >= m_baseNodes1.size() || newRelativeID2 >= m_baseNodes2.size() )
            {
                std::cerr << "new abs ID1: " << newAbsoluteID1 << ". new abs ID2: " << newAbsoluteID2 << std::endl;
                throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): new IDs dont match basenodes" );
            }

            correspondence[newRelativeID1] = newRelativeID2;
        }

        if( std::count( correspondence.begin(), correspondence.end(), nomatch) != 0 )
        {
            throw std::runtime_error( "ERROR @ treeCompare::simpleCorrespondence(): error in correspondence table" );
        }


    }
    else
    {
        correspondence = protoCorrespondence;
    }

    if( m_verbose )
    {
        std::cout << reportBaseNodes() << std::endl;
    }

    m_newCorrespondence = correspondence;

    m_newCorrespReverse.clear();
    m_newCorrespReverse.resize(m_newCorrespondence.size(),0);
    for( size_t i = 0; i < m_newCorrespondence.size(); ++i )
    {
        m_newCorrespReverse[m_newCorrespondence[i]]=i;
    }


    //////////////////////////7


//    // randomly match unmatched nodes
//    if(!leftNodes1.empty())
//    {
//        std::cout<<"there are "<<leftNodes1.size()<<" unmatched nodes due to physical distance, matching randomly"<<std::endl;
//        for (size_t i = 0; i < leftNodes1.size()-1; ++i )
//        {

//            size_t chosen( std::floor( u01Prng() * leftNodes2.size() ) );
//            if( chosen == leftNodes2.size() )
//                --chosen;

//            protoCorrespondence[leftNodes1[i]] = ( leftNodes2[chosen] );
//            leftNodes2.erase( leftNodes2.begin() + chosen );
//        }
//        if( leftNodes2.size() != 1 )
//        {
//            throw std::runtime_error( "ERROR @ treeCompare::randomCorrespondence(): error while matchign last rogue node" );
//        }
//        protoCorrespondence[leftNodes1.back()] = leftNodes2.front();
//    }


    return;
} // end treeComparer::randomCorrespondence() -------------------------------------------------------------------------------------



void treeComparer::writeFinalCorrespondence( std::string filename )
{

    if( m_baseNodes1.empty() || m_baseNodes1.size() != m_baseNodes2.size() )
    {
        std::cerr<< "ERROR @ treeCompare::writeCorrespondence(): base node vectors are empty or have different sizes"<<std::endl;
        return;
    }
    if( m_baseNodes1.size() != m_newCorrespondence.size() || m_baseNodes1.size() != m_correspDistances.size() )
    {
        std::cerr<< "ERROR @ treeCompare::writeCorrespondence(): correspondance vector size (" << m_newCorrespondence.size() << ")  or distances (" << m_correspDistances.size() << ") does not match basenodes vector size: "<< m_baseNodes1.size() <<std::endl;
        return;
    }
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR @ treeComparer::writeCorrespondence(): unable to open out file: \"" << outFile << "\"" << std::endl;
        return;
    }

    if( m_verbose )
    {
        std::cout<< "Writing down correspondence table in file: \""<< filename << "\"..." <<std::flush;
    }

    outFile<< "NodeTree1ID    NodeTree2ID    TractDistance    ClusterEuclideanDistance" << std::endl;
    outFile<< "#correspondence" << std::endl;
    for( size_t i = 0; i<m_newCorrespondence.size(); ++i )
    {
        outFile<< str( boost::format( "%06d" ) % m_baseNodes1[i] ) << " " << str( boost::format( "%06d" ) % m_baseNodes2[m_newCorrespondence[i]] ) << " " << m_correspDistances[i].first << " " << m_correspDistances[i].second << std::endl;
    }
    outFile<< "#endcorrespondence" << std::endl << std::endl;

    outFile<< "#relativecorresp" << std::endl;
    for( size_t i = 0; i<m_newCorrespondence.size(); ++i )
    {
        outFile<< str( boost::format( "%06d" ) % i ) << " " << str( boost::format( "%06d" ) % m_newCorrespondence[i] ) << std::endl;
    }
    outFile<< "#endrelativecorresp" << std::endl;

    if( m_verbose )
    {
        std::cout<< "Done" <<std::endl;
    }
    return;

} // end treeComparer::writeCorrespondence() -------------------------------------------------------------------------------------

void treeComparer::writeProtoCorrespondence( std::string filename )
{

    if( m_originalBaseNodes1.empty() ||m_originalBaseNodes2.empty() )
    {
        std::cerr<< "ERROR @ treeCompare::writeFullCorrespondence(): one or both of the base node vectors are empty"<<std::endl;
        return;
    }
    if( m_originalBaseNodes1.size() != m_fullCorrespondence.size() )
    {
        std::cerr<< "ERROR @ treeCompare::writeFullCorrespondence(): correspondance vector size (" << m_fullCorrespondence.size();
        std::cerr<< ") does not match basenodes 1 vector size (" << m_originalBaseNodes1.size() << ")"<<std::endl;
        return;
    }
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR @ treeComparer::writeFullCorrespondence(): unable to open out file: \"" << outFile << "\"" << std::endl;
        return;
    }

    if( m_verbose )
    {
        std::cout<< "Writing down correspondence table in file: \""<< filename << "\"..." <<std::flush;
    }

    outFile<< "#correspondence" << std::endl;
    for( size_t i = 0; i<m_fullCorrespondence.size(); ++i )
    {
        size_t abID1( m_originalBaseNodes1[i] );
        size_t nomatch( m_initialSizes.second ) ;

        long int abdID2;

        if( m_fullCorrespondence[i] == nomatch )
        {
            abdID2 = nomatch;
        }
        else
        {
            abdID2 = m_originalBaseNodes2[m_fullCorrespondence[i]];
        }
        outFile<< str( boost::format( "%06d" ) % abID1 ) << " " << str( boost::format( "%06d" ) %  abdID2 ) << std::endl;
    }
    outFile<< "#endcorrespondence" << std::endl << std::endl;

    outFile<< "#relativecorresp" << std::endl;
    for( size_t i = 0; i<m_fullCorrespondence.size(); ++i )
    {
        size_t abID1( i );
        size_t nomatch( m_initialSizes.second ) ;

        size_t abdID2( m_fullCorrespondence[i]);

        outFile<< str( boost::format( "%06d" ) % abID1 ) << " " << str( boost::format( "%06d" ) %  abdID2 ) << std::endl;
    }
    outFile<< "#endrelativecorresp" << std::endl;

    if( m_verbose )
    {
        std::cout<< "Done" <<std::endl;
    }
    return;

} // end treeComparer::writeFullCorrespondence() -------------------------------------------------------------------------------------



// PRIVATE member functionn



size_t treeComparer::findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t > &baseNodes ) const
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

size_t treeComparer::noiseBaseline( const bool treeCode, const float noiseAlpha )
{
    // initialize vectors
    WHtree* treePointer;
    std::vector< float >* noisePointer;
    std::vector< size_t >* basesPointer;
    std::vector< float > tractDists;

    if( treeCode == TREE1 )
    {
        treePointer = &m_tree1;
        noisePointer = &m_noiseLevels1;
        basesPointer = &m_baseNodes1;

        tractDists.reserve(m_correspDistances.size());
        for( size_t i=0; i < m_correspDistances.size(); ++i)
        {
            tractDists.push_back(m_correspDistances[i].first);
        }
    }
    else
    {
        treePointer = &m_tree2;
        noisePointer = &m_noiseLevels2;
        basesPointer = &m_baseNodes2;

        tractDists.reserve(m_newCorrespReverse.size());
        for( size_t i=0; i < m_newCorrespReverse.size(); ++i)
        {
            tractDists.push_back(m_correspDistances[m_newCorrespReverse[i]].first);
        }
    }



    WHtree& tree( *treePointer );
    std::vector< float >& noiseLevels( *noisePointer );
    std::vector< size_t >& baseNodes( *basesPointer );
    noiseLevels.clear();
    noiseLevels.resize( baseNodes.size(), 0 );


    if(noiseAlpha <= 0)
    {
        return baseNodes.size();
    }

    // start checking nodes recursively top-down

    std::list< size_t > worklist, flatSelection;
    worklist.push_back( tree.getRoot().getID() );
    size_t granCount(0);

    while( !worklist.empty() )
    {
        const WHnode& currentNode( tree.getNode( worklist.front() ) );
        worklist.pop_front();
        double noiseSum(0), sizeSum(0);


        // get contained base nodes for this cluster
        std::vector< size_t > currentBases( tree.getBaseNodes( currentNode.getID() ) );

        //obtaine noise level for this node
        for( size_t i=0; i < currentBases.size(); ++i)
        {
            // obtain relative base node ID
            size_t thisBaseID( findRelativeBasenodeID( currentBases[i], baseNodes ) );
            if( thisBaseID >= baseNodes.size() )
            {
                throw std::runtime_error(" basenode ID not found within bases ");
            }
            size_t thisSize( tree.getNode( currentBases[i] ).getSize() );
            noiseSum += tractDists[thisBaseID] * thisSize;
            sizeSum += thisSize;
        }
        double currentNoise( (noiseSum / sizeSum) * noiseAlpha );

        if( currentNode.getDistLevel() >= currentNoise )
        {
            // node is avobe the noise level
            if( currentNode.getHLevel() == 1 )
            {
                // if its a base node, we simply record its noise level
                size_t thisBaseID( findRelativeBasenodeID( currentNode.getID(), baseNodes ) );
                if( thisBaseID >= baseNodes.size() )
                {
                    throw std::runtime_error(" basenode ID not found within bases ");
                }
                noiseLevels[thisBaseID] = currentNoise;
                ++granCount;
            }
            else
            {
                // we continue checking its children
                std::vector<nodeID_t> kids(currentNode.getChildren());
                for( size_t i=0; i < kids.size(); ++i)
                {
                    //only include them if they are not real leaves
                    if( kids[i].first )
                    {
                        worklist.push_back( kids[i].second );
                    }
                }
            }
        }
        else
        {
            // node is below its current level, we flatten the inner structure of the node

            // if parent node is lower than the nose, equal the nose to the parent (parent had benn previously identified as avobe noise)
            float parentLevel( tree.getNode( currentNode.getParent() ).getDistLevel() );
            if( parentLevel < currentNoise )
            {
                currentNoise = parentLevel;
            }

            WHnode* currentModifiable( tree.fetchNode( currentNode.getID() ) );
            currentModifiable->setDistLevel( currentNoise );

            // mark that node for flattening
            flatSelection.push_back( currentNode.getID() );
            ++granCount;

            // record noise levels of contained leaves
            for( size_t i=0; i < currentBases.size(); ++i)
            {
                // obtain relative base node ID
                size_t thisBaseID( findRelativeBasenodeID( currentBases[i], baseNodes ) );
                if( thisBaseID >= baseNodes.size() )
                {
                    throw std::runtime_error(" basenode ID not found within bases ");
                }
                noiseLevels[thisBaseID] = currentNoise;
            }
        }

    } // end of while

    // eliminate structure below nose level
    WHtreeProcesser treeProc(&tree);
    size_t nodesLost( treeProc.flattenSelection( flatSelection, true ) );

    return granCount;
} // end treeComparer::noiseBaseline() -------------------------------------------------------------------------------------



