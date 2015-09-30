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
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <cmath>

#include "treeManager.h"

treeManager::treeManager( WHtree* const tree, bool verbose ) :
    m_tree( *tree ), m_logfile( 0 ), m_verbose( verbose ), m_floatFlag( true ), m_zipFlag( false )
{
    fileManagerFactory fMFtestFormat;
    m_niftiMode = fMFtestFormat.isNifti();
    if( m_tree.getDataGrid() == HC_NIFTI && !m_niftiMode )
    {
        m_tree.convert2grid( HC_VISTA );
    }
    else if( m_tree.getDataGrid() == HC_VISTA && m_niftiMode )
    {
        m_tree.convert2grid( HC_NIFTI );
    }
}

void treeManager::writeDebugTree() const
{
    if( m_outputFolder.empty() )
    {
        std::cerr << "Location of output folder has not been specified, please initialize with treeManager::setOutputFolder()"
                        << std::endl;
        return;
    }
    std::string debugTreeFilename( m_outputFolder + "/debugTree.txt" );
    m_tree.writeTreeDebug( debugTreeFilename );
} // end "writeDebugTree()" -----------------------------------------------------------------


float treeManager::doCpcc()
{
    if( m_distMatrixFolder.empty() )
    {
        std::cerr
                        << "Location of of distance matrix has not been specified, please initialize with treeManager::setDistMatrixFolder()"
                        << std::endl;
        return -1;
    }


    double sumT( 0 ), sumM( 0 ), sqT( 0 ), sqM( 0 ), sumProd( 0 );
    bool firstIteration( true );
    size_t doneCount( 0 );
    size_t K( m_tree.m_coordinates.size() * ( m_tree.m_coordinates.size() - 1 ) / 2. );

    // initialize tracking iterators
    std::pair< std::vector< WHcoord >::iterator, std::vector< WHcoord >::iterator > rangeRoiRow( m_tree.m_coordinates.begin(),
                    m_tree.m_coordinates.begin() );
    std::pair< std::vector< WHcoord >::iterator, std::vector< WHcoord >::iterator > rangeRoiColumn( m_tree.m_coordinates.begin(),
                    m_tree.m_coordinates.begin() );

    // load distance block index and first block block
    if( m_verbose )
        std::cout << "Reading distance matrix index..." << std::flush;
    distBlock dBlock( m_distMatrixFolder );
    if( !dBlock.indexReady() )
    {
        std::cerr << "ERROR: distance matrix index did not load" << std::endl;
        exit( -1 );
    }
    unsigned int topBlock( dBlock.topBlock() );
    unsigned int numBlocks( dBlock.numBlocks() );
    if( m_verbose )
        std::cout << "OK. Whole matrix roi is " << dBlock.matrixSize() << " elements. " << topBlock + 1 << "x" << topBlock + 1
                        << " blocks (real blocks: " << numBlocks << "). " << std::flush;

    time_t loopStartTime( time( NULL ) );

    // obtain sums for all elements
    while( rangeRoiRow.first != m_tree.m_coordinates.end() )
    {
        // load block containing the next leaves in the tree
        dBlock.loadBlock( *( rangeRoiRow.first ), *( rangeRoiColumn.first ) );
        //std::cout<<std::endl<<"loading "<<*(rangeRoiRow.first)<<" and "<<*(rangeRoiColumn.first)<<std::endl;
        std::pair< unsigned int, unsigned int > blockID( dBlock.blockID() );

        if( firstIteration )
        {
            std::cout << "Block size: " << dBlock.size() << "x" << dBlock.size() << std::endl;
            if( m_verbose )
                std::cout << "Computing block: " << blockID.first << "-" << blockID.second << "..." << std::flush;
            firstIteration = false;
        }
        else
        {
            if( m_verbose )
            {
                std::cout << "\rComputing block: " << blockID.first << "-" << blockID.second << "..." << std::flush;
                float progress = ( doneCount * 100. ) / ( K );
                time_t current_time( time( NULL ) );
                std::cout << ( int )progress << " % completed. Expected remaining time: ";
                if( progress > 0 )
                {
                    int expected_remain( difftime( current_time, loopStartTime ) * ( ( 100. - progress ) / progress ) );
                    std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                    % 3600 ) % 60 ) << "\"  " << std::flush;
                }
            }
        }

        // get range of block and find the position of the first leaf in the tree not contained in the block (for both rows and columns)
        std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > currentRange( dBlock.getBlockRange() );
        std::pair< WHcoord, WHcoord > rangeBlockRow( currentRange.first );
        std::pair< WHcoord, WHcoord > rangeBlockColumn( currentRange.second );

        rangeRoiRow.second = std::lower_bound( rangeRoiRow.first, m_tree.m_coordinates.end(), rangeBlockRow.second );
        if( *( rangeRoiRow.second ) == rangeBlockRow.second )
            ++( rangeRoiRow.second );

        if( rangeRoiColumn.first == rangeRoiRow.first )
        {
            rangeRoiColumn.second = rangeRoiRow.second;
        }
        else
        {
            rangeRoiColumn.second = std::lower_bound( rangeRoiColumn.first, m_tree.m_coordinates.end(), rangeBlockColumn.second );
            if( *( rangeRoiColumn.second ) == rangeBlockColumn.second )
                ++( rangeRoiColumn.second );
        }

        size_t beginPosRow( rangeRoiRow.first - m_tree.m_coordinates.begin() );
        size_t endPosRow( rangeRoiRow.second - m_tree.m_coordinates.begin() );
        size_t beginPosCol( rangeRoiColumn.first - m_tree.m_coordinates.begin() );
        size_t endPosCol( rangeRoiColumn.second - m_tree.m_coordinates.begin() );
        size_t sizeRow( rangeRoiRow.second - rangeRoiRow.first );
        size_t sizeCol( rangeRoiColumn.second - rangeRoiColumn.first );

        std::vector< std::vector< dist_t > > treeDistM;
        std::vector< std::vector< dist_t > > matrixDistM;
        {
            std::vector< dist_t > tempVector( sizeCol, 0 );
            treeDistM.resize( sizeRow, tempVector );
            matrixDistM = treeDistM;
        }

        if( beginPosRow == beginPosCol )
        {
            // its a block on the diagonal, only use values over the diagonal
#pragma omp parallel for schedule(guided)
            for( size_t i = beginPosRow; i < endPosRow; ++i )
            {
                for( size_t j = i + 1; j < endPosCol; ++j )
                {
                    matrixDistM[i - beginPosRow][j - beginPosCol] = ( dBlock.getDistance( m_tree.getCoordinate4leaf( i ),
                                    m_tree.getCoordinate4leaf( j ) ) );
                    treeDistM[i - beginPosRow][j - beginPosCol] = ( m_tree.getLeafDistance( i, j ) );
                }
            }
        }
        else
        {
            // its NOT a block on the diagonal, use all values
#pragma omp parallel for schedule(guided)
            for( size_t i = beginPosRow; i < endPosRow; ++i )
            {
                for( size_t j = beginPosCol; j < endPosCol; ++j )
                {
                    matrixDistM[i - beginPosRow][j - beginPosCol] = ( dBlock.getDistance( m_tree.getCoordinate4leaf( i ),
                                    m_tree.getCoordinate4leaf( j ) ) );
                    treeDistM[i - beginPosRow][j - beginPosCol] = ( m_tree.getLeafDistance( i, j ) );
                }
            }
        }

#pragma omp parallel sections
        {
#pragma omp section
            {
                double preSumM( 0 );
                for( size_t i = 0; i < sizeRow; ++i )
                    for( size_t j = 0; j < sizeCol; ++j )
                        preSumM += matrixDistM[i][j];
                sumM += preSumM;
            }
#pragma omp section
            {
                double preSqM( 0 );
                for( size_t i = 0; i < sizeRow; ++i )
                    for( size_t j = 0; j < sizeCol; ++j )
                        preSqM += matrixDistM[i][j] * matrixDistM[i][j];
                sqM += preSqM;
            }
#pragma omp section
            {
                double preSumT( 0 );
                for( size_t i = 0; i < sizeRow; ++i )
                    for( size_t j = 0; j < sizeCol; ++j )
                        preSumT += treeDistM[i][j];
                sumT += preSumT;
            }
#pragma omp section
            {
                double preSqT( 0 );
                for( size_t i = 0; i < sizeRow; ++i )
                    for( size_t j = 0; j < sizeCol; ++j )
                        preSqT += treeDistM[i][j] * treeDistM[i][j];
                sqT += preSqT;
            }
#pragma omp section
            {
                double preSumProd( 0 );
                for( size_t i = 0; i < sizeRow; ++i )
                    for( size_t j = 0; j < sizeCol; ++j )
                        preSumProd += treeDistM[i][j] * matrixDistM[i][j];
                sumProd += preSumProd;
            }
        }

        if( beginPosRow == beginPosCol )
        {
            doneCount += ( sizeRow * ( sizeRow - 1 ) / 2. );
        }
        else
        {
            doneCount += ( sizeRow * sizeCol );
        }

        if( rangeRoiColumn.second == m_tree.m_coordinates.end() )
        {
            // if we just finished the last column of the roi, move to the diagonal element of the next row
            rangeRoiRow.first = rangeRoiRow.second;
            rangeRoiColumn.first = rangeRoiRow.first;
        }
        else
        {
            // move to the next column
            rangeRoiColumn.first = rangeRoiColumn.second;
        }
    }

    if( m_verbose )
        std::cout << "\rAll " << topBlock + 1 << "x" << topBlock + 1 << " blocks processed, doing final caculations..."
                        << std::flush;

    // do the final computations
    double meanM( sumM / K );
    double meanT( sumT / K );
    double numerator( ( sumProd / K ) - ( meanM * meanT ) );
    double denominator1( ( sqM / K ) - ( meanM * meanM ) );
    double denominator2( ( sqT / K ) - ( meanT * meanT ) );
    float CPCC( numerator / sqrt( denominator1 * denominator2 ) );

    if( m_verbose )
        std::cout << "Done. CPCC: " << boost::lexical_cast< std::string >( CPCC ) << std::endl;

    m_tree.m_cpcc = CPCC;

    return CPCC;
} // end "cpcc()" -----------------------------------------------------------------


compactTract treeManager::getMeanTract( const size_t inNode ) const
{
    if( m_leavesTractFolder.empty() )
        throw std::runtime_error(
                        "Single tracts folder has not been specified, please initialize with treeManager::setSingleTractFolder()" );



    //create file reader class
    fileManagerFactory fileMF( m_leavesTractFolder );
    fileManager& fileSingle( fileMF.getFM() );
    fileSingle.readAsLog();
    fileSingle.readAsUnThres();

    //loop through nodes
    compactTract testTract;
    std::vector< size_t > nodeLeaves( m_tree.getLeaves4node( inNode) );

    fileSingle.readLeafTract( nodeLeaves.front(), m_tree.m_trackids, m_tree.m_coordinates, &testTract );

    std::vector< size_t >sumVector( testTract.size(), 0 );



    for( size_t j = 0; j < nodeLeaves.size(); ++j )
    {
        compactTract tempTract;
        size_t leafID = nodeLeaves[j];

        fileSingle.readLeafTract( leafID, m_tree.m_trackids, m_tree.m_coordinates, &tempTract );
        tempTract.unLog( m_tree.m_logFactor ); //to sum tractograms they must be in natural units

        std::vector< float > tempVector( tempTract.tract() );

        for( size_t i = 0; i <  tempVector.size(); ++i )
        {
            sumVector[i] += lroundf( tempVector[i] );
        }
    }
    std::vector< float >meanFloat( sumVector.size(), 0 );
    double divisor( nodeLeaves.size() );
    for( size_t i = 0; i <  sumVector.size(); ++i )
    {
        meanFloat[i] = sumVector[i]/divisor;
    }

    compactTract sumTract( meanFloat );
    sumTract.unLog( 0 ); // set as in natural units
    sumTract.doLog( m_tree.m_logFactor ); // we want to view it in logarithmic units

    return sumTract;
} // end "getMeanTract()" -----------------------------------------------------------------


void treeManager::writeFullTract( const nodeID_t input )
{
    std::vector< nodeID_t > inNodes( 1, input );
    writeFullTract( inNodes );
    return;
}
void treeManager::writeFullTract( std::vector< size_t > input )
{
    std::vector< nodeID_t > inNodes;
    inNodes.reserve(input.size());

    for( size_t i = 0; i < input.size(); ++i )
    {
        inNodes.push_back(std::make_pair(true, input[i]));
    }
    writeFullTract( inNodes );
    return;
}
void treeManager::writeFullTract( std::vector< nodeID_t > input )
{


    if( m_maskFilename.empty() )
    {
        std::cerr << "Location of mask file has not been specified, please initialize with treeManager::setMaskFilename()"
                        << std::endl;
        return;
    }
    if( m_fullTractFolder.empty() )
    {
        std::cerr
                        << "Location of full tracts folder has not been specified, please initialize with treeManager::setFullTractFolder()"
                        << std::endl;
        return;
    }
    if( m_leavesTractFolder.empty() && m_meanTractFolder.empty() )
    {
        std::cerr << "Neither single tracts nor mean tracts folder has been specified, please initialize with"
                  << " treeManager::setSingleTractFolder() or treeManager::setMeanTractFolder()" << std::endl;
        return;
    }

    std::vector< size_t > inNodes, inLeaves;
    for( size_t i = 0; i < input.size(); ++i )
    {
        if( input[i].first )
        {
            inNodes.push_back( input[i].second );
        }
        else
        {
            inLeaves.push_back( input[i].second );
        }
    }



    if( !inLeaves.empty() )
    {
        if( m_leavesTractFolder.empty() )
        {
            std::cerr
                            << "Single tracts folder has not been specified, please initialize with treeManager::setSingleTractFolder()"
                            << std::endl;
            return;
        }

        fileManagerFactory leavesFileMF( m_leavesTractFolder );
        fileManager& leavesFM( leavesFileMF.getFM() );
        leavesFM.readAsLog();
        leavesFM.readAsUnThres();


        fileManagerFactory fullTractFileMF( m_fullTractFolder );
        fileManager& fullTractFM( fullTractFileMF.getFM() );
        if ( m_floatFlag )
        {
            fullTractFM.writeInFloat();
        }
        else
        {
            fullTractFM.writeInChar();
        }
        if ( m_zipFlag )
        {
            fullTractFM.storeZipped();
        }
        else
        {
            fullTractFM.storeUnzipped();
        }

        fullTractFM.loadMaskImage( m_maskFilename ); //load full tract mask (for writing)

        //loop through leaves
        for( size_t i = 0; i < inLeaves.size(); ++i )
        {
            size_t thisLeafID ( inLeaves[i] );
            compactTract leaftract;
            leavesFM.readLeafTract(  thisLeafID, m_tree.m_trackids, m_tree.m_coordinates, &leaftract );
            fullTractFM.writeFullLeafTract( thisLeafID, m_tree.m_trackids, m_tree.m_coordinates, leaftract );
            std::string tractFilename( fullTractFM.getFullLeafTractFilename( thisLeafID, m_tree.m_trackids, m_tree.m_coordinates ) );
            if( m_verbose )
            {
                std::cout << " Full tract for leaf " << thisLeafID << "(" << m_tree.getCoordinate4leaf( thisLeafID )  << ") written in \"" << tractFilename
                                << "\"" << std::endl;
            }
            if( m_logfile != 0 )
            {
                ( *m_logfile ) << " Full tract for leaf " << thisLeafID << "(" << m_tree.getCoordinate4leaf( thisLeafID ) << ") written in \""
                                << tractFilename << "\"" << std::endl;
            }
        }
    }

    if( !inNodes.empty() )
    {
        bool reportProgress( inNodes.size() > 25 );


        if( m_verbose )
            std::cout << "loading leaves for each node..." << std::flush;
        m_tree.loadContainedLeaves();
        if( m_verbose )
            std::cout << "Done" << std::endl;

        // checking validity of input nodes
        for( std::vector< size_t >::iterator iter( inNodes.begin() ); iter != inNodes.end(); )
        {
            if( *iter >= m_tree.getNumNodes() )
            {
                std::cerr << "WARNING @ treeManager::writeFullTract(): input node " << *iter
                                << " does not correspond to a node on the tree. It will be ignored" << std::endl;
                iter = inNodes.erase( iter );
            }
            else
            {
                ++iter;
            }
        }

        if( !m_meanTractFolder.empty() )
        {
            // meanTract already exists in file

            if( m_verbose )
                std::cout << "Obtaining mean tractograms directly from file" << std::endl;
            //create writer class
            fileManagerFactory nodeFileMF( m_meanTractFolder );
            fileManager& nodeFM( nodeFileMF.getFM() );
            nodeFM.readAsLog();
            nodeFM.readAsUnThres();


            fileManagerFactory nodeFullFileMF( m_fullTractFolder );
            fileManager& nodeFullMF( nodeFullFileMF.getFM() );
            if ( m_floatFlag )
            {
                nodeFullMF.writeInFloat();
            }
            else
            {
                nodeFullMF.writeInChar();
            }
            if ( m_zipFlag )
            {
                nodeFullMF.storeZipped();
            }
            else
            {
                nodeFullMF.storeUnzipped();
            }

            nodeFullMF.loadMaskImage( m_maskFilename ); //load full tract mask (for writing)

            //loop through nodes
            for( size_t i = 0; i < inNodes.size(); ++i )
            {
                compactTract nodetract;
                nodeFM.readNodeTract( inNodes[i], &nodetract );
                nodeFullMF.writeFullNodeTract( inNodes[i], nodetract );
                std::string tractFilename( nodeFullMF.getFullNodeTractFilename( inNodes[i] ) );
                if( m_verbose )
                {
                    std::cout << " Full tract for node " << inNodes[i] << " with " << m_tree.getNode( inNodes[i] ).getSize()
                                    << " leaves written in \"" << tractFilename << "\"" << std::endl;
                }
                if( m_logfile != 0 )
                {
                    ( *m_logfile ) << "Full tract for node " << inNodes[i] << " with " << m_tree.getNode( inNodes[i] ).getSize()
                                    << " leaves written in \"" << tractFilename << "\"" << std::endl;
                }
            }
        }
        else
        {
            // we must obtain the mean tractograms from the single tracts

            if( m_verbose )
                std::cout << "Computing mean tractograms from single tracts" << std::endl;

            // if theres one more node to compute, precompute coordinates of all tree
            if( inNodes.size() > 1 )
                m_tree.loadContainedLeaves();

            //create  writer class
            fileManagerFactory fullTractFileMF( m_fullTractFolder );
            fileManager& fullTractFM( fullTractFileMF.getFM() );
            if ( m_floatFlag )
            {
                fullTractFM.writeInFloat();
            }
            else
            {
                fullTractFM.writeInChar();
            }
            if ( m_zipFlag )
            {
                fullTractFM.storeZipped(); //zip tract files after storing
            }
            else
            {
                fullTractFM.storeUnzipped(); //dont zip tract files after storing
            }

            fullTractFM.loadMaskImage( m_maskFilename ); //load full tract mask (for writing)

            time_t startTime( time( NULL ) ), lastTime( time( NULL ) );

            size_t progCount ( 0 );

            //loop through nodes
            #pragma omp parallel for schedule( dynamic, 1 )
            for( size_t i = 0; i < inNodes.size(); ++i )
            {
                compactTract meanTract( getMeanTract( inNodes[i] ) );

                fullTractFM.writeFullNodeTract( inNodes[i], meanTract );
                std::string tractFilename( fullTractFM.getFullNodeTractFilename( inNodes[i] ) );

                if( reportProgress )
                {
                    #pragma omp atomic
                    ++progCount;

                    #pragma omp single nowait // only one thread executes output
                    if( m_verbose )
                    {
                        time_t currentTime( time( NULL ) );
                        if( currentTime - lastTime > 1 )
                        {
                            // output progress every 2 seconds
                            lastTime = currentTime;
                            size_t currentCount( progCount );
                            float progress = ( currentCount ) * 100. / (inNodes.size() );
                            size_t expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                            std::cout << "\r" << ( int )progress << " % Completed (" << currentCount << " node tracts)"
                                            << ". Expected remaining time: ";
                            std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                            % 3600 ) % 60 ) << "\"  " << std::flush;
                        }
                    }
                }
                else
                {
                    if( m_verbose )
                        std::cout << " Mean tract of node " << inNodes[i] << " (" << m_tree.getNode( inNodes[i] ).getSize()
                                        << " tracts) written in \"" << tractFilename << "\"" << std::endl;

                    if( m_logfile != 0 )
                        ( *m_logfile ) << "Mean tract for node " << inNodes[i] << " with " << m_tree.getNode( inNodes[i] ).getSize()
                                        << " leaves written in \"" << tractFilename << "\"" << std::endl;
                }
            }

            if( reportProgress && m_verbose )
            {
                size_t timeTaken( difftime( time( NULL ), startTime ) );
                std::cout << "\r100% Completed (" << progCount << " node tracts). Time taken: "
                          << timeTaken / 3600 << "h " << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken
                                             % 3600 ) % 60 ) << "\"  " << std::endl;
            }

            if( inNodes.size() > 1 )
                m_tree.clearContainedLeaves();
        }
    }

    return;
} // end "writeFullTract()" -----------------------------------------------------------------


void treeManager::writeClusterMasks( const nodeID_t input )
{
    std::vector< size_t > inNodes( input.second );
    writeClusterMasks( inNodes );
    return;
}
void treeManager::writeClusterMasks( std::vector< nodeID_t > input )
{
    std::vector< size_t > inNodes;
    inNodes.reserve(input.size());

    for( size_t i = 0; i < input.size(); ++i )
    {
        inNodes.push_back(input[i].second);
    }
    writeClusterMasks( inNodes );
    return;
}
void treeManager::writeClusterMasks( std::vector< size_t > input )
{

    if( m_verbose )
        std::cout << "Writing cluster masks" << std::endl;

    fileManagerFactory maskFileMF( m_fullTractFolder );
    fileManager& maskManager( maskFileMF.getFM() );
    if( m_zipFlag )
    {
        maskManager.storeZipped();
    }
    else
    {
        maskManager.storeUnzipped();
    }
    if( !m_maskFilename.empty() )
    {
        maskManager.loadMaskImage( m_maskFilename );
    }

    time_t startTime( time( NULL ) ), lastTime( time( NULL ) );
    size_t progCount( 0 );


    // write cluster masks
    for( size_t i = 0; i < input.size(); ++i )
    {
        std::vector<std::vector< std::vector<bool> > > maskMatrix( m_tree.getDataSize().m_x,
                    std::vector< std::vector<bool> > ( m_tree.getDataSize().m_y,
                                 std::vector<bool>   ( m_tree.getDataSize().m_z, false ) ) );

        std::vector<WHcoord> nodeCoords( m_tree.getCoordinates4node( input[i] ) );
        for( size_t j = 0; j < nodeCoords.size(); ++j )
        {
            maskMatrix[nodeCoords[j].m_x][nodeCoords[j].m_y][nodeCoords[j].m_z] = true;
        }
        std::string maskFilename( maskManager.getClusterMaskFilename( input[i] ) );

        maskManager.writeMask( maskFilename, maskMatrix );
        ++progCount;

        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                // output progress every 2 seconds
                lastTime = currentTime;
                size_t currentCount( progCount );
                float progress = ( currentCount ) * 100. / (input.size() );
                size_t expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                std::cout << "\r" << ( int )progress << " % Completed (" << currentCount << " node masks)"
                                << ". Expected remaining time: ";
                std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                % 3600 ) % 60 ) << "\"  " << std::flush;
            }
        }
    }


    if( m_verbose )
    {
        size_t timeTaken( difftime( time( NULL ), startTime ) );
        std::cout << "\r100% Completed (" << progCount << " node masks). Time taken: "
                  << timeTaken / 3600 << "h " << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken
                                     % 3600 ) % 60 ) << "\"  " << std::endl;
    }

return;
} // end "writeClusterMasks()" -----------------------------------------------------------------



void treeManager::writeFullBaseNodeTracts( bool onlyMasks )
{
    if( m_fullTractFolder.empty() )
    {
        std::cerr
                        << "Location of full tracts folder has not been specified, please initialize with treeManager::setFullTractFolder()"
                        << std::endl;
        return;
    }

    if( m_maskFilename.empty() )
    {
        std::cerr << "Location of mask file has not been specified, please initialize with treeManager::setMaskFilename()"
                        << std::endl;
        return;
    }

    if( m_verbose )
        std::cout << "Getting tree base nodes" << std::endl;

    if( !m_tree.testRootBaseNodes() )
    {
        std::cerr << "ERROR: Base nodes are not all the the lowest level" << std::endl;
        return;
    }
    std::vector<size_t> baseNodes = m_tree.getRootBaseNodes();


    if( m_verbose )
        std::cout << "Writing cluster masks" << std::endl;

    fileManagerFactory maskFileMF( m_fullTractFolder );
    fileManager& maskManager( maskFileMF.getFM() );
    if( m_zipFlag )
    {
        maskManager.storeZipped();
    }
    else
    {
        maskManager.storeUnzipped();
    }

    time_t startTime( time( NULL ) ), lastTime( time( NULL ) );
    size_t progCount( 0 );



    // write cluster masks
    for( size_t i = 0; i < baseNodes.size(); ++i )
    {
        std::vector<std::vector<std::vector<bool> > > maskMatrix( m_tree.getDataSize().m_z,
                    std::vector<std::vector<bool> > ( m_tree.getDataSize().m_y, std::vector<bool> ( m_tree.getDataSize().m_x, false ) ) );

        std::vector<WHcoord> nodeCoords( m_tree.getCoordinates4node( baseNodes[i] ) );

        for( size_t j = 0; j < nodeCoords.size(); ++j )
        {
            maskMatrix[nodeCoords[j].m_z][nodeCoords[j].m_y][nodeCoords[j].m_x] = true;
        }

        std::string maskFilename( maskManager.getClusterMaskFilename( baseNodes[i] ) );
        maskManager.writeMask( maskFilename, maskMatrix );
        ++progCount;

        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                // output progress every 2 seconds
                lastTime = currentTime;
                size_t currentCount( progCount );
                float progress = ( currentCount ) * 100. / (baseNodes.size() );
                size_t expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                std::cout << "\r" << ( int )progress << " % Completed (" << currentCount << " node masks)"
                                << ". Expected remaining time: ";
                std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                % 3600 ) % 60 ) << "\"  " << std::flush;
            }
        }
    }


    if( m_verbose )
    {
        size_t timeTaken( difftime( time( NULL ), startTime ) );
        std::cout << "\r100% Completed (" << progCount << " node masks). Time taken: "
                  << timeTaken / 3600 << "h " << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken
                                     % 3600 ) % 60 ) << "\"  " << std::endl;
    }

    if( onlyMasks )
    {
        return;
    }

    if( m_verbose )
        std::cout << "Writing base nodes full mean tracts" << std::endl;

    writeFullTract( baseNodes );

return;
} // end "writeFullBaseNodeTracts()" -----------------------------------------------------------------



void treeManager::writeAllNodeTracts( const float memory ) const
{

    if( m_leavesTractFolder.empty() )
    {
        std::cerr
                        << "Location of single tracts folder has not been specified, please initialize with treeManager::setMeanTractFolder()"
                        << std::endl;
        return;
    }
    if( m_meanTractFolder.empty() )
    {
        std::cerr
                        << "Location of mean tracts folder has not been specified, please initialize with treeManager::setSingleTractFolder()"
                        << std::endl;
        return;
    }

    time_t lastTime( time( NULL ) );
    time_t loopStartTime( time( NULL ) );
    size_t tractProg( 0 );

    // compute cache size
    float tractMb( 0 );
    {
        fileManagerFactory leafFileMF( m_leavesTractFolder );
        fileManager& leafFM( leafFileMF.getFM() ); // object to read single tractograms
        leafFM.readAsUnThres();
        leafFM.readAsLog();
        compactTract tempTract;
        leafFM.readLeafTract( 0, m_tree.m_trackids, m_tree.m_coordinates, &tempTract );
        tractMb = tempTract.mBytes();
        if( m_verbose )
            std::cout << "Tractogram size is: " << tempTract.size() << " (" << tractMb << " MB)" << std::endl;
        if( m_logfile != 0 )
            ( *m_logfile ) << "Tractogram size:\t" << tempTract.size() << " (" << tractMb << " MB)" << std::endl;
    }
    size_t cacheSize( memory * 1024 / tractMb );
    if( m_verbose )
        std::cout << "Cache size is: " << cacheSize << " tracts" << std::endl;
    if( m_logfile != 0 )
        ( *m_logfile ) << "Cache size:\t" << cacheSize << " tracts" << std::endl;

    listedCache< compactTract > cache( m_tree.getNumNodes() );
    cache.setLimit( cacheSize );
    size_t maxLeaves( 2 * cacheSize ); // we can hold cacheSize nodes in cache, which means process 2*cacheSize leaves at a time

    // divide the three into branches with the most similar possible sizes
    size_t maxSize( maxLeaves );


    if( omp_get_max_threads() > 1 )
    maxSize = std::min( maxLeaves/omp_get_max_threads(), m_tree.getNumLeaves()/3*omp_get_max_threads() );


    // divide the three into branches with the maximum size equal to cache
    std::vector< nodeID_t > partition;
    WHtreePartition partitioner(&m_tree);

    size_t biggest = partitioner.partitionClassic( maxSize, &partition, HTP_SIZE, HTC_VALUE, false, m_tree.getRoot().getID() );
    m_tree.sortBySize( &partition );
    std::reverse( partition.begin(), partition.end() ); // put the biggest cluster to be processed first (makes better use of parallelism)

    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Branches:\t" << partition.size() << "+1" << std::endl;
        ( *m_logfile ) << "Branch max. size:\t" << biggest << std::endl;
    }

    if( m_verbose )
    {
        std::cout << "Dividing task into " << partition.size() << " branches of max. size " << biggest << std::endl;
        if( partition.size() > cacheSize )
            std::cout << "WARNING @ treeManager::writeAllNodeTracts(): partition (" << partition.size()
                            << ") is bigger than cache " << cacheSize << " memory problems might arise" << std::endl;
    }

#pragma omp parallel for schedule ( dynamic, 1 ) // use parallel threads
    for( size_t i = 0; i < partition.size(); ++i )
    {
        // If its a leaf, skip it
        if( !partition[i].first )
            continue;

        // get the vector of nodes from this branch of the tree
        std::list< size_t > worklist;
        std::vector< size_t > partNodes;
        partNodes.reserve( m_tree.getNode( partition[i] ).getSize() );
        worklist.push_back( partition[i].second );
        partNodes.push_back( partition[i].second );
        while( !worklist.empty() )
        {
            std::vector< nodeID_t > kids( m_tree.getNode( worklist.front() ).getChildren() );
            worklist.pop_front();
            for( std::vector< nodeID_t >::iterator iter( kids.begin() ); iter != kids.end(); ++iter )
            {
                if( iter->first )
                {
                    worklist.push_back( iter->second );
                    partNodes.push_back( iter->second );
                }
            }
        }
        std::sort( partNodes.begin(), partNodes.end() );
        writeNodeTracts( &partNodes, &cache, &tractProg, &lastTime, loopStartTime );
    }

    if( partition.size() > 1 )
    {
        if( m_verbose )
            std::cout << std::endl << " Final round" << std::endl;

        // get remainder treetop nodes
        std::vector< size_t > treetop;
        {
            std::list< size_t > worklist, outlist;
            for( size_t i = 0; i < partition.size(); ++i )
                worklist.push_back( ( m_tree.getNode( partition[i] ).getParent() ).second );
            worklist.sort();
            worklist.unique(); // eliminate duplicated nodes
            while( !worklist.empty() )
            {
                size_t current( worklist.front() );
                worklist.pop_front();
                outlist.push_back( current );
                if( m_tree.getNode( current ).isRoot() )
                    continue;
                worklist.push_back( ( m_tree.getNode( current ).getParent() ).second );
            }
            outlist.sort();
            outlist.unique();

            treetop.reserve( outlist.size() );
            for( std::list< size_t >::iterator iter( outlist.begin() ); iter != outlist.end(); ++iter )
                treetop.push_back( *iter );
        }
        writeNodeTracts( &treetop, &cache, &tractProg, &lastTime, loopStartTime );
    }
    if( m_verbose )
    {
        float progress = ( tractProg ) * 100. / ( m_tree.m_nodes.size() );
        std::cout << "\r" << ( int )progress << " % Completed (" << tractProg << " node tracts)" << std::endl;
    }
    if( m_logfile != 0 )
        ( *m_logfile ) << "Written tractograms:\t" << tractProg << std::endl;

    return;
} // end "writeAllNodeTracts()" -----------------------------------------------------------------

void treeManager::writeMeanTracts( std::vector< size_t > inNodes ) const
{
    size_t totalSize( 0 ), progSize( 0 );
    for( size_t i = 0; i < inNodes.size(); ++i )
    {
        totalSize += m_tree.getNode( inNodes[i] ).getSize();
    }

    if( m_verbose )
        std::cout << "Writing mean tracts for " << inNodes.size() << " nodes, containing a total of " << totalSize
                        << " single voxels..." << std::endl;

    if( m_leavesTractFolder.empty() )
    {
        std::cerr
                        << "Location of single tracts folder has not been specified, please initialize with treeManager::setMeanTractFolder()"
                        << std::endl;
        return;
    }
    if( m_meanTractFolder.empty() )
    {
        std::cerr
                        << "Location of mean tracts folder has not been specified, please initialize with treeManager::setSingleTractFolder()"
                        << std::endl;
        return;
    }

    fileManagerFactory meanTractFileMF( m_meanTractFolder );
    fileManager& meanTractFM( meanTractFileMF.getFM() );
    if ( m_floatFlag )
    {
        meanTractFM.writeInFloat();
    }
    else
    {
        meanTractFM.writeInChar();
    }
    if ( m_zipFlag )
    {
        meanTractFM.storeZipped();
    }
    else
    {
        meanTractFM.storeUnzipped();
    }

    m_tree.sortBySize( &inNodes );
    std::reverse( inNodes.begin(), inNodes.end() ); // put the biggest cluster to be processed first (makes better use of parallelism)


    time_t lastTime( time( NULL ) ), startTime( time( NULL ) );

#pragma omp parallel for schedule(guided)
    for( size_t i = 0; i < inNodes.size(); ++i )
    {
        compactTract meanTract( getMeanTract( inNodes[i] ) );
        meanTractFM.writeNodeTract( inNodes[i], meanTract );

#pragma omp critical
        progSize += m_tree.getNode( inNodes[i] ).getSize(); // increment thread counter

#pragma omp single nowait // only one thread executes output
        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                // output progress every 2 seconds
                lastTime = currentTime;
                size_t currentCount( progSize );
                float progress = ( currentCount ) * 100. / ( totalSize );
                int expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                std::cout << "\r" << ( int )progress << " % Completed (" << currentCount << " single voxels accounted for)"
                                << ". Expected remaining time: ";
                std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                % 3600 ) % 60 ) << "\"  " << std::flush;
            }
        }
    }

    std::cout << "\r100 % Completed (" << inNodes.size() << " mean tracts containing a total of " << totalSize
                    << " single voxels)    " << std::endl;
} // end "writeMeanTracts()" -----------------------------------------------------------------


void treeManager::writeMeanTracts( std::vector< nodeID_t > inNodes ) const
{
    std::vector< size_t > inVector;
    inVector.reserve( inNodes.size() );
    for( size_t i = 0; i < inNodes.size(); ++i )
    {
        if( inNodes[i].first )
        {
            inVector.push_back( inNodes[i].second );
        }
        else
        {
            continue;
        }
    }
    writeMeanTracts( inVector );
    return;
} // end "writeMeanTracts()" -----------------------------------------------------------------


void treeManager::flipX()
{
    if( m_verbose )
        std::cout << "Flipping tree around X axis ...";

    for( std::vector< WHcoord >::iterator iter( m_tree.m_coordinates.begin() ); iter != m_tree.m_coordinates.end(); ++iter )
        iter->m_x = m_tree.m_datasetSize.m_x - 1 - ( iter->m_x );
    for( std::list< WHcoord >::iterator iter( m_tree.m_discarded.begin() ); iter != m_tree.m_discarded.end(); ++iter )
        iter->m_x = m_tree.m_datasetSize.m_x - 1 - ( iter->m_x );

    m_tree.m_treeName += ( "_flipX" );

    if( m_verbose )
        std::cout << "Done" << std::endl;

    return;
} // end "flipXwrite()" -----------------------------------------------------------------


void treeManager::writeTree() const
{
    if( m_outputFolder.empty() )
    {
        std::cerr << "ERROR @ treeManager::writeTree(): Tree is not ready, or outputfolder is not set" << std::endl;
        return;
    }
    m_tree.writeTree( m_outputFolder + "/" + m_tree.m_treeName + ".txt", m_niftiMode );
    if( m_verbose )
    {
        std::cout << "Written standard tree file in: " << m_outputFolder << "/" << m_tree.m_treeName << ".txt" << std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Standard tree file in:\t" << m_outputFolder << "/" << m_tree.m_treeName << ".txt" << std::endl;
    }

    if( m_debug )
    {
        m_tree.writeTreeDebug( m_outputFolder + "/" + m_tree.m_treeName + "_debug.txt" );
        //    m_tree.writeTreeOldWalnut(m_outputFolder+"/"+m_tree.m_treeName+"_4ow.txt");
        if( m_verbose )
        {
            std::cout << "Written debug tree file in: " << m_outputFolder << "/" << m_tree.m_treeName << "_debug.txt" << std::endl;
            //        std::cout<<"Written walnut tree file in: "<< m_outputFolder <<"/"<<m_tree.m_treeName<<"_4ow.txt"<<std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Debug tree file in:\t" << m_outputFolder << "/" << m_tree.m_treeName << "_debug.txt" << std::endl;
            //        (*m_logfile)<<"Walnut tree file in:\t"<< m_outputFolder <<"/"<<m_tree.m_treeName<<"_4ow.txt"<<std::endl;
        }
    }
    return;
} // end treeManager::writeTree() -------------------------------------------------------------------------------------


void treeManager::writeNodeTracts( std::vector< size_t >* const nodeVector, listedCache< compactTract >* const cache, size_t* const tractProg,
                time_t* const lastTime, const time_t &startTime ) const
{
    std::vector< size_t >& nodeVectorRef( *nodeVector );
    time_t& lastTimeRef( *lastTime );

    fileManagerFactory leavesFileMF( m_leavesTractFolder );
    fileManager& leavesFM( leavesFileMF.getFM() );
    leavesFM.readAsLog();
    leavesFM.readAsUnThres();

    volatile size_t threadcounter( 0 ); // keeps track of unfinished generated threads

    // write down the tractograms for each node
    for( size_t i = 0; i < nodeVectorRef.size(); ++i )
    {
        const WHnode& currentNode( m_tree.getNode( nodeVectorRef[i] ) );

        std::vector< nodeID_t > kids( currentNode.getChildren() );

        // compute the mean tract

        // get tractogram of first element
        compactTract meanTract;
        if( kids[0].first )
        {
            compactTract* tractPointer;
            tractPointer = ( cache->getNoUpdate( kids[0].second ) );
            if( tractPointer == 0 )
            {
                throw std::runtime_error( "ERROR @ hTree::writeNodeTracts(): tractogram not found in memory" );
            }
            else
            {
                meanTract.steal( tractPointer );
#pragma omp critical( cache )
                cache->erase( kids[0].second );
            }
        }
        else
        {
            leavesFM.readLeafTract( kids[0].second, m_tree.m_trackids, m_tree.m_coordinates, &meanTract );
            meanTract.unLog( m_tree.m_logFactor );
        }
        size_t meanSize( ( m_tree.getNode( kids[0] ) ).getSize() );

        // iteratively add tractograms of following elements
        for( size_t j = 1; j < kids.size(); ++j )
        {
            compactTract addedTract;
            if( kids[j].first )
            {
                compactTract* tractPointer;
                tractPointer = ( cache->getNoUpdate( kids[j].second ) );
                if( tractPointer == 0 )
                {
                    throw std::runtime_error( "ERROR @ hTree::writeNodeTracts(): tractogram not found in memory" );
                }
                else
                {
                    addedTract.steal( tractPointer );
#pragma omp critical( cache )
                    cache->erase( kids[j].second );
                }
            }
            else
            {
                leavesFM.readLeafTract( kids[j].second, m_tree.m_trackids, m_tree.m_coordinates, &addedTract );
                addedTract.unLog( m_tree.m_logFactor );
            }
            size_t addedSize( ( m_tree.getNode( kids[j] ) ).getSize() );
            compactTract tempTract( addedTract, meanTract, addedSize, meanSize );
            meanTract.steal( &tempTract );
            meanSize += addedSize;
        }

#pragma omp critical( cache )
        cache->insert( m_tree.getNode( nodeVectorRef[i] ).getID(), meanTract );

#pragma omp atomic
        ++( *tractProg );

        // log and write the meantract, boost::bind creates a self contained function object, so ts ok if the function takes everything by reference
        //volatile size_t *tcpointer();
        logWriteTract( currentNode.getID(), &meanTract, m_meanTractFolder );

#pragma omp single nowait // only one thread executes output
        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTimeRef > 1 )
            {
                // output progress every 2 seconds
                lastTimeRef = currentTime;
                size_t currentCount( *tractProg );
                float progress = ( currentCount ) * 100. / ( m_tree.m_nodes.size() );
                int expected_remain( difftime( currentTime, startTime ) * ( ( 100. - progress ) / progress ) );
                std::cout << "\r" << ( int )progress << " % Completed (" << currentCount << " node tracts)"
                                << ". Expected remaining time: ";
                std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                % 3600 ) % 60 ) << "\"  " << std::flush;
            }
        }
    } //end for

    std::cout << "\r100 % Completed (" << *tractProg << " node tracts)" << std::endl;


} // end treeManager::writeNodeTracts() -------------------------------------------------------------------------------------

void treeManager::logWriteTract( const size_t &nodeID,  compactTract* tractPtr, const std::string &meanTractFolder ) const
{
    compactTract& tract( *tractPtr );
    tract.doLog( m_tree.m_logFactor );

    fileManagerFactory meanTractFileMF( meanTractFolder );
    fileManager& meanTractFM( meanTractFileMF.getFM() );
    if ( m_floatFlag )
    {
        meanTractFM.writeInFloat();
    }
    else
    {
        meanTractFM.writeInChar();
    }
    if ( m_zipFlag )
    {
        meanTractFM.storeZipped();
    }
    else
    {
        meanTractFM.storeUnzipped();
    }

    meanTractFM.writeNodeTract( nodeID, tract );

    return;
} // end treeManager::logWriteTract() -------------------------------------------------------------------------------------


