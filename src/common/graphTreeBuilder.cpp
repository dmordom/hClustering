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

#include "graphTreeBuilder.h"


graphTreeBuilder::graphTreeBuilder( std::string roiFilename, bool verbose ):
    m_roiLoaded( false ), m_treeReady( false ), m_logfile( 0 ), m_verbose( verbose )
{
    fileManagerFactory fMFtestFormat;
    m_niftiMode = fMFtestFormat.isNifti();
    RoiLoader roiLoader( m_niftiMode );
    m_roiLoaded = roiLoader.readRoi( roiFilename, &m_datasetGrid, &m_datasetSize, &m_numStreamlines, &m_roi, &m_trackids );
}


void graphTreeBuilder::buildGraph( const TG_GRAPHTYPE graphMethod )
{
    if( !m_roiLoaded )
    {
        std::cerr << "ERROR @ treeBuilder::buildCentroid(): voxel roi is not loaded" << std::endl;
        return;
    }
    if( m_inputFolder.empty() || m_outputFolder.empty() )
    {
        std::cerr << "ERROR @ treeBuilder::buildCentroid(): Location of single tracts or output folder has not been specified,"
                  << " please initialize with treeBuilder::setInputFolder() and treeBuilder::setOutputFolder()" << std::endl;
        return;
    }

    // load matrix
    std::vector< std::vector< float > > distMatrix;
    loadDistMatrix( &distMatrix );

    //initialize leaves vector
    //create lookup table (translate position in the list to leaf/node ID
    std::vector< WHnode > leaves, nodes;
    std::vector< nodeID_t > lookup;
    leaves.reserve( m_roi.size() );
    lookup.reserve( m_roi.size() );
    nodes.reserve( m_roi.size() - 1 );
    for( size_t i = 0; i < m_roi.size(); ++i )
    {
        WHnode leaf( std::make_pair( false, i ) );
        leaves.push_back( leaf );
        lookup.push_back( leaf.getFullID() );
    }

    time_t loopStart( time( NULL ) ), lastTime( time( NULL ) );

    std::vector< dist_t > lowestDistVector( distMatrix.size(), 2 );
    std::vector< std::pair< size_t, size_t > > lowestLocationVector( distMatrix.size(), std::make_pair( 0, 0 ) );

    // keep track of lowest distance per row (accelerates process)
#pragma omp parallel for schedule(guided)
    for( size_t i = 0; i < distMatrix.size(); ++i )
    {
        for( int j = 0; j < i; ++j )
        {
            if( distMatrix[i][j] < lowestDistVector[i] )
            {
                lowestDistVector[i] = distMatrix[i][j];
                lowestLocationVector[i] = std::make_pair( j, i ); //greater number is always in the second position
            }
        }
    }

    // repeat until all nodes are added
    while( nodes.size() < ( leaves.size() - 1 ) )
    {
        // find closest pair
        dist_t lowestDist( 999 ); // reset lower dist to greater than 1
        std::pair< size_t, size_t > lowestLocation( 0, 0 );
        for( size_t i = 1; i < lowestDistVector.size(); ++i )
        {
            if( lowestDistVector[i] < lowestDist )
            {
                lowestDist = lowestDistVector[i];
                lowestLocation = lowestLocationVector[i];
            }
        }

        // get children and new node IDs
        nodeID_t node2join1ID( lookup[lowestLocation.first] );
        nodeID_t node2join2ID( lookup[lowestLocation.second] );
        nodeID_t newID( std::make_pair( true, nodes.size() ) );
        WHnode* node2join1( fetchNode( node2join1ID, &leaves, &nodes ) );
        WHnode* node2join2( fetchNode( node2join2ID, &leaves, &nodes ) );
        node2join1->setParent( newID );
        node2join2->setParent( newID );

        // inroduce new node in node vector
        std::vector< nodeID_t > newKids( 1, node2join1ID );
        newKids.push_back( node2join2ID );
        size_t newSize( node2join1->getSize() + node2join2->getSize() );
        size_t newHLevel( std::max( node2join1->getHLevel(), node2join2->getHLevel() ) + 1 );
        WHnode newNode( newID, newKids, newSize, lowestDist, newHLevel );
        nodes.push_back( newNode );

        //update matrix, values corresponding to the 1st joining node are changed to the ones corresponding to the new node
        // values corresponding to the 2nd joining node are changed to discarded value (3)
#pragma omp parallel for
        for( size_t i = 0; i < distMatrix.size(); ++i )
        {
            if( i < lowestLocation.first )
            {
                distMatrix[lowestLocation.first][i] = newGraphDist( distMatrix[lowestLocation.first][i],
                                distMatrix[lowestLocation.second][i], node2join1->getSize(), node2join2->getSize(), graphMethod );
                distMatrix[lowestLocation.second][i] = 3;
            }
            else if( i == lowestLocation.first )
            {
                distMatrix[lowestLocation.second][i] = 3;
            }
            else if( i < lowestLocation.second )
            {
                distMatrix[i][lowestLocation.first] = newGraphDist( distMatrix[i][lowestLocation.first],
                                distMatrix[lowestLocation.second][i], node2join1->getSize(), node2join2->getSize(), graphMethod );
                distMatrix[lowestLocation.second][i] = 3;
            }
            else if( i == lowestLocation.second )
            {
            }
            else
            {
                distMatrix[i][lowestLocation.first] = newGraphDist( distMatrix[i][lowestLocation.first],
                                distMatrix[i][lowestLocation.second], node2join1->getSize(), node2join2->getSize(), graphMethod );
                distMatrix[i][lowestLocation.second] = 3;
            }
        }

        //update lookup table
        lookup[lowestLocation.first] = newID;
        lookup[lowestLocation.second] = std::make_pair( false, 0 );

        //update lowest distances
#pragma omp parallel for schedule( guided )
        for( int row = 1; row < lowestDistVector.size(); ++row )
        {
            if( lowestDistVector[row] != 3 )
            {
                // if row is not discarded
                if( row < lowestLocation.first )
                {
                    // if row is above first joining node theres nothing to be changed
                }
                else if( row == lowestLocation.second )
                {
                    // this row has been eliminated
                    lowestDistVector[row] = 3;
                    lowestLocationVector[row] = std::make_pair( 0, 0 );
                }
                else if( ( row == lowestLocation.first ) || ( lowestLocationVector[row].first == lowestLocation.first )
                                || ( lowestLocationVector[row].first == lowestLocation.second ) )
                {
                    // if the old smallest distance is no longer valid
                    lowestDistVector[row] = 2;
                    lowestLocationVector[row] = std::make_pair( 0, 0 );

                    for( int col = 0; col < row; ++col )
                    {
                        if( distMatrix[row][col] < lowestDistVector[row] )
                        {
                            lowestDistVector[row] = distMatrix[row][col];
                            lowestLocationVector[row] = std::make_pair( col, row ); //greater number is always in the second position
                        }
                    }
                }
                else
                { // if old distance is still valid
                    if( lowestDistVector[row] > distMatrix[row][lowestLocation.first] )
                    { // if new element is the smallest, change
                        lowestDistVector[row] = distMatrix[row][lowestLocation.first];
                        lowestLocationVector[row] = std::make_pair( lowestLocation.first, row );
                    }
                }
            } // end if
        } // end parallel for

        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                lastTime = currentTime;
                float progress = nodes.size() * 100. / ( leaves.size() - 1. );
                size_t elapsedTime( difftime( currentTime, loopStart ) );
                std::stringstream message;
                message << "\r" << static_cast<int>( progress ) << " % of tree built (";
                message << nodes.size() << " nodes). ";
                if( progress > 0 )
                {
                    size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                    message << "Expected remaining time: ";
                    message << expectedRemain / 3600 << "h ";
                    message << ( expectedRemain % 3600 ) / 60 << "' ";
                    message << ( expectedRemain % 3600 ) % 60 << "\". ";
                }
                message << "Elapsed time: ";
                message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                message << ( elapsedTime % 3600 ) % 60 << "\". ";
                std::cout << message.str() <<std::flush;
            }
        } // end verbose



    } // end big loop

    if( m_verbose )
    {
        int timeTaken = difftime( time( NULL ), loopStart );
        std::cout << "\r" << std::flush << "100% of of tree built. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                        % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }

    std::string graphName;
    if( graphMethod == TG_SINGLE )
    {
        graphName = "single";
    }
    else if( graphMethod == TG_COMPLETE )
    {
        graphName = "complete";
    }
    else if( graphMethod == TG_AVERAGE )
    {
        graphName = "average";
    }
    else if( graphMethod == TG_WEIGHTED )
    {
        graphName = "weighted";
    }
    else if( graphMethod == TG_WARD )
    {
        graphName = "ward";
    }
    else
    {
        throw std::runtime_error( "ERROR @ treeBuilder::buildGraph(): graph method not recognized" );
    }

    {
        std::list< WHcoord > discarded;
        WHtree thisTree( graphName, m_datasetGrid, m_datasetSize, m_numStreamlines, 0, leaves, nodes, m_trackids, m_roi, discarded );
        m_tree = thisTree;
        std::vector< WHnode > emptyL, emptyN;
        leaves.swap( emptyL );
        nodes.swap( emptyN );
    }

    if( !m_tree.check() )
    {
        m_tree.writeTreeDebug( m_outputFolder + "/treedebug.txt" );
        throw std::runtime_error( "ERROR @ treeBuilder::buildGraph(): resulting tree is not valid" );
    }

    m_treeReady = true;

    if( m_verbose )
        std::cout << m_tree.getReport() << std::endl;
    if( m_logfile != 0 )
        ( *m_logfile ) << m_tree.getReport() << std::endl;

    m_tree.m_treeName = ( graphName );
    writeTree();

    return;
} // end treeBuilder::buildGraph() -------------------------------------------------------------------------------------


void graphTreeBuilder::writeTree() const
{
    if( ( !m_treeReady ) || m_outputFolder.empty() )
    {
        std::cerr << "ERROR @ graphTreeBuilder::writeTree(): Tree is not ready, or outputfolder is not set" << std::endl;
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
} // end treeBuilder::writeTree() -------------------------------------------------------------------------------------


WHnode* graphTreeBuilder::fetchNode( const nodeID_t thisNode, std::vector< WHnode >* leavesPointer, std::vector< WHnode >* nodesPointer ) const
{
    std::vector< WHnode >& leaves = *leavesPointer;
    std::vector< WHnode >& nodes = *nodesPointer;

    if( thisNode.first )
    {
        return &( nodes[thisNode.second] );
    }
    else
    {
        return &( leaves[thisNode.second] );
    }
} // end treeBuilder::fetchNode() -------------------------------------------------------------------------------------


void graphTreeBuilder::loadDistMatrix( std::vector< std::vector< float > >* const distMatrix ) const
{
    std::vector< std::vector< float > >& distMatrixRef( *distMatrix );
    // load distance block index
    if( m_verbose )
        std::cout << "Reading distance matrix index..." << std::flush;
    distBlock dBlock( m_inputFolder );
    if( !dBlock.indexReady() )
    {
        std::cerr << "ERROR @ treeBuilder::loadDistMatrix(): distance matrix index did not load" << std::endl;
        exit( -1 );
    }
    unsigned int topBlock( dBlock.topBlock() );
    unsigned int numBlocks( dBlock.numBlocks() );
    if( m_verbose )
        std::cout << "OK. Whole matrix is " << topBlock + 1 << "x" << topBlock + 1 << " blocks (real: " << numBlocks << "). "
                        << std::flush;

    // initialize matrix
    float floatVar( 0 );
    float usedMem( ( ( m_roi.size() * m_roi.size() / 2. ) * ( sizeof( floatVar ) * CHAR_BIT / 8. ) ) / ( 1024 * 1024 * 1024 ) );
    if( m_verbose )
        std::cout << "WARNING: Initializing roi distance matrix. Expected memory comsumption at least " << usedMem
                        << " GBytes... " << std::flush;
    distMatrixRef.clear();
    distMatrixRef.resize( m_roi.size() );
    for( size_t i = 0; i < m_roi.size(); ++i )
        distMatrixRef[i].resize( i, 0 );
    if( m_verbose )
        std::cout << "Done" << std::endl;

    // initialize tracking positions
    size_t rowStart( 0 ), rowEnd( 0 ), colStart( 0 ), colEnd( 0 );

    time_t loopStartTime( time( NULL ) );
    bool firstIteration( true );
    size_t doneCount( 0 );

    // load all elements
    while( rowStart < m_roi.size() )
    {
        // load block containing the next leaves in the roi
        // we invert the matrix, from file is the superior triangular matrix, but for calculation we transform it into the
        // inferior triangular matrix, so that we may keep the indices valid, and only take half the space
        dBlock.loadBlock( m_roi[rowStart], m_roi[colStart] );
        std::pair< unsigned int, unsigned int > blockID( dBlock.blockID() );

        if( firstIteration )
        {
            std::cout << "Block size: " << dBlock.size() << "x" << dBlock.size() << std::endl;
            if( m_verbose )
                std::cout << "Loading block: " << blockID.first << "-" << blockID.second << "..." << std::flush;
            firstIteration = false;
        }
        else
        {
            if( m_verbose )
            {
                float progress = doneCount * 100. / ( ( m_roi.size() * ( m_roi.size() - 1 ) ) / 2 );
                std::cout << "\rLoading block: " << blockID.first << "-" << blockID.second << "..." << ( int )progress
                                << " % completed. Expected remaining time: ";
                if( progress > 0 )
                {
                    int expected_remain( difftime( time( NULL ), loopStartTime ) * ( ( 100. - progress ) / progress ) );
                    std::cout << expected_remain / 3600 << "h " << ( expected_remain % 3600 ) / 60 << "' " << ( ( expected_remain
                                    % 3600 ) % 60 ) << "\"  " << std::flush;
                }
            }
        }

        // get range of block and find the position of the first voxel in the roi not contained in the block (for both rows and columns)
        std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > currentRange( dBlock.getBlockRange() );
        std::pair< WHcoord, WHcoord > rangeBlockRow( currentRange.first );
        std::pair< WHcoord, WHcoord > rangeBlockColumn( currentRange.second );

        rowEnd = std::lower_bound( m_roi.begin() + rowStart, m_roi.end(), rangeBlockRow.second ) - m_roi.begin();
        if( m_roi[rowEnd] == rangeBlockRow.second )
            ++rowEnd;

        if( colStart == rowStart )
        {
            colEnd = rowEnd;
        }
        else
        {
            colEnd = std::lower_bound( m_roi.begin() + colStart, m_roi.end(), rangeBlockColumn.second ) - m_roi.begin();
            if( m_roi[colEnd] == rangeBlockColumn.second )
                ++colEnd;
        }

        if( rowStart == colStart )
        {
            // its a diagonal element

#pragma omp parallel for schedule( guided )
            for( size_t i = rowStart; i < rowEnd; ++i )
            {
                for( size_t j = i + 1; j < colEnd; ++j )
                {
                    distMatrixRef[j][i] = ( dBlock.getDistance( m_roi[i], m_roi[j] ) );
                }
            }
            doneCount += ( ( rowEnd - rowStart ) * ( rowEnd - rowStart - 1 ) ) / 2;
        }
        else
        {
            // its a non diagonal block

#pragma omp parallel for schedule ( guided )
            for( size_t i = rowStart; i < rowEnd; ++i )
            {
                for( size_t j = colStart; j < colEnd; ++j )
                {
                    distMatrixRef[j][i] = ( dBlock.getDistance( m_roi[i], m_roi[j] ) );
                }
            }
            doneCount += ( rowEnd - rowStart ) * ( colEnd - colStart );
        }

        if( colEnd == m_roi.size() )
        {
            // if we just finished the last column of the roi, move to the diagonal element of the next row
            rowStart = rowEnd;
            colStart = rowStart;
        }
        else
        {
            // move to the next column
            colStart = colEnd;
        }
    } // end while

    if( m_verbose )
    {
        int timeTaken = difftime( time( NULL ), loopStartTime );
        std::cout << "\r" << std::flush << "100 % of Matrix loaded. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                        % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }
    if( m_logfile != 0 )
        ( *m_logfile ) << "Distance matrix loaded" << std::endl;

    return;
} // end treeBuilder::loadDistMatrix() -------------------------------------------------------------------------------------

dist_t graphTreeBuilder::newGraphDist( const dist_t distance1, const dist_t distance2, const size_t size1, const size_t size2,
                const TG_GRAPHTYPE graphMethod ) const
{
    if( graphMethod == TG_SINGLE ) // Single graphMethodage: D(k,i+j) = min[D(i,k),D(j,k)]
    {
        return std::min( distance1, distance2 );
    }
    else if( graphMethod == TG_COMPLETE ) // Complete graphMethodage: D(k,i+j) = MAX[D(i,k),D(j,k)]
    {
        return std::max( distance1, distance2 );
    }
    else if( graphMethod == TG_AVERAGE ) // Average graphMethodage: D(k,i+j) = [D(i,k)*Size(i),D(j,k)*size(j)]/[size(i)+size(j)]
    {
        return ( ( size1 * distance1 ) + ( size2 * distance2 ) ) / ( size1 + size2 );
    }
    else if( graphMethod == TG_WEIGHTED ) // Weighted graphMethodage: D(k,i+j) = [D(i,k)+D(i,k)]/2
    {
        return ( distance1 + distance2 ) / 2;
    }
    else if ( graphMethod == TG_WARD )
    {
        dist_t avrg( ( ( size1 * distance1 ) + ( size2 * distance2 ) ) / ( size1 + size2 ) );
        dist_t ward( (size1 * size2) * ( avrg - ( distance1 / 2.) - ( distance2 / 2.) ) / ( size1 + size2 )  );
        return ward;
    }
    else
    {
        throw std::runtime_error( "ERROR @ treeBuilder::loadDistMatrix(): graphMethodage option has an invalid value" );
    }
} // end treeBuilder::newGraphDist() -------------------------------------------------------------------------------------
