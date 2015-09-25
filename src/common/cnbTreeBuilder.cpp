//---------------------------------------------------------------------------
//
// Project: hClustering
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
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------


// std library
#include <vector>
#include <list>
#include <string>
#include <map>
#include <utility>
#include <algorithm>

#include "WStringUtils.h"

#include "cnbTreeBuilder.h"
#include "WHtreeProcesser.h"


CnbTreeBuilder::CnbTreeBuilder( const std::string& roiFilename, const float thresholdRatio, const float maxNbDist, bool verbose,  bool noLog )
{
    m_verbose = verbose;
    m_roiLoaded = false;
    m_treeReady = false;
    m_debug = false;
    m_logfile = 0;
    m_numComps = 0;
    m_ncHits = 0;
    m_ncMiss = 0;
    m_lcHits = 0;
    m_lcMiss = 0;

    fileManagerFactory fMFtestFormat;
    m_niftiMode = fMFtestFormat.isNifti();
    RoiLoader roiLoader( m_niftiMode );
    m_roiLoaded = roiLoader.readRoi( roiFilename, &m_datasetGrid, &m_datasetSize, &m_numStreamlines, &m_roi, &m_trackids );
    std::cout << "Roi loaded, " << m_roi.size() << " seed voxels" << std::endl;

    if( maxNbDist <= 0 || maxNbDist > 1 )
    {
        std::cerr << "WARNING @ cnbTreeBuilder::cnbTreeBuilder(): maximum closest neighbor distance provided (" << maxNbDist << ") is out of bounds (0,1], using a value of 1.0" << std::endl;
        m_maxNbDist = 1;
    }
    else
    {
        m_maxNbDist = maxNbDist;
    }

    if( noLog )
    {
        m_logFactor = 0;
    }
    else
    {
        if( m_numStreamlines == 0 )
        {
            std::cerr << "WARNING @ cnbTreeBuilder::cnbTreeBuilder(): provided a number of stramlines per voxel of 0, interpreting it as requesting no logarithmic normalization of tracts (input tracts must also be in natural units)" << std::endl;
            m_logFactor = 0;
        }
        else
        {
            m_logFactor = log10( m_numStreamlines );
        }
    }

    if( thresholdRatio <= 0 || thresholdRatio >= 1 )
    {
        if( thresholdRatio != 0 )
        {
            std::cerr << "WARNING @ cnbTreeBuilder::cnbTreeBuilder(): threshold ratio provided (" << thresholdRatio << ") is out of bounds [0,1), using a value of 0.0 (no thresholding)" << std::endl;
        }
        m_tractThreshold = 0;
    }
    else if ( m_logFactor == 0 )
    {
        m_tractThreshold = thresholdRatio; // if using natural units the normalized tract threshold is the same as the threshold ratio
    }
    else
    {
        m_tractThreshold = log10( m_numStreamlines * thresholdRatio ) / m_logFactor; // if using log units the threshold must be calculated
    }


}

void CnbTreeBuilder::buildCentroid( const unsigned int nbLevel, const float memory, const TC_GROWTYPE growType,
                                    const size_t baseSize, const bool keepDiscarded )
{
    m_numComps = 0;

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

    if( m_verbose )
    {
        std::cout << "Farthest nearest neighbour distance allowed: " << m_maxNbDist << std::endl;
        std::cout << "Tractogram threshold (in log units): " << m_tractThreshold << std::endl;
        std::cout << "Tractogram log factor: " << m_logFactor << std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Farthest nearest neighbour distance allowed: " << m_maxNbDist << std::endl;
        ( *m_logfile ) << "Tractogram threshold (in log units): " << m_tractThreshold << std::endl;
        ( *m_logfile ) << "Tractogram log factor: " << m_logFactor << std::endl;
    }

    // file io classes
    fileManagerFactory fileSingleMF(m_inputFolder);
    fileManager& fileSingle(fileSingleMF.getFM()); // file object to read single tractograms
    fileSingle.readAsUnThres();
    fileSingle.readAsLog();
    fileSingle.storeUnzipped();

    fileManagerFactory fileNatMeanMF(m_tempFolder);
    fileManager& fileNatMean(fileNatMeanMF.getFM());  // to read/write tractograms in float, natural units and unthresholded  (to be merged later)
    fileNatMean.writeInFloat();
    fileNatMean.readAsUnThres();
    fileNatMean.readAsNat();
    fileNatMean.storeUnzipped();

    // vectors for hierarchical and neighborhood information
    std::vector< protoNode > protoLeaves, protoNodes;
    std::vector< WHnode > leaves, nodes;

    // compute cache size
    float tractMb( 0 ), leafTractMb( 0 );
    size_t cacheSize( 0 );
    float leafCacheRatio( 1 );
    {
        compactTract tempTract;
        compactTractChar tempTractChar;
        fileSingle.readLeafTract( 0, m_trackids, m_roi, &tempTract );
        tractMb = tempTract.mBytes();
        fileSingle.readLeafTract( 0, m_trackids, m_roi, &tempTractChar );
        leafTractMb = tempTractChar.mBytes();
        if( m_verbose )
        {
            std::cout << "Tractogram size is: " << tempTract.size() << " (" << tractMb << " MB)" << std::endl;
            std::cout << "Leaf tractogram size is: " << leafTractMb << " MB" << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Tractogram size:\t" << tempTract.size() << " (" << tractMb << " MB)" << std::endl;
            ( *m_logfile ) << "Leaf tractogram size is: " << leafTractMb << " MB" << std::endl;
        }
        cacheSize = ( memory * 1024 / ( tractMb * 2 ) );
        leafCacheRatio = ( tractMb / leafTractMb );
        if( m_verbose )
        {
            std::cout << "Cache size is: " << cacheSize << " tracts. ("<< static_cast<size_t>(cacheSize*leafCacheRatio) <<" leaf tracts)" << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Cache size:\t" << cacheSize << " tracts. ("<< static_cast<size_t>(cacheSize*leafCacheRatio) <<" leaf tracts)" << std::endl;
        }
    }

    // precompute seed voxel norms
    computeNorms();

    // initialize neighborhood info for all seed voxels
    std::list< WHcoord > discarded = initialize( nbLevel, cacheSize*leafCacheRatio, &protoLeaves );
    std::list< size_t > baseNodes;


    { // ------- Tree build up ----------
        std::cout << "Starting tree build-up" << std::endl;
        // supernode, keeps track of the current top nodes in the hierarchy
        std::multimap< dist_t, nodeID_t > priorityNodes;
        std::set< size_t > currentNodes;

        size_t activeSize( 1 ), prioritySize( 1 );
        bool growingStage( true );
        bool clust10k( true );
        bool lastIsNewline ( true );

        if( growType == TC_GROWOFF || baseSize <= 1 )
        {
            growingStage = false;
            activeSize = protoLeaves.size();
            prioritySize = protoLeaves.size();
        }

        leaves.reserve( protoLeaves.size() );
        nodes.reserve( protoLeaves.size() );

        m_nodeNorms.clear();
        m_nodeNorms.reserve( protoLeaves.size() );
        size_t doneLeavesCounter( 0 );
        WHnode rootNode( std::make_pair( false, 0 ) ); // node where to keep all isolated clusters
        rootNode.setSize( 0 );

        std::vector< std::multimap< dist_t, nodeID_t >::iterator > priorityLeafIndex( protoLeaves.size(), priorityNodes.end() );
        std::vector< std::multimap< dist_t, nodeID_t >::iterator > priorityNodeIndex( protoLeaves.size(), priorityNodes.end() );
        for( size_t i = 0; i < protoLeaves.size(); ++i )
        {
            priorityLeafIndex[i] = priorityNodes.insert( std::make_pair( protoLeaves[i].nearDist(), std::make_pair( false, i ) ) );
            WHnode newLeaf( std::make_pair( false, i ) );
            leaves.push_back( newLeaf );
        }

        listedCache< compactTractChar > leavesCache( protoLeaves.size(), cacheSize*leafCacheRatio );
        listedCache< compactTract > nodesCache( protoLeaves.size(), cacheSize );

        time_t lastTime( time( NULL ) ), loopStart( time( NULL ) ); // time object
        size_t maxNbs( 0 ); // to keep track of maximum number of neighbours in an iteration during the program
        std::stringstream eventStream;
        volatile size_t threadCount( 0 ); // to keep track of running children threads

        m_ncHits = 0;
        m_ncMiss = 0;
        m_lcHits = 0;
        m_lcMiss = 0;


#if DEBUG
        if( verbose )
        {
            std::cout<< "P Size: "<< prioritySize << std::endl;
            std::cout<< "A Size: "<< activeSize << std::endl;
            std::cout<< "Pnumber: "<< priorityNodes.size() << std::endl;
            std::cout<< "Cnumber: "<< currentNodes.size() << std::endl;
        }
#endif


        while( !priorityNodes.empty() || currentNodes.size() > 1 )
        {

            while( !priorityNodes.empty() )
            {

                // get nodes to join
                WHnode* node2join1( fetchNode( priorityNodes.begin()->second, &leaves, &nodes ) );
                protoNode* protoNode2join1( fetchProtoNode( node2join1->getFullID(), &protoLeaves, &protoNodes ) );
                WHnode* node2join2( fetchNode( protoNode2join1->nearNb(), &leaves, &nodes ) );
                protoNode* protoNode2join2( fetchProtoNode( node2join2->getFullID(), &protoLeaves, &protoNodes ) );
                dist_t newDist( priorityNodes.begin()->first );
                size_t newID( nodes.size() );
                size_t newSize( node2join1->getSize() + node2join2->getSize() );
                size_t newHLevel( std::max( node2join1->getHLevel(), node2join2->getHLevel() ) + 1 );

                //std::cout<<std::endl<<node2join1->printAllData()<<std::endl<<*protoNode2join1<<std::endl;
                //std::cout<<node2join2->printAllData()<<std::endl<<*protoNode2join2<<std::endl;

                // if no priority node has an active neighbour, go to the next phase
                if( newDist == noNbDist )
                {
                    break;
                }

#if DEBUG
                // test if information is consistent
                {
                    bool thereIsError( false );
                    // data from priorityNodes vector and protonode1 dont match, there has been an error
                    if( ( newDist != protoNode2join1->nearDist() ) || ( protoNode2join1->nearNb() != node2join2->getFullID() ) || ( node2join1->getFullID() == node2join2->getFullID() ) )
                    {
                        thereIsError = true;
                    }
                    // if data doest match with protonode 2:
                    else if( ( newDist != protoNode2join2->nearDist() ) || ( protoNode2join2->nearNb() != node2join1->getFullID() )  )
                    {
                        // if not anymore on growing stage, or priority and active sizes are equal, it means an error
                        if ( !growingStage || prioritySize == activeSize )
                        {
                            thereIsError = true;
                        }
                        // otherwise its only an error if it is not on the neighborhood
                        else if( protoNode2join2->m_nbNodes.find( node2join1->getFullID() ) == protoNode2join2->m_nbNodes.end() )
                        {
                            thereIsError = true;
                        }
                    }

                    if (thereIsError)
                    {
                        std::cerr << "NewDist: " << newDist << std::endl;
                        std::cerr << "Priority nodes: " << priorityNodes.size() << std::endl;
                        std::cerr << "Current nodes: " << currentNodes.size() << std::endl;
                        std::cerr << "Done nodes size: " << nodes.size() << std::endl;
                        std::cerr << "protoNode2join1: " << *protoNode2join1 << std::endl;
                        std::cerr << "Node2join1: " << node2join1->printAllData() << std::endl;
                        std::cerr << "protoNode2join2: " << *protoNode2join2 << std::endl;
                        std::cerr << "Node2join2: " << node2join2->printAllData() << std::endl;
                        m_tree.writeTreeDebug( m_outputFolder + "/treeErrorDebug.txt" );
                        throw std::runtime_error( "ERROR @ treeBuilder::buildCentroid(): closest distance in prioritynodes does not agree with protoNode inner data" );
                    }
                }
#endif

                //get stored unloged tractograms
                compactTract tract1, tract2;

                if( node2join1->isNode() )
                {

                    fileNatMean.readNodeTract( node2join1->getID(), &tract1 );
                    // delete meantract, boost::bind creates a self contained function object,
                    // so ts ok if the function takes everything by reference

                    boost::thread delThread( boost::bind( &CnbTreeBuilder::deleteTract, this, node2join1->getID(),
                                                         &fileNatMean, &threadCount ) );
                }
                else
                {
                    fileSingle.readLeafTract( node2join1->getID(), m_trackids, m_roi, &tract1 );
                    tract1.unLog( m_logFactor );
                }

                if( node2join2->isNode() )
                {
                    fileNatMean.readNodeTract( node2join2->getID(), &tract2 );
                    // delete meantract, boost::bind creates a self contained function object,
                    // so ts ok if the function takes everything by reference

                    boost::thread delThread( boost::bind( &CnbTreeBuilder::deleteTract, this, node2join2->getID(),
                                                         &fileNatMean, &threadCount ) );
                }
                else
                {
                    fileSingle.readLeafTract( node2join2->getID(), m_trackids, m_roi, &tract2 );
                    tract2.unLog( m_logFactor );
                }

                if( node2join1->isNode() )
                {
#pragma omp critical( nodesCache )
                    nodesCache.erase( node2join1->getID() );
                }
                else
                {
#pragma omp atomic
                    ++doneLeavesCounter;
#pragma omp critical( leavesCache )
                    leavesCache.erase( node2join1->getID() );
                }

                if( node2join2->isNode() )
                {
#pragma omp critical( nodesCache )
                    nodesCache.erase( node2join2->getID() );
                }
                else
                {
#pragma omp atomic
                    ++doneLeavesCounter;
#pragma omp critical( leavesCache )
                    leavesCache.erase( node2join2->getID() );
                }


                // initialize data members of new node object
                std::pair< nodeID_t, dist_t > newNearNb( noNbID, noNbDist );
                std::map< nodeID_t, dist_t > newNbNodes;
                bool newIsActive( newSize <= activeSize );


                // eliminate children entries from current and priority vectors
                priorityNodes.erase( priorityNodes.begin() );
                if( node2join2->isNode() )
                {
                    if (node2join2->getSize() > prioritySize )
                    {
                        currentNodes.erase( node2join2->getID() );
                    }
                    else
                    {
                        priorityNodes.erase( priorityNodeIndex[node2join2->getID()] );
                    }
                }
                else
                {
                    priorityNodes.erase( priorityLeafIndex[node2join2->getID()] );
                }

                // update parent of joining nodes
                node2join1->setParent( std::make_pair( true, newID ) );
                node2join2->setParent( std::make_pair( true, newID ) );

                // start new protonode (merging nbhood tables)
                newNbNodes.insert( protoNode2join1->m_nbNodes.begin(), protoNode2join1->m_nbNodes.end() );
                newNbNodes.insert( protoNode2join2->m_nbNodes.begin(), protoNode2join2->m_nbNodes.end() );
                newNbNodes.erase( node2join1->getFullID() );
                newNbNodes.erase( node2join2->getFullID() );
                protoNode2join1->clearNbhood();
                protoNode2join1->inactivate();
                protoNode2join2->clearNbhood();
                protoNode2join2->inactivate();
                maxNbs = std::max( maxNbs, newNbNodes.size() );


                // get mean tractogram, write it to file in natural units, log it, threshold it, compute norm and put in cache
                compactTract* newTract;
                compactTract tempTract( tract1, tract2, node2join1->getSize(), node2join2->getSize() );

                fileNatMean.writeNodeTract( newID, tempTract );


                tempTract.doLog( m_logFactor );
                tempTract.threshold( m_tractThreshold );
                m_nodeNorms.push_back( tempTract.computeNorm() );

#pragma omp critical( nodesCache )
                newTract = nodesCache.insert( newID, compactTract() );
                newTract->steal( &tempTract );



                // load tracts from all neighbors
                std::vector< void* >nbTractVect;
                nbTractVect.reserve( newNbNodes.size() );

                for( std::map< nodeID_t, dist_t >::iterator nbIter( newNbNodes.begin() ); nbIter != newNbNodes.end(); ++nbIter )
                {
                    nbTractVect.push_back( loadTract( nbIter->first, &fileSingle, &fileNatMean, &leavesCache, &nodesCache ) );
                }

                // get distances to all neighbours
#pragma omp parallel for schedule( static )
                for( size_t i = 0; i < newNbNodes.size(); ++i )
                {
                    std::map< nodeID_t, dist_t >::iterator nbIter( newNbNodes.begin() );
                    for( int j = 0; j < i; ++j )
                        ++nbIter;

                    bool nbIsNode( nbIter->first.first );
                    size_t nbId( nbIter->first.second );
                    dist_t newNbDist( 0 );
                    bool isNbActive( false );

                    // if we are in the first go of the homogeneus building we set new distances to 1 and recompute all at the end


                    if( nbIsNode )
                    { //nb is a node
                        if( protoNodes[nbId].isActive() )
                        {
                            isNbActive = true;
                        }

                        // update distance
                        newNbDist=( newTract->tractDistance( * static_cast< compactTract* >( nbTractVect[i] ) ) );
                    }
                    else
                    { //its a leaf
                        isNbActive = true;

                        // update distance
                        newNbDist=( newTract->tractDistance( * static_cast< compactTractChar* >( nbTractVect[i] ) ) );
                    }
#pragma omp atomic
                    m_numComps++;

                    nbIter->second = newNbDist;
#pragma omp critical
                    if( isNbActive && newNbDist < newNearNb.second )
                    {
                        newNearNb.second = newNbDist;
                        newNearNb.first = nbIter->first;
                    }

                    // update neighbourhood in neighbour node object
                    protoNode* newProtoNb( fetchProtoNode( nbIter->first, &protoLeaves, &protoNodes ) );

                    bool nbhoodChanged( false );
                    nbhoodChanged = ( newProtoNb->updateActivhood( node2join1->getFullID(), node2join2->getFullID(),
                                                                  std::make_pair( true, newID ),newNbDist, newIsActive, protoNodes  ) );

                    if( nbhoodChanged )
                    {
                        //if nearest distance has changed and it is a node in the priority vector update priority vector entry
#pragma omp critical( supernode )
                        {
                            if ( !nbIsNode )
                            {
                                priorityNodes.erase( priorityLeafIndex[nbId] );
                                priorityLeafIndex[nbId] = priorityNodes.insert( std::make_pair( newProtoNb->nearDist(), nbIter->first ) );
                            }
                            else if( nodes[nbId].getSize() <= prioritySize )
                            {
                                priorityNodes.erase( priorityNodeIndex[nbId] );
                                priorityNodeIndex[nbId] = priorityNodes.insert( std::make_pair( newProtoNb->nearDist(), nbIter->first ) );
                            }

                        }
                    }
                } // end parallel for

                // set cache sizes, and clean up if overflowed, (50% of the memory for leaves and 50% for nodes unless all leaves are done)
                if( leavesCache.limit() != 0 )
                {
                    size_t leavesCacheSize(1);
                    if( growingStage )
                    {
                        leavesCacheSize = std::min( leaves.size() - doneLeavesCounter, static_cast<size_t>(leafCacheRatio * cacheSize) );

                    }
                    else
                    {
                        leavesCacheSize = std::min( leaves.size() - doneLeavesCounter, static_cast<size_t>(leafCacheRatio * cacheSize / 2) );

                    }
                    leavesCache.setLimit( leavesCacheSize );
                    if( leavesCacheSize == 0 )
                    {
                        leavesCache.shutdown();
                    }
                    else
                    {
                        leavesCache.cleanup(); // clean leaves cache
                    }
                    nodesCache.setLimit( cacheSize - (leavesCacheSize/leafCacheRatio) + 1 );
                }
                nodesCache.cleanup();


                // insert new node object
                std::vector< nodeID_t > newKids( 1, node2join1->getFullID() );
                newKids.push_back( node2join2->getFullID() );
                WHnode newNode( std::make_pair( true, newID ), newKids, newSize, newDist, newHLevel );
                nodes.push_back( newNode );

                // insert new protoNode object
                protoNode newProtoNode( newNearNb, newNbNodes, newIsActive );
                protoNodes.push_back( newProtoNode );

                // if new node is isolated
                if( ( newNbNodes.empty() ) )
                {
                    if( m_verbose && ( newSize != m_roi.size() ) )
                    {
                        if( !lastIsNewline )
                        {
                            std::cout << std::endl;
                        }
                        std::cout << "Node (1-" << newID << ") with " << newSize
                                  << " leaves has no more neighbours it wont be further considered for clustering."
                                  << std::endl;
                        lastIsNewline = true;
                    }

                    eventStream << "Node (1-" << newID << ") with " << newSize << " leaves is isolated" << std::endl;

                    // update top node
                    rootNode.setID( std::make_pair( true, newID + 1 ) );
                    rootNode.setHLevel( std::max( newHLevel + 1, rootNode.getHLevel() ) );
                    rootNode.setSize( rootNode.getSize() + newSize );
                    std::vector< nodeID_t > topKids( rootNode.getChildren() );
                    topKids.push_back( std::make_pair( true, newID ) );
                    rootNode.setChildren( topKids );

                    // if it is the biggest root node
                    if( newSize > ( m_roi.size() / 2 ) )
                    {
                        if( m_verbose && ( newSize != m_roi.size() ) )
                        {
                            std::cout << "This node contains " << newSize * 100. / m_roi.size() << "% of the total leaves,"
                                      << " it will be kept as the root of the tree, remaining isolated nodes will be eliminated"
                                      << std::endl;
                        }

                        if( m_debug )
                        {
                            // store last tract in output folder
                            compactTract rootTract;
                            fileNatMean.readNodeTract( newID, &rootTract );
                            rootTract.doLog( m_logFactor );

                            fileManagerFactory fileLastMF(m_outputFolder);
                            fileManager& fileLast(fileLastMF.getFM());
                            fileLast.writeInFloat();
                            fileLast.storeUnzipped();
                            fileLast.writeNodeTract( newID, rootTract );
                        }

                        // delete saved tracts
                        boost::thread delThread( boost::bind( &CnbTreeBuilder::deleteTract, this, newID, &fileNatMean, &threadCount ) );
                    }
                    else
                    {
                        if( m_verbose && ( newSize > m_roi.size() / 20 ) )
                        {
                            std::cout << "WARNING: " << newSize * 100 / m_roi.size() << "% of the total leaves are on this isolated node"
                                      << " that cant be further integrated in the tree, the corresponding branch will be eliminated from results"
                                      << std::endl;
                        }

                        // delete saved tracts
                        boost::thread delThread( boost::bind( &CnbTreeBuilder::deleteTract, this, newID, &fileNatMean, &threadCount ) );

                        std::list< nodeID_t > worklist;
                        worklist.push_back( std::make_pair( true, newID ) );
                        while( !worklist.empty() )
                        {
                            nodeID_t currentID( worklist.front() );
                            worklist.pop_front();
                            WHnode* currentNode = fetchNode( currentID, &leaves, &nodes );
                            currentNode->setFlag( true );
                            std::vector< nodeID_t > currentKids( currentNode->getChildren() );
                            worklist.insert( worklist.end(), currentKids.begin(), currentKids.end() );
                        }
                    }
                }
                else
                {
                    //add new node to corresponding vector
                    if( newSize > prioritySize )
                    {
                        currentNodes.insert( newID );
                    }
                    else
                    {
                        priorityNodeIndex[newID] = priorityNodes.insert( std::make_pair( newNearNb.second, std::make_pair( true, newID ) ) );
                    }
                }

                if( m_verbose )
                {
                    time_t currentTime( time( NULL ) );
                    if( currentTime - lastTime > 1 )
                    {
                        lastTime = currentTime;
                        float progress = nodes.size() * 100. / ( leaves.size() - 1. );
                        int progressInt( progress );
                        size_t elapsedTime( difftime( currentTime, loopStart ) );
                        std::stringstream message;
                        message << "\r" << progressInt << " % of tree built (";
                        message << nodes.size() << " nodes built. " << priorityNodes.size() + currentNodes.size() <<" current";
                        if( growingStage )
                        {
                            message<< ". P: "<< prioritySize <<". A: " << activeSize;
                        }
                        message <<"). ";
                        message << "Elapsed: ";
                        message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                        message << ( elapsedTime % 3600 ) % 60 << "\". ";
                        if( progress > 0 )
                        {
                            size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                            message << "Remaining: ";
                            message << expectedRemain / 3600 << "h ";
                            message << ( expectedRemain % 3600 ) / 60 << "' ";
                            message << ( expectedRemain % 3600 ) % 60 << "\". ";
                        }
                        std::cout << message.str() <<std::flush;
                        lastIsNewline = false;
                    }
                } // end verbose

                if( growingStage && clust10k && ( currentNodes.size() + priorityNodes.size() <= 10000 ) )
                {
                    clust10k = false;
                    std::vector< size_t > bases10k;
                    for( std::multimap< dist_t, nodeID_t >::iterator priorityIter( priorityNodes.begin() ); priorityIter != priorityNodes.end(); ++priorityIter )
                    {
                        if(priorityIter->second.first)
                        {
                            bases10k.push_back(priorityIter->second.second);
                        }
                    }
                    bases10k.insert( bases10k.end(), currentNodes.begin(), currentNodes.end() );
                    std::sort(bases10k.begin(),bases10k.end());
                    if( m_debug )
                    {
                        WHtree dummyTree;
                        WHtreeProcesser treeProcesser( &dummyTree );
                        treeProcesser.writeBases( bases10k, m_outputFolder + "/baselist_bin_10k.txt" );
                    }
                }

                if( growingStage && ( growType == TC_GROWNUM ) && ( currentNodes.size() + priorityNodes.size() <= baseSize ) )
                {
                    growingStage = false;
                    activeSize = protoLeaves.size();
                    prioritySize = protoLeaves.size();
                    baseNodes.clear();
                    for( std::multimap< dist_t, nodeID_t >::iterator priorityIter( priorityNodes.begin() ); priorityIter != priorityNodes.end(); ++priorityIter )
                    {
                        if(priorityIter->second.first)
                        {
                            baseNodes.push_back(priorityIter->second.second);
                        }
                    }
                    baseNodes.insert( baseNodes.end(), currentNodes.begin(), currentNodes.end() );
                    break;
                }
            } // end inner big loop (priority size)

            if( growingStage )
            {
                if ( !priorityNodes.empty()  )
                {
                    ++activeSize;
                }
                else if( !currentNodes.empty() )
                {
                    ++prioritySize;
                    if( ( growType == TC_GROWSIZE ) && ( prioritySize >= baseSize ) )
                    {
                        growingStage = false;
                        prioritySize = protoLeaves.size();
                        activeSize = protoLeaves.size();
                        for( std::multimap< dist_t, nodeID_t >::iterator priorityIter( priorityNodes.begin() ); priorityIter != priorityNodes.end(); )
                        {
                            if(priorityIter->second.first)
                            {
                                baseNodes.push_back(priorityIter->second.second);
                            }
                        }
                        baseNodes.insert( baseNodes.end(), currentNodes.begin(), currentNodes.end() );
                    }
                    else
                    {
                        activeSize = prioritySize;
                    }
                }

#if DEBUG
                if( verbose )
                {
                    std::cout<< "P Size: "<< prioritySize << std::endl;
                    std::cout<< "A Size: "<< activeSize << std::endl;
                }
#endif

            }

            if ( growingStage || !currentNodes.empty() )
            {

                // activate or deactivate clusters given new active size
                for( std::set< size_t >::const_iterator currentIter( currentNodes.begin() ); currentIter != currentNodes.end(); ++currentIter )
                {
                    size_t thisSize( nodes[*currentIter].getSize() );
                    if( thisSize <= activeSize )
                    {
                        protoNodes[*currentIter].reactivate();
                    }
                    else
                    {
                        protoNodes[*currentIter].inactivate();
                    }
                }
                //update nearest neighbors for nodes already in the priority list, save changed entries in a temporal list
                std::list< std::pair< dist_t, nodeID_t > > tempPnodes;
                for( std::multimap< dist_t, nodeID_t >::iterator priorityIter( priorityNodes.begin() ); priorityIter != priorityNodes.end(); )
                {
                    bool elementChanged( false );
                    bool isNode( priorityIter->second.first );
                    size_t thisNodeID(priorityIter->second.second );
                    if( isNode )
                    {
                        elementChanged = protoNodes[thisNodeID].updateActive( protoNodes );
                        if( elementChanged )
                        {
                            tempPnodes.push_back( std::make_pair( protoNodes[thisNodeID].nearDist(), priorityIter->second ) );
                        }
                    }
                    else
                    {
                        elementChanged = protoLeaves[thisNodeID].updateActive( protoNodes );
                        if( elementChanged )
                        {
                            tempPnodes.push_back( std::make_pair( protoLeaves[thisNodeID].nearDist(), priorityIter->second ) );
                        }
                    }
                    // if the entry changed delete it
                    if( elementChanged )
                    {
                        priorityNodes.erase( priorityIter++ );
                    }
                    else
                    {
                        ++priorityIter;
                    }
                }
                //insert updated elements on the priority map
                for( std::list< std::pair< dist_t, nodeID_t > >::const_iterator tempIter( tempPnodes.begin() ); tempIter != tempPnodes.end(); ++tempIter )
                {
                    bool isNode( tempIter->second.first );
                    size_t thisNodeID(tempIter->second.second );
                    if( isNode )
                    {
                        priorityNodeIndex[thisNodeID] = priorityNodes.insert( *tempIter );
                    }
                    else
                    {
                        priorityLeafIndex[thisNodeID] = priorityNodes.insert( *tempIter );
                    }
                }
                tempPnodes.clear();
                //update nearest neighbors for nodes in the current list and move into the priority list if necessary
                for( std::set< size_t >::const_iterator currentIter( currentNodes.begin() ); currentIter != currentNodes.end(); )
                {
                    size_t thisNodeID( *currentIter );
                    size_t thisSize( nodes[*currentIter].getSize() );
                    protoNodes[*currentIter].updateActive( protoNodes );
                    if( thisSize <= prioritySize )
                    {
                        priorityNodeIndex[*currentIter] = priorityNodes.insert( std::make_pair( protoNodes[*currentIter].nearDist(), std::make_pair( true, *currentIter) ) );
                        currentNodes.erase( currentIter++ );
                    }
                    else
                    {
                        ++currentIter;
                    }
                }

#if DEBUG
                if( verbose )
                {
                    std::cout<< "Pnumber: "<< priorityNodes.size() << std::endl;
                    std::cout<< "Cnumber: "<< currentNodes.size() << std::endl;
                }
#endif
            }


        } // end upper big loop

        while( threadCount != 0 ) // wait until all threads have finished
        {
            boost::this_thread::sleep( boost::posix_time::microseconds( 100 ) );
        }

        if( !priorityNodes.empty() )
        {
            std::cerr << "WARNING @ treeBuilder::buildCentroid(): after finish, supernode is not empty" << std::endl;
            WHnode* leftNode( fetchNode( priorityNodes.begin()->second, &leaves, &nodes ) );
            std::cerr << "Node info: " << leftNode << std::endl;
            protoNode* leftProtoNode( fetchProtoNode( leftNode->getFullID(), &protoLeaves, &protoNodes ) );
            std::cerr << "Protonode info: " << leftProtoNode << std::endl;
            m_tree.writeTreeDebug( m_outputFolder + "/treeWarningDebug.txt" );
        }

        nodesCache.shutdown();

        // fix last node
        rootNode.setDistLevel( 1 );
        std::vector< nodeID_t > topNodes( rootNode.getChildren() );
        if( topNodes.size() > 1 )
        {
            size_t numValidTopNodes( 0 );
            for( std::vector< nodeID_t >::iterator iter( topNodes.begin() ); iter != topNodes.end(); ++iter )
            {
                WHnode* thisTopNode( fetchNode( *iter, &leaves, &nodes ) );
                thisTopNode->setParent( rootNode.getFullID() );
                if( !thisTopNode->isFlagged() )
                {
                    rootNode.setDistLevel( thisTopNode->getDistLevel() );
                    ++numValidTopNodes;
                }
            }
            if( numValidTopNodes != 1 )
            {
                std::cerr << "WARNING @ treeBuilder::buildCentroid(): more than one valid top node" << std::endl;
                std::cerr << "Root node info: " << rootNode << std::endl;
                m_tree.writeTreeDebug( m_outputFolder + "/treeWarningDebug.txt" );
            }
            nodes.push_back( rootNode );
        }
        else
        {
            fetchNode( topNodes.front(), &leaves, &nodes )->setParent( std::make_pair( false, 0 ) );
        }

        // empty protoNode vectors
        {
            std::vector< protoNode > emptyL, emptyN;
            protoLeaves.swap( emptyL );
            protoNodes.swap( emptyN );
        }

        if( m_verbose )
        {
            int timeTaken = difftime( time( NULL ), loopStart );
            std::cout << "\r" << std::flush << "100% of of tree built. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                                                                                                                     % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
            std::cout << "maximum number of neighbours in one iteration: " << maxNbs << std::endl;
            std::cout << "Node cache. Hits: " << m_ncHits << ". Misses: " << m_ncMiss << std::endl;
            std::cout << "Leaf cache. Hits: " << m_lcHits << ". Misses: " << m_lcMiss << std::endl;
            std::cout << "Total Hits: " << m_lcHits + m_ncHits << ". Total Misses: " << m_lcMiss
                         + m_ncMiss << std::endl;
            std::cout << "Total correlations: " << m_numComps << std::endl;
        }



        if( m_logfile != 0 )
        {
            ( *m_logfile ) << eventStream.str();
            ( *m_logfile ) << "Max #Nbs during construction: " << maxNbs << std::endl;
            ( *m_logfile ) << "Node cache hits: " << m_ncHits << std::endl;
            ( *m_logfile ) << "Node cache misses: " << m_ncMiss << std::endl;
            ( *m_logfile ) << "Leaf cache hits: " << m_lcHits << std::endl;
            ( *m_logfile ) << "Leaf cache misses: " << m_lcMiss << std::endl;
            ( *m_logfile ) << "Total hits: " << m_lcHits + m_ncHits << std::endl;
            ( *m_logfile ) << "Total misses: " << m_lcMiss + m_ncMiss << std::endl;
            ( *m_logfile ) << "Total correlations: " << m_numComps << std::endl;
        }
    } // end tree build up -------------

    time_t procStart( time( NULL ) ); // time object

    std::cout << "Setting up and cleaning tree..." << std::endl;
    {
        std::string treeName( "centroid" + string_utils::toString( nbLevel ) );
        WHtree thisTree( treeName, m_datasetGrid, m_datasetSize, m_numStreamlines, m_logFactor, leaves, nodes, m_trackids, m_roi, discarded );
        m_tree = thisTree;
        std::vector< WHnode > emptyL, emptyN;
        leaves.swap( emptyL );
        nodes.swap( emptyN );
    }

    if( !m_tree.check() )
    {
        m_tree.writeTreeDebug( m_outputFolder + "/treeErrorDebug.txt" );
        throw std::runtime_error( "ERROR @ treeBuilder::buildCentroid(): resulting tree is not valid" );
    }

    if ( baseNodes.empty() )
    {
        std::pair< size_t, size_t > numPruned( m_tree.cleanup() );
        if( m_verbose )
        {
            std::cout << "Done. An additional " << numPruned.first << " leaves and " << numPruned.second
                      << " nodes were discarded" << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Pruned nodes:\t" << numPruned.second << std::endl;
            ( *m_logfile ) << "Total discarded leaves:\t" << m_tree.m_discarded.size() << std::endl;
        }
        if( !keepDiscarded )
        {
            m_tree.m_discarded.clear();
        }

        m_treeReady = true;

        if( m_verbose )
        {
            std::cout << m_tree.getReport() << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << m_tree.getReport() << std::endl;
        }

        if( m_debug )
        {
            m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) + "_bin_nmt" );
            writeTree();
        }

        m_tree.forceMonotonicity();

        if( m_verbose )
        {
            std::cout << "Monotonicity forced, " << m_tree.getReport( false ) << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Monotonicity forced, " << m_tree.getReport( false ) << std::endl;
        }

        if( m_debug )
        {
            m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) + "_bin" );
            writeTree();
        }

        m_tree.debinarize( false );

        if( m_verbose )
        {
            std::cout << "Debinarized, " << m_tree.getReport( false ) << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Debinarized, " << m_tree.getReport( false ) << std::endl;
        }


        m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) );
        writeTree();
    }
    else
    {
        m_treeReady = true;

        baseNodes.sort();
        std::vector< size_t > baseVector( baseNodes.begin(), baseNodes.end() );


        if( m_debug )
        {
            WHtreeProcesser treeProcesser(&m_tree);
            treeProcesser.writeBases( baseVector, m_outputFolder + "/baselist_bin.txt" );
            if( m_verbose )
            {
                std::cout << "Binary tree base list written in: "<< m_outputFolder << "/baselist_bin.txt" << std::endl;
            }
            if( m_logfile != 0 )
            {
                ( *m_logfile ) << "Binary tree base list written in: "<< m_outputFolder << "/baselist_bin.txt" << std::endl;
            }
            if( m_verbose )
            {
                std::cout << m_tree.getReport() << std::endl;
            }
            if( m_logfile != 0 )
            {
                ( *m_logfile ) << m_tree.getReport() << std::endl;
            }
            m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) + "_bin_nmt" );
            writeTree();

            WHtree treeUp(m_tree);
            WHtree treeDown(m_tree);

            treeUp.forceMonotonicityUp();
            WHtreeProcesser processerUp( &treeUp );
            processerUp.flattenSelection(baseNodes, false);
            treeUp.debinarize( true );
            treeUp.m_treeName = ( "c" + string_utils::toString( nbLevel ) +"_Up" );
            treeUp.writeTree( m_outputFolder + "/" + treeUp.m_treeName + ".txt", m_niftiMode );

            treeDown.forceMonotonicityDown();
            WHtreeProcesser processerDown( &treeDown );
            processerDown.flattenSelection(baseNodes, false);
            treeDown.debinarize( true );
            treeDown.m_treeName = ( "c" + string_utils::toString( nbLevel ) +"_Down" );
            treeDown.writeTree( m_outputFolder + "/" + treeDown.m_treeName + ".txt", m_niftiMode );
        }


        m_tree.forceMonotonicity();

        if( m_verbose )
        {
            std::cout << "Monotonicity forced, " << m_tree.getReport( false ) << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Monotonicity forced, " << m_tree.getReport( false ) << std::endl;
        }

        if( m_debug )
        {
            m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) + "_bin" );
            writeTree();
        }

        WHtreeProcesser processer( &m_tree );
        processer.flattenSelection(baseNodes, false);


        if( m_verbose )
        {
            std::cout << "BaseNodes flattened, and tree pruned" << m_tree.getReport( false ) << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "BaseNodes flattened,  and tree pruned" << m_tree.getReport( false ) << std::endl;
        }

        if( m_debug )
        {
            m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) + "_bin_granlimit" );
            writeTree();
        }

        if( !keepDiscarded )
        {
            m_tree.m_discarded.clear();
        }

        m_tree.debinarize( true );


        if( m_verbose )
        {
            std::cout << "Tree Debinarized, " << m_tree.getReport( false ) << std::endl;
        }
        if( m_logfile != 0 )
        {
            ( *m_logfile ) << "Tree Debinarized, " << m_tree.getReport( false ) << std::endl;
        }

        m_tree.m_treeName = ( "c" + string_utils::toString( nbLevel ) );
        writeTree();

        if( m_tree.testRootBaseNodes() )
        {
            baseVector = m_tree.getRootBaseNodes();
            std::sort( baseVector.begin(), baseVector.end() );
            WHtreeProcesser treeProcesser( &m_tree );
            treeProcesser.writeBases( baseVector, m_outputFolder + "/baselist.txt" );

            if( m_verbose )
            {
                std::cout << "Final base list written in: "<< m_outputFolder << "/baselist.txt" << std::endl;
            }
            if( m_logfile != 0 )
            {
                ( *m_logfile ) << "Final base list written in: "<< m_outputFolder << "/baselist.txt" << std::endl;
            }
        }
        else
        {
            if( m_verbose )
            {
                std::cout << "Final tree is not a pure basenode tree" << std::endl;
            }
            if( m_logfile != 0 )
            {
                ( *m_logfile ) << "Final tree is not a pure basenode tree" << std::endl;
            }
        }
    }

    int timeTaken = difftime( time( NULL ), procStart );
    if( m_verbose )
    {
        std::cout << "Tree processed. time taken: " << timeTaken / 3600 << "h "<<std::flush;
        std::cout << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Tree processed. time taken: " << timeTaken / 3600 << "h "<<std::flush;
        ( *m_logfile ) << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }

    return;
} // end cnbTreeBuilder::buildcentroid() -------------------------------------------------------------------------------------


void CnbTreeBuilder::writeTree() const
{
    if( ( !m_treeReady ) || m_outputFolder.empty() )
    {
        std::cerr << "ERROR @ CnbTreeBuilder::writeTree(): Tree is not ready, or outputfolder is not set" << std::endl;
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


protoNode* CnbTreeBuilder::fetchProtoNode( const nodeID_t& thisNode, std::vector< protoNode >* protoLeavesPointer, std::vector< protoNode >* protoNodesPointer ) const
{
    std::vector< protoNode >& protoLeaves = *protoLeavesPointer;
    std::vector< protoNode >& protoNodes = *protoNodesPointer;

    if( thisNode.first )
    {
        return &( protoNodes[thisNode.second] );
    }
    else
    {
        return &( protoLeaves[thisNode.second] );
    }
}

WHnode* CnbTreeBuilder::fetchNode( const nodeID_t& thisNode, std::vector< WHnode >* leavesPointer, std::vector< WHnode >* nodesPointer ) const
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
}

void CnbTreeBuilder::computeNorms()
{
    // loop  through all the seed voxels and compute tractogram norms
    std::cout << "Precomputing tractogram norms" << std::endl;
    time_t loopStart( time( NULL ) ), lastTime( time( NULL ) );
    m_leafNorms.assign( m_roi.size(), 0 );
    size_t progCount( 0 );


    fileManagerFactory fileSingleMF(m_inputFolder);
    fileManager& fileSingle(fileSingleMF.getFM());
    fileSingle.readAsUnThres();
    fileSingle.readAsLog();

#pragma omp parallel for schedule( static )
    for( size_t i = 0; i < m_roi.size(); ++i )
    {
        compactTractChar thisTract;

        fileSingle.readLeafTract( i, m_trackids, m_roi, &thisTract );
        thisTract.threshold( m_tractThreshold );
        m_leafNorms[i] = thisTract.computeNorm();

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
                float progress( currentCount * 100. / m_roi.size() );
                size_t intProgress( progress );
                size_t elapsedTime( difftime( currentTime, loopStart ) );
                std::stringstream message;
                message << "\r" << intProgress << " % of norms computed (" << currentCount << " tracts). ";
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
    } // end  for


    int timeTaken = difftime( time( NULL ), loopStart );
    if( m_verbose )
    {
        std::cout << "\r" << std::flush << "100 % of norms computed. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                        % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }
    if( m_logfile != 0 )
    {
        ( *m_logfile ) << "Norms computed. Time taken: " << timeTaken / 3600 << "h ";
        ( *m_logfile ) << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }

    if( m_verbose )
    {
        for( size_t i = 0; i < m_leafNorms.size(); ++i )
        {
            if( m_leafNorms[i] == 0)
            {
                std::cerr << "track " << i << " has norm 0" <<std::endl;
            }
        }
    }
} // end treeBuilder::computeNorms() -------------------------------------------------------------------------------------

std::list< WHcoord > CnbTreeBuilder::initialize( const unsigned int nbLevel, const size_t cacheSize, std::vector< protoNode >* protoLeavesPointer )
{
    std::cout << "Initializing seed neighbor dissimilarity information" << std::endl;
    std::vector< protoNode >& protoLeaves = *protoLeavesPointer;

    // Translate neighborhood level
    unsigned int nbLevel1( 0 ), nbLevel2( 0 ); // nbhood level indicators
    switch( nbLevel )
    {
        case 6:
        case 18:
        case 26:
        case 32:
            nbLevel1 = nbLevel;
            nbLevel2 = 0;
            break;
        case 92:
            nbLevel1 = 18;
            nbLevel2 = 18;
            break;
        case 124:
            nbLevel1 = 26;
            nbLevel2 = 26;
            break;
        default:
            throw std::runtime_error( "ERROR @ treeBuilder::buildCentroid(): invalid neighbourhood level value" );
    }

    //create a matrix with a mask of the seed voxels and a roi map
    std::list< WHcoord > discarded;
    std::map< WHcoord, size_t > roimap;
    std::vector< std::vector< std::vector< bool > > > roimatrix; // mask matrix indicating seed voxel presence
    {
        std::vector< bool > emptyRow( m_datasetSize.m_z, false );
        std::vector< std::vector< bool > > emptySlice( m_datasetSize.m_y, emptyRow );
        roimatrix.assign( m_datasetSize.m_x, emptySlice );
    }
    for( size_t i = 0; i < m_roi.size(); ++i )
    {
        roimatrix[m_roi[i].m_x][m_roi[i].m_y][m_roi[i].m_z] = true;
        roimap[m_roi[i]] = i;
    }

    protoLeaves.clear();
    protoLeaves.reserve( m_roi.size() );

    fileManagerFactory fileSingleMF(m_inputFolder);
    fileManager& fileSingle(fileSingleMF.getFM());
    fileSingle.readAsUnThres();
    fileSingle.readAsLog();
    listedCache< compactTractChar > cache( m_roi.size(), cacheSize );
    volatile size_t threadCount( 0 );

    //initialize proto-leaves
    time_t loopStart( time( NULL ) ), lastTime( time( NULL ) );
    for( size_t roiID = 0; roiID < m_roi.size(); ++roiID )
    {
        bool discard( true ); // discard voxel flag
        compactTractChar* thisTract( cache.get( roiID ) );

        if( thisTract == 0 )
        { // if tract was not in cche, get it from file
            compactTractChar tempTract;
            fileSingle.readLeafTract( roiID, m_trackids, m_roi, &tempTract );
            tempTract.threshold( m_tractThreshold );
            tempTract.setNorm( m_leafNorms[roiID] );
            thisTract = cache.insert( roiID, compactTractChar() );
            thisTract->steal( &tempTract );
        }
        // get coordinates of neighbouring voxels
        std::vector< WHcoord > nbCoords( m_roi[roiID].getPhysNbs( m_datasetSize, nbLevel1 ) );
        // dicard coordinates that are not part of the roi
        {
            std::vector< WHcoord >::iterator cleanIter( nbCoords.begin() );
            while( cleanIter != nbCoords.end() )
            {
                if( !roimatrix[cleanIter->m_x][cleanIter->m_y][cleanIter->m_z] )
                {
                    cleanIter = nbCoords.erase( cleanIter );
                }
                else
                {
                    ++cleanIter;
                }
            }
        }


        if( ( nbLevel2 != 0 ) )
        { // theres a second nbhood search phase
            std::list< WHcoord > nbAllCoordslist( nbCoords.begin(), nbCoords.end() );


            // loop through level 1 neighbours
            for( std::vector< WHcoord >::const_iterator level1NbsIter = nbCoords.begin(); level1NbsIter != nbCoords.end(); ++level1NbsIter )
            {
                // get level 2 neighbors
                std::vector< WHcoord > nbCoords2part( level1NbsIter->getPhysNbs( m_datasetSize, nbLevel2 ) );
                nbAllCoordslist.insert(nbAllCoordslist.end(),nbCoords2part.begin(),nbCoords2part.end());
            }
            nbAllCoordslist.sort();
            nbAllCoordslist.unique();

            // discard neighbors that are not seeds or are already level 1 ,
            // or that have already been scanned in a previous iteration
            std::list< WHcoord >::iterator cleanIter( nbAllCoordslist.begin() );
            while( cleanIter != nbAllCoordslist.end() )
            {
                if( !roimatrix[cleanIter->m_x][cleanIter->m_y][cleanIter->m_z]  ) // it is not a seed voxel
                {
                    cleanIter = nbAllCoordslist.erase( cleanIter );
                }
                else if( *cleanIter == m_roi[roiID] ) // nb seed is the current seed
                {
                    cleanIter = nbAllCoordslist.erase( cleanIter );
                }
                else
                {
                    ++cleanIter;
                }
            }

            nbCoords.assign(nbAllCoordslist.begin(),nbAllCoordslist.end());
        }

        //convert vector of neighbor coordinates to vector of neighbor ids
        std::vector< size_t > nbIDs;
        nbIDs.reserve( nbCoords.size() );
        for( size_t i = 0; i < nbCoords.size(); ++i )
        {
            nbIDs.push_back( roimap[nbCoords[i]] );
        }
        // get neighborhood information
        std::map< size_t, dist_t > nbLeaves; // variable to keep the current seed voxel neighbours
        discard = scanNbs( roiID, thisTract, protoLeaves, nbIDs, &nbLeaves, &cache );

        if( !discard )
        { // it is a valid seed voxel


            //create a valid proto leaf
            std::pair< nodeID_t, dist_t > nearNb( std::make_pair( std::make_pair( false, 0 ), 999 ) );
            std::map< nodeID_t, dist_t > nbNodes;
            for( std::map< size_t, dist_t >::const_iterator fillIter = nbLeaves.begin(); fillIter != nbLeaves.end(); ++fillIter )
            {
                nbNodes.insert( std::make_pair( std::make_pair( false, fillIter->first ), fillIter->second ) );
                if( fillIter->second < nearNb.second )
                    nearNb = ( std::make_pair( std::make_pair( false, fillIter->first ), fillIter->second ) );
            }
            protoNode thisProtoLeaf( nearNb, nbNodes );
            protoLeaves.push_back( thisProtoLeaf );
        }
        else
        {
            //create discarded proto leaf
            std::pair< nodeID_t, dist_t > nearNbEmpty( std::make_pair( std::make_pair( false, 0 ), 1 ) );
            std::map< nodeID_t, dist_t > nbNodesEmpty;
            protoNode thisProtoLeaf( nearNbEmpty, nbNodesEmpty );
            thisProtoLeaf.discard();
            protoLeaves.push_back( thisProtoLeaf );
        }

        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                lastTime = currentTime;
                float progress( roiID * 100. / m_roi.size() );
                size_t elapsedTime( difftime( currentTime, loopStart ) );
                std::stringstream message;
                message << "\r" << static_cast<int>( progress ) << " % of leaves initialized (" << roiID << "). ";
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

        cache.erase( roiID ); //we wont need this tract anymore
        cache.cleanup(); // if cache is full, delete overflowed tracts
    } // end big loop

    while( threadCount != 0 ) // wait until all threads have finished
        boost::this_thread::sleep( boost::posix_time::microseconds( 100 ) );

    if( m_verbose )
    {
        int timeTaken = difftime( time( NULL ), loopStart );
        std::cout << "\r" << std::flush << "100 % of leaves initialized. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                        % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
        std::cout << "Cleaning up discarded voxels..." << std::endl;
    }


    // cleanup proto leaves: first create alookup table for oldID->newID
    const size_t INVALID_PROTONODE( m_roi.size() );
    size_t validCounter( 0 );
    std::vector< size_t > lookuptable( m_roi.size(), INVALID_PROTONODE );
    for( size_t i = 0; i < protoLeaves.size(); ++i )
    {
        if( !protoLeaves[i].isDiscarded() )
            lookuptable[i] = validCounter++;
    }

    // re-assign names
    for( size_t i = 0; i < protoLeaves.size(); ++i )
    {
        if( !protoLeaves[i].isDiscarded() )
        {
            protoLeaves[i].m_nearNb.first.second = lookuptable[protoLeaves[i].m_nearNb.first.second]; // near nb id
            std::map< nodeID_t, dist_t > newNbs;
            for( std::map< nodeID_t, dist_t >::iterator nbIter( protoLeaves[i].m_nbNodes.begin() ); nbIter
                            != protoLeaves[i].m_nbNodes.end(); ++nbIter )
            {
                size_t nbNewID = ( lookuptable[nbIter->first.second] );
                if( nbNewID != INVALID_PROTONODE )
                    newNbs[std::make_pair( 0, nbNewID )] = nbIter->second;
            }
            protoLeaves[i].m_nbNodes.swap( newNbs );
        }
    }

    // eliminate the discarded coordinates and corresponding proto leaves form the vector
    {
        std::vector< WHcoord >::iterator roiIter( m_roi.begin() );
        std::vector< size_t >::iterator trackidIter( m_trackids.begin() );
        std::vector< protoNode >::iterator protoIter( protoLeaves.begin() );
        std::vector< double >::iterator normIter( m_leafNorms.begin() );

        while( roiIter != m_roi.end() )
        {
            if( protoIter->isDiscarded() )
            {
                discarded.push_back( *roiIter );
                roiIter = m_roi.erase( roiIter );
                protoIter = protoLeaves.erase( protoIter );
                normIter = m_leafNorms.erase( normIter );
                trackidIter = m_trackids.erase( trackidIter );
            }
            else
            {
                ++roiIter;
                ++protoIter;
                ++normIter;
                ++trackidIter;
            }
        }
    }

    float meanNbs( 0 );
    for( size_t i = 0; i < protoLeaves.size(); ++i )
        meanNbs += protoLeaves[i].m_nbNodes.size();
    meanNbs = meanNbs / ( protoLeaves.size() );
    if( m_verbose )
    {
        std::cout << "Done. Mean number of neighbors: " << meanNbs << ". Discarded " << discarded.size() << " seeds" << std::endl;
    }
    if( m_logfile != 0 )
    {
        int timeTaken = difftime( time( NULL ), loopStart );
        ( *m_logfile ) << "Leaves initialized. Time taken: " << timeTaken / 3600 << "h ";
        ( *m_logfile ) << ( timeTaken % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"" << std::endl;
        ( *m_logfile ) << "Mean # of nbs:\t" << meanNbs << std::endl;
        ( *m_logfile ) << "Seeds discarded on Init.:\t" << discarded.size() << std::endl;
    }
    discarded.sort();
    return discarded;
} // end treeBuilder::initialize() -------------------------------------------------------------------------------------


bool CnbTreeBuilder::scanNbs( const size_t currentSeedID,
                              const compactTractChar* const currentTract,
                              const std::vector< protoNode >& protoLeaves,
                              const std::vector< size_t >& nbIDs,
                              std::map< size_t, dist_t >* nbLeavesPointer,
                              listedCache< compactTractChar >* cachePointer )
{
    std::map< size_t, dist_t >& nbLeaves = *nbLeavesPointer;
    listedCache< compactTractChar >& cache = *cachePointer;

    fileManagerFactory fileSingleMF(m_inputFolder);
    fileManager& fileSingle(fileSingleMF.getFM());
    fileSingle.readAsUnThres();
    fileSingle.readAsLog();

    bool discard( true );

    std::vector< compactTractChar* > tractVect;
    tractVect.reserve( nbIDs.size() );
    std::vector< std::pair< size_t, dist_t > > distPairVect;
    distPairVect.reserve( nbIDs.size() );


    // first loop obtains tractograms for non computed distances and recovers the ones already computed
    for( int j = 0; j < nbIDs.size(); ++j )
    {
        size_t thisNbID( nbIDs[j] );

        if( currentSeedID < thisNbID )
        { // neighbour voxel has not yet been processed
            compactTractChar* tractPointer;
            //#pragma omp critical( cache )
            tractPointer = cache.get( thisNbID );

            if( tractPointer == 0 )
            { // if tract was not in cache, get it from file
                compactTractChar nbTract; // get neighbour tract
                fileSingle.readLeafTract( thisNbID, m_trackids, m_roi, &nbTract );
                nbTract.threshold( m_tractThreshold );
                nbTract.setNorm( m_leafNorms[thisNbID] );
                //#pragma omp critical( cache ) // we insert it this way so that a copy of the tractogram doesnt have to be made
                tractPointer = cache.insert( thisNbID, compactTractChar() );
                tractPointer->steal( &nbTract );
            }
            tractVect.push_back( tractPointer );
            distPairVect.push_back( std::make_pair( thisNbID, 2 ) );
        }
        else
        { // neighbour was already processed as seed voxel
            if( thisNbID >= protoLeaves.size() )
            {
                //#pragma omp critical( error )
                {
                    std::cerr<< "Seed: " <<  m_roi[currentSeedID] << ". Nb: " << m_roi[thisNbID] << std::endl;
                    std::cerr<< "nbID: " <<  thisNbID << ". protoLeaves.size(): " << protoLeaves.size() << std::endl;
                    throw std::runtime_error( "ERROR @ treeBuilder::buildCentroid(): neighbor is not in protoLeaves vector" );
                }
            }
            else if( protoLeaves[thisNbID].isDiscarded() )
            { // seed voxel was discarded, skip
                continue;
            }

            std::map< nodeID_t, dist_t >::const_iterator searchIter( protoLeaves[thisNbID].m_nbNodes.find( std::make_pair( false, currentSeedID ) ) );
            if( searchIter == protoLeaves[thisNbID].m_nbNodes.end() )
            {
                //#pragma omp critical( error )
                {
                    std::cerr<< "nb was supposedly already processed but seed is not found in nb data " << std::endl;
                    std::cerr<< "Seed: " <<  m_roi[currentSeedID] << ". Nb: " << m_roi[thisNbID] << std::endl;
                    std::cerr<< "SeedID: " <<  currentSeedID << ". nbID: " <<  thisNbID << ". protoLeaves.size(): " << protoLeaves.size() << std::endl;
                    std::cerr<< "nbInfo: " <<   protoLeaves[thisNbID] << std::endl;

                    throw std::runtime_error( "ERROR @ treeBuilder::scanNbs(): neighborhood data not found" );
                }
            }
            //#pragma omp critical( leaves )
            {
                nbLeaves.insert( std::make_pair( thisNbID, searchIter->second ) ); // save the pair neighbour-distance
                if( searchIter->second <= m_maxNbDist )
                { // neighbour is valid -> seed is valid
                    discard = false;
                }
            }
        }
    }


    // compute distances from tracts in  parallel
    #pragma omp parallel for schedule( static ) // loop through neighbours
    for( int j = 0; j < distPairVect.size(); ++j )
    {
        distPairVect[j].second = currentTract->tractDistance( *tractVect[j] );
    } // end parallel for

    m_numComps += distPairVect.size();

    for( int j = 0; j < distPairVect.size(); ++j )
    {
        nbLeaves.insert( distPairVect[j] ); // save the pair neighbour-distance
        if( distPairVect[j].second <= m_maxNbDist )
        { // neighbour is valid -> seed is valid
            discard = false;
        }
        else if( distPairVect[j].second == 2 )
        {
            throw std::runtime_error( "ERROR @ treeBuilder::scanNbs(): dist value is still 2" );
        }
    }
    return discard;
} // end treeBuilder::scanNbs() -------------------------------------------------------------------------------------


void CnbTreeBuilder::deleteTract( const size_t &nodeID, const fileManager* const fileMngrPointer, volatile size_t* threadcounterPointer ) const
{
#pragma omp atomic
    ++( *threadcounterPointer );

    fileMngrPointer->deleteTractFile( nodeID );

#pragma omp atomic
    --( *threadcounterPointer );
    return;
} // end cnbTreeBuilder::deleteTract() -------------------------------------------------------------------------------------


compactTract* CnbTreeBuilder::loadNodeTract( const size_t nodeID,
                                             const fileManager* const nodeMngrPointer,
                                             listedCache< compactTract >* nodesCachePointer)
{

        compactTract* newNbTract( 0 );

        #pragma omp critical( nodesCache )
        newNbTract = nodesCachePointer->get( nodeID );
        if( newNbTract == 0 )
        { //tractogram is not in cache, it must be in file then
            compactTract nbTractogram;
            nodeMngrPointer->readNodeTract( nodeID, &nbTractogram );
            nbTractogram.doLog( m_logFactor );
            nbTractogram.threshold( m_tractThreshold );
            nbTractogram.setNorm( m_nodeNorms[nodeID] );
            #pragma omp critical( nodesCache )
            newNbTract = nodesCachePointer->insert( nodeID, compactTract() );
            newNbTract->steal( &nbTractogram );
            #pragma omp atomic
            ++m_ncMiss;
        }
        else
        {
            #pragma omp atomic
            ++m_ncHits;
        }
        return newNbTract;

} // end cnbTreeBuilder::loadNodeTract() -------------------------------------------------------------------------------------


compactTractChar* CnbTreeBuilder::loadLeafTract( const size_t leafID,
                                 const fileManager* const leafMngrPointer,
                                 listedCache< compactTractChar >* leavesCachePointer)
{
        compactTractChar* newNbTract( 0 );
        #pragma omp critical( leavesCache )
        newNbTract = leavesCachePointer->get( leafID );
        if( newNbTract == 0 )
        { //tractogram is not in cache, it must be in file then
            compactTractChar nbTractogram;
            leafMngrPointer->readLeafTract( leafID, m_trackids, m_roi, &nbTractogram );
            nbTractogram.threshold( m_tractThreshold );
            nbTractogram.setNorm( m_leafNorms[leafID] );
            #pragma omp critical( leavesCache )
            newNbTract = leavesCachePointer->insert( leafID, compactTractChar() );
            newNbTract->steal( &nbTractogram );
            #pragma omp atomic
            ++m_lcMiss;
        }
        else
        {
            #pragma omp atomic
            ++m_lcHits;
        }
        return newNbTract;
} // end cnbTreeBuilder::loadLeafTract() -------------------------------------------------------------------------------------

void* CnbTreeBuilder::loadTract( const nodeID_t nodeID,
                                 const fileManager* const leafMngrPointer,
                                 const fileManager* const nodeMngrPointer,
                                 listedCache< compactTractChar >* leavesCachePointer,
                                 listedCache< compactTract >* nodesCachePointer )
{
    if( nodeID.first )
    {
        return loadNodeTract( nodeID.second, nodeMngrPointer, nodesCachePointer );
    }
    else
    {
        return loadLeafTract( nodeID.second, leafMngrPointer, leavesCachePointer );
    }

} // end cnbTreeBuilder::loadTract() -------------------------------------------------------------------------------------


