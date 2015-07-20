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

#include <vector>
#include <list>
#include <string>
#include <utility>
#include <map>

#include "WStringUtils.h"

#include "WHtreeProcesser.h"


WHtreeProcesser::WHtreeProcesser( WHtree* const tree ) : m_tree( *tree )
{
}

WHtreeProcesser::~WHtreeProcesser()
{
    //Cleanup
}

// === PUBLIC MEMBER FUNCTIONS ===


std::pair< size_t, size_t > WHtreeProcesser::pruneTree( float condition, size_t safeSize, const HTPROC_MODE pruneType )
{
    if( safeSize == 0 )
    {
        safeSize = m_tree.getNumLeaves(); // any cluster no matter what size may be pruned if he meets the conditions
    }

    if( pruneType == HTPR_SIZERATIO )
    {
        if( condition < 2 )
            throw std::runtime_error( "ERROR @ WHtreeProcesser::pruneTree(): size pruning ratio must be equal or greater than 2" );
    }
    else if( pruneType == HTPR_JOINSIZE )
    {
        condition = std::floor( condition );
        if( condition < safeSize )
            throw std::runtime_error(
                            "ERROR @ WHtreeProcesser::pruneTree(): size pruning lower boundary must be smaller than greater boundary" );
    }
    else if( pruneType == HTPR_JOINLEVEL )
    {
        if( condition <= 0 || condition >= 1 )
        {
            throw std::runtime_error( "ERROR @ WHtreeProcesser::pruneTree(): condition is out of boundaries" );
        }
        if( safeSize >= m_tree.getNumLeaves() )
        {
            throw std::runtime_error(
                     "ERROR @ WHtreeProcesser::pruneTree(): when pruning by distance level a safe size smaller than the roi size must be entered" );
        }
    }

    size_t prunedLeaves( 0 ), prunedNodes( 0 );

    // loop through all leaves and set them to prune if they match discardidng conditions
    for( std::vector< WHnode >::iterator leavesIter( m_tree.m_leaves.begin() ); leavesIter != m_tree.m_leaves.end(); ++leavesIter )
    {
        size_t parentID( leavesIter->getParent().second );
        size_t parentLevel( m_tree.getNode( parentID ).getDistLevel() );

        if( ( pruneType == HTPR_JOINLEVEL ) && ( parentLevel > condition ) )
        {
            leavesIter->setFlag( true );
        }
        else if( ( pruneType == HTPR_SIZERATIO ) || ( pruneType == HTPR_JOINSIZE ) )
        {
            size_t biggerSize( 0 );
            std::vector< nodeID_t > kids( m_tree.getNode( parentID ).getChildren() );
            for( size_t i = 0; i < kids.size(); ++i )
            {
                size_t brotherSize( m_tree.getNode( kids[i] ).getSize() );
                if ( brotherSize > biggerSize )
                {
                    biggerSize = brotherSize;
                }
            }
            if(  biggerSize > condition )
            {
                leavesIter->setFlag( true );
            }
        }
    }

    // loop through all nodes and set them to prune if they match discarding conditions
    for( std::vector< WHnode >::iterator nodesIter( m_tree.m_nodes.begin() ); nodesIter != m_tree.m_nodes.end() - 1; ++nodesIter )
    { // dont check last node
        size_t parentID( nodesIter->getParent().second );
        size_t nodeSize( nodesIter->getSize() );
        size_t parentLevel( m_tree.getNode( parentID ).getDistLevel() );
        bool pruneBranch( false );

        if( nodeSize < safeSize )
        {
            if( ( pruneType == HTPR_JOINLEVEL ) && ( parentLevel > condition ) )
            {
                pruneBranch = true;
            }
            else if( ( pruneType == HTPR_SIZERATIO ) || ( pruneType == HTPR_JOINSIZE ) )
            {
                size_t biggerSize( 0 );
                std::vector< nodeID_t > kids( m_tree.getNode( parentID ).getChildren() );
                for( size_t i = 0; i < kids.size(); ++i )
                {
                    if ( kids[i] == nodesIter->getFullID() )
                    {
                        continue;
                    }
                    size_t brotherSize( m_tree.getNode( kids[i] ).getSize() );
                    if ( brotherSize > biggerSize )
                    {
                        biggerSize = brotherSize;
                    }
                }

                if( ( pruneType == HTPR_SIZERATIO ) && ( biggerSize > ( nodeSize * condition ) ) )
                {
                    pruneBranch = true;
                }

                if( ( pruneType == HTPR_JOINSIZE ) && ( biggerSize >=  condition ) )
                {
                    pruneBranch = true;
                }
            }
        }

        if ( pruneBranch )
        {
            std::list< nodeID_t > worklist;
            worklist.push_back( nodesIter->getFullID() );
            while( !worklist.empty() )
            {
                WHnode* currentNode( m_tree.fetchNode( worklist.front() ) );
                worklist.pop_front();
                // if current node has already been pruned we continue with the next iteration
                if( currentNode->isFlagged() )
                {
                    continue;
                }
                currentNode->setFlag( true );
                std::vector< nodeID_t > currentKids( currentNode->getChildren() );
                worklist.insert( worklist.end(), currentKids.begin(), currentKids.end() );
            }
        }
    }

    // count total pruned leaves
    for( std::vector< WHnode >::iterator leavesIter( m_tree.m_leaves.begin() ); leavesIter != m_tree.m_leaves.end(); ++leavesIter )
    {
        if( leavesIter->isFlagged() )
        {
            ++prunedLeaves;
        }
    }

    // count total pruned nodes
    for( std::vector< WHnode >::iterator nodesIter( m_tree.m_nodes.begin() ); nodesIter != m_tree.m_nodes.end() - 1; ++nodesIter )
    {
        if( nodesIter->isFlagged() )
            ++prunedNodes;
    }

    if( pruneType == HTPR_SIZERATIO )
    {
        m_tree.m_treeName += ( "_prunedR" + string_utils::toString( safeSize ) + ":" + string_utils::toString( condition ) );
    }
    else if( pruneType == HTPR_JOINSIZE )
    {
        m_tree.m_treeName += ( "_prunedS" + string_utils::toString( condition ) + ":" + string_utils::toString( safeSize ) );
    }
    else if( pruneType == HTPR_JOINLEVEL )
    {
        m_tree.m_treeName += ( "_prunedL" + string_utils::toString( safeSize ) + ":" + string_utils::toString( condition ) );
    }

    std::pair< size_t, size_t > pruned( m_tree.cleanup() );

    pruned.second += m_tree.debinarize();

    return pruned;
} // end "pruneTree()" -----------------------------------------------------------------


std::pair< size_t, size_t > WHtreeProcesser::pruneRandom( const size_t numberPruned, unsigned int seed )
{
    if( numberPruned >= m_tree.getNumLeaves() )
    {
        throw std::runtime_error(
                        "ERROR @ WHtreeProcesser::pruneRandom(): too many seeds to be pruned! (more than leaves in the tree)" );
    }

    std::vector< size_t > prunedIDs( m_tree.getNumLeaves(), 0 );
    for( size_t i = 0; i < prunedIDs.size(); ++i )
    {
        prunedIDs[i] = i;
    }
    size_t prunedLeaves( 0 );

    while( prunedLeaves < numberPruned )
    {
        size_t prunedPosition = ( rand_r( &seed ) % prunedIDs.size() );
        if( m_tree.fetchLeaf( prunedIDs[prunedPosition] )->isFlagged() )
        {
            continue;
        }
        m_tree.fetchLeaf( prunedIDs[prunedPosition] )->setFlag( true );
        ++prunedLeaves;
        prunedIDs.erase( prunedIDs.begin() + prunedPosition );
    }
    unsigned int perthousand( numberPruned * 1000. / m_tree.getNumLeaves() );
    float perone( perthousand / 1000. );
    m_tree.m_treeName += ( "_randpruned" + string_utils::toString( perone ) );

    std::pair< size_t, size_t > pruned( m_tree.cleanup() );

    pruned.second += m_tree.debinarize();

    return pruned;
} // end "pruneRandom()" -----------------------------------------------------------------


size_t WHtreeProcesser::collapseTree( const dist_t flatGap, const dist_t distLevelLimit, const bool keepBaseNodes )
{
    if( flatGap >= 1 )
    {
        std::cerr << "ERROR @ WHtreeProcesser::collapseTree(): flattening interval must be smaller than 1" << std::endl;
        return 0;
    }

    if( flatGap == 0 )
    {
        m_tree.forceMonotonicityDown();
    }
    else
    {
        // flatten all levels
        for( int32_t i = m_tree.getRoot().getID(); i >= 0; --i )
        {
            if( m_tree.getNode( i ).getDistLevel() > distLevelLimit )
            {
                collapseNode( i, flatGap );
            }
        }
        m_tree.m_treeName += ( "_flat" + str( boost::format( "%1.3f" ) % flatGap ) );
    }

    size_t collapsed( m_tree.debinarize( keepBaseNodes ) );

    return collapsed;
} // end "collapseTree()" -----------------------------------------------------------------

size_t WHtreeProcesser::collapseTreeLinear( const dist_t coefficient, const bool keepBaseNodes )
{
    if( coefficient <= 0 || coefficient >= 1 )
    {
        std::cerr << "ERROR @ WHtreeProcesser::collapseTreeLinear(): coefficient must be in the range (0,1)" << std::endl;
        return 0;
    }

    else
    {
        // flatten all levels
        for( int32_t i = m_tree.getRoot().getID(); i >= 0; --i )
        {
            collapseNode( i, coefficient, HTPR_C_LINEAR );
        }
        m_tree.m_treeName += ( "_flatL" + str( boost::format( "%1.3f" ) % coefficient ) );
    }

    size_t collapsed( m_tree.debinarize( keepBaseNodes ) );

    return collapsed;
} // end "collapseTreeLinear()" -----------------------------------------------------------------

size_t WHtreeProcesser::collapseTreeSquare( const dist_t coefficient, const bool keepBaseNodes )
{
    if( coefficient <= 0 || coefficient >= 1 )
    {
        std::cerr << "ERROR @ WHtreeProcesser::collapseTreeLinear(): coefficient must be in the range (0,1)" << std::endl;
        return 0;
    }

    else
    {
        // flatten all levels
        for( int32_t i = m_tree.getRoot().getID(); i >= 0; --i )
        {
            collapseNode( i, coefficient, HTPR_C_SQ );
        }
        m_tree.m_treeName += ( "_flatSQ" + str( boost::format( "%1.3f" ) % coefficient ) );
    }

    size_t collapsed( m_tree.debinarize( keepBaseNodes ) );

    return collapsed;
} // end "collapseTreeSquare()" -----------------------------------------------------------------


size_t WHtreeProcesser::collapseBranch( const dist_t flatGap, const dist_t distLevelLimit, size_t root, const bool keepBaseNodes )
{
    if( flatGap >= 1 )
    {
        std::cerr << "ERROR @ WHtreeProcesser::flattenTree(): flattening interval must be smaller than 1" << std::endl;
        return 0;
    }

    std::vector< size_t > branchNodes( m_tree.getBranchNodes( root ) );

    //put higher hierarchy first:
    std::reverse( branchNodes.begin(), branchNodes.end() );
    // flatten levels
    for( std::vector< size_t >::iterator nodesIter( branchNodes.begin() ); nodesIter != branchNodes.end(); ++nodesIter )
    {
        if( m_tree.getNode( *nodesIter ).getDistLevel() > distLevelLimit )
        {
            collapseNode( *nodesIter, flatGap );
        }
    }

    size_t collapsed( m_tree.debinarize( keepBaseNodes ) );

    return collapsed;
} // end "collapseBranch()" -----------------------------------------------------------------


//std::pair< size_t, size_t > WHtreeProcesser::smoothTree( const size_t safeSize, const size_t smoothSize, const dist_t smoothGap )
//{
//    std::pair< size_t, size_t > firstPruned, secondPruned, finalPruned;
//    size_t firstCollapsed( 0 ), secondCollapsed( 0 ), biggestSize( 0 );

//    //1st stage: prune clusters smaller than safeSize that join clusters bigger than smoothSize
//    firstPruned = pruneTree( smoothSize, safeSize, HTPR_JOINSIZE );

//    //2nd stage: find smooth clusters with gaps smaller than smoothGap and flatten them
//    std::vector< nodeID_t > partition;
//    m_tree.partitionSmooth( smoothGap, &partition, true, m_tree.getRoot().getID() );
//    firstCollapsed = flattenSelection( partition, false );

//    //3rd stage: prune clusters smaller than half of smoothSize that join clusters bigger than double smoothSize
//    secondPruned = pruneTree( ( 2 * smoothSize ), ( smoothSize / 2 ), HTPR_JOINSIZE );

//    //4th stage: find a partition of clusters not bigger than smoothSize (unless undivisible) and flatten them
//    partition.clear();
//    biggestSize = m_tree.partitionClassic( smoothSize, &partition, HTP_SIZE, HTC_VALUE, true, m_tree.getRoot().getID() );
//    secondCollapsed = flattenSelection( partition, false );

//    // output value (eliminated leaves and nodes)
//    finalPruned.first = firstPruned.first + firstCollapsed + secondPruned.first + secondCollapsed;
//    finalPruned.second = firstPruned.second + secondPruned.second;

//    // change name
//    m_tree.m_treeName += ( "_smooth" + string_utils::toString( safeSize ) + ":" + string_utils::toString(
//                    smoothSize ) + ":" + str( boost::format( "%1.3f" ) % smoothGap ) );

//    return finalPruned;
//} // end "smoothTree()" -----------------------------------------------------------------


void WHtreeProcesser::coarseTree( const unsigned int coarseRatio )
{
    if( coarseRatio < 2 )
    {
        return;
    }

    std::map< WHcoord, size_t > roimap;
    size_t count( 0 );
    //create a matrix with seed voxels positions
    std::vector< std::vector< std::vector< bool > > > roimatrix;
    roimatrix.resize( m_tree.m_datasetSize.m_x );
    for( size_t i = 0; i < roimatrix.size(); ++i )
    {
        roimatrix[i].resize( m_tree.m_datasetSize.m_y );
        for( size_t j = 0; j < roimatrix[i].size(); ++j )
        {
            roimatrix[i][j].resize( m_tree.m_datasetSize.m_z, false );
        }
    }

    for( std::vector< WHcoord >::const_iterator coordIter = m_tree.m_coordinates.begin(); coordIter != m_tree.m_coordinates.end(); ++coordIter )
    {
        //std::cout<<*coord_iter<<std::endl;
        roimatrix[coordIter->m_x][coordIter->m_y][coordIter->m_z] = true;
        roimap.insert( std::make_pair( *coordIter, count++ ) );
    }

    // convert previously discarded elements to new grid
    std::list< WHcoord > newDiscarded( m_tree.m_discarded );
    for( std::list< WHcoord >::iterator iter( newDiscarded.begin() ); iter != newDiscarded.end(); ++iter )
    {
        iter->m_x = iter->m_x / coarseRatio;
        iter->m_y = iter->m_y / coarseRatio;
        iter->m_z = iter->m_z / coarseRatio;
    }
    newDiscarded.sort();
    newDiscarded.unique();
    WHcoord newDataSetSize( m_tree.m_datasetSize.m_x / coarseRatio, m_tree.m_datasetSize.m_y / coarseRatio,
                    m_tree.m_datasetSize.m_z / coarseRatio );

    //loop through coarser grid
    for( coord_t k = 0; k < m_tree.m_datasetSize.m_z; k += coarseRatio )
    {
        for( coord_t j = 0; j < m_tree.m_datasetSize.m_y; j += coarseRatio )
        {
            for( coord_t i = 0; i < m_tree.m_datasetSize.m_x; i += coarseRatio )
            {
                //loop through finer grid
                std::list< WHcoord > bigGridVoxel;

                for( unsigned int n = 0; n < coarseRatio; ++n )
                {
                    for( unsigned int m = 0; m < coarseRatio; ++m )
                    {
                        for( unsigned int l = 0; l < coarseRatio; ++l )
                        {
                            if( roimatrix[i + l][j + m][k + n] )
                            {
                                WHcoord tempCoord( i + l, j + m, k + n );
                                bigGridVoxel.push_back( tempCoord );
                            }
                        }
                    }
                }
                if( !bigGridVoxel.empty() )
                {
                    WHcoord keptCoord( bigGridVoxel.front() );
                    bigGridVoxel.pop_front();
                    size_t keptID( roimap[keptCoord] );
                    WHcoord newCoord( keptCoord.m_x / coarseRatio, keptCoord.m_y / coarseRatio, keptCoord.m_z / coarseRatio );
                    m_tree.m_coordinates[keptID] = newCoord;

                    for( std::list< WHcoord >::iterator iter( bigGridVoxel.begin() ); iter != bigGridVoxel.end(); ++iter )
                    {
                        size_t decimatedID( roimap[*iter] );
                        m_tree.fetchLeaf( decimatedID )->setFlag( true );
                        WHcoord emptyCoord;
                        m_tree.m_coordinates[decimatedID] = emptyCoord;
                    }
                }
            }
        }
    }

    m_tree.cleanup();
    m_tree.debinarize();
    m_tree.m_discarded = newDiscarded;
    m_tree.m_datasetSize = newDataSetSize;
    m_tree.m_treeName += ( "_coarse" + string_utils::toString( coarseRatio ) );

    return;
} // end "coarseTree()" -----------------------------------------------------------------


 size_t WHtreeProcesser::baseNodes2Leaves()
 {
     if ( !m_tree.testRootBaseNodes() )
     {
         std::cerr<< "ERROR @ WHtreeProcesser::baseNodes2Leaves(): base nodes have both leaves and other nodes as children, tree wont be processed"<<std::endl;
         return m_tree.getNumLeaves();
     }

     std::vector <size_t> bases( m_tree.getRootBaseNodes() );
     for ( size_t i = 0; i < bases.size(); ++i )
     {
         std::vector <size_t> leaves4node( m_tree.getLeaves4node( bases[i] ) );
         // start in j=1 so that we always leave one leaf not pruned
         for ( size_t j = 1; j < leaves4node.size(); ++j )
         {
             WHnode* currentLeaf( m_tree.fetchLeaf( leaves4node[j] ) );
             currentLeaf->setFlag( true );
         }
     }
     m_tree.cleanup();
     m_tree.m_treeName += ( "_bases" );

     return m_tree.getNumLeaves();
 } // end "baseNodes2Leaves()" -----------------------------------------------------------------

 void WHtreeProcesser::forceMonotonicityUp()
 {
     m_tree.forceMonotonicityUp();
     return;
 }

 void WHtreeProcesser::forceMonotonicityDown()
 {
     m_tree.forceMonotonicityDown();
     return;
 }

 void WHtreeProcesser::forceMonotonicity()
 {
     m_tree.forceMonotonicity();
     return;
 }

 size_t WHtreeProcesser::debinarize( const bool keepBaseNodes )
 {
     return m_tree.debinarize( keepBaseNodes );
 }

// === PRIVATE MEMBER FUNCTIONS ===


size_t WHtreeProcesser::flattenBranch( size_t root, const bool keepBaseNodes )
{
    collapseNode( root, 1 );

    size_t collapsed( m_tree.debinarize( keepBaseNodes ) );

    return collapsed;
} // end "flattenBranch()" -----------------------------------------------------------------


size_t WHtreeProcesser::flattenSelection( std::list< size_t > selection, bool keepBaseNodes )
{
    if( keepBaseNodes && !m_tree.testRootBaseNodes()  )
    {
        std::cerr << "WARNING@ flattenSelection: base nodes have mixed nodes and leaves, flattening will be standard "<< std::endl;
        keepBaseNodes = false;
    }
    while( !selection.empty() )
    {
        size_t thisNode(selection.front());
        selection.pop_front();
        std::vector< nodeID_t > kids( m_tree.getNode(thisNode).getChildren() );
        for( size_t i=0; i<kids.size(); ++i )
        {
            if(kids[i].first)
            {
                if ( keepBaseNodes &&  (m_tree.getNode(kids[i]).getHLevel() == 1 ) )
                {
                    continue;
                }
                else
                {
                    m_tree.fetchNode(kids[i].second)->setFlag(true);
                    selection.push_back(kids[i].second);
                }
            }
        }
    }

    std::pair< size_t, size_t > pruned(m_tree.cleanup() );

    return pruned.second;
} // end "flattenSelection()" -----------------------------------------------------------------
size_t WHtreeProcesser::flattenSelection( const std::vector< size_t > &selection, const bool keepBaseNodes )
{
    std::list< size_t > worklist(selection.begin(),selection.end());
    return flattenSelection(worklist, keepBaseNodes);
} // end "flattenSelection()" -----------------------------------------------------------------
size_t WHtreeProcesser::flattenSelection( const std::vector< nodeID_t > &selection, const bool keepBaseNodes )
{
    std::list< size_t > worklist;
    for( std::vector< nodeID_t >::const_iterator iter( selection.begin() ); iter != selection.end(); ++iter )
    {
        if( iter->first )
        {
            worklist.push_back( iter->second );
        }
    }
    return flattenSelection(worklist, keepBaseNodes);
} // end "flattenSelection()" -----------------------------------------------------------------


std::pair< size_t, size_t > WHtreeProcesser::pruneSelection( const std::vector< size_t > &selection )
{
    flagSelection(selection);
    std::pair< size_t, size_t > pruned( m_tree.cleanup() );
    return pruned;
} // end "pruneSelection()" -----------------------------------------------------------------
std::pair< size_t, size_t > WHtreeProcesser::pruneSelection( const std::vector< nodeID_t > &selection )
{
    flagSelection(selection);
    std::pair< size_t, size_t > pruned( m_tree.cleanup() );
    return pruned;
} // end "pruneSelection()" -----------------------------------------------------------------

void WHtreeProcesser::flagSelection( const std::vector< size_t > &selection )
{
    for( std::vector< size_t >::const_iterator iter( selection.begin() ); iter != selection.end(); ++iter )
    {
        std::vector< size_t > pruneLeaves( m_tree.getLeaves4node( *iter ) );
        flagLeaves( pruneLeaves );
    }
    return;
} // end "flagSelection()" -----------------------------------------------------------------
void WHtreeProcesser::flagSelection( const std::vector< nodeID_t > &selection )
{
    for( std::vector< nodeID_t >::const_iterator iter( selection.begin() ); iter != selection.end(); ++iter )
    {
        std::vector< size_t > pruneLeaves( m_tree.getLeaves4node( *iter ) );
        flagLeaves( pruneLeaves );
    }
    return;
} // end "flagSelection()" -----------------------------------------------------------------
void WHtreeProcesser::flagLeaves( const std::vector< size_t > &selection )
{
    for( std::vector< size_t >::const_iterator iter( selection.begin() ); iter != selection.end(); ++iter )
    {
        m_tree.fetchLeaf( *iter )->setFlag( true );
    }
    return;
} // end "flagSelection()" -----------------------------------------------------------------



void WHtreeProcesser::collapseNode( const nodeID_t &thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode )
{
    if( !thisNodeID.first )
    {
        return;
    }
    else
    {
        collapseNode( thisNodeID.second, coefficient, collapseMode );
    }
    return;
} // end collapseNode() -------------------------------------------------------------------------------------

void WHtreeProcesser::collapseNode( const size_t thisNodeID, const dist_t coefficient, HTPROC_COLLAPSE collapseMode )
{
    if( thisNodeID > m_tree.getRoot().getID() )
    {
        throw std::runtime_error( "ERROR @ collapseNode::flattenBranch(): nodeID is out of boundaries" );
    }

    const WHnode& thisNode( m_tree.getNode( thisNodeID ) );
    dist_t upperLevel( thisNode.getDistLevel() );
    std::list< nodeID_t > worklist;
    worklist.push_back( thisNode.getFullID() );
    while( !worklist.empty() )
    {
        WHnode* currentNode( m_tree.fetchNode( worklist.front() ) );
        worklist.pop_front();
        std::vector< nodeID_t > currentKids( currentNode->getChildren() );
        dist_t parentLevel( currentNode->getDistLevel() );
        for( std::vector< nodeID_t >::iterator kidIter( currentKids.begin() ); kidIter != currentKids.end(); ++kidIter )
        {
            if( m_tree.getNode( *kidIter ).isNode() )
            {
                dist_t kidLevel( m_tree.getNode( *kidIter ).getDistLevel() );

                bool doCollapse(false);

                switch( collapseMode )
                {
                case HTPR_C_CONSTANT:
                    doCollapse = ( ( parentLevel - kidLevel ) < coefficient );
                    break;
                case HTPR_C_LINEAR:
                    doCollapse = ( ( parentLevel - kidLevel ) < ( kidLevel * coefficient ) );
                    break;
                case HTPR_C_SQ:
                    doCollapse = ( ( parentLevel - kidLevel ) < ( kidLevel * kidLevel * coefficient ) );
                    break;
                default:
                    break;
                }

                if( doCollapse )
                {
                    worklist.push_back( *kidIter );
                }
            }
        }
        currentNode->setDistLevel( upperLevel );
    }
    return;
} // end collapseNodeLinear() -------------------------------------------------------------------------------------
