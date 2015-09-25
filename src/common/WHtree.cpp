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

// std library
#include <utility>
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <fstream>

#include "WStringUtils.h"

#include "WHtree.h"



WHtree::WHtree(): m_loadStatus( false ), m_cpcc( 0 )
{
}

WHtree::WHtree( std::string filename ): m_loadStatus( false ), m_cpcc( 0 )
{
    readTree( filename );
}

WHtree::WHtree( std::string treeName, HC_GRID datasetGridInit, WHcoord datasetSizeInit, size_t numStreamlinesInit, float logFactorInit,
                std::vector< WHnode > leavesInit, std::vector< WHnode > nodesInit, std::vector< size_t > trackidInit,
                std::vector< WHcoord > coordInit, std::list< WHcoord > discardInit, float cpccInit )
    : m_loadStatus( false ),
      m_datasetSize( datasetSizeInit ),
      m_datasetGrid( datasetGridInit ),
      m_numStreamlines( numStreamlinesInit ),
      m_logFactor( logFactorInit ),
      m_cpcc( cpccInit ),
      m_treeName( treeName ),
      m_leaves( leavesInit ),
      m_nodes( nodesInit ),
      m_coordinates( coordInit ),
      m_trackids( trackidInit ),
      m_discarded( discardInit )
{
    if( check() )
    {
        m_loadStatus = true;
    }
}

WHtree::WHtree( const WHtree &object )
    : m_loadStatus( object.m_loadStatus ),
      m_datasetSize( object.m_datasetSize ),
      m_datasetGrid( object.m_datasetGrid ),
      m_numStreamlines( object.m_numStreamlines ),
      m_logFactor( object.m_logFactor ),
      m_cpcc( object.m_cpcc ),
      m_treeName( object.m_treeName ),
      m_leaves( object.m_leaves ),
      m_nodes( object.m_nodes ),
      m_coordinates( object.m_coordinates ),
      m_trackids( object.m_trackids ),
      m_discarded( object.m_discarded ),
      m_containedLeaves( object.m_containedLeaves )
{
}

WHtree::~WHtree()
{
    //Cleanup
}



// === PUBLIC MEMBER FUNCTIONS ===



std::string WHtree::getReport( bool longMsg ) const
{
    if( !m_loadStatus)
    {
        return "tree not loaded";
    }
    std::string reportMessage( "Tree has " + string_utils::toString( m_leaves.size() ) + " leaves and " +
                               string_utils::toString( m_nodes.size() ) + " nodes" );
    if( longMsg )
    {
        reportMessage += "\nDataset size is: " + m_datasetSize.getNameString() + " in " + getGridString( m_datasetGrid ) + " format";
        if( m_cpcc != 0 )
        {
            reportMessage += ". CPCC: " + string_utils::toString( m_cpcc );
        }
    }
    return reportMessage;
} // end getReport() -------------------------------------------------------------------------------------


bool WHtree::check() const
{
    if( m_nodes.empty() || m_leaves.size() < 2 )
    {
        std::cerr << "ERROR @ WHtree::check(): only 0-1 leaf / no nodes" << std::endl;
        return false;
    }
    if( m_nodes.size() >= m_leaves.size() )
    {
        std::cerr << "ERROR @ WHtree::check(): same number of nodes as leaves" << std::endl;
        return false;
    }

    std::vector<size_t> sumLeafParents( m_leaves.size(), 0 );
    std::vector<size_t> sumNodeParents( m_nodes.size(), 0 );
    std::vector<size_t> sumNodeKids( m_nodes.size(), 0 );

    // loop through leaves
    for( std::vector<WHnode>::const_iterator leafIter( m_leaves.begin() ); leafIter != m_leaves.end(); ++leafIter )
    {
        nodeID_t parentID( leafIter->getParent() );
        if( !parentID.first )
        {
            std::cerr << "ERROR @ WHtree::check(): leaf has a leaf as parent" << std::endl;
            return false;
        }
        std::vector<nodeID_t> kids = getNode( parentID ).getChildren();
        if( find( kids.begin(), kids.end(), leafIter->getFullID() ) == kids.end() )
        {
            std::cerr << "ERROR @ WHtree::check(): leaf parent doesnt have leaf ID among its children" << std::endl;
            return false;
        }
        if( leafIter->getSize() != 1 )
        {
            std::cerr << "ERROR @ WHtree::check(): leaf has a size other than 1" << std::endl;
            return false;
        }
        if( leafIter->getHLevel() != 0 )
        {
            std::cerr << "ERROR @ WHtree::check(): leaf has  hLevel other than 0" << std::endl;
            return false;
        }
        ++sumNodeKids[parentID.second];
    }

    // loop through nodes
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ); nodeIter != m_nodes.end(); ++nodeIter)
    {
        std::vector<nodeID_t> kids = nodeIter->getChildren();
        size_t currentHLevel( 0 ), currentSize( 0 );

        for( std::vector<nodeID_t>::const_iterator kidIter( kids.begin() ); kidIter != kids.end(); ++kidIter )
        {
            currentHLevel = std::max( currentHLevel, getNode( *kidIter ).getHLevel()+1 );
            currentSize += getNode( *kidIter ).getSize();
            nodeID_t kidParent = getNode( *kidIter ).getParent();
            if( kidParent != nodeIter->getFullID() )
            {
                std::cerr << "ERROR @ WHtree::check(): node child (" <<  kidIter->first << "-" <<  kidIter->second << ") doesnt have node ID (";
                std::cerr << nodeIter->isNode() << "-" <<  nodeIter->getID() << ") as its parent but instead has (" <<  kidParent.first << "-"
                          <<  kidParent.second << ")" << std::flush;
                return false;
            }
            if( kidIter->first )
            {
                ++sumNodeParents[kidIter->second];
            }
            else
            {
                ++sumLeafParents[kidIter->second];
            }
        }

        nodeID_t parentID( nodeIter->getParent() );
        if( ( !parentID.first) && ( ( nodeIter+1 ) != m_nodes.end() ) )
        {
            std::cerr << "ERROR @ WHtree::check(): node has a leaf as parent" << std::endl;
            return false;
        }
        if( ( !nodeIter->isRoot() ) && ( ( nodeIter+1 ) == m_nodes.end() ) )
        {
            std::cerr << "ERROR @ WHtree::check(): last node does not have 0-0 as parent" << std::endl;
            return false;
        }
        if( !nodeIter->isRoot() )
        {
            std::vector<nodeID_t> kids = getNode( parentID ).getChildren();
            if( find( kids.begin(), kids.end(), nodeIter->getFullID() ) == kids.end() )
            {
                std::cerr << "ERROR @ WHtree::check(): node parent doesnt have node ID among its children" << std::endl;
                return false;
            }
            ++sumNodeKids[parentID.second];
        }
        if( nodeIter->getSize() != currentSize )
        {
            std::cerr << "ERROR @ WHtree::check(): node " <<  nodeIter->getID() << " size (" << nodeIter->getSize()
                      << ") is not sum of its children sizes (" <<  currentSize << ")" << std::endl;
            return false;
        }
        if( nodeIter->getHLevel() != currentHLevel )
        {
            std::cerr << "ERROR @ WHtree::check(): node hLevel is not one more than its highest child" << std::endl;
            return false;
        }
    }

    //check consistency of counters
    for( std::vector<size_t>::const_iterator iter( sumLeafParents.begin() ); iter != sumLeafParents.end(); ++iter )
    {
        if( *iter != 1 )
        {
            std::cerr << "ERROR @ WHtree::check(): more than one node has the same leaf as child" << std::endl;
            return false;
        }
    }
    for( std::vector<size_t>::const_iterator iter( sumNodeParents.begin() ); iter != sumNodeParents.end()-1; ++iter )
    {
        if( *iter != 1 )
        {
            std::cerr << "ERROR @ WHtree::check(): more than one node has the same node as child" << std::endl;
            return false;
        }
    }
    if( sumNodeParents.back() != 0 )
    {
        std::cerr << "ERROR @ WHtree::check(): at least one node has the root node as child" << std::endl;
        return false;
    }
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ); nodeIter != m_nodes.end(); ++nodeIter)
    {
        std::vector<nodeID_t> kids = nodeIter->getChildren();
        if( kids.size() != sumNodeKids[nodeIter->getID()] )
        {
            std::cerr << "ERROR @ WHtree::check(): node children vector size does not match the number of nodes/leafs that have it as parent";
            std::cerr << std::endl << nodeIter->printAllData() << std::endl;
            return false;
        }
    }
    return true; // everything checked out fine
} // end check() -------------------------------------------------------------------------------------


const WHnode& WHtree::getNode( const nodeID_t &thisNode ) const
{
    if( thisNode.first )
    {
        return getNode( thisNode.second );
    }
    else
    {
        return getLeaf( thisNode.second );
    }
} // end getNode() -------------------------------------------------------------------------------------
const WHnode& WHtree::getNode( const size_t thisNode ) const
{
    if( thisNode >= m_nodes.size() )
    {
        std::cerr << "ERROR @ WHtree::getNode: index is out of boundaries(" << thisNode << ". total nodes: "
                  << m_nodes.size() << "), returning first node" << std::endl;
        return m_nodes.front();
    }
    else
    {
        return ( m_nodes[thisNode] );
    }
} // end getNode() -------------------------------------------------------------------------------------


const WHnode& WHtree::getLeaf( const size_t thisLeaf ) const
{
    if( thisLeaf >= m_leaves.size() )
    {
        std::cerr << "ERROR @ WHtree::getLeaf: index is out of boundaries (" << thisLeaf << ". total leaves: "
                  << m_leaves.size() << "), returning first leaf" << std::endl;
        return m_leaves.front();
    }
    else
    {
        return ( m_leaves[thisLeaf] );
    }
} // end getLeaf() -------------------------------------------------------------------------------------


const WHnode& WHtree::getRoot() const
{
    return ( m_nodes.back() );
} // end getRoot() -------------------------------------------------------------------------------------


size_t WHtree::getLeafID( const WHcoord &thisCoord ) const
{
    std::vector<WHcoord>::const_iterator findIter( std::find( m_coordinates.begin(), m_coordinates.end(), thisCoord ) );
    if( findIter == m_coordinates.end() )
    {
        throw std::runtime_error( "ERROR @ WHtree::getLeafID(): coordinate is not in the tree" );
    }
    else
    {
        return( findIter-m_coordinates.begin() );
    }
} // end getLeafID() -------------------------------------------------------------------------------------


std::vector<size_t> WHtree::getLeaves4node( const size_t nodeID ) const
{
    std::vector<size_t> returnVector;

    if( nodeID >= m_nodes.size() )
    {
        std::cerr << "ERROR @ WHtree::coordinate4leaf(): leafID is out of boundaries" << std::endl;
        return returnVector;
    }
    else
    {
        if( !m_containedLeaves.empty() ) // if contained leaves for each node have already been calculated, just return them
        {
            return m_containedLeaves[nodeID];
        }
        else  // if not, calculate them for this node
        {
            std::list<size_t> worklist;
            worklist.push_back( nodeID );
            while( !worklist.empty() )
            {
                size_t currentNode( worklist.front() );
                worklist.pop_front();
                std::vector<nodeID_t> kids( getNode( currentNode ).getChildren() );
                for( std::vector<nodeID_t>::const_iterator iter( kids.begin() ); iter != kids.end(); ++iter )
                {
                    if( iter->first) // is node
                    {
                        worklist.push_back( iter->second );
                    }
                    else    // is leaf
                    {
                        returnVector.push_back( iter->second );
                    }
                }
            }
            std::sort( returnVector.begin(), returnVector.end() );
            return returnVector;
         }
    }
} // end getLeaves4node() -------------------------------------------------------------------------------------
std::vector<size_t> WHtree::getLeaves4node( const nodeID_t &nodeID ) const
{
    if( nodeID.first )
    {
        return getLeaves4node( nodeID.second );
    }
    else
    {
        return std::vector<size_t>( 1, nodeID.second );
    }
} // end getLeaves4node() -------------------------------------------------------------------------------------


std::vector<size_t> WHtree::getBranchNodes( const size_t nodeID ) const
{
    std::vector<size_t> returnVector;
    if( nodeID >= getNumNodes() )
    {
        std::cerr << "ERROR @ WHtree::getBranchNodes(): nodeID is out of boundaries" << std::endl;
        return returnVector;
    }
    else
    {
        std::list<size_t> worklist;
        worklist.push_back( nodeID );
        while( !worklist.empty() )
        {
            size_t currentNode( worklist.front() );
            worklist.pop_front();
            returnVector.push_back( currentNode );
            std::vector<nodeID_t> kids( getNode( currentNode ).getChildren() );
            for( std::vector<nodeID_t>::const_iterator iter( kids.begin() ); iter != kids.end(); ++iter )
            {
                if( iter->first )
                {
                    worklist.push_back( iter->second ); // is node
                }
            }
        }
        std::sort( returnVector.begin(), returnVector.end() );
        return returnVector;
    }
} // end getBranchNodes() -------------------------------------------------------------------------------------


WHcoord WHtree::getCoordinate4leaf( const size_t leafID ) const
{
    if( leafID >= m_coordinates.size() )
    {
        std::cerr << "ERROR @ WHtree::coordinate4leaf(): leafID is out of boundaries" << std::endl;
        return WHcoord( 0, 0, 0 );
    }
    else
    {
        return m_coordinates[leafID];
    }
} // end getCoordinate4leaf() -------------------------------------------------------------------------------------


std::vector<WHcoord> WHtree::getCoordinates4node( const size_t nodeID ) const
{
    std::vector<WHcoord> returnVector;
    std::vector<size_t> containedLeaves( getLeaves4node( nodeID ) );

    for( std::vector<size_t>::const_iterator iter( containedLeaves.begin() ); iter != containedLeaves.end(); ++iter )
    {
        returnVector.push_back( getCoordinate4leaf( *iter ) );
    }
    return returnVector;
} // end getCoordinates4node() -------------------------------------------------------------------------------------
std::vector<WHcoord> WHtree::getCoordinates4node( const nodeID_t &nodeID ) const
{
    if( nodeID.first )
    {
        return getCoordinates4node( nodeID.second );
    }
    else
    {
        return std::vector<WHcoord>( 1, getCoordinate4leaf( nodeID.second ) );
    }
} // end getCoordinates4node() -------------------------------------------------------------------------------------


WHcoord WHtree::getMeanCoordinate4node( const size_t nodeID ) const
{
    WHcoord baseCoord;
    {
        std::vector<WHcoord> bnCoords( getCoordinates4node( nodeID ) );
        size_t sumX( 0 ), sumY( 0 ), sumZ( 0 );
        for( std::vector<WHcoord>::const_iterator coordIter( bnCoords.begin() ); coordIter != bnCoords.end(); ++coordIter )
        {
            sumX += coordIter->m_x;
            sumY += coordIter->m_y;
            sumZ += coordIter->m_z;
        }
        baseCoord.m_x = sumX/bnCoords.size();
        baseCoord.m_y = sumY/bnCoords.size();
        baseCoord.m_z = sumZ/bnCoords.size();
    }
    return baseCoord;
} // end getMeanCoordinate4node() -------------------------------------------------------------------------------------
WHcoord WHtree::getMeanCoordinate4node( const nodeID_t &nodeID ) const
{
    if( nodeID.first )
    {
        return getMeanCoordinate4node( nodeID.second );
    }
    else
    {
        return getCoordinate4leaf( nodeID.second );
    }
} // end getMeanCoordinate4node() -------------------------------------------------------------------------------------


size_t WHtree::getCommonAncestor( const size_t nodeID1, const size_t nodeID2 ) const
{
    if( nodeID1 == nodeID2 )
    {
        return nodeID1;
    }
    else
    {
        size_t tempID1( getNode( nodeID1 ).getID() ), tempID2( getNode( nodeID2 ).getID() );
        while( tempID1 != tempID2 )
        {
            if( tempID1 < tempID2 ) // node 1 is lower in the hierarchy
            {
                tempID1 = ( getNode( tempID1 ).getParent() ).second;
            }
            else  // seed 2 is is lower in the hierarchy
            {
                tempID2 = ( getNode( tempID2 ).getParent() ).second;
            }
        }
        return tempID1;
    }
} // end hTree::getCommonAncestor() -------------------------------------------------------------------------------------
nodeID_t WHtree::getCommonAncestor( const nodeID_t &nodeID1, const nodeID_t &nodeID2 ) const
{
    if( nodeID1 == nodeID2 )
    {
        return nodeID1;
    }
    else
    {
        size_t node1, node2;
        if( !nodeID1.first ) // its a leaf
        {
            node1 = ( getLeaf( nodeID1.second ).getParent() ).second;
        }
        else
        {
            node1 = ( nodeID1.second );
        }

        if( !nodeID2.first ) // its a leaf
        {
            node2 = ( getLeaf( nodeID2.second ).getParent() ).second;
        }
        else
        {
            node2 = ( nodeID2.second );
        }
        return std::make_pair( true, getCommonAncestor( node1, node2 ) );
    }
} // end hTree::getCommonAncestor() -------------------------------------------------------------------------------------


std::vector<nodeID_t> WHtree::getRoute2Root( const nodeID_t &nodeID ) const
{
    std::vector<nodeID_t> returnVector;
    if( ( nodeID.first ) && ( nodeID.second >= m_nodes.size() ) )
    {
        std::cerr << "ERROR @ WHtree::route2Root(): leafID is out of boundaries" << std::endl;
        return returnVector;
    }
    const WHnode& root( getRoot() );
    const WHnode* current( &getNode( nodeID ) );
    returnVector.reserve( root.getHLevel() - current->getHLevel() );
    returnVector.push_back( nodeID );

    while( !current->isRoot() )
    {
        current = &getNode( current->getParent() );
        returnVector.push_back( current->getFullID() );
    }

    return returnVector;
} // end hTree::getRoute2Root() -------------------------------------------------------------------------------------


unsigned int WHtree::getTripletOrder( const nodeID_t &nodeIDa, const nodeID_t &nodeIDb, const nodeID_t &nodeIDc ) const
{
    //0="non resolved" (a,b,c join same parent), 1="ab before c", 2="ac before b", 3="bc before a"

    nodeID_t ancestorAB( getCommonAncestor( nodeIDa, nodeIDb ) );
    nodeID_t ancestorAC( getCommonAncestor( nodeIDa, nodeIDc ) );
    nodeID_t ancestorBC( getCommonAncestor( nodeIDb, nodeIDc ) );

    if( ancestorAB == ancestorAC )
    {
        if( ancestorAB == ancestorBC )
        {
            return 0; // all join toghether
        }
        else
        {
            return 3; // a is the last to join
        }
    }
    else if( ancestorAB == ancestorBC )
    {
        return 2; // b is the last to join
    }
    else
    {
        return 1; // c is the last to join
    }
} // end hTree::getTripletOrder() -------------------------------------------------------------------------------------


std::vector<size_t> WHtree::getBaseNodes( const size_t root ) const
{
    if( root > getRoot().getID() )
    {
        std::cerr << "ERROR @ WHtree::getBaseNodes(): branch root ID is out of boundaries (ID: "
                  <<  root << ", # nodes: " << getNumNodes() << ")." << std::endl;
        return std::vector<size_t> ();
    }

    std::list<size_t> baseList;
    for( std::vector<WHnode>::const_iterator leafIter( m_leaves.begin() ); leafIter != m_leaves.end(); ++leafIter )
    {
        baseList.push_back( ( leafIter->getParent() ).second );
    }

    baseList.sort();
    baseList.unique();
    std::vector<size_t> returnVector;

    if( root != getRoot().getID() )
    {
        for( std::list<size_t>::iterator listIter( baseList.begin() ); listIter != baseList.end(); ++listIter )
        {
            const WHnode* current( &getNode( *listIter ) );
            while( !current->isRoot() )
            {
                if( current->getID() ==root )
                {
                    returnVector.push_back( getNode( *listIter ).getID() );
                    break;
                }
                current = &getNode( current->getParent() );
            }
        }
    }
    else
    {
        for( std::list<size_t>::const_iterator baseIter( baseList.begin() ); baseIter != baseList.end(); ++baseIter )
        {
            returnVector.push_back( *baseIter );
        }
    }
    return returnVector;
} // end getBaseNodes() -------------------------------------------------------------------------------------
std::vector<nodeID_t> WHtree::getBaseNodes( const nodeID_t &root ) const
{
    std::vector<nodeID_t> returnVector;

    if( !root.first ) // if its a leaf
    {
        return std::vector<nodeID_t>();
    }
    else
    {
        std::vector<size_t> bases( getBaseNodes( root.second ) );
        for( size_t i = 0; i < bases.size(); ++i )
        {
            returnVector.push_back( std::make_pair( true, bases[i] ) );
        }
        return returnVector;
    }
} // end getBaseNodes() -------------------------------------------------------------------------------------


std::vector<size_t> WHtree::getRootBaseNodes() const
{
    return getBaseNodes( getRoot().getID() );
} // end getRootBaseNodes() -------------------------------------------------------------------------------------


bool WHtree::testRootBaseNodes() const
{
    std::vector<size_t> bases( getRootBaseNodes() );
    if( bases.empty() )
    {
        return false;
    }
    for( size_t i = 0; i <  bases.size(); ++i )
    {
        if( getNode( bases[i] ).getHLevel() > 1 )
        {
            return false;
        }
    }
    return true;
} // end testRootBaseNodes() -------------------------------------------------------------------------------------


dist_t WHtree::getDistance( const size_t nodeID1, const size_t nodeID2 ) const
{
    size_t ancestor( getCommonAncestor( nodeID1, nodeID2 ) );
    return getNode( ancestor ).getDistLevel();
} // end hTree::getNodeDistance() -------------------------------------------------------------------------------------
dist_t WHtree::getDistance( const nodeID_t &nodeID1, const nodeID_t &nodeID2 ) const
{
    nodeID_t ancestor( getCommonAncestor( nodeID1, nodeID2 ) );
    return getNode( ancestor ).getDistLevel();
} // end hTree::getDistance() -------------------------------------------------------------------------------------
dist_t WHtree::getDistance( const WHcoord &coord1, const WHcoord &coord2 ) const
{
    nodeID_t nodeID1( std::make_pair( false, getLeafID( coord1 ) ) );
    nodeID_t nodeID2( std::make_pair( false, getLeafID( coord2 ) ) );
    return getDistance( nodeID1, nodeID2 );
} // end hTree::distance() -------------------------------------------------------------------------------------


dist_t WHtree::getLeafDistance( const size_t leafID1, const size_t leafID2 ) const
{
    nodeID_t ancestor( getCommonAncestor( std::make_pair( false, leafID1 ), std::make_pair( false, leafID2 ) ) );
    return getNode( ancestor ).getDistLevel();
} // end hTree::getLeafDistance() -------------------------------------------------------------------------------------


void WHtree::sortBySize( std::vector<size_t>* nodeVector ) const
{
    std::vector<size_t> &nodeVectorRef( *nodeVector );
    for( size_t i = 0; i < nodeVectorRef.size(); ++i )
    {
        if( nodeVectorRef[i] >= getNumNodes() )
        {
            std::cerr << "ERROR @ WHtree::sortBySize(): indices out of bounds (" << nodeVectorRef[i];
            std::cerr << ". total nodes: " << getNumNodes() << ")" << std::endl;
            return;
        }
    }
    std::sort( nodeVectorRef.begin(), nodeVectorRef.end(), compSize( this ) );
    return;
} // end sortBySize() -------------------------------------------------------------------------------------
void WHtree::sortBySize( std::list<size_t>* nodeList ) const
{
    std::list<size_t> &nodeListRef( *nodeList );
    for( std::list<size_t>::const_iterator iter = nodeListRef.begin(); iter != nodeListRef.end(); ++iter )
    {
        if( *iter >= getNumNodes() )
        {
            std::cerr << "ERROR @ WHtree::sortBySize(): indices out of bounds (" << *iter << ". total nodes: " << getNumNodes() << ")" << std::endl;
            return;
        }
    }
    nodeListRef.sort( compSize( this ) );
    return;
} // end sortBySize() -------------------------------------------------------------------------------------
void WHtree::sortBySize( std::vector<nodeID_t>* nodeVector ) const
{
    std::vector<nodeID_t> &nodeVectorRef( *nodeVector );
    for( size_t i = 0; i < nodeVectorRef.size(); ++i )
    {
        if( ( nodeVectorRef[i].first && nodeVectorRef[i].second >= getNumNodes() ) ||
            ( !nodeVectorRef[i].first && nodeVectorRef[i].second >= getNumLeaves() ) )
        {
            std::cerr << "ERROR @ WHtree::sortBySize(): indices out of bounds(" << nodeVectorRef[i].first << "-" << nodeVectorRef[i].second
                      << ". total leaves: " << getNumLeaves() << ". total nodes: " << getNumNodes() << ")" << std::endl;
            return;
        }
    }
    std::sort( nodeVectorRef.begin(), nodeVectorRef.end(), compSize( this ) );
    return;
} // end sortBySize() -------------------------------------------------------------------------------------
void WHtree::sortBySize( std::list<nodeID_t>* nodeList ) const
{
    std::list<nodeID_t> &nodeListRef( *nodeList );
    for( std::list<nodeID_t>::const_iterator iter = nodeListRef.begin(); iter != nodeListRef.end(); ++iter )
    {
        if( ( iter->first && iter->second >= getNumNodes() ) ||
            ( !iter->first && iter->second >= getNumLeaves() ) )
        {
            std::cerr << "ERROR @ WHtree::sortBySize(): indices out of bounds(" << iter->first << "-" << iter->second
                      << ". total leaves: " << getNumLeaves() << ". total nodes: " << getNumNodes() << ")" << std::endl;
            return;
        }
    }
    nodeListRef.sort( compSize( this ) );
    return;
} // end sortBySize() -------------------------------------------------------------------------------------

void WHtree::sortByHLevel( std::vector<size_t>* nodeVector ) const
{
    std::sort( nodeVector->begin(), nodeVector->end(), compHLevel( this ) );
    return;
} // end sortByHLevel() -------------------------------------------------------------------------------------
void WHtree::sortByHLevel( std::vector<nodeID_t>* nodeVector ) const
{
    std::sort( nodeVector->begin(), nodeVector->end(), compHLevel( this ) );
    return;
} // end sortByHLevel() -------------------------------------------------------------------------------------
void WHtree::sortByHLevel( std::list<size_t>* nodeList ) const
{
    nodeList->sort( compHLevel( this ) );
    return;
} // end sortByHLevel() -------------------------------------------------------------------------------------
void WHtree::sortByHLevel( std::list<nodeID_t>* nodeList ) const
{
    nodeList->sort( compHLevel( this ) );
    return;
} // end sortByHLevel() -------------------------------------------------------------------------------------

void WHtree::loadContainedLeaves()
{
    std::vector<size_t> emptyVect;
    m_containedLeaves.resize( m_nodes.size(), emptyVect );

    for( std::vector<WHnode>::const_iterator leafIter( m_leaves.begin() ); leafIter != m_leaves.end(); ++leafIter )
    {
        nodeID_t thisParent( leafIter->getParent() );
        m_containedLeaves[thisParent.second].push_back( leafIter->getID() );
    }
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ); nodeIter != m_nodes.end(); ++nodeIter )
    {
        std::sort( m_containedLeaves[nodeIter->getID()].begin(), m_containedLeaves[nodeIter->getID()].end() );
        if( !nodeIter->isRoot() )
        {
            m_containedLeaves[nodeIter->getParent().second].insert( m_containedLeaves[nodeIter->getParent().second].end(),
                                                                    m_containedLeaves[nodeIter->getID()].begin(),
                                                                    m_containedLeaves[nodeIter->getID()].end() );
        }
    }
} // end loadContainedLeaves() -------------------------------------------------------------------------------------

void WHtree::clearContainedLeaves()
{
    std::vector<std::vector<size_t> > emptyThing;
    m_containedLeaves.swap( emptyThing );
} // end clearContainedLeaves() -------------------------------------------------------------------------------------


bool WHtree::convert2grid( const HC_GRID newGrid )
{
    if( m_datasetGrid == newGrid )
    {
        return false;
    }
    else if( m_datasetGrid == HC_VISTA && newGrid == HC_NIFTI )
    {
        for( std::vector<WHcoord>::iterator iter( m_coordinates.begin() ); iter != m_coordinates.end(); ++iter )
        {
            WHcoord tempCoord = *iter;
            *iter = tempCoord.vista2nifti( m_datasetSize );
        }
        for( std::list<WHcoord>::iterator iter( m_discarded.begin() ); iter != m_discarded.end(); ++iter )
        {
            WHcoord tempCoord = *iter;
            *iter = tempCoord.vista2nifti( m_datasetSize );
        }
        m_datasetGrid = HC_NIFTI;
        return true;
    }
    else if( m_datasetGrid == HC_NIFTI && newGrid == HC_VISTA )
    {
        for( std::vector<WHcoord>::iterator iter( m_coordinates.begin() ); iter != m_coordinates.end(); ++iter )
        {
            WHcoord tempCoord = *iter;
            *iter = tempCoord.nifti2vista( m_datasetSize );
        }
        for( std::list<WHcoord>::iterator iter( m_discarded.begin() ); iter != m_discarded.end(); ++iter )
        {
            WHcoord tempCoord=*iter;
            *iter = tempCoord.nifti2vista( m_datasetSize );
        }
        m_datasetGrid = HC_VISTA;
        return true;
    }
    else
    {
        std::cerr << "ERROR @ WHtree::convert2grid(): coordinate grid not recognized" << std::endl;
        return false;
    }
} // end convert2grid() -------------------------------------------------------------------------------------


bool WHtree::readTree( const std::string &filename )
{
    m_nodes.clear();
    m_leaves.clear();
    m_coordinates.clear();
    m_trackids.clear();
    m_discarded.clear();

    WFileParser parser( filename );
    if( !parser.readFile() )
    {
        std::cerr << "ERROR @ WHtree::readTree(): Parser error" << std::endl;
        return false;
    }

    std::vector<std::string> lines = parser.getRawLines();
    if( lines.size() == 0 )
    {
        std::cerr << "ERROR @ WHtree::readTree(): File is empty" << std::endl;
        return false;
    }

    {
        std::vector< std::vector< std::string> >datasetStrings = parser.getLinesForTagSeparated( "imagesize" );
        if( datasetStrings.size() == 0 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): Dataset size was not found in tree file" << std::endl;
            return false;
        }
        if( datasetStrings.size() > 1 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): Dataset attribute had multiple lines" << std::endl;
            return false;
        }
        WHcoord datasetSize( string_utils::fromString<coord_t>( datasetStrings[0][0] ),
                                string_utils::fromString<coord_t>( datasetStrings[0][1] ),
                                string_utils::fromString<coord_t>( datasetStrings[0][2] ) );
        std::string gridString( datasetStrings[0][3] );
        if( gridString == getGridString( HC_VISTA ) )
        {
            m_datasetGrid = HC_VISTA;
        }
        else if( gridString == getGridString( HC_NIFTI ) )
        {
            m_datasetGrid = HC_NIFTI;
        }
        else
        {
            std::cerr << "ERROR @ WHtree::readTree(): Dataset grid type string \"" << gridString << "\" could not be identified" << std::endl;
            return false;
        }
        m_datasetSize = datasetSize;
    }

    {
        std::vector< std::vector< std::string > > streamNumberStrings = parser.getLinesForTagSeparated( "streams" );
        if( streamNumberStrings.size() == 0 )
        {
            std::cerr << "WARNING @ WHtree::readTree(): tracking streams number was not found in tree file,";
            std::cerr << " assuming streams=0 for compatibility" << std::endl;
            m_numStreamlines = 0;
        }
        if( streamNumberStrings.size() > 1 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): tracking streams number attribute has multiple lines" << std::endl;
            return false;
        }
        if( streamNumberStrings[0].size() > 1 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): tracking streams number attribute has multiple elements" << std::endl;
            return false;
        }
        m_numStreamlines = string_utils::fromString< size_t >( streamNumberStrings[0][0] );
    }

    {
        std::vector< std::vector< std::string > >logFactorStrings = parser.getLinesForTagSeparated( "logfactor" );
        if( logFactorStrings.size() == 0 )
        {
            std::cerr << "WARNING @ WHtree::readTree(): logarithmic normalization factor was not found in tree file,";
            std::cerr << " assuming logFactor=0 for compatibility" << std::endl;
            m_logFactor = 0;
        }
        if( logFactorStrings.size() > 1 )
        {
            std::cerr << "ERROR @ WHtree::readTree():";
            std::cerr << "logarithmic normalization factor attribute has multiple lines" << std::endl;
            return false;
        }
        if( logFactorStrings[0].size() > 1 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): logarithmic normalization factor attribute has multiple elements" << std::endl;
            return false;
        }
        m_logFactor = string_utils::fromString< float >( logFactorStrings[0][0] );

        if( m_logFactor != 0 && m_numStreamlines != 0 &&  std::fabs( m_logFactor - log10( m_numStreamlines ) ) > 0.00001 )
        {
            std::cerr << "ERROR @ WHtree::readTree(): tracking streams number (";
            std::cerr << m_numStreamlines << ") and logarithmic normalization factor (";
            std::cerr << m_logFactor << ") are a missmatch . Log factor should be: "<< log10( m_numStreamlines) << std::endl;
            return false;
        }
    }

    {
        std::vector< std::vector< std::string> >coordStrings = parser.getLinesForTagSeparated( "coordinates" );
        m_coordinates.reserve( coordStrings.size() );
        m_leaves.reserve( coordStrings.size() );
        size_t leafCount( 0 );
        for( size_t i = 0; i < coordStrings.size(); ++i )
        {
            WHcoord tempCoord( string_utils::fromString<coord_t>( coordStrings[i][0] ),
                                  string_utils::fromString<coord_t>( coordStrings[i][1] ),
                                  string_utils::fromString<coord_t>( coordStrings[i][2] ) );
            WHnode tempNode( std::make_pair( 0, leafCount ) );
            m_leaves.push_back( tempNode );
            m_coordinates.push_back( tempCoord );
            ++leafCount;
        }
    }
    {
        std::vector< std::vector< std::string > > indexStrings = parser.getLinesForTagSeparated( "trackindex" );
        if( indexStrings.empty() )
        {
            if( m_datasetGrid == HC_NIFTI )
            {
                std::cerr << "ERROR @ WHtree::readTree(): no tract ids in roi file, necessary to work on nifti mode" << std::endl;
                return false;
            }
            else
            {
                m_trackids.reserve( m_coordinates.size() );
                for( size_t i = 0; i < m_coordinates.size(); ++i )
                {
                    m_trackids.push_back( i );
                }
            }
        }
        else
        {
            m_trackids.reserve( indexStrings.size() );
            for( size_t i = 0; i < indexStrings.size(); ++i )
            {
                size_t tempIndex( string_utils::fromString< size_t >( indexStrings[i][0] ) );
                m_trackids.push_back( tempIndex );
            }
        }
    }
    {
        std::vector< std::vector< std::string> >clusterStrings = parser.getLinesForTagSeparated( "clusters" );
        m_nodes.reserve( clusterStrings.size() );
        size_t nodeCount( 0 );
        for( size_t i = 0; i < clusterStrings.size(); ++i )
        {
            dist_t distance( string_utils::fromString< dist_t > ( clusterStrings[i][0] ) );
            std::vector<nodeID_t>joiningNodes;
            for( size_t j = 1; j < clusterStrings[i].size(); j += 2 )
            {
                joiningNodes.push_back( std::make_pair( string_utils::fromString<bool>( clusterStrings[i][j] ),
                                                        string_utils::fromString<size_t>( clusterStrings[i][j+1] ) ) );
            }

            size_t tempSize( 0 ), tempHLevel( 0 );
            nodeID_t tempID( std::make_pair( 1, nodeCount ) );
            for( std::vector<nodeID_t>::const_iterator iter( joiningNodes.begin() );  iter != joiningNodes.end(); ++iter)
            {
                WHnode* kid( fetchNode( *iter ) );
                if( kid == 0 )
                {
                    std::cerr << "ERROR @ WHtree::readTree(): kid id (" <<  iter->first << "-" <<  iter->second << ") was out of boundaries. Nodes: "
                              << m_nodes.size() << std::endl;
                    return false;
                }
                tempSize += kid->getSize();
                tempHLevel = std::max( tempHLevel, kid->getHLevel() );
                kid->setParent( tempID );
            }
            ++tempHLevel;

            WHnode tempNode( tempID, joiningNodes, tempSize, distance, tempHLevel );
            //std::cout << tempNode.printAllData() << std::endl;
            m_nodes.push_back( tempNode );
            ++nodeCount;
        }
    }

    {
        std::vector< std::vector< std::string> >discardedStrings = parser.getLinesForTagSeparated( "discarded" );
        for( size_t i = 0; i < discardedStrings.size(); ++i )
        {
            WHcoord tempCoord( string_utils::fromString<coord_t>( discardedStrings[i][0] ),
                                  string_utils::fromString<coord_t>( discardedStrings[i][1] ),
                                  string_utils::fromString<coord_t>( discardedStrings[i][2] ) );
            m_discarded.push_back( tempCoord );
        }
        m_discarded.sort();
    }

    {
        std::vector< std::vector< std::string> >cpccStrings = parser.getLinesForTagSeparated( "cpcc" );
        if( !cpccStrings.empty() )
        {
            if( cpccStrings.size() > 1 || cpccStrings[0].size() > 1 )
            {
                std::cerr << "ERROR @ WHtree::readTree(): multiple objects on cpcc attribute" << std::endl;
                return false;
            }
            m_cpcc = string_utils::fromString<float>( cpccStrings[0][0] );
        }
    }

    {
        m_selectedValues.clear();
        m_selectedPartitions.clear();
        std::vector< std::vector< std::string> >partValueStrings = parser.getLinesForTagSeparated( "partvalues" );
        if( !partValueStrings.empty() )
        {
            m_selectedValues.reserve( partValueStrings.size() );
            for( size_t i = 0; i < partValueStrings.size(); ++i )
            {
                float value( string_utils::fromString< float > ( partValueStrings[i][0] ) );
                m_selectedValues.push_back( value );
            }
        }

        std::vector< std::vector< std::string> >partitionStrings = parser.getLinesForTagSeparated( "partitions" );
        if( !partitionStrings.empty() )
        {
            m_selectedPartitions.reserve( partitionStrings.size() );
            for( size_t i = 0; i < partitionStrings.size(); ++i )
            {
                std::vector< size_t > thisPartition;
                thisPartition.reserve( partitionStrings[i].size() );
                for( size_t j = 0; j < partitionStrings[i].size(); ++j )
                {
                    size_t partCluster( string_utils::fromString< size_t > ( partitionStrings[i][j] ) );
                    thisPartition.push_back( partCluster );
                }
                m_selectedPartitions.push_back( thisPartition );
            }
        }

        std::vector< std::vector< std::string > >partcolorStrings = parser.getLinesForTagSeparated( "partcolors" );
        if( !partcolorStrings.empty() )
        {
            m_selectedColors.reserve( partcolorStrings.size() );
            for( size_t i = 0; i < partcolorStrings.size(); ++i )
            {
                std::vector< WHcoord > thisPartColors;
                thisPartColors.reserve( partcolorStrings[i].size() );

                for( size_t j = 0; j < partcolorStrings[i].size(); ++j )
                {
                    std::string thisCoordstring( partcolorStrings[i][j] );
                    if( thisCoordstring.size() != 11 )
                    {
                        std::cerr << "ERROR @ WHtree::readTree(): partition colors have wrong size (" << thisCoordstring.size();
                        std::cerr << ") while it should be 11. string: " << thisCoordstring << std::endl;
                        m_selectedColors.clear();
                        break;
                    }

                        std::string colorR( thisCoordstring.substr( 0, 3 ) );
                        std::string colorG( thisCoordstring.substr( 4, 3 ) );
                        std::string colorB( thisCoordstring.substr( 8, 3 ) );


                        WHcoord thisColor( string_utils::fromString< coord_t > ( colorR ),
                                           string_utils::fromString< coord_t > ( colorG ),
                                           string_utils::fromString< coord_t > ( colorB ) );
                        thisPartColors.push_back( thisColor );
                }
                m_selectedColors.push_back( thisPartColors );
            }
        }

        if( m_selectedColors.size() != 0 )
        {
            if( m_selectedColors.size() != m_selectedPartitions.size() )
            {
                std::cerr << "ERROR @ WHtree::readTree(): partition and colors dimensions dont match. Color field will be left empty" << std::endl;
                m_selectedColors.clear();
            }
            else
            {
                for( size_t i = 0; i < m_selectedColors.size(); ++i )
                {
                    if( m_selectedColors[i].size() != m_selectedPartitions[i].size() )
                    {
                        std::cerr << "ERROR @ WHtree::readTree(): partition and colors dimensions dont match.";
                        std::cerr << " Color field will be left empty" << std::endl;
                        m_selectedColors.clear();
                        break;
                    }
                }
            }
        }
        if( m_selectedPartitions.size() != m_selectedValues.size() )
        {
            std::cerr << "ERROR @ WHtree::readTree(): partition and value dimensions dont match. Fields will be left empty" << std::endl;
            clearPartitions();
        }
    }



    if( !check() )
    {
        std::cerr << "ERROR @ WHtree::readTree(): loaded tree is not consistent" << std::endl;
        return false;
    }

    m_treeName = boost::filesystem::path( filename ).stem().string();


    m_loadStatus = true;
    return true;
} // end readTree() -------------------------------------------------------------------------------------



bool WHtree::writeTree( const std::string &filename, const bool niftiMode ) const
{
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR: unable to open out file: \"" << outFile << "\"" << std::endl;
        exit( -1 );
    }

    std::string gridString;
    if( niftiMode )
    {
        gridString = getGridString( HC_NIFTI );
    }
    else
    {
        gridString = getGridString( HC_VISTA );
    }

    outFile << "#imagesize" << std::endl << m_datasetSize << " " << gridString << std::endl << "#endimagesize" << std::endl;
    outFile << std::endl;
    if( m_cpcc != 0 )
    {
        outFile << "#cpcc" << std::endl << string_utils::toString( m_cpcc ) << std::endl << "#endcpcc" << std::endl << std::endl;
    }

    outFile << "#streams" << std::endl << m_numStreamlines << std::endl << "#endstreams" << std::endl;

    outFile << "#logfactor" << std::endl << m_logFactor << std::endl << "#endlogfactor" << std::endl;

    outFile << "#coordinates" << std::endl;
    for( std::vector<WHcoord>::const_iterator coordIter( m_coordinates.begin() ) ; coordIter != m_coordinates.end() ; ++coordIter )
    {
        WHcoord currentCoord( *coordIter );
        if( niftiMode )
        {
            if( m_datasetGrid == HC_VISTA )
            {
                currentCoord = currentCoord.vista2nifti( m_datasetSize );
            }
        }
        else
        {
            if( m_datasetGrid == HC_NIFTI )
            {
                currentCoord = currentCoord.nifti2vista( m_datasetSize );
            }
        }
        outFile << currentCoord << std::endl;
    }
    outFile << "#endcoordinates" << std::endl << std::endl;

    outFile << "#trackindex" << std::endl;
    for( std::vector<size_t>::const_iterator indexIter( m_trackids.begin() ) ; indexIter != m_trackids.end() ; ++indexIter )
    {
        outFile << *indexIter << std::endl;
    }
    outFile << "#endtrackindex" << std::endl << std::endl;

    outFile << "#clusters" << std::endl;
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ) ; nodeIter != m_nodes.end() ; ++nodeIter )
    {
        outFile << *nodeIter << std::endl;
    }
    outFile << "#endclusters" << std::endl << std::endl;

    outFile << "#discarded" << std::endl;
    for( std::list<WHcoord>::const_iterator coordIter( m_discarded.begin() ) ; coordIter != m_discarded.end() ; ++coordIter )
    {
        WHcoord currentCoord( *coordIter );
        if( niftiMode )
        {
            if( m_datasetGrid == HC_VISTA )
            {
                currentCoord = currentCoord.vista2nifti( m_datasetSize );
            }
        }
        else
        {
            if( m_datasetGrid == HC_NIFTI )
            {
                currentCoord = currentCoord.nifti2vista( m_datasetSize );
            }
        }
        outFile << currentCoord << std::endl;
    }
    outFile << "#enddiscarded" << std::endl;

    if( m_selectedValues.size() !=0 )
    {
        outFile << std::endl << "#partvalues" << std::endl;
        for( size_t i = 0; i < m_selectedValues.size(); ++i )
        {
            outFile << m_selectedValues[i] << std::endl;
        }
        outFile << "#endpartvalues" << std::endl;

        outFile << std::endl << "#partitions" << std::endl;
        for( size_t i = 0; i < m_selectedPartitions.size(); ++i )
        {
            for( size_t j = 0; j < m_selectedPartitions[i].size(); ++j )
            {
                outFile <<m_selectedPartitions[i][j]<< " ";
            }
            outFile << std::endl;
        }
        outFile << "#endpartitions" << std::endl;

        if( m_selectedColors.size() !=0 )
        {
            outFile << std::endl << "#partcolors" << std::endl;
            for( size_t i = 0; i < m_selectedColors.size(); ++i )
            {
                for( size_t j = 0; j < m_selectedColors[i].size(); ++j )
                {
                    WHcoord thisColor( m_selectedColors[i][j] );
                    size_t xcolor( thisColor.m_x );
                    size_t ycolor( thisColor.m_y );
                    size_t zcolor( thisColor.m_z );

                    outFile << boost::format( "%03d;%03d;%03d " ) % xcolor  % ycolor  % zcolor;
                }
                outFile << std::endl;
            }
            outFile << "#endpartcolors" << std::endl;
        }
    }



    return true;
} // end writeTree() -------------------------------------------------------------------------------------


bool WHtree::writeTreeDebug( const std::string &filename ) const
{
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR: unable to open out file: \"" << outFile << "\"" << std::endl;
        exit( -1 );
    }

    outFile << "Dataset size: " << m_datasetSize << " " << getGridString( m_datasetGrid ) << std::endl;

    if( m_cpcc != 0 )
    {
        outFile << "CPCC: " <<  string_utils::toString( m_cpcc ) << std::endl << std::endl;
    }

    outFile << "Streamlines per seed voxel: " << m_numStreamlines << std::endl;

    outFile << "Logarithmic normalization factor: " << m_logFactor << std::endl << std::endl;

    outFile << "============LEAVES============" << std::endl << std::endl;
    for( std::vector<WHnode>::const_iterator leafIter( m_leaves.begin() ); leafIter != m_leaves.end(); ++leafIter )
    {
        WHcoord currentCoord( getCoordinate4leaf( leafIter->getID() ) );
        outFile << "Coord: " << currentCoord << " " <<  leafIter->printAllData() << std::endl;
    }
    outFile << std::endl << std::endl << "============NODES============" << std::endl << std::endl;
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ) ; nodeIter != m_nodes.end() ; ++nodeIter )
    {
        outFile << nodeIter->printAllData() << std::endl;
    }

    return true;
} // end writeTreeDebug() -------------------------------------------------------------------------------------


bool WHtree::writeTreeOldWalnut( const std::string &filename ) const
{
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR: unable to open out file: \"" << outFile << "\"" << std::endl;
        exit( -1 );
    }


    outFile << "#coordinates" << std::endl;
    for( std::vector<WHcoord>::const_iterator coordIter( m_coordinates.begin() ) ; coordIter != m_coordinates.end() ; ++coordIter )
    {
        WHcoord currentCoord( *coordIter );
        if( m_datasetGrid == HC_VISTA )
        {
            currentCoord = currentCoord.vista2nifti( m_datasetSize );
        }
        outFile << boost::format( "%03d,%03d,%03d\n" ) % currentCoord.m_x  % currentCoord.m_y  % currentCoord.m_z;
    }
    outFile << "#endcoordinates" << std::endl << std::endl;

    outFile << "#clusters" << std::endl;
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ); nodeIter != m_nodes.end(); ++nodeIter)
    {
        std::vector<nodeID_t> currentKids( nodeIter->getChildren() );
        for( size_t i = 0; i < currentKids.size(); ++i )
        {
            size_t currentID( currentKids[i].second );
            if( currentKids[i].first )
            {
                currentID += getNumLeaves();
            }
            outFile << boost::format( "%06d" ) % currentID << "," << std::flush;
        }
        outFile << string_utils::toString( nodeIter->getDistLevel() ) << std::endl;
    }
    outFile << "#endclusters" << std::endl << std::endl;

return true;
} // end writeTreeOldWalnut() -------------------------------------------------------------------------------------

bool WHtree::writeTreeSimple( const std::string &filename ) const
{
    std::ofstream outFile( filename.c_str() );
    if( !outFile )
    {
        std::cerr << "ERROR: unable to open out file: \"" << outFile << "\"" << std::endl;
        exit( -1 );
    }

    outFile << getNumLeaves() << std::endl;
    for( std::vector<WHnode>::const_iterator nodeIter( m_nodes.begin() ) ; nodeIter != m_nodes.end() ; ++nodeIter )
    {
        outFile << *nodeIter << std::endl;
    }
    return true;
} // end writeTreeSimple() -------------------------------------------------------------------------------------

void WHtree::insertPartitions( const std::vector<std::vector<size_t> > &selectedPartitions,  const std::vector< float > &selectedValues,
                               const std::vector<std::vector<WHcoord> > &selectedColors )
{
    clearPartitions();
    if( selectedPartitions.size() != selectedValues.size() )
    {
        std::cerr << "ERROR @ WHtree::insertPartitions(): inserted partition set and partition value vector have different dimensions" << std::endl;
    }
    else
    {
        m_selectedPartitions = selectedPartitions;
        m_selectedValues = selectedValues;
    }

    if( !selectedColors.empty() )
    {
        if( selectedColors.size() != selectedPartitions.size() )
        {
            std::cerr << "ERROR @ WHtree::insertPartitions(): inserted partition color set and partition set have different dimensions" << std::endl;
        }
        else
        {
            for( size_t i = 0; i < selectedColors.size(); ++i )
            {
                if( selectedColors[i].size() != selectedPartitions[i].size() )
                {
                    std::cerr << "ERROR @ WHtree::insertPartitions(): partition and colors dimensions dont match";
                    std::cerr << " (" << selectedPartitions[i].size() << "-" << selectedColors[i].size();
                    std::cerr << ") Color field will be left empty" << std::endl;
                    return;
                }
            }
            m_selectedColors = selectedColors;
        }
    }
    return;
} // end insertPartitions() -------------------------------------------------------------------------------------

void WHtree::insertPartColors( const std::vector<std::vector<WHcoord> > &selectedColors )
{
    if( selectedColors.size() != m_selectedPartitions.size() )
    {
        std::cerr << "ERROR @ WHtree::insertPartColors(): inserted partition color set and partition set have different dimensions" << std::endl;
    }
    else
    {
        for( size_t i = 0; i < selectedColors.size(); ++i )
        {
            if( selectedColors[i].size() != m_selectedPartitions[i].size() )
            {
                std::cerr << "ERROR @ WHtree::insertPartColors(): partition and colors dimensions dont match.";
                std::cerr << " Color field will be left empty" << std::endl;
                clearPartColors();
                return;
            }
        }
        m_selectedColors = selectedColors;
    }
    return;
} // end insertPartColors() -------------------------------------------------------------------------------------

void WHtree::clearPartitions()
{
    std::vector<std::vector<size_t> > empty1;
    std::vector< float > empty2;
    std::vector<std::vector<WHcoord> > empty3;
    m_selectedPartitions.swap( empty1 );
    m_selectedValues.swap( empty2 );
    m_selectedColors.swap( empty3 );
    return;
} // end clearPartitions() -------------------------------------------------------------------------------------

void WHtree::clearPartColors()
{
    std::vector<std::vector<WHcoord> > empty;
    m_selectedColors.swap( empty );
    return;
} // end clearPartColors() -------------------------------------------------------------------------------------

std::vector< std::vector< unsigned int > > WHtree::getBranching( const std::vector < nodeID_t > &thisPartition,
                                                                 size_t depthLevel,
                                                                 std::vector< std::vector < nodeID_t > > *partitionSet,
                                                                 const bool excludeLeaves )
{
    if( depthLevel == 0 )
    {
        return std::vector< std::vector< unsigned int > >();
    }

    std::vector< std::vector < nodeID_t > > addedPartitionSet;
    std::vector< std::vector< unsigned int > > addedIndexTable;
    addedIndexTable.reserve( thisPartition.size() );
    addedPartitionSet.reserve( thisPartition.size() );

    for( size_t i = 0; i < thisPartition.size(); ++i)
    {
        // if its a base node and flag is set not to divide them, skip
        if( getNode( thisPartition[i] ).getHLevel() == 1 && excludeLeaves)
        {
            continue;
        }
        // get branched sub/partition for every cluster
        std::vector<nodeID_t> branch;
        {
            std::vector<nodeID_t> kids( getNode( thisPartition[i] ).getChildren() );
            branch.reserve( kids.size() );
            for( size_t j = 0; j < kids.size(); ++j)
            {
                branch.push_back( kids[j] );
            }
        }
        // insert branched partition
        {
            std::vector < nodeID_t > newPartition( thisPartition );
            newPartition.erase( newPartition.begin()+i );
            newPartition.insert( newPartition.begin()+i, branch.begin(), branch.end() );

            addedPartitionSet.push_back( newPartition );
            addedIndexTable.push_back( std::vector<unsigned int>( 1, i ) );
        }
        // obtain further sub-partitions if depth level continues
        if( depthLevel > 1 )
        {
            std::vector< std::vector < nodeID_t > > subPartitionSet;
            std::vector< std::vector< unsigned int > > subIndexTable( getBranching( branch, depthLevel-1, &subPartitionSet, excludeLeaves ) );

            addedPartitionSet.reserve( addedPartitionSet.size() + subPartitionSet.size() );
            addedIndexTable.reserve( addedPartitionSet.size() + subPartitionSet.size() );

            if( subPartitionSet.size() != subIndexTable.size() )
            {
                throw std::runtime_error( "ERROR @ WHtree::getBranching(): dimension error on obtained vectors" );
            }

            for( size_t j = 0; j < subIndexTable.size(); ++j)
            {
                std::vector < nodeID_t > newPartition( thisPartition );
                newPartition.erase( newPartition.begin()+i );
                newPartition.insert( newPartition.begin()+i, subPartitionSet[j].begin(), subPartitionSet[j].end() );
                addedPartitionSet.push_back( newPartition );

                std::vector< unsigned int > newIndexEntry( 1, i );
                newIndexEntry.insert( newIndexEntry.end(), subIndexTable[j].begin(), subIndexTable[j].end() );
                addedIndexTable.push_back( newIndexEntry );
            }
        }
    } //end for loop

    // add the obtained partitions to the vector
    partitionSet->insert( partitionSet->end(), addedPartitionSet.begin(), addedPartitionSet.end() );

    return addedIndexTable;
} // end getBranching() -------------------------------------------------------------------------------------

std::vector< std::vector< unsigned int > > WHtree::getBranching( const std::vector < size_t > &thisPartition,
                                                                 size_t depthLevel,
                                                                 std::vector< std::vector < size_t > > *partitionSet )
{
    if( !partitionSet->empty() )
    {
        throw std::runtime_error( "ERROR @ WHtree::getBranching(): partition set wasnt empty" );
    }

    std::vector<nodeID_t> partFullId;
    std::vector< std::vector < nodeID_t > > partitionFullIdSet;

    partFullId.reserve( thisPartition.size() );
    for( size_t i = 0; i < thisPartition.size(); ++i )
    {
        partFullId.push_back( std::make_pair( true, thisPartition[i] ) );
    }

    std::vector< std::vector< unsigned int > > indexTable( getBranching( partFullId, depthLevel, &partitionFullIdSet, true ) );

    partitionSet->reserve( partitionFullIdSet.size() );
    for( size_t i = 0; i < partitionFullIdSet.size(); ++i )
    {
        std::vector<size_t> thisSet;
        thisSet.reserve( partitionFullIdSet[i].size() );
        for( size_t j = 0; j < partitionFullIdSet[i].size(); ++j )
        {
            if( partitionFullIdSet[i][j].first)
            {
                thisSet.push_back( partitionFullIdSet[i][j].second );
            }
            else
            {
                std::cerr << "WARNING @  WHtree::getBranching(), leaves were returned" << std::endl;
            }
        }
        partitionSet->push_back( thisSet );
    }

    return indexTable;
}

// === PRIVATE MEMBER FUNCTIONS ===



WHnode* WHtree::fetchNode( const size_t thisNode )
{
    if( thisNode >= m_nodes.size() )
    {
        return 0;
    }
    else
    {
        return &( m_nodes[thisNode] );
    }
} // end fetchNode() -------------------------------------------------------------------------------------
WHnode* WHtree::fetchNode( const nodeID_t &thisNode )
{
    if( thisNode.first )
    {
        return fetchNode( thisNode.second );
    }
    else
    {
        return fetchLeaf( thisNode.second );
    }
} // end fetchNode() -------------------------------------------------------------------------------------


WHnode* WHtree::fetchLeaf( const size_t thisLeaf )
{
    if( thisLeaf >= m_leaves.size() )
    {
        return 0;
    }
    else
    {
        return &( m_leaves[thisLeaf] );
    }
} // end fetchLeaf() -------------------------------------------------------------------------------------


WHnode* WHtree::fetchRoot()
{
    return fetchNode( getNumNodes()-1 );
} // end fetchRoot() -------------------------------------------------------------------------------------


std::pair<size_t, size_t> WHtree::cleanup( std::vector<size_t> *outLookup )
{
    //finish setting pruning flags on unnecessary nodes

    //reset nodes size and
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end(); ++nodesIter )
    {
        nodesIter->setSize( 0 );
        nodesIter->setHLevel( 0 );
    }
    //initialize base nodes size
    for( std::vector<WHnode>::const_iterator leavesIter( m_leaves.begin() ); leavesIter != m_leaves.end(); ++leavesIter )
    {
        WHnode *parentNode( fetchNode( leavesIter->getParent() ) );
        size_t currentSize( parentNode->getSize() );
        if( !( leavesIter->isFlagged() ) )
        {
            // leaf will not be pruned
            parentNode->setSize( currentSize+1 );
            parentNode->setHLevel( 1 );
        }
    }
    //initialize remaining nodes
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end()-1; ++nodesIter )
    {
        WHnode *parentNode( fetchNode( nodesIter->getParent() ) );
        size_t currentNodeSize( nodesIter->getSize() );
        size_t currentPapaSize( parentNode->getSize() );
        parentNode->setSize( currentPapaSize + currentNodeSize );
        if( currentNodeSize < 2 )
        {
            // if a node has less than 2 leaves it is unnecesary in a hierarchical tree, so in that case the prune bit must be set too
            nodesIter->setFlag( true );
            if( currentNodeSize > 0 )
            {
                parentNode->setHLevel( nodesIter->getHLevel() );
            }
        }
        else
        {
            parentNode->setHLevel( std::max( parentNode->getHLevel(), nodesIter->getHLevel()+1 ) );
        }
    }
    // loop through nodes and check that there are no hanging nodes ( nodes with one or no children) without structural function
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end(); ++nodesIter )
    {
        size_t numNewKids( 0 );
        std::vector<nodeID_t> kids( nodesIter->getChildren() );

        for( std::vector<nodeID_t>::iterator kidsIter( kids.begin() ); kidsIter != kids.end(); ++kidsIter )
        {
            if( fetchNode( *kidsIter )->isLeaf() )
            {
                if( !( getNode( *kidsIter ).isFlagged() ) )
                {
                    ++numNewKids;
                }
            }
            else
            {
                if( ( getNode( *kidsIter ).getHLevel() ) != 0 )
                {
                    ++numNewKids;
                }
            }
        }

        if( numNewKids <= 1 )
        {
            nodesIter->setFlag( true );
        }
    }

    // create new IDs lookup tables
    const size_t INVALID( getNumLeaves()+1 );
    std::vector<size_t> lookupLeafID( getNumLeaves(), INVALID ), lookupNodeID( getNumNodes(), INVALID ), lookupParentID( getNumNodes(), INVALID );
    size_t leafCounter( 0 ), nodeCounter( 0 ); // new ID counters
    // loop through leaves and create lookup table for new leaves IDs
    for( size_t i = 0; i < m_leaves.size(); ++i )
    {
        if( !m_leaves[i].isFlagged() )
        {
            lookupLeafID[i] = leafCounter++;
        }
    }
    // loop through nodes and create lookup table for new nodes IDs
    for( size_t i = 0; i < m_nodes.size(); ++i )
    {
        if( !m_nodes[i].isFlagged() )
        {
            lookupNodeID[i] = nodeCounter++;
        }
    }
    // loop through nodes and create lookup table for new nodes parent_IDs
    for( size_t i = 0; i < m_nodes.size(); ++i )
    {
        if( m_nodes[i].isFlagged() )
        {
            const WHnode* searchNode( &getNode( i ) );
            while( searchNode->isFlagged() )  // while parents are set to prune, keep going up the tree until a valid node is reached
            {
                if( searchNode->isRoot() )
                {
                    break; // this was the top node of the tree
                }
                searchNode = &getNode( searchNode->getParent() );
            }
            if( searchNode->isRoot() && searchNode->isFlagged() )
            {
                lookupParentID[i] = 0;
            }
            else
            {
                lookupParentID[i] = lookupNodeID[searchNode->getID()];
            }

            if( lookupParentID[i] == INVALID )
            {
                throw std::runtime_error( "ERROR @ WHtree::cleanup(): error filling new parent ID lookup table" );
            }
        }
        else
        {
            lookupParentID[i] = lookupNodeID[m_nodes[i].getID()];
        }
    }

    //  delete discarded elements from containers
    // eliminate discarded leaves from the leaves vector and the coordinates vector
    size_t discardedLeaves( 0 ), discardedNodes( 0 );
    std::vector<WHnode>::iterator leavesDelIter( m_leaves.begin() );
    std::vector<WHcoord>::iterator coordDelIter( m_coordinates.begin() );
    while( leavesDelIter != m_leaves.end() )
    {
        if( leavesDelIter->isFlagged() )
        {
            ++discardedLeaves;
            m_discarded.push_back( *coordDelIter );
            leavesDelIter = m_leaves.erase( leavesDelIter );
            coordDelIter = m_coordinates.erase( coordDelIter );
        }
        else
        {
            ++leavesDelIter;
            ++coordDelIter;
        }
    }
    // eliminate discarded nodes from the vector
    std::vector<WHnode>::iterator nodesDelIter( m_nodes.begin() );
    while( nodesDelIter != m_nodes.end() )
    {
        if( nodesDelIter->isFlagged() )
        {
            ++discardedNodes;
            nodesDelIter = m_nodes.erase( nodesDelIter );
        }
        else
        {
            ++nodesDelIter;
        }
    }

    // update names of non-discarded elements
    // loop through leaves and change leaves IDs and Parents IDs
    for( std::vector<WHnode>::iterator leavesIter( m_leaves.begin() ); leavesIter != m_leaves.end(); ++leavesIter )
    {
        size_t newID( lookupLeafID[leavesIter->getID()] );
        size_t newParentID( lookupParentID[leavesIter->getParent().second] );

        if( ( newID == INVALID ) || ( newParentID == INVALID ) )
        {
            std::cerr << "Discarded " << discardedLeaves<< " and " << discardedNodes<< " nodes" << std::endl;
            std::cerr << "Old ID: " << leavesIter->getID() << ". New ID: " << newID<< ". Old parent ID: " << leavesIter->getParent().second
                      << ". New parent ID: " << newParentID << std::endl;

            throw std::runtime_error( "ERROR @ WHtree::cleanup(): error updating leaf IDs, invalid lookup table value" );
        }

        leavesIter->setID( std::make_pair( false, newID ) );
        leavesIter->setParent( std::make_pair( true, newParentID ) );
    }
    // loop through nodes and change names of nodes IDs and parents
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end(); ++nodesIter )
    {
        size_t newID( lookupNodeID[nodesIter->getID()] );

        size_t newParentID( 0 );
        if( ( nodesIter+1 ) != m_nodes.end() )
        {
            newParentID = ( lookupParentID[nodesIter->getParent().second] );
        }

        std::vector<nodeID_t> emptyKids;

        if( ( newID == INVALID) || ( newParentID == INVALID ) )
        {
            throw std::runtime_error( "ERROR @ WHtree::cleanup(): error updating node IDs, invalid lookup table value" );
        }

        if( newParentID == 0 ) // if node is top of the tree parent must be set to 0-0
        {
            if( ( nodesIter+1 ) != m_nodes.end() ) // this should only happen at the last node
            {
                std::cerr << std::endl << "Node says its root: " << nodesIter->printAllData() << std::endl;
                std::cerr << "New ID: " << newID << ". New parent ID: " << newParentID << std::endl;
                std::cerr << "But last node is: " << ( m_nodes.end()-1 )->printAllData() << std::endl;
                std::cerr << "New ID: " << lookupNodeID[( m_nodes.end()-1 )->getID()] << std::endl;
                throw std::runtime_error( "ERROR @ WHtree::cleanup(): pruning failed, top of tree is not last node in vector" );
            }
            nodesIter->setParent( std::make_pair( false, 0 ) );
        }
        else
        {
            nodesIter->setParent( std::make_pair( true, newParentID ) );
        }

        nodesIter->setID( std::make_pair( true, newID ) );
        nodesIter->setChildren( emptyKids );
        nodesIter->setHLevel( 0 );
    }

    // fill up children and hierarchical level data
    for( std::vector<WHnode>::iterator leavesIter( m_leaves.begin() ); leavesIter != m_leaves.end(); ++leavesIter )
    {
        WHnode *parentNode( fetchNode( leavesIter->getParent() ) );
        std::vector<nodeID_t> currentKids( parentNode->getChildren() );
        currentKids.push_back( leavesIter->getFullID() );
        parentNode->setChildren( currentKids );
        parentNode->setHLevel( 1 );
    }
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end(); ++nodesIter )
    {
        if( nodesIter->isRoot() )
        {
            continue;
        }
        WHnode *parentNode( fetchNode( nodesIter->getParent() ) );
        std::vector<nodeID_t> currentKids( parentNode->getChildren() );
        currentKids.push_back( nodesIter->getFullID() );
        parentNode->setChildren( currentKids );
        parentNode->setHLevel( std::max( parentNode->getHLevel(), ( nodesIter->getHLevel() )+1 ) );
    }

    // check if final tree is consistent
    if( !check() )
    {
        throw std::runtime_error( "ERROR @ WHtree::cleanup(): resulting tree is not consistent" );
    }
    if( discardedLeaves != 0 || discardedNodes != 0 )
    {
        m_cpcc = 0;
        clearPartitions();
    }

    if( outLookup != 0 ) // if a pointer was passed as argument
    {
        *outLookup = lookupNodeID;
    }

    clearPartitions();
    clearPartColors();
    return std::make_pair( discardedLeaves, discardedNodes );
} // end cleanup() -------------------------------------------------------------------------------------


size_t WHtree::debinarize( bool keepBaseNodes )
{
    if( keepBaseNodes && !testRootBaseNodes() )
    {
        std::cerr << "WARNING@ Debinarize: base nodes have mixed nodes and leaves, debinarize will be standard " << std::endl;
        keepBaseNodes = false;
    }
    const size_t origNumNodes( getNumNodes() );

    std::vector<bool> validNode( getNumNodes(), true );
    std::vector<std::vector<nodeID_t> > realChildren;
    realChildren.resize( getNumNodes() );
    std::vector<size_t> realParentsforLeaves( getNumLeaves(), 0 );
    std::vector<size_t> realParentsforNodes( getNumNodes(), 0 );


    if( !keepBaseNodes )
    {
        // first loop through leaves
        for( unsigned int id = 0; id < getNumLeaves(); ++id )
        {
            size_t currentNode( fetchLeaf( id )->getParent().second );
            dist_t currentDist( fetchNode( currentNode )->getDistLevel() );

            if( fetchNode( currentNode )->isRoot() ) // its parent is the last node
            {
                realParentsforLeaves[id] = currentNode;
                realChildren[currentNode].push_back( std::make_pair( false, id ) );
                continue;
            }
            size_t nextParent( fetchNode( currentNode )->getParent().second );
            dist_t nextDist( fetchNode( nextParent )->getDistLevel() );

            while( currentDist == nextDist )
            {
                validNode[currentNode] = false;
                currentNode = nextParent;
                currentDist = nextDist;

                if( fetchNode( nextParent )->isRoot() ) // if we reach the last node we stop
                {
                    break;
                }

                nextParent = fetchNode( currentNode )->getParent().second;
                nextDist = fetchNode( nextParent )->getDistLevel();
            }
            realParentsforLeaves[id] = currentNode;
            realChildren[currentNode].push_back( std::make_pair( false, id ) );
        }
    }
    else
    {
        for( unsigned int id = 0; id < getNumLeaves(); ++id )
        {
            size_t currentNode( fetchLeaf( id )->getParent().second );

            if( fetchNode( currentNode )->isRoot() ) // its parent is the last node
            {
                realParentsforLeaves[id] = currentNode;
                realChildren[currentNode].push_back( std::make_pair( false, id ) );
                continue;
            }

            realParentsforLeaves[id] = currentNode;
            realChildren[currentNode].push_back( std::make_pair( false, id ) );
        }
    }

    // then loop through nodes
    for( unsigned int id = 0; id < getNumNodes()-1; ++id ) //we never process the last node (has no parent)
    {
        size_t currentNode( fetchNode( id )->getParent().second );
        dist_t currentDist( fetchNode( currentNode )->getDistLevel() );

        if( fetchNode( currentNode )->isRoot() ) // its parent is the last node
        {
            realParentsforNodes[id] = currentNode;
            if( validNode[id] )
            {
                realChildren[currentNode].push_back( std::make_pair( true, id ) );
            }
            continue;
        }
        size_t nextParent( fetchNode( currentNode )->getParent().second );
        dist_t nextDist( fetchNode( nextParent )->getDistLevel() );

        while( currentDist == nextDist )
        {
            validNode[currentNode] = false;

            currentNode = nextParent;
            currentDist = nextDist;

            if( fetchNode( nextParent )->isRoot() ) // its parent is the last node
            {
                break;
            }

            nextParent = fetchNode( currentNode )->getParent().second;
            nextDist = fetchNode( nextParent )->getDistLevel();
        }
        realParentsforNodes[id] = currentNode;
        if( validNode[id] )
            realChildren[currentNode].push_back( std::make_pair( true, id ) );
    }
    realParentsforNodes[getNumNodes()-1 ]= 0; // parent of root node


    // check validity
    for( unsigned int i = 0; i < getNumNodes(); ++i ) //we never process the last node (has no parent)
    {
            if( validNode[i] && ( realChildren[i].empty() ) )
            {
                std::cerr << "node (1-" << i << ") has no real children" << std::endl;
                throw std::runtime_error( "ERROR @ WHtree::debinarize(): node has no real children" );
            }
    }

    // lookuptable to account for id change due to invalid nodes
    size_t nbCount( 0 );
    const size_t invalid( validNode.size()+1 );
    std::vector<size_t> changeLookup( validNode.size(), invalid );
    for( size_t i = 0; i < validNode.size(); ++i )
    {
        if( validNode[i] )
        {
            changeLookup[i] = nbCount;
            ++nbCount;
        }
    }

    //rename real children
    for( size_t i = 0; i < realChildren.size(); ++i )
    {
        for( size_t j = 0; j < realChildren[i].size(); ++j )
        {
            if( realChildren[i][j].first )
            {
                size_t newName( changeLookup[realChildren[i][j].second] );
                if( newName == invalid )
                {
                    throw std::runtime_error( "ERROR @ WHtree::debinarize(): error renaming node children" );
                }
                realChildren[i][j].second = newName;
            }
        }
    }

    // change leaves info
    for( unsigned int id = 0; id < getNumLeaves(); ++id )
    {
        size_t realDad( changeLookup[realParentsforLeaves[id]] );
        if( realDad == invalid )
        {
            throw std::runtime_error( "ERROR @ WHtree::debinarize(): error renaming nb leaf parents" );
        }
        fetchLeaf( id )->setParent( std::make_pair( true, realDad ) );
    }

    {
        // create new nodes vector
        std::vector<WHnode> nbNodes;
        nbNodes.reserve( getNumNodes() );
        for( unsigned int id = 0; id < getNumNodes()-1; ++id )
        {
            if( validNode[id] )
            {
                size_t realDad( changeLookup[realParentsforNodes[id]] );
                if( realDad == invalid )
                {
                    std::cerr << "node (1-" << id << ") is valid but has invalid dad, ( preDad was 1-" << realParentsforNodes[id] << ")" << std::endl;
                    throw std::runtime_error( "ERROR @ WHtree::debinarize(): error renaming nb node parents" );
                }
                WHnode thisnode( std::make_pair( true, changeLookup[id] ), realChildren[id],
                                 fetchNode( id )->getSize(), fetchNode( id )->getDistLevel(), 0 );
                thisnode.setParent( std::make_pair( true, realDad ) );
                nbNodes.push_back( thisnode );
            }
        }
        WHnode rootnode( std::make_pair( true, changeLookup[getNumNodes()-1] ), realChildren[getNumNodes()-1],
                         fetchNode( getNumNodes()-1 )->getSize(), fetchNode( getNumNodes()-1 )->getDistLevel(), 0 );
        nbNodes.push_back( rootnode );


        // replace old node vector with new one
        m_nodes = nbNodes;
    }
    clearContainedLeaves();
    clearPartitions();
    // refill hLevel data
    for( unsigned int id = 0; id < getNumNodes(); ++id )
    {
        WHnode* currentNode( fetchNode( id ) );
        std::vector<nodeID_t> currentKids( currentNode->getChildren() );
        size_t currentHLevel( 0 );
        for( size_t j = 0; j < currentKids.size(); ++j )
        {
            currentHLevel = std::max( currentHLevel, fetchNode( currentKids[j] )->getHLevel()+1 );
        }
        currentNode->setHLevel( currentHLevel );
    }

    if( !check() )
    {
        throw std::runtime_error( "ERROR @ WHtree::debinarize(): resutling tree is not consistent" );
    }

    size_t discardedNodes( origNumNodes - getNumNodes() );

    if( discardedNodes != 0 )
    {
        m_cpcc = 0;
        clearPartitions();
    }

    return discardedNodes;
} // end debinarize() -------------------------------------------------------------------------------------

void WHtree::forceMonotonicityUp()
{
    // loop through nodes and force monotonicity
    for( std::vector<WHnode>::iterator nodesIter( m_nodes.begin() ); nodesIter != m_nodes.end(); ++nodesIter )
    {
        WHnode* currentNode( fetchNode( nodesIter->getID() ) );
        std::vector<nodeID_t> currentKids( currentNode->getChildren() );

        for( size_t j = 0; j < currentKids.size(); ++j )
        {
            WHnode* kid( fetchNode( currentKids[j] ) );
            // if its kid has a higher level than the node, there was a non-monotonic step, and the highest level is kept
            if( kid->getDistLevel() > currentNode->getDistLevel() )
            {
                currentNode->setDistLevel( kid->getDistLevel() );
            }
        }
    }
    return;
} // end forceMonotonicityUp() -------------------------------------------------------------------------------------

void WHtree::forceMonotonicityDown()
{
    // loop through nodes and force monotonicity
    for( std::vector<WHnode>::reverse_iterator nodesIter( m_nodes.rbegin() ); nodesIter != m_nodes.rend(); ++nodesIter )
    {
        WHnode* currentNode( fetchNode( nodesIter->getID() ) );
        std::vector<nodeID_t> currentKids( currentNode->getChildren() );

        for( size_t j = 0; j < currentKids.size(); ++j )
        {
            WHnode* kid( fetchNode( currentKids[j] ) );
            // if its kid has a higher level than the node, there was a non-monotonic step, and the lowest level is kept
            if( kid->getDistLevel() > currentNode->getDistLevel() )
            {
                kid->setDistLevel( currentNode->getDistLevel() );
            }
        }
    }
    return;
} // end forceMonotonicityDown() -------------------------------------------------------------------------------------

void WHtree::forceMonotonicity( double errorMult )
{
    if( errorMult > 100 || errorMult <= 0 )
    {
        errorMult = 1;
    }
    double errorMargin( errorMult * 1E-5 );

    for( int64_t i = m_nodes.size()-1; i >= 0; )
    {
        WHnode* currentNode( fetchNode( m_nodes[i].getID() ) );
        size_t currentSize( currentNode->getSize() );
        dist_t currentLevel( currentNode->getDistLevel() );

        std::vector<nodeID_t> currentKids( currentNode->getChildren() );

        double newLevelSum( 0 );
        size_t remainingSize( currentSize );
        bool doCorrect( false );

        for( size_t j = 0; j < currentKids.size(); ++j )
        {
            WHnode* kid( fetchNode( currentKids[j] ) );
            // if its kid has a higher level than the node, there was a non-monotonic step
            if( kid->getDistLevel() > currentLevel + errorMargin )
            {
                newLevelSum += kid->getDistLevel() * kid->getSize();
                remainingSize -= kid->getSize();
                doCorrect = true;
            }
        }
        if( doCorrect )
        {
//            std::cout << i << " -> " << std::flush;

            double correctedLevel( ( newLevelSum + ( remainingSize * currentLevel )  ) / currentSize );
            // correct level of current node
            currentNode->setDistLevel( correctedLevel );
            // correct level of non-monotonic children nodes
            for( size_t j = 0; j < currentKids.size(); ++j )
            {
                WHnode* kid( fetchNode( currentKids[j] ) );
                // if its kid has a higher level than the node, there was a non-monotonic step, and the highest level must be kept to avoid info loss
                if( kid->getDistLevel() > correctedLevel )
                {
                    kid->setDistLevel( correctedLevel );
                }
            }

            // if node is not root check if corrected level is higher than parent node
            if( currentNode->isRoot() )
            {
                --i;
            }
            else
            {
                WHnode* parentNode( fetchNode( currentNode->getParent() ) );
                double parentLevel( parentNode->getDistLevel() );
                // if new level is higher than parent level, go back to parent node id and continue from there
                if( correctedLevel > parentLevel + errorMargin )
                {
                    i = parentNode->getID();
                }
                else
                {
                    --i;
                }
            }
//DEBUG            std::cout << i << " : " <<currentLevel<< " -> " << correctedLevel << std::endl;
        }
        else
        {
            --i;
        }
    }

    forceMonotonicityDown();

    return;
} // end forceMonotonicity() -------------------------------------------------------------------------------------



