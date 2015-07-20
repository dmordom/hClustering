#include <map>

#include "protoNode.h"

bool protoNode::updateNbhood( nodeID_t oldNode1, nodeID_t oldNode2, nodeID_t newNode, dist_t newDist )
{
    // update nbhood table
    m_nbNodes.erase( oldNode1 );
    m_nbNodes.erase( oldNode2 );
    m_nbNodes.insert( std::make_pair( newNode, newDist ) );
    // update nearest neighbour
    if( ( m_nearNb.first == oldNode1 ) || ( m_nearNb.first == oldNode2 ) )
    { //if one of the deleted nb was the nearest, check all of them again
        m_nearNb = *( m_nbNodes.begin() );
        for( std::map< nodeID_t, dist_t >::const_iterator iter = m_nbNodes.begin(); iter != m_nbNodes.end(); ++iter )
            if( iter->second < m_nearNb.second )
                m_nearNb = *iter;
        return true;
    }
    else
    { // if deleted nb was not the nearest, simply check if distance to the new one is smaller
        if( newDist < m_nearNb.second )
        {
            m_nearNb.first = newNode;
            m_nearNb.second = newDist;
            return true;
        }
        else
        {
            return false; // nearest nb has not changed
        }
    }
    return false;
} // end "updateNbhood()" -----------------------------------------------------------------

bool protoNode::updateActivhood( const nodeID_t oldNode1, const nodeID_t oldNode2, const nodeID_t newNode, const dist_t newDist,
                                 const bool isActive, const std::vector< protoNode > &protoNodes )
{
    bool changed( false );
    // update nbhood table
    m_nbNodes.erase( oldNode1 );
    m_nbNodes.erase( oldNode2 );
    //if one of the deleted nb was the nearest, check all of them again
    if( ( m_nearNb.first == oldNode1 ) || ( m_nearNb.first == oldNode2 ) )
    {
        m_nearNb.first = noNbID;
        m_nearNb.second = noNbDist;
        updateActive( protoNodes );
        changed = true;
    }
    m_nbNodes.insert( std::make_pair( newNode, newDist ) );
    // if deleted nb was not the nearest, simply check if distance to the new one is smaller
    if( isActive && ( newDist < m_nearNb.second ) )
    {
        m_nearNb.first = newNode;
        m_nearNb.second = newDist;
        changed = true;
    }

    return changed;
} // end "updateActivhood()" -----------------------------------------------------------------

void protoNode::updateDist( const nodeID_t updatedNode, const dist_t updatedDist )
{
    m_nbNodes[updatedNode]=updatedDist;
} // end "updateDist()" -----------------------------------------------------------------

bool protoNode::updateActive( const std::vector< protoNode > &protoNodes )
{
    bool changed( false );
    if( m_nearNb.first.first && ( m_nearNb.second != noNbDist ) )
    {
        if( !protoNodes[m_nearNb.first.second].isActive() )
        {
            m_nearNb = std::make_pair( noNbID, noNbDist );
            changed = true;
        }
    }
    for( std::map< nodeID_t, dist_t >::const_iterator iter = m_nbNodes.begin(); iter != m_nbNodes.end(); ++iter )
    {
        bool isNode( iter->first.first );
        size_t thisNodeID( iter->first.second );
        bool isactive( false );

        if ( isNode )
        {
            isactive = protoNodes[thisNodeID].isActive();
        }
        else
        {
            isactive = true;
        }

        if( isactive && ( iter->second < m_nearNb.second ) )
        {
            m_nearNb = *iter;
            changed = true;
        }
    }
    return changed;
} // end "updateActivehood()" -----------------------------------------------------------------




// non-member operators

std::ostream& operator <<( std::ostream& os, const protoNode& object )
{
    os << "Near Nb: " << object.m_nearNb.first.first << "-" << object.m_nearNb.first.second << "|" << object.m_nearNb.second
                    << std::flush;
    os << ". Nbs: " << std::flush;
    for( std::map< nodeID_t, dist_t >::const_iterator iter( object.m_nbNodes.begin() ); iter != object.m_nbNodes.end(); ++iter )
    {
        os << "(" << iter->first.first << "-" << iter->first.second << "|" << iter->second << ") " << std::flush;
    }
    return os;
} // end operator << -----------------------------------------------------------------
