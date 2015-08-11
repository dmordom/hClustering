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


#include <map>

#include "protoNode.h"

bool protoNode::updateNbhood( const nodeID_t& oldNode1, const nodeID_t& oldNode2, const nodeID_t& newNode, const dist_t newDist )
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

bool protoNode::updateActivhood( const nodeID_t& oldNode1, const nodeID_t& oldNode2, const nodeID_t& newNode, const dist_t newDist,
                                 const bool isActive, const std::vector< protoNode >& protoNodes )
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

void protoNode::updateDist( const nodeID_t& updatedNode, const dist_t updatedDist )
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
