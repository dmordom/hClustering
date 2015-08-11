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
#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include "WStringUtils.h"

#include "WHnode.h"

// === PUBLIC MEMBER FUNCTIONS ===

WHnode::WHnode( nodeID_t idInit ): m_fullID( idInit ), m_parent( std::make_pair( false, 0 ) ),
                                   m_nodeSize( 1 ), m_distanceLevel( 0 ), m_hLevel( 0 ), m_flag( false )
{
}

WHnode::WHnode( nodeID_t idInit, std::vector<nodeID_t>childrenInit, size_t nodeSizeInit, dist_t distanceLevelInit, size_t hLevelInit ):
    m_fullID( idInit ), m_parent( std::make_pair( false, 0 ) ), m_children( childrenInit ), m_nodeSize( nodeSizeInit ),
    m_distanceLevel( distanceLevelInit ), m_hLevel( hLevelInit ), m_flag( false )
{
}

WHnode::~WHnode()
{
    // Cleanup
}

bool WHnode::isRoot() const
{
    return ( ( m_fullID.first )&&( !m_parent.first )&&( !m_parent.second ) );
}

void WHnode::setID( const nodeID_t newID )
{
    m_fullID = newID;
}

void WHnode::setParent( const nodeID_t newDad )
{
    m_parent = newDad;
}

void WHnode::setSize( const size_t newSize )
{
    m_nodeSize = newSize;
}

void WHnode::setHLevel( const size_t newLevel )
{
    m_hLevel = newLevel;
}

void WHnode::setDistLevel( const dist_t newLevel )
{
    m_distanceLevel = newLevel;
}

void WHnode::setChildren( std::vector<nodeID_t> newKids )
{
    m_children = newKids;
}

void WHnode::setFlag( const bool newFlag )
{
    m_flag = newFlag;
}

std::string WHnode::printAllData() const
{
    std::string oString;
    oString += ( "ID: " + string_utils::toString( m_fullID.first )+ "-" + str( boost::format( "%06d" ) % m_fullID.second ) );
    oString += ".  Dad: " + string_utils::toString( m_parent.first ) + "-" + str( boost::format( "%06d" ) % m_parent.second );
    oString += ".  Size: " + str( boost::format( "%06d" ) % m_nodeSize );
    oString += ".  HLevel: " + string_utils::toString( m_hLevel );
    oString += ".  DLevel: " + string_utils::toString( m_distanceLevel );
    oString += ".  Kids: (";
    for( size_t i = 0; i < m_children.size(); ++i )
    {
        oString += " " + string_utils::toString( m_children[i].first ) + "-" + str( boost::format( "%06d" ) % m_children[i].second )+ " ";
        if( i < ( m_children.size()-1 ) )
        {
            oString += ",";
        }
    }
    oString += ")";
    if( m_flag )
    {
        oString+=" F";
    }
    return oString;
} // end hNode::printAllData() -----------------------------------------------------------------

std::string WHnode::printJointData() const
{
    std::string oString;
    oString += string_utils::toString( m_distanceLevel );
    for( size_t i = 0; i < m_children.size(); ++i )
    {
        oString += " " + string_utils::toString( m_children[i].first ) + " " + str( boost::format( "%06d" ) % m_children[i].second );
    }
    return oString;
} // end hNode::printJointData() -----------------------------------------------------------------



// === NON-MEMBER OPERATORS ===



std::ostream& operator <<( std::ostream& os, const WHnode& object )
{
    os << object.printJointData();
    return os;
} // end operator << -----------------------------------------------------------------
