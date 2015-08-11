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
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "WHcoord.h"

WHcoord::WHcoord(): m_x( 0 ), m_y( 0 ), m_z( 0 )
{
}

WHcoord::WHcoord( coord_t x_init, coord_t y_init, coord_t z_init ): m_x( x_init ), m_y( y_init ), m_z( z_init )
{
}

WHcoord::WHcoord( const WHcoord &object ): m_x( object.m_x ), m_y( object.m_y ), m_z( object.m_z )
{
}

WHcoord::~WHcoord()
{
    //Cleanup
}


// === PUBLIC MEMBER FUNCTIONS ===


// "phys_dist()": returns the euclidean distance between this voxel and the input voxel
float WHcoord::getPhysDist( const WHcoord voxel ) const
{
    float xdif( m_x-voxel.m_x ), ydif( m_y-voxel.m_y ), zdif( m_z-voxel.m_z );
    return sqrt( ( xdif*xdif )+( ydif*ydif )+( zdif*zdif ) );
} // end "phys_dist()" -----------------------------------------------------------------


// "get_phys_neighbours()": returns a vector containing the coordinates of the physical neighbours adjacent to the voxel,
// the neighbourhood level is defined by nb_level
std::vector<WHcoord> WHcoord::getPhysNbs( const WHcoord dataSize, const unsigned int nbLevel ) const
{
    std::vector<WHcoord> physNeighbours;

    int condition( 0 );
    int range( 0 );

    switch( nbLevel )
    {
    case 6:
        range = 1;
        condition = 1;
        break;
    case 18:
        range = 1;
        condition = 2;
        break;
    case 26:
    case 32:
        range = 1;
        condition = 3;
        break;
    case 56:
        range = 2;
        condition = 3;
        break;
    case 92:
        range = 2;
        condition = 4;
        break;
    case 116:
        range = 2;
        condition = 5;
        break;
    case 124:
        range = 2;
        condition = 6;
        break;
    default:
        std::cerr << "ERROR @ coordinate::getPhysNbs(): Unrecognized nbhood level: " << nbLevel << std::endl;
        return physNeighbours;
        break;
    }

    for( coord_t i = -range; i <= range ; ++i )
    {
        for( coord_t j = -range; j <= range ; ++j )
        {
            for( coord_t k = -range; k <= range ; ++k )
            {
                const int test( abs( i )+abs( j )+abs( k ) );
                if( test == 0 )
                {
                    continue; // dont use the voxel itself as neighbour
                }
                else if( test > condition )
                {
                    continue;
                }
                else
                {
                    WHcoord temp;
                    temp.m_z = this->m_z + i;
                    temp.m_y = this->m_y + j;
                    temp.m_x = this->m_x + k;
                    if( ( temp.m_z < 0 ) || ( temp.m_y < 0 ) || ( temp.m_x < 0 ) || ( temp.m_z > dataSize.m_z ) ||
                         ( temp.m_y > dataSize.m_y ) ||( temp.m_x > dataSize.m_x ) )
                    {
                        continue; // voxel is out of mask boundaries
                    }
                    physNeighbours.push_back( temp );
                }
            }
        }
    }

    // extra nbs for nb level 32
    if( nbLevel == 32 )
    {
        for( coord_t i =- 2; i <= 2; i += 2 )
        {
            for( coord_t j =- 2; j <= 2; j += 2 )
            {
                for( coord_t k =- 2; k <= 2; k += 2 )
                {
                    const int test( abs( i )+abs( j )+abs( k ) );
                    if( test != 2 )
                    {
                        continue;
                    }
                    else
                    {
                        WHcoord temp;
                        temp.m_z = this->m_z + i;
                        temp.m_y = this->m_y + j;
                        temp.m_x = this->m_x + k;
                        if( ( temp.m_z < 0 ) || ( temp.m_y < 0 ) || ( temp.m_x < 0 ) || ( temp.m_z > dataSize.m_z ) ||
                             ( temp.m_y > dataSize.m_y ) || ( temp.m_x > dataSize.m_x ) )
                        {
                            continue; // voxel is out of mask boundaries
                        }
                        physNeighbours.push_back( temp );
                    }
                }
            }
        }
    }

    std::sort( physNeighbours.begin(), physNeighbours.end() );
    return physNeighbours;
} // end "phys_dist()" -----------------------------------------------------------------


// "getNameString()": returns a string with the coordinates of the voxel in the form "xxx_yyy_zzz"
std::string WHcoord::getNameString() const
{
    std::stringstream stream;
    std::string name;
    stream << m_x << "_" << m_y << "_" << m_z;
    stream >> name;
    return name;
} // end "getNameString()" -----------------------------------------------------------------


WHcoord WHcoord::nifti2vista( const WHcoord dataSize ) const
{
        WHcoord returnVista;
        returnVista.m_x = this->m_x;
        returnVista.m_y = ( dataSize.m_y-1 )-this->m_y;
        returnVista.m_z = ( dataSize.m_z-1 )-this->m_z;
        return returnVista;
} // end "nifti2Vista()" -----------------------------------------------------------------



WHcoord WHcoord::vista2nifti( const WHcoord dataSize ) const
{
        WHcoord returnNifti;
        returnNifti.m_x = this->m_x;
        returnNifti.m_y = ( dataSize.m_y-1 )-this->m_y;
        returnNifti.m_z = ( dataSize.m_z-1 )-this->m_z;
        return returnNifti;
} // end "vista2Nifti()" -----------------------------------------------------------------

WHcoord WHcoord::surf2vista( const WHcoord dataSize ) const
{
        WHcoord returnVista;
        returnVista.m_x = this->m_x+( ( dataSize.m_x-1.0 )/2.0 );
        returnVista.m_y = ( ( dataSize.m_y-1.0 )/2.0 )-this->m_y;
        returnVista.m_z = ( ( dataSize.m_z-1.0 )/2.0 )-this->m_z;
        return returnVista;
} // end "surf2vista()" -----------------------------------------------------------------

WHcoord WHcoord::surf2nifti( const WHcoord dataSize ) const
{
        WHcoord returnVista = this->surf2vista( dataSize );
        return returnVista.vista2nifti( dataSize );
} // end "surf2nifti()" -----------------------------------------------------------------


// === MEMBER OPERATORS ===



WHcoord& WHcoord::operator =( const WHcoord &rhs )
{
    m_x = rhs.m_x;
    m_y = rhs.m_y;
    m_z = rhs.m_z;
    return *this;
} // end "operator =" -----------------------------------------------------------------



// === NON-MEMBER OPERATORS ===



std::ostream& operator <<( std::ostream& os, const WHcoord& object )
{
    os << boost::format( "%03d %03d %03d" ) % object.m_x  % object.m_y  % object.m_z;
    return os;
} // end "operator << " -----------------------------------------------------------------


bool operator <( const WHcoord &lhs, const WHcoord &rhs )
{
    if( lhs.m_z < rhs.m_z ) // z1<z2
    {
        return true;
    }
    else if( lhs.m_z > rhs.m_z ) // z1>z2
    {
        return false;
    }
    else if( lhs.m_y < rhs.m_y ) // z1=z2, y1<y2
    {
        return true;
    }
    else if( lhs.m_y > rhs.m_y ) // z1=z2, y1>y2
    {
        return false;
    }
    else if( lhs.m_x < rhs.m_x ) // z1=z2, y1=y2, x1<x2
    {
        return true;
    }
    else if( lhs.m_x > rhs.m_x ) // z1=z2, y1=y2, x1>x2
    {
        return false;
    }
    else  // coord1 = coord2
    {
        return false;
    }
} // end "operator <" -----------------------------------------------------------------

bool operator ==( const WHcoord &lhs, const WHcoord &rhs )
{
    if( lhs.m_z == rhs.m_z )
    {
        if( lhs.m_y == rhs.m_y )
        {
            if( lhs.m_x == rhs.m_x )
            {
                return true;
            }
        }
    }

    return false;
} // end "operator ==" -----------------------------------------------------------------

bool operator !=( const WHcoord &lhs, const WHcoord &rhs )
{
    return !( lhs == rhs );
} // end "operator !=" -----------------------------------------------------------------

std::string getGridString( const HC_GRID gridType )
{
    if( gridType == HC_VISTA )
    {
        return "vista";
    }
    else if( gridType == HC_NIFTI )
    {
        return "nifti";
    }
    else if( gridType == HC_SURF )
    {
        return "surf";
    }
    else
    {
        std::string emptyString;
        return emptyString;
    }
}
