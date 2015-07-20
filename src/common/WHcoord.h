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

#ifndef WHCOORD_H
#define WHCOORD_H

// std library
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

// boost library
#include <boost/format.hpp>

// typedefs
typedef int16_t coord_t;
typedef enum
{
    HC_VISTA, HC_NIFTI
} HC_GRID;

/**
 * Class to contain a seed voxel coordinate. consisting on x,y,z also implements coordinate gird changes and physical neighbor search
 */
class WHcoord
{
public:
    //! Constructor
    WHcoord();

    //! Constructor
    //! \param x_init x coordinate initializer
    //! \param y_init y coordinate initializer
    //! \param z_init z coordinate initializer
    WHcoord( coord_t x_init, coord_t y_init, coord_t z_init );

    //! Constructor
    //! \param object coordinate initializer
    WHcoord( const WHcoord &object );

    //! Destructor
    ~WHcoord();

    //! = member operator
    //! \param rhs object to copy
    //! \return reference to copied object
    WHcoord& operator =( const WHcoord &rhs );

    // === PUBLIC MEMBER FUNCTIONS ===


    //! returns the euclidean distance between this voxel and the input voxel
    //! \param voxel coordinate to calculate distance to
    //! \return euclidean distance value
    float getPhysDist( const WHcoord voxel ) const;

    //! returns a vector containing the coordinates of the physical neighbours adjacent to the voxel, the neighbourhood level is defined by nb_level
    //! \param dataSize dataset limits
    //! \param nbLevel neighborhood value
    //! \return coordinate vector
    std::vector< WHcoord > getPhysNbs( const WHcoord dataSize, const unsigned int nbLevel ) const;

    //! returns a string with the coordinates of the voxel in the form "xxx_yyy_zzz"
    //! \return output string
    std::string getNameString() const;

    //! transform coordinates to vista format
    //! \param dataSize dataset size
    //! \return converted coordinate
    WHcoord nifti2vista( const WHcoord dataSize ) const;

    //! transform coordinates to vista format
    //! \param dataSize dataset size
    //! \return converted coordinate
    WHcoord vista2nifti( const WHcoord dataSize ) const;

    // === PUBLIC DATA MEMBERS ===

    coord_t m_x; //!< x coordinate
    coord_t m_y; //!< y coordinate
    coord_t m_z; //!< z coordinate

protected:
private:
};

// === NON-MEMBER OPERATORS ===


//! << operator for the coordinate class
//! \param os output stream
//! \param object coordinate to print out
//! \return ostream with the coordinate information
std::ostream& operator <<( std::ostream& os, const WHcoord& object );

//! < operator for the coordinate class
//! \param lhs left coordinate operator
//! \param rhs right coordinate operator
//! \return result of the inequality
bool operator <( const WHcoord &lhs, const WHcoord &rhs );

//! == operator for the coordinate class
//! \param lhs left coordinate operator
//! \param rhs right coordinate operator
//! \return result of the equality
bool operator ==( const WHcoord &lhs, const WHcoord &rhs );

//! != operator for the coordinate class
//! \param lhs left coordinate operator
//! \param rhs right coordinate operator
//! \return result of the inequality
bool operator !=( const WHcoord &lhs, const WHcoord &rhs );

//! helper function to print a string with the coordinate grid name
//! \param gridType type of coordinate grid to get the string from
//! \return identifier string
std::string getGridString( const HC_GRID gridType );
#endif  // WHCOORD_H
