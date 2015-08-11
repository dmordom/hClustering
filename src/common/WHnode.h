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

#ifndef WHNODE_H
#define WHNODE_H

// std library
#include <vector>
#include <string>
#include <utility>
#include <iostream>

// boost library
#include <boost/format.hpp>

// typedefs
typedef float dist_t;
typedef std::pair<bool, size_t> nodeID_t;

/**
 * this class implements a hierarchical tree node with several relevant attributes
 */
class WHnode
{
public:
    /**
     * Constructor
     * \param idInit node id initializer
     */
    explicit WHnode( nodeID_t idInit );

    /**
     * Constructor
     * \param idInit node id initializer
     * \param childrenInit children id vector initializer
     * \param nodeSizeInit node size initializer
     * \param distanceLevelInit node distance level initializer
     * \param hLevelInit node hierarchical level initializer
     */
    WHnode( nodeID_t idInit, std::vector<nodeID_t>childrenInit, size_t nodeSizeInit, dist_t distanceLevelInit, size_t hLevelInit );

    //! Destructor
    ~WHnode();

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * returns true if object is a node
     * \return node indicator
     */
    bool isNode() const;

    /**
     * returns true if object is a leaf
     * \return leaf indicator
     */
    bool isLeaf() const;

    /**
     * returns true if object is the root of the tree
     * \return root indicator
     */
    bool isRoot() const;

    /**
     * returns true if object is flagged for pruning
     * \return flag indicator
     */
    bool isFlagged() const;

    /**
     * returns node/leaf ID
     * \return id indicator
     */
    size_t getID()   const;

    /**
     * returns full ID
     * \return full id indicator
     */
    nodeID_t getFullID()     const;

    /**
     * returns parent ID
     * \return parent id indicator
     */
    nodeID_t getParent()     const;

    /**
     * returns number of elements contained by the node
     * \return size
     */
    size_t getSize()       const;

    /**
     * returns distance level at which the element was formed
     * \return distance
     */
    dist_t getDistLevel()  const;

    /**
     * returns maximum number of nodes between the node and a leaf element
     * \return hierarchical level
     */
    size_t getHLevel()     const;

    /**
     * returns a vector with the ids of the children of that node
     * \return children vector
     */
    std::vector<nodeID_t> getChildren() const;

    /**
     * sets the ID field to the specified value
     * \param newID new ID
     */
    void setID( const nodeID_t newID );

    /**
     * sets the parent id field to the specified value
     * \param newDad new parent id
     */
    void setParent( const nodeID_t newDad );

    /**
     * sets the size field to the specified value
     * \param newSize new size
     */
    void setSize( const size_t newSize );

    /**
     * sets the hierarchical level field to the specified value
     * \param newLevel new hierarchical level
     */
    void setHLevel( const size_t newLevel );

    /**
     * sets the distance level field to the specified value
     * \param newLevel new distance level
     */
    void setDistLevel( const dist_t newLevel );

    /**
     * sets the children vector to the input values
     * \param newKids vector with the new children ids
     */
    void setChildren( std::vector<nodeID_t> newKids );

    /**
     * sets the prune flag to the indicated value
     * \param newFlag new value for the prune flag
     */
    void setFlag( const bool newFlag );

    /**
     * prints out a string with the node data, in order to implement the operator <<
     * \return a string with the node information
     */
    std::string printAllData() const;

    /**
     * prints out a string with the node join data, (reduced version of the node data)
     * \return a string with the node join information
     */
    std::string printJointData() const;

private:
    // === PRIVATE DATA MEMBERS ===

    nodeID_t                m_fullID; //!< node ID
    nodeID_t                m_parent; //!< parent node ID
    std::vector<nodeID_t>   m_children; //!< chlidren nodes IDs
    size_t                  m_nodeSize; //!< number of leaves contained on the node
    dist_t                  m_distanceLevel; //!< distance level where the node was created
    size_t                  m_hLevel; //!< level in the hierarchy (max number of nodes until reaching a leaf)
    bool                    m_flag; //!< flag bit
};


// === IN-LINE MEMBER FUNCTIONS ===

inline bool WHnode::isNode() const
{
    return m_fullID.first;
}

inline bool WHnode::isLeaf() const
{
    return !( m_fullID.first );
}

inline bool WHnode::isFlagged() const
{
    return m_flag;
}

inline size_t WHnode::getID() const
{
    return m_fullID.second;
}

inline nodeID_t WHnode::getFullID() const
{
    return m_fullID;
}

inline nodeID_t WHnode::getParent() const
{
    return m_parent;
}

inline size_t WHnode::getSize() const
{
    return m_nodeSize;
}

inline dist_t WHnode::getDistLevel() const
{
    return m_distanceLevel;
}

inline size_t WHnode::getHLevel() const
{
    return m_hLevel;
}

inline std::vector<nodeID_t> WHnode::getChildren() const
{
    return m_children;
}


// === NON-MEMBER OPERATORS ===

/**
 * << operator for the node class
 * \param os output stream
 * \param object node to print out
 * \return ostream with the node information
 */
std::ostream&   operator <<( std::ostream& os, const WHnode& object );
#endif  // WHNODE_H
