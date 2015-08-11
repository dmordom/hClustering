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


#ifndef PROTONODE_H
#define PROTONODE_H

// std library
#include <utility>
#include <map>


// hClustering
#include "WHnode.h"

const dist_t noNbDist( 999 );
const nodeID_t noNbID( std::make_pair( false, 0 ) );

/**
 * This class implements a tree node information container without hierarchical relationships to use while tree building
 */
class protoNode
{
public:
    /**
     * Constructor
     * \param nearNb information on the nearest neighbor of this protonode
     * \param nbNodes a map wit  the information on all the neighbor of this protonode
     * \param isactive a boolean flag indicating whther this protonode is currently active (if its allowed to be merged with another protonode)
     */
    protoNode( std::pair< nodeID_t, dist_t > nearNb, std::map< nodeID_t, dist_t > nbNodes, bool isactive = true ) :
               m_nearNb( nearNb ), m_nbNodes( nbNodes ), m_discarded( false ),  m_active( isactive ) {}

    //! Destructor
    ~protoNode() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * returns true if protonode is active (allowed to merge with other protonode)
     * \return active flag
     */
    inline bool isActive() const { return m_active; }

    /**
     * returns true if protonode has been discarded as outlier
     * \return discarded flag
     */
    inline bool isDiscarded() const { return m_discarded; }

    /**
     * returns the distance (dissimilarity) to the closest neighbor of the protonode
     * \return nearest neighbor distance value
     */
    inline dist_t nearDist() const { return m_nearNb.second; }

    /**
     * returns the ID of the closest neighbor of the protonode
     * \return nearest neighbor full-ID
     */
    inline nodeID_t nearNb() const { return m_nearNb.first; }

    /**
     * erases the stored neighbor information
     */
    inline void clearNbhood() { m_nbNodes.clear(); }

    /**
     * discards the protonode as outlier
     */
    inline void discard() { m_nbNodes.clear(); m_discarded = true; }

    /**
     * inactivates the protonode so that it may not be merged until reactivated
     */
    inline void inactivate() { m_active = false; }

    /**
     * reactivates the protonode to allow for it being merged
     */
    inline void reactivate() { m_active = true; }


    // === PUBLIC MEMBER FUNCTIONS ===


    /**
     * replaces old neighbourhood data and substitutes it with updated data when neighbour voxels join in the tree
     * \param oldNode1 ID of the first protonode of the pair that was merged
     * \param oldNode2 ID of the second protonode of the pair that was merged
     * \param newNode ID of the newly created protonode
     * \param newDist distance to the newly created protonode
     * \return if true the new nearest neighbor for this protonode has changed
     */
    bool updateNbhood( const nodeID_t& oldNode1,const  nodeID_t& oldNode2, const nodeID_t& newNode, const dist_t newDist );

    /**
     * replaces old neighbourhood data and substitutes it with updated data when neighbour voxels join in the tree
     * takes into account active flag information so that nearest neighbor will be the nearest active neighbor
     * \param oldNode1 ID of the first protonode of the pair that was merged
     * \param oldNode2 ID of the second protonode of the pair that was merged
     * \param newNode ID of the newly created protonode
     * \param newDist distance to the newly created protonode
     * \param isActive active flag of the newly created protonode
     * \param protoNodes the vector containing all the protonodes
     * \return if true the new nearest neighbor for this protonode has changed
     */
    bool updateActivhood( const nodeID_t& oldNode1, const nodeID_t& oldNode2,
                          const nodeID_t& newNode, const dist_t newDist,
                          const bool isActive, const std::vector< protoNode >& protoNodes );

    /**
     * update the stored distance value for the indicated neighbor
     * \param updatedNode neighbor ID which distance value needs updating
     * \param updatedDist new distance value to the indicated neighbor
     */
    void updateDist( const nodeID_t& updatedNode, const dist_t updatedDist);

   /**
     * scans the neighbourhood information and updates the nearest neighbor information to that of the closest active neighbor protonode
     * \param protoNodes the vector containing all the protonodes
     * \return if true the new nearest neighbor for this protonode has changed
     */
    bool updateActive( const std::vector< protoNode >& protoNodes );


    // === PUBLIC DATA MEMBERS ===

    std::pair< nodeID_t, dist_t > m_nearNb; //!< current nearest neighbour data
    std::map< nodeID_t, dist_t > m_nbNodes; //!< list of current neighbours data

private:
    // === PRIVATE DATA MEMBERS ===

    bool m_active;    //!< active flag
    bool m_discarded; //!< discarded flag
};


// === NON-MEMBER OPERATORS ===

/**
 * << operator for the protonode class
 * \param os output stream
 * \param object protonode to print out
 * \return ostream with the protonode information
 */
std::ostream& operator <<( std::ostream& os, const protoNode& object );

#endif  // PROTONODE_H
