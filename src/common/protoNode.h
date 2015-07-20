#ifndef PROTONODE_H
#define PROTONODE_H

// std library
#include <utility>
#include <map>


// hClustering
#include "WHnode.h"

const dist_t noNbDist( 999 );
const nodeID_t noNbID( std::make_pair( false, 0 ) );

class protoNode
{
public:
    protoNode( std::pair< nodeID_t, dist_t > nearNb, std::map< nodeID_t, dist_t > nbNodes, bool isactive = true ) :
        m_nearNb( nearNb ), m_nbNodes( nbNodes ), m_discarded( false ),  m_active( isactive )
    {
    }
    ~protoNode()
    {
    }

    bool isActive() const
    {
        return m_active;
    }
    bool isDiscarded() const
    {
        return m_discarded;
    }
    dist_t nearDist() const
    {
        return m_nearNb.second;
    }
    nodeID_t nearNb() const
    {
        return m_nearNb.first;
    }

    void clearNbhood()
    {
        m_nbNodes.clear();
    }

    void discard()
    {
        m_nbNodes.clear();
        m_discarded = true;
    }

    void inactivate()
    {
        m_active = false;
    }

    void reactivate()
    {
        m_active = true;
    }

    // "updateNbhood(): replaces old neighbourhood data and substitutes it with updated data when neighbour voxels join in the tree
    // return value indicates if the new nearest neighbor has changed
    bool updateNbhood( nodeID_t oldNode1, nodeID_t oldNode2, nodeID_t newNode, dist_t newDist );


    bool updateActivhood( const nodeID_t oldNode1, const nodeID_t oldNode2, const nodeID_t newNode, const dist_t newDist,
                                     const bool isActive, const std::vector< protoNode > &protoNodes );

    void updateDist( const nodeID_t updatedNode, const dist_t updatedDist);

    bool updateActive( const std::vector< protoNode > &protoNodes );


    // member variables

    std::pair< nodeID_t, dist_t > m_nearNb; // current nearest neighbour data
    std::map< nodeID_t, dist_t > m_nbNodes; // list of current neighbours data
    bool m_active;
    bool m_discarded;
};

// non-member operators
std::ostream& operator <<( std::ostream& os, const protoNode& object );

#endif  // PROTONODE_H
