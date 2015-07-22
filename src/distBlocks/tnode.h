#ifndef TNODE_H
#define TNODE_H

// std library
#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>

// pre-declare classes to use as arguments in member functions
class tnode;

// typedefs
typedef float dist_t;
typedef std::vector<tnode>::size_type nodes_size_t;
typedef std::pair<bool, nodes_size_t> nodeID;

class tnode
{
public:
    tnode (nodeID ID_init, nodeID parent_init): ID(ID_init), parent(parent_init), children(std::make_pair(std::make_pair(false,0),std::make_pair(false,0))), nleaves(1), level(0), prune(false)  {}
    tnode (nodeID ID_init, nodeID parent_init, std::pair<nodeID,nodeID> children_init, nodes_size_t nleaves_init, dist_t level_init): ID(ID_init), parent(parent_init), children(children_init), nleaves(nleaves_init), level(level_init), prune(false)  {}

    // implicit member functions
    nodeID                      get_ID() const {return ID;}
    nodeID                      get_parent() const {return parent;}
    std::pair<nodeID,nodeID>    get_children() const {return children;}
    nodes_size_t                get_nleaves() const {return nleaves;}
    dist_t                      get_level() const {return level;}
    bool                        get_prune() const {return prune;}

    void                        set_prune() {prune=true;}
    void                        change_ID(nodeID newID) {ID = newID;}
    void                        change_parent(nodeID newpapa) {parent = newpapa;}
    void                        change_children(std::pair<nodeID,nodeID> new_children) {children = new_children;}
    void                        change_nleaves(nodes_size_t newleaves) {nleaves = newleaves;}
    void                        change_level(dist_t newlevel) {level = newlevel;}



    // explicit member functions

    // "printdata(): prints out a string with the node data, in order to implement the operator <<"
    void printdata(std::ostream& os) const;

protected:

    // data members
    nodeID ID;                          // node ID
    nodeID parent;                      // parent node ID
    std::pair<nodeID,nodeID> children;  // chlidren nodes IDs
    nodes_size_t nleaves;               // number of leaves contained on the node
    dist_t level;                       // distance level where the node was created
    bool prune;                         // prune bit
};


// non-member operators
std::ostream&   operator <<(std::ostream& os, const tnode& object);


#endif // TNODE_H
