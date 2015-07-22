#include "tnode.h"

// member functions



// "printdata(): prints out a string with the node data, in order to implement the operator <<"
void tnode::printdata(std::ostream& os) const {

    os <<"ID: "<< ID.first <<"-"<< ID.second;
    os <<".  Dad: "<< parent.first <<"-"<< parent.second;
    os <<".  Kids: ( "<< children.first.first <<"-"<<children.first.second <<" , "<< children.second.first <<"-"<<children.second.second <<" )";
    os <<".  #Leaves: "<< nleaves;
    os <<".  Level: "<< level;
    os <<".  Prune: "<< prune;

    return;
}// end "printdata()" -----------------------------------------------------------------



// non-member operators

std::ostream& operator <<(std::ostream& os, const tnode& object) {

    tnode temp_node(object);
    temp_node.printdata(os);
    return os;
}// end operator << -----------------------------------------------------------------
