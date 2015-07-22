#include "get_names.h"

// "get_Vname()": returns the filename of the tractogram corresponding to the given seed voxel coordinates
std::string get_Vname (std::string path, coordinate seed) {

    std::string filename(path + "/connect_" + seed.get_name_string() + ".v");
    return filename;

}// end "get_Vname()" -----------------------------------------------------------------


// "get_Vname()": returns the filename of the tractogram corresponding to the given leaf identifier
std::string get_Vname (std::string path, nodeID tree_node, std::vector<coordinate> &roivect) {

    std::string filename;
    if(!tree_node.first)
        filename = get_Vname(path,roivect[tree_node.second]);
    else
        throw std::runtime_error ("ERROR [get_Vname()]: calling leaf filename retrieval for a node");
    return filename;

}// end "get_Vname()" -----------------------------------------------------------------


