#ifndef GET_NAMES_H
#define GET_NAMES_H

// std library
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>

// classes
#include "coordinate.h"
#include "tnode.h"

// "get_Vname()": returns the filename of the tractogram corresponding to the given seed voxel coordinates
std::string get_Vname (std::string path, coordinate seed);

// "get_Vname()": returns the filename of the tractogram corresponding to the given leaf identifier
std::string get_Vname (std::string path, nodeID tree_node, std::vector<coordinate> &roivect);


#endif // GET_NAMES_H
