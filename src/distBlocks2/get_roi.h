#ifndef READ_TREE_H
#define READ_TREE_H

// std library
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

// classes
#include "coordinate.h"


void read_tree(std::string roi_file,
               std::vector<coordinate> &roivect,
               std::map<coordinate,size_t> &roimap);

#endif // READ_TREE_H
