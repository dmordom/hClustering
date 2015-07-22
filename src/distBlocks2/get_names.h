#ifndef GET_NAMES_H
#define GET_NAMES_H

// std library
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdio>


// "get_Vname()": returns the filename of the tractogram corresponding to the given seed voxel coordinates
std::string get_Vname (std::string path, size_t seed);



#endif // GET_NAMES_H
