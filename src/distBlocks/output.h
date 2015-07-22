#ifndef OUTPUT_H
#define OUTPUT_H

// std library
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <stdexcept>

// classes
#include "coordinate.h"

// "write_output():" Writes the output files (tree, pruned roi, leaves, nodes, discarded seeds)
void write_output(std::string path,
                  std::vector<coordinate> &roivect,
                  std::vector<std::pair<size_t,size_t> > &roi_block_index);

#endif // OUTPUT_H
