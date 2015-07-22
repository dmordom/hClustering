#ifndef DO_VISTA_H
#define DO_VISTA_H

// std library
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

// vista
#include <viaio/Vlib.h>
#include <viaio/VImage.h>

// parallel execution
//#include <omp.h>


// "get_Vtract_th()": extracts tractogram data from a vista file and returns it thresholded in float array form
void get_Vtract_th (std::string filename, unsigned char *tractogram, size_t tract_length, float tractThreshold);

// "get_Vtract()": extracts tractogram data from a vista file and returns it non-thresholded in std vector form
std::vector<char> get_Vtract (std::string tract_filename);

// "write_dist_block()": writes down distance block to file
void write_dist_block (std::string filename, std::vector<std::vector<float> > &dist_block);

// "write_Vtract()": write compact tractogram data into vista format file
void write_Vtract (std::string filename, std::vector<float> &tractogram);

// "WriteVImage()": write Vista image
int WriteVImage( const VString Name, VImage &Image );


#endif // DO_VISTA_H
