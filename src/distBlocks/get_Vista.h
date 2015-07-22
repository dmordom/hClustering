#ifndef GET_VISTA_H
#define GET_VISTA_H

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

// classes
#include "coordinate.h"

// "get_Vtract_th()": extracts tractogram data from a vista file and returns it thresholded in float array form
void get_Vtract_th (std::string filename, float *tractogram, unsigned int tract_length, float threshold = 0.4);
unsigned int get_Vtract_th (std::string filename, float **tractogram, float threshold = 0.4);

// "get_Vtract()": extracts tractogram data from a vista file and returns it in std vector form
std::vector<float> get_Vtract (std::string tract_filename);

void write_dist_block (std::string filename, std::vector<std::vector<float> > &dist_block);

// "write_Vtract()": write compact tractogram data into vista format file
void write_Vtract (std::string filename, std::vector<float> &tractogram);

// "store_Vtract_Image()": write full image tractogram into vista format file and compress
void store_Vtract_Image ( std::string tract_filename, std::string mask_filename, std::vector<float> &tractogram);

// "WriteVImage()": write Vista image
int WriteVImage( const VString Name, VImage &Image );

// "ReadImage()": read Vista image file
int ReadImage( const VString Name, VImage &Image );

#endif // GET_VISTA_H
