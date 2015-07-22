#ifndef TRACT_DIST_H
#define TRACT_DIST_H

// std library
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

// classes
#include "tnode.h"

double tract_distance(std::vector<float> &tractogram1, std::vector<float> &tractogram2);

// "vectprod()": returns  the vecotial product between two tractograms (tractograms must be in logarithmic units and thresholded).
double vectprod(std::vector<float> &tractogram1, std::vector<float> &tractogram2);

// "correlate()": returns  the correlation coefficient between two tractograms (tractograms must be in logarithmic units and thresholded).
float correlate(std::vector<float> &tractogram1, std::vector<float> &tractogram2);

// "join_tracts()": returns a vector containing the weighted mean tractogram of the two input tracts (tractograms must be in natural units)
std::vector<float> join_tracts(std::vector<float> tractogram1, std::vector<float> tractogram2, nodes_size_t nleaves1, nodes_size_t nleaves2);

// "unLog()": transforms the input tractogram doing a 10^x exponential
void unLog(std::vector<float> &tractogram);

// "doLog()": transforms the input tractogram doing a base-10 logarithm
void doLog(std::vector<float> &tractogram);

// "threshold()": thresholds the input tractogram, if the value of a point is less than the given threshold, it is set to 0
//void threshold(std::vector<float> &tractogram, float threshold= 0.4);


#endif // TRACT_DIST_H
