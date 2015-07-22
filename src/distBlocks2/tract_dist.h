#ifndef TRACT_DIST_H
#define TRACT_DIST_H

// std library
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>


double tract_distance(std::vector<float> &tractogram1, std::vector<float> &tractogram2);

// "vectprod()": returns  the vecotial product between two tractograms (tractograms must be in logarithmic units and thresholded).
double vectprod(std::vector<float> &tractogram1, std::vector<float> &tractogram2);

// "correlate()": returns  the correlation coefficient between two tractograms (tractograms must be in logarithmic units and thresholded).
float correlate(std::vector<float> &tractogram1, std::vector<float> &tractogram2);



#endif // TRACT_DIST_H
