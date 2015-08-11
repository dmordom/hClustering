//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
//
// For more reference on the underlying algorithm and research they have been used for refer to:
// - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
//   A hierarchical method for whole-brain connectivity-based parcellation.
//   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528
// - Moreno-Dominguez, D. (2014).
//   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
//   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
//   ISBN 978-3-941504-45-5
//
// hClustering is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------


#ifndef COMPACTTRACTCHAR_H
#define COMPACTTRACTCHAR_H

// parallel execution
#include <omp.h>

// std library
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <climits>
#include <ctime>

// boost library
#include "boost/date_time/posix_time/posix_time.hpp"

#include "compactTract.h"

class compactTract;

/**
 * This class stores the data from a vector-compacted probabilistic tractogram in 8-bit (unsigned char) precision
 * it also keeps track of the thresholded status of the data and implements the tractogram dissimilarity measures
 */
class compactTractChar
{
public:
    /**
     * Constructor
     */
    compactTractChar() :
        m_norm( 0 ), m_thresholded( false ), m_normReady( false ) {}
    /**
     * \overload
     * \param tractInit a vector with the tractogram data
     */
    explicit compactTractChar( std::vector< unsigned char > tractInit ) :
        m_tract( tractInit ), m_norm( 0 ), m_thresholded( false ), m_normReady( false ) {}
    /**
     * \overload
     * \param object a compactTractChar object to copy the member data from
     */
    compactTractChar( const compactTractChar &object ) :
        m_tract( object.m_tract ), m_norm( object.m_norm ), m_thresholded( object.m_thresholded ),
        m_normReady( object.m_normReady ) {}

    //! Destructor
    ~compactTractChar() {}

    // === MEMBER OPERATORS ===

    /**
     * = Operator
     * \param rhs object to be copied
     */
    compactTractChar& operator =( const compactTractChar &rhs );

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * returns the size of the compact tract vector
     * \return tractogram size (number of voxels in the white matter)
     */
    inline size_t size() const { return m_tract.size(); }

    /**
     * returns true if the tractogram vector norm has been precomputed and saved in the tract object
     * \return norm ready boolean flag
     */
    inline bool normReady() const { return m_normReady; }

    /**
     * returns true if the tractogram vector data has been thresholded
     * \return thresholded boolean flag
     */
    inline bool thresholded() const { return m_thresholded; }

    /**
     * returns a copy of the tractogram data stored
     * \return a vector containing the tractogram data
     */
    inline std::vector< unsigned char > tract() const { return m_tract; }

    /**
     * saves a precomputed vector norm value in the tractogram object
     * \param norm the norm data in double precision
     */
    inline void setNorm( double norm ) { m_norm = norm; m_normReady = true; }

    // === PUBLIC MEMBER FUNCTIONS ===


    /**
     * returns the total size of the compactTractChar object in bytes (including the data vector)
     * \return size in bytes of a compactTractChar object
     */
    size_t bytes() const;

    /**
     * returns the total size of the compactTractChar object in megaBytes (including the data vector)
     * \return size in megaBytes of a compactTractChar object
     */
    float mBytes() const;

    /**
     * swaps the tractogram memory from another tractogram object into this one and copies its data members
     * implemented as a time and memory efficient alternative of the copy constructor if a new object is needed but the old object is not required anymore
     * \param stolen pointer to the tractogram object with the data to be "stolen"
     */
    void steal( compactTractChar* const stolen );

    /**
     * computes the distance (dissimilarity) between this tract and a tract defined by the parameter
     * by default calls to normDotProduct(), but implemented to easily change the distance measure used with minimum changes to the source code
     * if planning to use this code/software for functional connectivity analysis/processing, call correlation() instead
     * \param tractogram compactTract object to compute the distance to
     * \return the distance value between the two tracts in double precision
     */
    double tractDistance( const compactTractChar &tractogram ) const;
    /**
     * \overload
     */
    double tractDistance( const compactTract &tractogram ) const;

    /**
     * computes and returns the norm (rooted square-sum) of the tractogram and saves it in the tarct object data member
     * if the norm was previously computed will re-compute and overwrite
     * \return the tractogram norm value
     */
    double computeNorm();

    /**
     * returns a compactTract in float precision with the data values transformed doing a 10^(x*f) exponential where x is the tract data and f is the normalization-related logFactor
     * \param logFactor normalization-related parameter to properly switch between logarithmic and natural units
     * \return the transformed tractogram in float precision
     */
    compactTract unLog( float logFactor ) const;

    /**
     * thresholds the tractogram data, if the value of a point is less than the given threshold, it is set to 0
     * \param threshold threshold value
     */
    void threshold( float threshold );


    friend class fileManager;
    friend class vistaManager;
    friend class niftiManager;
    friend class compactTract;


protected:
    // === PROTECTED DATA MEMBERS ===

    std::vector< unsigned char > m_tract; //!< unsigned char vector with the compact tractogram data
    double  m_norm;         //!< norm of the data vector
    bool    m_thresholded;  //!< thresholded flag (if true, data vector has been thresholded)
    bool    m_normReady;    //!< norm ready flag (if true, norm has been precomputed and saved)

private:
    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * computes the normalized dot product between this tract and a tract defined by the parameter (tractograms must be thresholded).
     * \param tractogram compactTract object to compute the norm. dot product with
     * \return the norm. dot product value between the two tracts in double precision
     */
    double normDotProduct( const compactTractChar &tractogram ) const;

    /**
     * computes pearsons correlation coefficient between this tract and a tract defined by the parameter
     * DEPRECATED IN FAVOR OF NORMALIZED DOT PRODUCT VALUE (for tractograms/anatomical connectivity)
     * if planning to use this code for functional correlation use this measure over nor. dot product
     * \param tractogram compactTract object to compute the correlation to
     * \return the correlation value between the two tracts in double precision
     */
    double correlation( const compactTractChar &tractogram ) const;
};

// === NON-MEMBER OPERATORS ===


/**
 * << operator for the compact tract class
 * \param os output stream
 * \param object tract to print out
 * \return ostream with the tract information
 */
std::ostream& operator <<( std::ostream& os, const compactTractChar& object );

#endif  // COMPACTTRACTCHAR_H
