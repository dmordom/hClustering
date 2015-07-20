#ifndef COMPACTTRACT_H
#define COMPACTTRACT_H

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

#include "compactTractChar.h"

class compactTractChar;

class compactTract
{
public:
    //! Constructor
    compactTract();
    explicit compactTract( std::vector< float > tractInit );
    compactTract( const compactTract &object );
    compactTract( const compactTractChar &charTract );
    compactTract( const compactTract &tract1, const compactTract &tract2, const size_t size1, const size_t size2 );




    //! Destructor
    ~compactTract()
    {
    }

    //! = Operator
    compactTract& operator =( const compactTract &rhs );
    compactTract& operator =( const compactTractChar &rhs );

    //! swaps the tractogram memory to the other tract class and copies its properties
    void steal( compactTract* const stolen );

    size_t size() const
    {
        return m_tract.size();
    }



    bool normReady() const
    {
        return m_normReady;
    }

    bool thresholded() const
    {
        return m_thresholded;
    }

    std::vector< float > tract() const
    {
        return m_tract;
    }

    double tractDistance( const compactTract &tractogram ) const;
    double tractDistance( const compactTractChar &tractogram ) const;


    // "normInnerProduct()": returns  the nprmalized inner product between two tractograms (tractograms must be thresholded).
    double normDotProduct( const compactTract &tractogram ) const;
    double normDotProduct( const compactTractChar &tractogram ) const;


    // "getNorm()": returns  the rooted square-sum the tractogram.
    double getNorm();

    void setNorm( double norm )
    {
        m_norm = norm;
        m_normReady = true;
    }

    size_t bytes() const;

    float mBytes() const;

    // "unLog()": transforms the input tractogram doing a 10^x exponential
    void unLog( float logFactor );

    // "doLog()": transforms the input tractogram doing a base-10 logarithm
    void doLog( float logFactor );

    // "threshold()": thresholds the input tractogram, if the value of a point is less than the given threshold, it is set to 0
    void threshold( float threshold );

    // "add()": sum input tractogram top the tractogram object
    void add( const compactTract &tractogram );

    // "divide()": divide tractogram by a given value
    void divide( float divisor );

    friend class vistaManager;
    friend class randCnbTreeBuilder;


protected:
    std::vector< float > m_tract;
    double m_norm;
    bool m_thresholded;
    bool m_normReady;
    bool m_inLogUnits;

private:
    // "correlation()": returns  the correlation coefficient between two tractograms (tractograms must be thresholded).
    double correlation( const compactTract &tractogram ) const;
};

#endif  // COMPACTTRACT_H
