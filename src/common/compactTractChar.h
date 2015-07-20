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

class compactTractChar
{
public:
    //! Constructor
    compactTractChar() :
        m_norm( 0 ), m_thresholded( false ), m_normReady( false )
    {
    }
    explicit compactTractChar( std::vector< unsigned char > tractInit ) :
        m_tract( tractInit ), m_norm( 0 ), m_thresholded( false ), m_normReady( false )
    {
    }

    compactTractChar( const compactTractChar &object ) :
        m_tract( object.m_tract ), m_norm( object.m_norm ), m_thresholded( object.m_thresholded ),
        m_normReady( object.m_normReady )
    {
    }

    //! Destructor
    ~compactTractChar()
    {
    }

    //! = Operator
    compactTractChar& operator =( const compactTractChar &rhs );

    //! swaps the tractogram memory to the other tract class and copies its properties
    void steal( compactTractChar* const stolen );

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

    std::vector< unsigned char > tract() const
    {
        return m_tract;
    }

    double tractDistance( const compactTractChar &tractogram ) const;
    double tractDistance( const compactTract &tractogram ) const;


    // "normInnerProduct()": returns  the nprmalized inner product between two tractograms (tractograms must be thresholded).
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
    compactTract unLog( float logFactor ) const;


    // "threshold()": thresholds the input tractogram, if the value of a point is less than the given threshold, it is set to 0
    void threshold( float threshold );


    friend class vistaManager;
    friend class compactTract;


protected:
    std::vector< unsigned char > m_tract;
    double m_norm;
    bool m_thresholded;
    bool m_normReady;

private:
    // "correlation()": returns  the correlation coefficient between two tractograms (tractograms must be thresholded).
    double correlation( const compactTractChar &tractogram ) const;
};

#endif  // COMPACTTRACTCHAR_H
