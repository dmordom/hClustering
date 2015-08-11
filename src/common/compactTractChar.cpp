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


#include <vector>

#include "compactTractChar.h"



/*
 This program allows you to calculate the similarity between tractograms.
 Other possible similarity indices would be:

 * Correlation   "signal product"  co = 1/N*sum(x*y)
 * Cross-correlation coefficient   cc = 1/N*sum{(x-mean_x)*(y-mean_y)}/(stddev_x*stddev_y)
 * Euclidian Distance              eu = sqrt[sum{(x-y)*(x-y)}]
 * Mutual Information              mi = sum [ P(x,y)*log2{P(x,y)/(P(x)*P(y))} ]
 * Covariance                      covar = 1/N*sum{(x-mean_x)*(y-mean_y)}= 1/N*sum(x*y) -mean_x*mean_y
 */
double compactTractChar::tractDistance( const compactTractChar &tractogram ) const
{
    return 1 - normDotProduct( tractogram );
}
double compactTractChar::tractDistance( const compactTract &tractogram ) const
{
    return tractogram.tractDistance(*this);
} // end "tractDistance()" -----------------------------------------------------------------


/* "normInnerProduct()":
 Calculates the correlation coefficient between two tractograms. (tractograms must be in logarithmic units and thresholded)
 The correlation is calculated and normalized as:

 considering X = tractogram1 ; Y = tractogram2

 X.*Y                          Sum(Xi*Yi)
 NIP = --------------------- = -----------------------------
 ||X||*||Y||            sqrt[Sum(Xi2)] * sqrt[Sum(Yi2)]

 */
double compactTractChar::normDotProduct( const compactTractChar &tractogram ) const
{
    if( m_tract.size() != tractogram.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTractChar::normInnerProduct(): Tractograms are not of the same size" );
    }
    else if( ( !this->m_normReady ) || ( !tractogram.m_normReady ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTractChar::normInnerProduct(): one (or both) of the tracts has no available precomputed norm" );
    }
    else if( ( !this->m_thresholded ) || ( !tractogram.m_thresholded ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTractChar::normInnerProduct(): one (or both) of the tracts has not been thresholded" );
    }
    else if( ( this->m_norm == 0. ) || ( tractogram.m_norm == 0. ) )
    {
        std::cerr << "WARNING @ compactTractChar::normInnerProduct(): At least one of the tractograms is a zero vector, inner product will be set to 0"
                  << std::endl;
        return 0.;
    }


        // Initialize variables
    double inProd( 0 ), dotprod_sum( 0 );


    const unsigned char* piter1, *piter2;
    piter1 = &m_tract.front();
    piter2 = &tractogram.m_tract.front();
    size_t tractSize( m_tract.size() );

    for( size_t i = 0; i < tractSize; ++i )
    {
        dotprod_sum += ( *piter1++ ) * ( *piter2++ );
    }

    inProd = dotprod_sum / ( this->m_norm * tractogram.m_norm );

    if( inProd < 0 )
    {
        if( inProd < -0.0001 )
            std::cerr << std::endl << "WARNING @ compactTractChar::normInnerProduct(): Negative inner product (" << inProd << ")"
                            << std::endl;
        inProd = 0;
    }
    else if( inProd > 1 )
    {
        if( inProd > 1.0001 )
        {
            std::cerr << std::endl << "WARNING @ compactTract::normInnerProduct(): Bad inner product (" << inProd << ")" << std::endl;
        }
        inProd = 1;
    }

    return inProd;
} // end "normInnerProduct()" -----------------------------------------------------------------


double compactTractChar::computeNorm()
{
    if( !this->m_thresholded )
    {
        std::cerr << "ERROR @ compactTractChar::correlation(): tract has not been thresholded" << std::endl;
        return 0;
    }
    else
    {
        m_norm = 0;
        for( std::vector< unsigned char >::const_iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
            m_norm += ( *iter ) * ( *iter );

        m_norm = sqrt( m_norm );
        m_normReady = true;
        return m_norm;
    }
} // end "getNorm()" -----------------------------------------------------------------


size_t compactTractChar::bytes() const
{
    size_t elementSize( 0 );
    if( m_tract.empty() )
    {
        unsigned char el0( 0 );
        elementSize = sizeof( el0 );
    }
    else
    {
        elementSize = sizeof( m_tract.front() );
    }
    return ( sizeof( *this ) + ( elementSize * m_tract.size() ) ) * CHAR_BIT / 8.;
}

float compactTractChar::mBytes() const
{
    return bytes() / ( 1024. * 1024. );
}


// "threshold()": thresholds the input tractogram, if the value of a point is less than the given threshold, it is set to 0
void compactTractChar::threshold( const float threshold )
{

    if( this->m_thresholded )
    {
        std::cerr << "WARNING @ compactTractChar::threshold(): tract has already been thresholded" << std::endl;
    }
    else
    {
        if( threshold != 0 )
        {
            unsigned char charThreshold( 255 * threshold);
            for( std::vector< unsigned char >::iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
                if( *iter < charThreshold )
                    *iter = 0.;
        }
        m_thresholded = true;
    }
    return;
} // end "threshold()" -----------------------------------------------------------------




void compactTractChar::steal( compactTractChar* const stolen )
{
    m_tract.swap( stolen->m_tract );
    m_norm = stolen->m_norm;
    m_thresholded = stolen->m_thresholded;
    m_normReady = stolen->m_normReady;
} // end "steal()" -----------------------------------------------------------------

// member operators
compactTractChar& compactTractChar::operator =( const compactTractChar &rhs )
{
    m_tract = rhs.m_tract;
    m_norm = rhs.m_norm;
    m_thresholded = rhs.m_thresholded;
    m_normReady = rhs.m_normReady;
    return *this;
} // end "operator =" -----------------------------------------------------------------

// PRIVATE MEMBERS


/* "correlation()":
 Calculates the correlation coefficient between two tractograms. (tractograms must be in logarithmic units and thresholded)
 The correlation is calculated and normalized as:

 considering X = tractogram1 ; Y = tractogram2

 Cov(X,Y)                E( [X-E(X)]*[Y-E(Y)] )                sum( [X-E(X)]*[Y-E(Y)] )
 CC = --------------------- = ----------------------------- = --------------------------------------------
 stddev(X)*stddev(Y)     sqrt[Var(X)] * sqrt[Var(Y)]     N * sqrt[E(X^2)-E(X)^2] * sqrt[E(Y^2)-E(Y)^2]


 This program allows you to calculate the similarity between tractograms.
 */
double compactTractChar::correlation( const compactTractChar &tractogram ) const
{
    if( m_tract.size() != tractogram.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTractChar::correlation(): Tractograms are not of the same size" );
    }
    else if( ( !this->m_normReady ) || ( !tractogram.m_normReady ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTractChar::correlation(): one (or both) of the tracts has no available precomputed norm" );
    }
    else if( ( !this->m_thresholded ) || ( !tractogram.m_thresholded ) )
    {
        throw std::runtime_error( "ERROR @ compactTractChar::correlation(): one (or both) of the tracts has not been thresholded" );
    }
    else if( ( this->m_norm == 0. ) || ( tractogram.m_norm == 0. ) )
    {
        std::cerr
                        << "WARNING @ compactTractChar::correlation(): At least one of the tractograms is a zero vector, correlation will be set to 0"
                        << std::endl;
        return 0.;
    }
    else
    {
        // Initialize variables
        double sum1( 0 ), avr1( 0 ), var1( 0 ), stddev1( 0 );
        double sum2( 0 ), avr2( 0 ), var2( 0 ), stddev2( 0 );
        double corr( 0 ), cov( 0 );

        // Compute average and stddev;  E(X) = sum(x)/N;  Sttdev = sqrt[E(X)-E(X)]
#pragma omp parallel sections
        {
#pragma omp section
            {
                for( std::vector< unsigned char >::const_iterator iter( this->m_tract.begin() ); iter != this->m_tract.end(); ++iter )
                    sum1 += *iter; // sum(X)
            }
#pragma omp section
            {
                for( std::vector< unsigned char >::const_iterator iter( tractogram.m_tract.begin() ); iter != tractogram.m_tract.end(); ++iter )
                    sum2 += *iter; // sum(Y)
            }
        }

        avr1 = sum1 / this->m_tract.size(); // E(X)
        avr2 = sum2 / tractogram.m_tract.size(); // E(Y)
        var1 = ( ( this->m_norm * this->m_norm ) / this->m_tract.size() ) - ( avr1 * avr1 ); // Var(X) = E(X2)-E2(X)  [E(X2) = ||X||2/N ]
        var2 = ( ( tractogram.m_norm * tractogram.m_norm ) / tractogram.m_tract.size() ) - ( avr2 * avr2 ); // Var(Y) = E(Y2)-E2(Y)

        if( ( var1 == 0. ) || ( var2 == 0 ) )
        {
            std::cerr << "WARNING @ compactTractChar::correlation(): One (or both) of the tractograms is a non-zero constant vector,"
                      << " correlation will be set to 0" << std::endl;
            return 0.;
        }

        stddev1 = sqrt( var1 ); // sttdev(X) = sqrt(Var(X))
        stddev2 = sqrt( var2 ); // sttdev(Y) = sqrt(Var(Y))

        // Compute Covariance Cov(X,Y) = E( [X-E(X)]*[Y-E(Y)] )
        for( std::vector< unsigned char >::const_iterator iter1( this->m_tract.begin() ), iter2( tractogram.m_tract.begin() ); iter1
                        != this->m_tract.end(); ++iter1, ++iter2 )
        {
            cov += ( *iter1 - avr1 ) * ( *iter2 - avr2 );
        }
        corr = cov / ( this->m_tract.size() * stddev1 * stddev2 );

        if( corr < 0 )
        {
            std::cerr << std::endl << "WARNING @ compactTractChar::correlation(): Negative correlation (" << corr
                            << ") saving it as 0" << std::endl;
            corr = 0;
        }
        else if( corr > 1.0001 )
        {
            std::cerr << std::endl << "WARNING @ compactTractChar::correlation(): Bad correlation (" << corr << ")" << std::endl;
            corr = 1;
        }

        return corr;
    }
} // end "correlation()" -----------------------------------------------------------------

std::ostream& operator <<( std::ostream& os, const compactTractChar& object )
{
    std::vector< unsigned char > tract( object.tract() );
    size_t counter0 ( 0 );
    for( size_t i = 0; i < tract.size() && counter0 < 15; ++i )
    {
        size_t datapoint(tract[i]);
        if(datapoint)
        {
            os << datapoint << " ";
            ++counter0;
        }
    }
    return os;
} // end "operator << " -----------------------------------------------------------------
