
#include <vector>

#include "compactTract.h"

compactTract::compactTract() :
    m_norm( 0 ), m_thresholded( false ), m_normReady( false ), m_inLogUnits( true )
{
}
compactTract::compactTract( std::vector< float > tractInit ) :
    m_tract( tractInit ), m_norm( 0 ), m_thresholded( false ), m_normReady( false ), m_inLogUnits( true )
{
}

compactTract::compactTract( const compactTract &object ) :
    m_tract( object.m_tract ), m_norm( object.m_norm ), m_thresholded( object.m_thresholded ),
    m_normReady( object.m_normReady ), m_inLogUnits( object.m_inLogUnits )
{
}

compactTract::compactTract( const compactTractChar &charTract ) :
    m_norm( charTract.m_norm / 255. ), m_thresholded( charTract.m_thresholded ), m_normReady( charTract.m_normReady ),
    m_inLogUnits( true )
{
    m_tract.clear();
    m_tract.reserve( charTract.m_tract.size() );
    for( size_t i = 0; i < charTract.m_tract.size(); ++i )
    {
        m_tract.push_back( charTract.m_tract[i] / 255.0 );
    }
}
// creates a tract containing the weighted mean tractogram of the two input tracts (tractograms must be in natural units)
compactTract::compactTract( const compactTract &tract1, const compactTract &tract2, const size_t size1, const size_t size2 ) :
    m_norm( 0 ), m_thresholded( false ), m_normReady( false ), m_inLogUnits( false )
{
    if( tract1.m_tract.size() != tract2.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTract::joinTracts(): Tractograms are not of the same size" );
    }
    else if( ( tract1.m_thresholded ) || ( tract2.m_thresholded ) )
    {
        throw std::runtime_error( "ERROR @ compactTract::joinTracts(): one (or both) of the tracts has been thresholded" );
    }
    else if( ( tract1.m_inLogUnits ) || ( tract2.m_inLogUnits ) )
    {
        throw std::runtime_error( "ERROR @ compactTract::joinTracts(): one (or both) of the tracts is in logarithmic units" );
    }

    m_tract.clear();
    m_tract.reserve( tract1.m_tract.size() );
    for( size_t i = 0; i< tract1.m_tract.size(); ++i )
    {
        m_tract.push_back( ( ( tract1.m_tract[i] * size1 ) + ( tract2.m_tract[i] * size2 ) ) / ( size1 + size2 ) );
    }
}


/*
 This program allows you to calculate the similarity between tractograms.
 Other possible similarity indices would be:

 * Correlation   "signal product"  co = 1/N*sum(x*y)
 * Cross-correlation coefficient   cc = 1/N*sum{(x-mean_x)*(y-mean_y)}/(stddev_x*stddev_y)
 * Euclidian Distance              eu = sqrt[sum{(x-y)*(x-y)}]
 * Mutual Information              mi = sum [ P(x,y)*log2{P(x,y)/(P(x)*P(y))} ]
 * Covariance                      covar = 1/N*sum{(x-mean_x)*(y-mean_y)}= 1/N*sum(x*y) -mean_x*mean_y
 */
double compactTract::tractDistance( const compactTract &tractogram ) const
{
    return 1 - normDotProduct( tractogram );
}
double compactTract::tractDistance( const compactTractChar &tractogram ) const
{
    return 1 - normDotProduct( tractogram );
} // end "tractDistance()" -----------------------------------------------------------------

/* "normInnerProduct()":
 Calculates the correlation coefficient between two tractograms. (tractograms must be in logarithmic units and thresholded)
 The correlation is calculated and normalized as:

 considering X = tractogram1 ; Y = tractogram2

 X.*Y                          Sum(Xi*Yi)
 NIP = --------------------- = -----------------------------
 ||X||*||Y||            sqrt[Sum(Xi2)] * sqrt[Sum(Yi2)]

 */
double compactTract::normDotProduct( const compactTract &tractogram ) const
{
    if( m_tract.size() != tractogram.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTract::normInnerProduct(): Tractograms are not of the same size" );
    }
    else if( ( !this->m_normReady ) || ( !tractogram.m_normReady ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::normInnerProduct(): one (or both) of the tracts has no available precomputed norm" );
    }
    else if( ( !this->m_thresholded ) || ( !tractogram.m_thresholded ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::normInnerProduct(): one (or both) of the tracts has not been thresholded" );
    }
    else if( ( this->m_norm == 0. ) || ( tractogram.m_norm == 0. ) )
    {
        std::cerr << "WARNING @ compactTract::normInnerProduct(): At least one of the tractograms is a zero vector, inner product will be set to 0"
                  << std::endl;
        return 0.;
    }
    else if( ( !this->m_inLogUnits ) || ( !tractogram.m_inLogUnits ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::normInnerProduct(): one (or both) of the tracts is not in logaritmic units" );
    }


        // Initialize variables
    double inProd( 0 ), dotprod_sum( 0 );

#if 0

    for( std::vector<float>::const_iterator iter1( this->m_tract.begin() ), iter2( tractogram.m_tract.begin() );
                    iter1 != this->m_tract.end(); ++iter1, ++iter2 )
    {
        dotprod_sum += ( *iter1 )*( *iter2 );
    }

#else

    const float* piter1, *piter2;
    piter1 = &m_tract.front();
    piter2 = &tractogram.m_tract.front();
    size_t tractSize( m_tract.size() );

    for( size_t i = 0; i < tractSize; ++i )
    {
        dotprod_sum += ( *piter1++ ) * ( *piter2++ );
    }

#endif

    inProd = dotprod_sum / ( this->m_norm * tractogram.m_norm );

    if( inProd < 0 )
    {
        if( inProd < -0.0001 )
            std::cerr << std::endl << "WARNING @ compactTract::normInnerProduct(): Negative inner product (" << inProd << ")"
                            << std::endl;
        inProd = 0;
    }
    else if( inProd > 1 )
    {
        if( inProd > 1.0001 )
            std::cerr << std::endl << "WARNING @ compactTract::normInnerProduct(): Bad inner product (" << inProd << ")"
                            << std::endl;
        inProd = 1;
    }

    return inProd;
}
double compactTract::normDotProduct( const compactTractChar &charTract ) const
{
    if( m_tract.size() != charTract.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTract::normInnerProduct(): Tractograms are not of the same size" );
    }
    else if( ( !this->m_normReady ) || ( !charTract.m_normReady ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::normInnerProduct(): one (or both) of the tracts has no available precomputed norm" );
    }
    else if( ( !this->m_thresholded ) || ( !charTract.m_thresholded ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::normInnerProduct(): one (or both) of the tracts has not been thresholded" );
    }
    else if( ( this->m_norm == 0. ) || ( charTract.m_norm == 0. ) )
    {
        std::cerr << "WARNING @ compactTract::normInnerProduct(): At least one of the tractograms is a zero vector, inner product will be set to 0"
                  << std::endl;
        return 0.;
    }
    else if( ( !this->m_inLogUnits ) )
    {
        throw std::runtime_error("ERROR @ compactTract::normInnerProduct(): float tract is not in logaritmic units" );
    }


        // Initialize variables
    double inProd( 0 ), dotprod_sum( 0 );



    const float* piterFloat;
    const unsigned char* piterChar;

    piterFloat = &m_tract.front();
    piterChar = &charTract.m_tract.front();
    size_t tractSize( m_tract.size() );

    for( size_t i = 0; i < tractSize; ++i )
    {
        dotprod_sum += ( *piterFloat++ ) * ( *piterChar++ / 255. );
    }


    inProd = dotprod_sum / ( this->m_norm * charTract.m_norm / 255. );

    if( inProd < 0 )
    {
        if( inProd < -0.0001 )
            std::cerr << std::endl << "WARNING @ compactTract::normInnerProduct(): Negative inner product (" << inProd << ")"
                            << std::endl;
        inProd = 0;
    }
    else if( inProd > 1 )
    {
        if( inProd > 1.0001 )
            std::cerr << std::endl << "WARNING @ compactTract::normInnerProduct(): Bad inner product (" << inProd << ")"
                            << std::endl;
        inProd = 1;
    }

    return inProd;
} // end "normDotProduct()" -----------------------------------------------------------------


double compactTract::getNorm()
{
    if( !this->m_thresholded )
    {
        std::cerr << "ERROR @ compactTract::correlation(): tract has not been thresholded" << std::endl;
        return 0;
    }
    else if( !this->m_inLogUnits )
    {
        std::cerr << "ERROR @ compactTract::correlation(): tract has is not in logarithmic units" << std::endl;
        return 0;
    }
    else
    {
        m_norm = 0;
        for( std::vector< float >::const_iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
            m_norm += ( *iter ) * ( *iter );

        m_norm = sqrt( m_norm );
        m_normReady = true;
        return m_norm;
    }
} // end "getNorm()" -----------------------------------------------------------------


size_t compactTract::bytes() const
{
    size_t elementSize( 0 );
    if( m_tract.empty() )
    {
        float el0( 0 );
        elementSize = sizeof( el0 );
    }
    else
    {
        elementSize = sizeof( m_tract.front() );
    }
    return ( sizeof( *this ) + ( elementSize * m_tract.size() ) ) * CHAR_BIT / 8.;
}

float compactTract::mBytes() const
{
    return bytes() / ( 1024. * 1024. );
}


// "unLog()": transforms the tractogram doing a 10^x exponential
void compactTract::unLog( float logFactor )
{
    if( logFactor == 0 )
    {
        m_inLogUnits = false;
        return;
    }

    if( this->m_thresholded )
    {
        std::cerr << "ERROR @ compactTract::unLog(): tract has been thresholded" << std::endl;
    }
    else if( !this->m_inLogUnits )
    {
        std::cerr << "ERROR @ compactTract::unLog(): tract is not in logarithmic units" << std::endl;
    }
    else
    {
        for( std::vector< float >::iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
            *iter = ( pow( 10., ( ( *iter ) * logFactor ) ) );
        m_inLogUnits = false;
    }
    return;
} // end "unLog()" -----------------------------------------------------------------


// "doLog()": transforms the tractogram doing a base-10 logarithm
void compactTract::doLog( float logFactor )
{
    if( logFactor == 0 )
    {
        m_inLogUnits = true;
        return;
    }

    if( this->m_thresholded )
    {
        std::cerr << "ERROR @ compactTract::doLog(): tract has been thresholded" << std::endl;
    }
    else if( this->m_inLogUnits )
    {
        std::cerr << "ERROR @ compactTract::doLog(): tract is already in logarithmic units" << std::endl;
    }
    else
    {
        for( std::vector< float >::iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
            *iter = ( ( log10( *iter ) ) / logFactor );
        m_inLogUnits = true;
    }
    return;
} // end "doLog()" -----------------------------------------------------------------


// "threshold()": thresholds the input tractogram, if the value of a point is less than the given threshold, it is set to 0
void compactTract::threshold( const float threshold )
{
    if( !this->m_inLogUnits )
    {
        std::cerr << "ERROR @ compactTract::threshold(): tract is not in logarithmic units" << std::endl;
    }
    else if( this->m_thresholded )
    {
        std::cerr << "WARNING @ compactTract::threshold(): tract has already been thresholded" << std::endl;
    }
    else
    {
        if( threshold != 0 )
        {
            for( std::vector< float >::iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
                if( *iter < threshold )
                    *iter = 0.;
        }
        m_thresholded = true;
    }
    return;
} // end "threshold()" -----------------------------------------------------------------


// "add()": sum input tractogram top the tractogram object
void compactTract::add( const compactTract &tractogram )
{
    if( m_tract.size() != tractogram.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTract::add(): Tractograms are not of the same size" );
    }
    else if( this->m_thresholded || tractogram.thresholded() )
    {
        throw std::runtime_error( "WARNING @ compactTract::add(): one (or both) of the tracts has been thresholded" );
    }
    else if( this->m_inLogUnits || tractogram.m_inLogUnits )
    {
        throw std::runtime_error(
         "WARNING @ compactTract::add(): one (or both) of the tracts is in logarithmic units, natural units are necessary for summing operations" );
    }
    else
    {
        for( int i = 0; i < this->m_tract.size(); ++i )
            this->m_tract[i] += tractogram.m_tract[i];
    }
    return;
} // end "add()" -----------------------------------------------------------------

// "divide()": divide tractogram by a given value
void compactTract::divide( const float divisor )
{
    for( std::vector< float >::iterator iter = m_tract.begin(); iter != m_tract.end(); ++iter )
        *iter = *iter / divisor;
    return;
} // end "divide()" -----------------------------------------------------------------


void compactTract::steal( compactTract* const stolen )
{
    m_tract.swap( stolen->m_tract );
    m_norm = stolen->m_norm;
    m_thresholded = stolen->m_thresholded;
    m_normReady = stolen->m_normReady;
    m_inLogUnits = stolen->m_inLogUnits;
} // end "steal()" -----------------------------------------------------------------

// member operators
compactTract& compactTract::operator =( const compactTract &rhs )
{
    m_tract = rhs.m_tract;
    m_norm = rhs.m_norm;
    m_thresholded = rhs.m_thresholded;
    m_normReady = rhs.m_normReady;
    m_inLogUnits = rhs.m_inLogUnits;
    return *this;
}
compactTract& compactTract::operator =( const compactTractChar &rhs )
{
    m_norm = rhs.m_norm / 255.;
    m_thresholded = rhs.m_thresholded;
    m_normReady = rhs.m_normReady;
    m_inLogUnits = true;
    m_tract.clear();
    m_tract.reserve( rhs.m_tract.size() );
    for( size_t i = 0; i < rhs.m_tract.size(); ++i )
    {
        m_tract.push_back( rhs.m_tract[i] / 255.0 );
    }
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
double compactTract::correlation( const compactTract &tractogram ) const
{
    if( m_tract.size() != tractogram.m_tract.size() )
    {
        throw std::runtime_error( "ERROR @ compactTract::correlation(): Tractograms are not of the same size" );
    }
    else if( ( !this->m_normReady ) || ( !tractogram.m_normReady ) )
    {
        throw std::runtime_error(
                        "ERROR @ compactTract::correlation(): one (or both) of the tracts has no available precomputed norm" );
    }
    else if( ( !this->m_thresholded ) || ( !tractogram.m_thresholded ) )
    {
        throw std::runtime_error( "ERROR @ compactTract::correlation(): one (or both) of the tracts has not been thresholded" );
    }
    else if( ( this->m_norm == 0. ) || ( tractogram.m_norm == 0. ) )
    {
        std::cerr
                        << "WARNING @ compactTract::correlation(): At least one of the tractograms is a zero vector, correlation will be set to 0"
                        << std::endl;
        return 0.;
    }
    else if( ( !this->m_inLogUnits ) || ( !tractogram.m_inLogUnits ) )
    {
        throw std::runtime_error( "ERROR @ compactTract::correlation(): one (or both) of the tracts is not in logaritmic units" );
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
                for( std::vector< float >::const_iterator iter( this->m_tract.begin() ); iter != this->m_tract.end(); ++iter )
                    sum1 += *iter; // sum(X)
            }
#pragma omp section
            {
                for( std::vector< float >::const_iterator iter( tractogram.m_tract.begin() ); iter != tractogram.m_tract.end(); ++iter )
                    sum2 += *iter; // sum(Y)
            }
        }

        avr1 = sum1 / this->m_tract.size(); // E(X)
        avr2 = sum2 / tractogram.m_tract.size(); // E(Y)
        var1 = ( ( this->m_norm * this->m_norm ) / this->m_tract.size() ) - ( avr1 * avr1 ); // Var(X) = E(X2)-E2(X)  [E(X2) = ||X||2/N ]
        var2 = ( ( tractogram.m_norm * tractogram.m_norm ) / tractogram.m_tract.size() ) - ( avr2 * avr2 ); // Var(Y) = E(Y2)-E2(Y)

        if( ( var1 == 0. ) || ( var2 == 0 ) )
        {
            std::cerr << "WARNING @ compactTract::correlation(): One (or both) of the tractograms is a non-zero constant vector,"
                      << " correlation will be set to 0" << std::endl;
            return 0.;
        }

        stddev1 = sqrt( var1 ); // sttdev(X) = sqrt(Var(X))
        stddev2 = sqrt( var2 ); // sttdev(Y) = sqrt(Var(Y))

        // Compute Covariance Cov(X,Y) = E( [X-E(X)]*[Y-E(Y)] )
        for( std::vector< float >::const_iterator iter1( this->m_tract.begin() ), iter2( tractogram.m_tract.begin() ); iter1
                        != this->m_tract.end(); ++iter1, ++iter2 )
        {
            cov += ( *iter1 - avr1 ) * ( *iter2 - avr2 );
        }
        corr = cov / ( this->m_tract.size() * stddev1 * stddev2 );

        if( corr < 0 )
        {
            std::cerr << std::endl << "WARNING @ compactTract::correlation(): Negative correlation (" << corr
                            << ") saving it as 0" << std::endl;
            corr = 0;
        }
        else if( corr > 1.0001 )
        {
            std::cerr << std::endl << "WARNING @ compactTract::correlation(): Bad correlation (" << corr << ")" << std::endl;
            corr = 1;
        }

        return corr;
    }
} // end "correlation()" -----------------------------------------------------------------
