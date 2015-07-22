#include "tract_dist.h"

/* "correlate()":
Calculates the correlation coefficient between two tractograms. (tractograms must be in logarithmic units and thresholded)
The correlation is calculated and normalized as:

 considering X = tractogram1 ; Y = tractogram2

          Cov(X,Y)                E( [X-E(X)]*[Y-E(Y)] )                sum( [X-E(X)]*[Y-E(Y)] )
 CC = --------------------- = ----------------------------- = --------------------------------------------
       stddev(X)*stddev(Y)     sqrt[Var(X)] * sqrt[Var(Y)]     N * sqrt[E(X)-E(X)] * sqrt[E(Y)-E(Y)]


This program allows you to calculate the similarity between tractograms.
Other possible similarity indices would be:

* Correlation   "signal product"  co = 1/N*sum(x*y)
* Cross-correlation coefficient   cc = 1/N*sum{(x-mean_x)*(y-mean_y)}/(stddev_x*stddev_y)
* Euclidian Distance              eu = sqrt[sum{(x-y)*(x-y)}]
* Mutual Information              mi = sum [ P(x,y)*log2{P(x,y)/(P(x)*P(y))} ]
* Covariance                      covar = 1/N*sum{(x-mean_x)*(y-mean_y)}= 1/N*sum(x*y) -mean_x*mean_y
*/


double tract_distance(std::vector<float> &tractogram1, std::vector<float> &tractogram2) {
    double distance(1-vectprod(tractogram1, tractogram2));
    return (distance);
}



double vectprod(std::vector<float> &tractogram1, std::vector<float> &tractogram2) {


    // Initialize variables
    double sqr_sum1(0), sqr_sum2(0), corr(0), cov(0);

    std::vector<float>::size_type tract_size(tractogram1.size());

    if ( tract_size != tractogram2.size() )
        throw std::runtime_error ("ERROR: Tractograms are not of the same size");

    // Compute average and stddev;  E(X) = sum(x)/N;  Sttdev = sqrt[E(X)-E(X)]
    // Compute Covariance Cov(X,Y) = E( [X-E(X)]*[Y-E(Y)] )


    // Compute average and stddev;  E(X) = sum(x)/N;  Sttdev = sqrt[E(X)-E(X)]

    for(std::vector<float>::const_iterator iter1=tractogram1.begin(), iter2=tractogram2.begin() ; iter1!= tractogram1.end(); ++iter1, ++iter2){
        sqr_sum1 += (*iter1) * (*iter1);
        sqr_sum2 += (*iter2) * (*iter2);
        cov += (*iter1)*(*iter2);
    }


    if ( (sqr_sum1==0.)  || (sqr_sum2==0. ) ) {
        std::cerr <<"WARNING [correlate()]: At least one of the tractograms is a zero vector, correlation will be set to 0"<< std::endl;
        return 0.;
        }


    corr = cov / sqrt(sqr_sum1*sqr_sum2);

    if (corr<0)
        std::cerr <<std::endl<<"Negative correlation: "<< corr << std::endl;
    if (corr>1)
        std::cerr <<std::endl<<"Bad correlation: "<< corr << std::endl;


     //Debug
//**    std::cout <<"Tractogram1 average: "<< avr1 <<" Std dev: "<< stddev1 << std::endl;
//**    std::cout <<"Tractogram2 average: "<< avr2 <<" Std dev: "<< stddev2 << std::endl;
//**    std::cout <<"correlation = "<<corr << std::endl;

    return corr;

}// end "vectprod()" -----------------------------------------------------------------

float correlate(std::vector<float> &tractogram1, std::vector<float> &tractogram2) {


    // Initialize variables
    double sum1(0), sqr_sum1(0), avr1(0), var1(0), stddev1(0);
    double sum2(0), sqr_sum2(0), avr2(0), var2(0), stddev2(0);
    double corr(0), cov(0);

    std::vector<float>::const_iterator iter1, iter2;
    std::vector<float>::size_type tract_size(tractogram1.size());

    if ( tract_size != tractogram2.size() )
        throw std::runtime_error ("ERROR: Tractograms are not of the same size");

    // Compute average and stddev;  E(X) = sum(x)/N;  Sttdev = sqrt[E(X)-E(X)]
    for(iter1=tractogram1.begin(), iter2=tractogram2.begin() ; iter1!= tractogram1.end(); ++iter1, ++iter2){
        sum1 += *iter1; // sum(X)
        sum2 += *iter2; // sum(Y)
        sqr_sum1 += (*iter1) * (*iter1); // sum(X)
        sqr_sum2 += (*iter2) * (*iter2) ; // sum(Y)
    }

    avr1 = sum1 / tract_size; // E(X)
    avr2 = sum2 / tract_size; // E(Y)
    var1 = (sqr_sum1/tract_size) - (avr1 * avr1) ; // Var(X) = E(X2)-E(X)
    var2 = (sqr_sum2/tract_size) - (avr2 * avr2) ; // Var(Y) = E(Y2)-E(Y)

    if ((var1==0. && avr1==0.) || (var2==0. && avr2==0.)) {
        std::cerr <<"WARNING [correlate()]: One of the tractograms is a zero vector, correlation will be set to 0"<< std::endl;
        return 0.;
        }
    if ((var1==0. && avr1!=0.) && (var2==0. && avr2!=0.)){
        std::cerr <<"WARNING [correlate()]: Both tractograms are a non-zero constant vector, correlation will be set to 1"<< std::endl;
        return 1.;
    }
    if ((var1==0. && avr1!=0.) || (var2==0. && avr2!=0.)){
        std::cerr <<"WARNING [correlate()]: One of the tractograms is a non-zero constant vector, correlation will be set to 0"<< std::endl;
        return 0.;
    }

    stddev1= sqrt(var1); // sttdev(X) = sqrt(Var(X))
    stddev2= sqrt(var2); // sttdev(Y) = sqrt(Var(Y))

    // Compute Covariance Cov(X,Y) = E( [X-E(X)]*[Y-E(Y)] )
    for(iter1=tractogram1.begin(), iter2=tractogram2.begin() ; iter1!= tractogram1.end(); ++iter1, ++iter2){
        cov += (*iter1-avr1)*(*iter2-avr2);
    }
    corr = cov / (tract_size * stddev1 * stddev2);

     //Debug
//**    std::cout <<"Tractogram1 average: "<< avr1 <<" Std dev: "<< stddev1 << std::endl;
//**    std::cout <<"Tractogram2 average: "<< avr2 <<" Std dev: "<< stddev2 << std::endl;
//**    std::cout <<"correlation = "<<corr << std::endl;

    // FIX FOR CONSTANT VECTORS


    if (corr<0) {
        //std::cerr <<std::endl<<"Negative correlation: "<< corr <<"saving it as 0"<< std::endl;
        corr = 0;
    }
    if (corr>1)
        std::cerr <<std::endl<<"Bad correlation (should not happen): "<< corr << std::endl;


     //Debug
//**    std::cout <<"Tractogram1 average: "<< avr1 <<" Std dev: "<< stddev1 << std::endl;
//**    std::cout <<"Tractogram2 average: "<< avr2 <<" Std dev: "<< stddev2 << std::endl;
//**    std::cout <<"correlation = "<<corr << std::endl;

    return corr;

}// end "old_correlate()" -----------------------------------------------------------------


