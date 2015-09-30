

// std library
#include <vector>
#include <list>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>


// parallel execution
#include <omp.h>


#include "surfProjecter.h"
#include "WStringUtils.h"

surfProjecter::surfProjecter( std::string roiFilename, std::string surfFfilename, bool verbose ):
    m_loaded( false ), m_gaussKernel(false), m_kernelRadius(0), m_sigma(1), m_verbose( verbose )
{
    fileManagerFactory fMFtestFormat;
    m_niftiMode = fMFtestFormat.isNifti();

    RoiLoader loader( m_niftiMode, true );
    bool roiLoaded( loader.readRoi( roiFilename, &m_roiDatasetGrid, &m_roiDatasetSize, &m_numStreamlines, &m_roiCoords, &m_trackids ) );
    size_t dummySurfStreamlines;
    std::vector< size_t > dummySurfTrackids;
    bool surfLoaded( loader.readRoi( surfFfilename, &m_surfDatasetGrid, &m_surfDatasetSize, &dummySurfStreamlines, &m_surfCoords, &dummySurfTrackids ) );
    m_loaded = ( roiLoaded && surfLoaded );

    if( m_roiDatasetSize != m_surfDatasetSize )
    {
        throw std::runtime_error( "Roi and surface dataset sizes do not match" );
    }
    if( m_roiDatasetGrid != m_surfDatasetGrid )
    {
        throw std::runtime_error( "Roi and surface dataset grids do not match" );
    }
    if( m_roiDatasetGrid != HC_VISTA && m_roiDatasetGrid != HC_NIFTI )
    {
        throw std::runtime_error( "Roi dataset grid is neither nifti nor vista" );
    }

    if( verbose )
    {
        std::cout <<"Roi file has " << m_roiCoords.size() << " seeds." << std::endl;
        std::cout <<"Surf file has " << m_surfCoords.size() << " vertices." << std::endl;
    }
} // end constructor



void surfProjecter::matchCoordsKernel()
{
    m_coordMatch.clear();
    m_matchDists.clear();

    if ( !m_loaded )
    {
        std::cerr << "ERROR @ seedMatcher::matchCoords(): coordinates are not loaded" << std::endl;
        return;
    }

    std::cout << "Using a kernel size of: " << m_kernelRadius << std::endl;

    m_coordMatch.resize( m_surfCoords.size() );
    m_matchDists.resize( m_surfCoords.size() );

    time_t lastTime( time( NULL ) ), loopStart( time( NULL ) ); // time object
    size_t doneCount(0);


    // get for every surface vertex seed the roi seeds within a distance
    for ( size_t i = 0; i < m_surfCoords.size(); ++i )
    {
        WHcoord thisSurfCoord( m_surfCoords[i] );

        for ( size_t j = 0; j < m_roiCoords.size(); ++j )
        {
            // center the roi coord so that we calculate physical distances to the center of the voxel
            WHcoord thisCenteredRoiCoord( m_roiCoords[j] );
//            thisCenteredRoiCoord.m_x += 0.5;
//            thisCenteredRoiCoord.m_y += 0.5;
//            thisCenteredRoiCoord.m_z += 0.5;

            float thisDist( thisSurfCoord.getPhysDist( thisCenteredRoiCoord ) );
            if ( thisDist <= m_kernelRadius )
            {
                m_coordMatch[i].push_back( m_roiCoords[j] );
                m_matchDists[i].push_back( std::make_pair( thisDist, j ) );
            }

        } // end inner for

        #pragma omp atomic
        ++doneCount;

        if( m_verbose )
        {
            #pragma omp single
            {
                time_t currentTime( time( NULL ) );
                if( currentTime - lastTime > 1 )
                {
                    lastTime = currentTime;
                    float progress =  doneCount * 100. / ( m_surfCoords.size() );
                    size_t elapsedTime( difftime( currentTime, loopStart ) );
                    std::stringstream message;
                    message << "\r" << static_cast<int>( progress ) << " % of distances computed (";
                    message << doneCount << " surf vertex assigned). ";
                    message << "Elapsed: ";
                    message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                    message << ( elapsedTime % 3600 ) % 60 << "\". ";
                    if( progress > 0 )
                    {
                        size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                        message << "Remaining: ";
                        message << expectedRemain / 3600 << "h ";
                        message << ( expectedRemain % 3600 ) / 60 << "' ";
                        message << ( expectedRemain % 3600 ) % 60 << "\". ";
                    }
                    std::cout << message.str() <<std::flush;
                }
            } // end pragma single
        } // end verbose
    } // end parallel for

    std::cout << "\r 100 % of distances computed (";
    std::cout << doneCount << " surf vertex asssigned).      " << std::endl;

    // order the resulting vectors, so that the closest rois will be at the top of the list for each vertex
    for ( size_t i = 0; i < m_surfCoords.size(); ++i )
    {
        std::sort(m_coordMatch[i].begin(), m_coordMatch[i].end() );
        std::sort(m_matchDists[i].begin(), m_matchDists[i].end() );
    }

    return;
} // end seedMatcher::matchCoords() -------------------------------------------------------------------------------------




void surfProjecter::matchCoordsNearestNb()
{
    m_coordMatch.clear();
    m_matchDists.clear();

    if ( !m_loaded )
    {
        std::cerr << "ERROR @ seedMatcher::matchCoordsSimple(): coordinates are not loaded" << std::endl;
        return;
    }

    std::vector< WHcoord > singleRoiMatchCoords;
    singleRoiMatchCoords.resize( m_roiCoords.size() );

    std::vector< float > singleRoiMatchDists;
    singleRoiMatchDists.resize( m_roiCoords.size(), -1 );

    time_t lastTime( time( NULL ) ), loopStart( time( NULL ) ); // time object
    size_t doneCount(0);

    // get for every difussion roi seed the closest surface vertex
    for ( size_t i = 0; i < m_roiCoords.size(); ++i )
    {
        // center the roi coord so that we calculate physical distances to the center of the voxel
        WHcoord thisCenteredRoiCoord( m_roiCoords[i] );
//        thisCenteredRoiCoord.m_x += 0.5;
//        thisCenteredRoiCoord.m_y += 0.5;
//        thisCenteredRoiCoord.m_z += 0.5;

        float closestDist( 999 );
        WHcoord closestSurfPoint;

        for ( size_t j = 0; j < m_surfCoords.size(); ++j )
        {
            WHcoord thisSurfCoord( m_surfCoords[j] );

            float thisDist( thisCenteredRoiCoord.getPhysDist( thisSurfCoord ) );
            if ( thisDist < closestDist )
            {
                closestSurfPoint = thisSurfCoord;
                closestDist = thisDist;
            }
        } // end inner for
        singleRoiMatchCoords[i] = closestSurfPoint;
        singleRoiMatchDists[i] = closestDist;

        #pragma omp atomic
        ++doneCount;

        if( m_verbose )
        {
            #pragma omp single
            {
                time_t currentTime( time( NULL ) );
                if( currentTime - lastTime > 1 )
                {
                    lastTime = currentTime;
                    float progress =  doneCount * 100. / ( m_roiCoords.size() );
                    size_t elapsedTime( difftime( currentTime, loopStart ) );
                    std::stringstream message;
                    message << "\r" << static_cast<int>( progress ) << " % of distances computed (";
                    message << doneCount << " roi seeds asigned). ";
                    message << "Elapsed: ";
                    message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                    message << ( elapsedTime % 3600 ) % 60 << "\". ";
                    if( progress > 0 )
                    {
                        size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                        message << "Remaining: ";
                        message << expectedRemain / 3600 << "h ";
                        message << ( expectedRemain % 3600 ) / 60 << "' ";
                        message << ( expectedRemain % 3600 ) % 60 << "\". ";
                    }
                    std::cout << message.str() <<std::flush;
                }
            } // end pragma single
        }
    } // end parallel for
    std::cout << "\r 100 % of distances computed (";
    std::cout << doneCount << " roi seeds asigned).      " << std::endl;

    m_coordMatch.resize( m_surfCoords.size() );
    m_matchDists.resize( m_surfCoords.size() );

    // create a vector of lists indicating the corresponding roi seeds for every vertex
    for ( size_t i = 0; i < m_roiCoords.size(); ++i )
    {
        std::vector< WHcoord >::const_iterator surfFindIter( std::find( m_surfCoords.begin(), m_surfCoords.end(), singleRoiMatchCoords[i] ) );
        if ( surfFindIter == m_surfCoords.end() )
        {
            std::cerr << "ERROR @ seedMatcher::matchCoordsSimple(): surf coord not found" << std::endl;
            return;
        }
        size_t surfPos( surfFindIter - m_surfCoords.begin() );
        m_coordMatch[surfPos].push_back( m_roiCoords[i] );
        m_matchDists[surfPos].push_back( std::make_pair( singleRoiMatchDists[i], i ) );
    }

    // order the resulting vectors, so that the closest rois will be at the top of the list for each vertex
    for ( size_t i = 0; i < m_surfCoords.size(); ++i )
    {
        std::sort(m_coordMatch[i].begin(), m_coordMatch[i].end() );
        std::sort(m_matchDists[i].begin(), m_matchDists[i].end() );
    }
    return;
} // end seedMatcher::matchCoordsSimple() -------------------------------------------------------------------------------------

void surfProjecter::writeMeanTracts( std::string singleTractFolder, std::string outputTractFolder, bool useFloat, bool doZip )
{
    if ( m_matchDists.empty() )
    {
        std::cerr << "ERROR @ seedMatcher::writeMeanTracts(): matching was not previously performed" << std::endl;
        return;
    }
    std::cout << "Kernel radius: " << m_kernelRadius << std::endl;


    std::cout << "Kernel type: ";
    if( !m_gaussKernel )
    {
        if( m_kernelRadius == 0 )
        {
            std::cout << "Nearest Neighbor" << std::endl;
        }
        else
        {
            std::cout << "Mean (square function)" << std::endl;
        }
    }
    else
    {
        std::cout << "Gaussian. sigma = "<< m_sigma << std::endl;
    }



    m_leafTractFolder = singleTractFolder;
    fileManagerFactory meanFileMF( outputTractFolder );
    fileManager& fileMean( meanFileMF.getFM() );
    if( useFloat )
    {
        fileMean.writeInFloat();
    }
    else
    {
        fileMean.writeInChar();
    }

    if( doZip )
    {
        fileMean.storeZipped();
    }
    else
    {
        fileMean.storeUnzipped();
    }

    time_t lastTime( time( NULL ) ), loopStart( time( NULL ) ); // time object
    size_t doneCount(0), emptyCount(0);

//    #pragma omp parallel for
    for ( size_t i = 0; i < m_surfCoords.size(); ++i )
    {
        #pragma omp atomic
        ++doneCount;

        if(m_matchDists[i].empty() )
        {
            std::pair< float, size_t> surfNearest( matchSurfNearest( m_surfCoords[i] ) );

            m_matchDists[i].push_back( surfNearest );
            m_coordMatch[i].push_back( m_roiCoords[surfNearest.second] );

//            std::cerr << std::endl <<"WARNING, surface point "<< i <<" did not get any seeds. assigned closest at dist: "<< surfNearest.first <<std::endl;
//            std::cerr << "\r" << static_cast<int>( doneCount * 100. / ( m_surfCoordList.size() ) ) << " % of mean tracts computed (" << doneCount << "). ";
            #pragma omp atomic
            ++emptyCount;
        }
        else
        {
//            std::cout << "surf voxel "<< i << " ("<< m_surfCoordList[i] << ") got "<< m_matchDists[i].size() <<" roi seeds"<< std::endl;
        }
        compactTract thisMeanTract( getMeanTract( m_matchDists[i] ) );

        fileMean.writeNodeTract( i, thisMeanTract );

        if( m_verbose )
        {
            #pragma omp single
            {
                time_t currentTime( time( NULL ) );
                if( currentTime - lastTime > 1 )
                {
                    lastTime = currentTime;
                    float progress =  doneCount * 100. / ( m_surfCoords.size() );
                    size_t elapsedTime( difftime( currentTime, loopStart ) );
                    std::stringstream message;
                    message << "\r" << static_cast<int>( progress ) << " % of mean tracts computed (";
                    message << doneCount << "). ";
                    message << "Elapsed: ";
                    message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                    message << ( elapsedTime % 3600 ) % 60 << "\". ";
                    if( progress > 0 )
                    {
                        size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                        message << "Remaining: ";
                        message << expectedRemain / 3600 << "h ";
                        message << ( expectedRemain % 3600 ) / 60 << "' ";
                        message << ( expectedRemain % 3600 ) % 60 << "\". ";
                    }
                    std::cout << message.str() <<std::flush;
                }
            }
        }
    }// end outer for
    std::cout << "\r 100 % of mean tracts computed (";
    std::cout << doneCount << ").      " << std::endl;
    std::cout << "WARNING, "<< emptyCount << " surface points did not get any seeds"<<std::endl;

} // end seedMatcher::writeMeanTracts() -------------------------------------------------------------------------------------


compactTract surfProjecter::getMeanTract( std::vector< std::pair< float, size_t > > iIDs ) const
{
    if( m_leafTractFolder.empty() )
    {
        throw std::runtime_error( "Leaf tracts folder has not been specified" );
    }

    compactTract sumTract;
    if( iIDs.empty() )
    {
        std::cerr <<  "ERROR @ seedMatcher::writeMeanTracts() ID list is empty" << std::endl;
        return sumTract;
    }


    //create file reader class
    fileManagerFactory fileMF( m_leafTractFolder );
    fileManager& fileLeaves( fileMF.getFM() );
    fileLeaves.readAsLog();
    fileLeaves.readAsUnThres();

    fileLeaves.readLeafTract( iIDs[0].second , m_trackids, m_roiCoords, &sumTract );

    if( iIDs.size() == 1 )
    {
        return sumTract;
    }

    float tempDist( iIDs.front().first );
    float tempWeight(0);


    double weightsum(0);
    float gcoeff(0), gexp(0);

    if( !m_gaussKernel )
    {
        tempWeight = 1;
    }
    else
    {
        gcoeff = 1.0 / (m_sigma * std::sqrt( 2.0 * 3.1415926536 ) );
        gexp = -1.0 / ( 2 * m_sigma * m_sigma );
        tempWeight = gcoeff * std::exp( gexp * tempDist * tempDist );
    }




    sumTract.unLog( m_logFactor ); //to sum tractograms they must be in natrual units
    sumTract.mult( tempWeight );
    weightsum += tempWeight;

    for( size_t j = 1; j < iIDs.size(); ++j )
    {
        compactTract tempTract;
        size_t tempID = iIDs[j].second;
        fileLeaves.readLeafTract( tempID, m_trackids, m_roiCoords, &tempTract );
        tempTract.unLog( m_logFactor ); //to sum tractograms they must be in natrual units
        float weight(0);

        if( !m_gaussKernel )
        {
            weight = 1;
        }
        else
        {
            float dist( iIDs[j].first );
            weight = gcoeff * std::exp( gexp * dist * dist );
        }
        tempTract.mult( weight );
        sumTract.add( tempTract );
        weightsum += weight;


    }
    sumTract.divide( weightsum ); // we must divide the sum tract by the the sum of weights in order to normalize
    sumTract.doLog( m_logFactor ); // we want to view it in logarithmic units

    return sumTract;
} // end "getMeanTract()" -----------------------------------------------------------------

std::pair< float, size_t > surfProjecter::matchSurfNearest( WHcoord surfCoord )
{
    float closestDist( 999 );
    size_t closestRoiPointID;

    // get for every difussion roi seed the closest surface vertex
    for ( size_t i = 0; i < m_roiCoords.size(); ++i )
    {
        // center the roi coord so that we calculate physical distances to the center of the voxel
        WHcoord thisCenteredRoiCoord( m_roiCoords[i] );
//        thisCenteredRoiCoord.m_x += 0.5;
//        thisCenteredRoiCoord.m_y += 0.5;
//        thisCenteredRoiCoord.m_z += 0.5;

        float thisDist( thisCenteredRoiCoord.getPhysDist( surfCoord ) );
        if ( thisDist < closestDist )
        {
            closestRoiPointID = i;
            closestDist = thisDist;
        }
    }
    return std::make_pair( closestDist, closestRoiPointID );
} // end "matchSurfNearest()" -----------------------------------------------------------------
