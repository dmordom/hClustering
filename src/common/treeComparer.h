#ifndef TREECOMPARER_H
#define TREECOMPARER_H

// parallel execution
#include <omp.h>

// std library
#include <vector>
#include <list>
#include <string>
#include <utility>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <sstream>

// boost library
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

// hClustering
#include "compactTract.h"
#include "WHcoord.h"
#include "distBlock.h"
#include "WHtree.h"
#include "listedCache.hpp"
#include "vistaManager.h"
#include "treeManager.h"
#include "WHtreeProcesser.h"

const bool TREE1 = true;
const bool TREE2 = false;

class treeComparer
{
public:
    treeComparer( WHtree* const tree1, WHtree* const tree2, bool verbose = true) :
        m_tree1( *tree1 ), m_tree2( *tree2 ), m_maxPhysDist( 0 ), m_tractThreshold( 0.4 ), m_logFactor( 5.0 ), m_logfile( 0 ),
        m_realBaseNodes( false ), m_coordsFromFile( false ), m_meanTractsFromFile( false ), m_verbose( verbose )
    {
        fetchBaseNodes( false );
        m_initialSizes = std::make_pair(m_tree1.getNumLeaves(), m_tree2.getNumLeaves());
    }
    treeComparer( WHtree* const tree1, WHtree* const tree2, treeComparer &comparer ) :
      m_verbose( comparer.m_verbose ),
      m_tree1(*tree1),
      m_tree2(*tree2),
      m_singleTractFolder1(comparer.m_singleTractFolder1),
      m_singleTractFolder2(comparer.m_singleTractFolder2),
      m_meanTractFolder1(comparer.m_meanTractFolder1),
      m_meanTractFolder2(comparer.m_meanTractFolder2),
      m_maxPhysDist(comparer.m_maxPhysDist),
      m_tractThreshold(comparer.m_tractThreshold),
      m_logFactor(comparer.m_logFactor),
      m_logfile(comparer.m_logfile),
      m_baseNodes1(comparer.m_baseNodes1),
      m_originalBaseNodes1(comparer.m_originalBaseNodes1),
      m_baseCoords1(comparer.m_baseCoords1),
      m_noiseLevels1(comparer.m_noiseLevels1),
      m_baseNodes2(comparer.m_baseNodes2),
      m_originalBaseNodes2(comparer.m_originalBaseNodes2),
      m_baseCoords2(comparer.m_baseCoords2),
      m_noiseLevels2(comparer.m_noiseLevels2),
      m_initialSizes(comparer.m_initialSizes),
      m_baseDistMatrix(comparer.m_baseDistMatrix),
      m_realBaseNodes(comparer.m_realBaseNodes),
      m_coordsFromFile(comparer.m_coordsFromFile),
      m_meanTractsFromFile(comparer.m_meanTractsFromFile),
      m_fullCorrespondence(comparer.m_fullCorrespondence),
      m_newCorrespondence(comparer.m_newCorrespondence),
      m_newCorrespReverse(comparer.m_newCorrespReverse),
      m_correspDistances(comparer.m_correspDistances)
    {
    }

    ~treeComparer()
    {
    }

    bool setSingleTractFolder1( std::string singleTractFolder1 )
    {
        m_singleTractFolder1 = singleTractFolder1;
    }

    bool setSingleTractFolder2( std::string singleTractFolder2 )
    {
        m_singleTractFolder2 = singleTractFolder2;
    }

    bool setMeanTractFolder1( std::string meanTractFolder1 )
    {
        m_meanTractFolder1 = meanTractFolder1;
    }

    bool setMeanTractFolder2( std::string meanTractFolder2 )
    {
        m_meanTractFolder2 = meanTractFolder2;
    }

    void setMaxPhysDist( size_t maxPhysDist )
    {
        m_maxPhysDist = maxPhysDist;
    }
    void log( std::ofstream* const logfile )
    {
        m_logfile = logfile;
    }
    void setCoordsFromFile( bool coordsFromFile )
    {
        m_coordsFromFile = coordsFromFile;
    }
    void setMeanTractsFromFile( bool meanTractsFromFile )
    {
        m_meanTractsFromFile = meanTractsFromFile;
    }

    void getBaseDistMatrix();
    void readBaseDistMatrix( std::string matrixFilename );
    void writeBaseDistMatrix( std::string matrixFilename );


    std::pair< std::pair< float, float >, std::pair< float, float > > doTcpcc() const;
    std::pair< float, float > simpleTriplets( size_t sampleFreq = 1 ) const;

    bool leafCorrespondence();
    void simpleCorrespondence( float threshold, bool redoCoords = true );
    void mergedUpCorrespondence();
    void randomCorrespondence();
    void writeCorrespondence( std::string outFilename );
    void writeFullCorrespondence( std::string outFilename );


    std::vector<float> rateCorrespondence();
    std::pair< float, float > applyNoiseBaseline( const float addedNoise = 0 );


    std::string reportBaseNodes() const;
    bool areRealBaseNodes() const
    {
        return m_realBaseNodes;
    }

    void reduceBaseNodes( size_t number );

private:
    bool m_verbose;
    WHtree &m_tree1;
    WHtree &m_tree2;

    std::string m_singleTractFolder1;
    std::string m_singleTractFolder2;
    std::string m_meanTractFolder1;
    std::string m_meanTractFolder2;
    float m_maxPhysDist;
    float m_tractThreshold;
    float m_logFactor;
    std::ofstream *m_logfile;

    std::vector< size_t > m_baseNodes1;
    std::vector< size_t > m_originalBaseNodes1;
    std::vector< WHcoord > m_baseCoords1;
    std::vector< float > m_noiseLevels1;

    std::vector< size_t > m_baseNodes2;
    std::vector< size_t > m_originalBaseNodes2;
    std::vector< WHcoord > m_baseCoords2;
    std::vector< float > m_noiseLevels2;

    std::pair< size_t, size_t > m_initialSizes;
    std::vector< std::vector< dist_t > > m_baseDistMatrix;
    bool m_realBaseNodes;
    bool m_coordsFromFile;
    bool m_meanTractsFromFile;


    std::vector< size_t > m_fullCorrespondence;
    std::vector< size_t > m_newCorrespondence;
    std::vector< size_t > m_newCorrespReverse;
    std::vector< std::pair< float, float > > m_correspDistances; //correspondence distances (tractogram distance, euclidean distance of cluster centres)


    bool fetchBaseNodes(bool doGetCoords = true);
    void getBaseCoords();

    bool mergeNodes( const WHtree &tree1, std::vector< size_t > &merged1, std::vector< std::vector< size_t > > &containedBnodes1,
                    std::vector< std::vector< size_t > > &targets1, std::vector< size_t > &merged2,
                    std::vector< std::vector< size_t > > &containedBnodes2, std::vector< std::vector< size_t > > &targets2 );

    size_t noiseBaseline( const bool treeCode, const float addedNoise );
    size_t findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t > &baseNodes ) const;



    void getTreeRouteCorresp3G( const WHtree &tree1, const std::string &meanTractFolder1, std::vector< nodeID_t > &baseNodes1,
                    const WHtree &tree2, const std::string &meanTractFolder2, std::vector< nodeID_t > &baseNodes2, std::vector<
                                    std::vector< std::pair< nodeID_t, dist_t > > > &correspVector );
    void getTreeRouteCorresp1G( const WHtree &tree1, const std::string &meanTractFolder1, std::vector< nodeID_t > &baseNodes1,
                    const WHtree &tree2, const std::string &meanTractFolder2, std::vector< nodeID_t > &baseNodes2, std::vector<
                                    std::vector< std::pair< nodeID_t, dist_t > > > &correspVector );
    long int distance2kids( const WHcoord &thatCoord, const compactTract &thatTract, const nodeID_t &thisNode,
                    const WHtree &thisTree, const vistaManager &nodeVista, std::vector< nodeID_t > &baseNodes,
                    std::vector< std::pair< nodeID_t, dist_t > > &kidDists );
    float minPhysDist( const WHcoord thisCoord, std::vector< WHcoord > nodeCoords );
};

#endif  // TREECOMPARER_H
