#ifndef CNBTREEBUILDER_H
#define CNBTREEBUILDER_H

// parallel execution
#include <omp.h>

// std library
#include <vector>
#include <list>
#include <string>
#include <map>
#include <fstream>
#include <ctime>
#include <climits>
#include <sstream>

// boost library
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

// hClustering
#include "compactTract.h"
#include "compactTractChar.h"
#include "WHcoord.h"
#include "WFileParser.h"
#include "WHtree.h"
#include "protoNode.h"
#include "listedCache.hpp"
#include "vistaManager.h"

typedef enum
{
    TC_GROWOFF, TC_GROWNUM, TC_GROWSIZE
} TC_GROWTYPE;

// typical values:
// m_maxNbDist = 0.2
// m_tractThreshold = 0.4
// m_logFactor = 5.0

// for rand tracts:
// m_maxNbDist = 1.0
// m_tractThreshold = 0.0
// m_logFactor = 0.0


class cnbTreeBuilder
{
public:
#if 1
    explicit cnbTreeBuilder( std::string roiFilename, bool verbose = true):
        m_roiLoaded( false ), m_treeReady( false ), m_logfile( 0 ), m_maxNbDist( 0.2 ), m_tractThreshold( 0.4 ),
        m_logFactor( 5.0 ), m_numComps( 0 ), m_ncHits( 0 ), m_ncMiss( 0 ), m_lcHits( 0 ), m_lcMiss( 0 ), m_verbose( verbose )
    {
        readRoi( roiFilename );
    }
#else
    explicit treeBuilder( std::string roiFilename ):
        m_roiLoaded( false ), m_treeReady( false ), m_logfile( 0 ), m_maxNbDist( 1.0 ), m_tractThreshold( 0.0 ),
                        m_logFactor( 0.0 ), m_numComps( 0 ), m_ncHits( 0 ), m_ncMiss( 0 ), m_lcHits( 0 ), m_lcMiss( 0 )
    {
        readRoi( roiFilename );
    }
#endif
    ~cnbTreeBuilder()
    {
    }

    void log( std::ofstream* const logfile )
    {
        m_logfile = logfile;
    }
    bool ready() const
    {
        return m_roiLoaded;
    }
    void setInputFolder( std::string inputFolder )
    {
        m_inputFolder = inputFolder;
    }
    void setOutputFolder( std::string outputFolder )
    {
        m_outputFolder = outputFolder;
    }
    void setMaxNbDist( dist_t maxNbDist )
    {
        m_maxNbDist = maxNbDist;
    }


    size_t roiSize() const
    {
        return m_roi.size();
    }

    bool readRoi( const std::string roiFilename );

    void buildCentroid( const unsigned int nbLevel, const float memory, const std::string meanTractFolder, bool keepDiscarded, TC_GROWTYPE growType, size_t baseSize = 0 );

    void buildC2( const unsigned int nbLevel, const float memory, const std::string meanTractFolder, bool keepDiscarded, TC_GROWTYPE growType, size_t baseSize = 0 );


    void writeTree() const;

private:
    bool m_verbose;
    bool m_roiLoaded;
    bool m_treeReady;
    WHtree m_tree;
    std::vector< WHcoord > m_roi;
    std::vector< double > m_leafNorms;
    std::vector< double > m_nodeNorms;
    size_t m_numComps;

    WHcoord m_datasetSize;
    HC_GRID m_datasetGrid;
    std::ofstream *m_logfile;
    dist_t m_maxNbDist;
    float m_tractThreshold;
    float m_logFactor;
    std::string m_inputFolder;
    std::string m_outputFolder;

    volatile size_t m_ncHits;
    volatile size_t m_ncMiss;
    volatile size_t m_lcHits;
    volatile size_t m_lcMiss;


    protoNode* getProtoNode( const nodeID_t &thisNode, std::vector< protoNode > &protoLeaves,
                    std::vector< protoNode > &protoNodes ) const;
    WHnode* fetchNode( const nodeID_t &thisNode, std::vector< WHnode > &leaves, std::vector< WHnode > &nodes ) const;

    void computeNorms();

    std::list< WHcoord > initialize( const unsigned int &nbLevel, const size_t &cacheSize, std::vector< protoNode > &protoLeaves );

    bool scanNbs( const size_t &currentSeedID, const compactTractChar* const currentTract, std::map< size_t, dist_t > *nbLeaves,
                    std::vector< WHcoord > &nbCoords, std::vector< protoNode > &protoLeaves, std::map< WHcoord, size_t > &roimap,
                    listedCache< compactTractChar > *cache );

    compactTract* loadNodeTract( const size_t nodeID,
                                 const vistaManager* const nodeMngr,
                                 listedCache< compactTract > *nodesCache,
                                 const volatile size_t &natCounter );

    compactTractChar* loadLeafTract( const size_t leafID,
                                     const vistaManager* const leafMngr,
                                     listedCache< compactTractChar > *leavesCache);

    void* loadTract( const nodeID_t nodeID,
                     const vistaManager* const leafMngr,
                     const vistaManager* const nodeMngr,
                     listedCache< compactTractChar > *leavesCache,
                     listedCache< compactTract > *nodesCache,
                     const volatile size_t &natCounter );


    void writeTract( const size_t &nodeID, const compactTract &tract, const vistaManager* const vistaMngr,
                    volatile size_t* threadcounter ) const;
    void deleteTract( const size_t &nodeID, const vistaManager* const vistaMngr, volatile size_t* threadcounter ) const;

    void writeBases( const std::vector< size_t > &baseNodes, const std::string filename ) const;

};

#endif  // CNBTREEBUILDER_H
