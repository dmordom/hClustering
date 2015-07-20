#ifndef RANDCNBTREEBUILDER_H
#define RANDCNBTREEBUILDER_H

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
#include "WHcoord.h"
#include "distBlock.h"
#include "WFileParser.h"
#include "WHtree.h"
#include "protoNode.h"
#include "vistaManager.h"
#include "cnbTreeBuilder.h"


// typical values:
// m_maxNbDist = 0.2
// m_tractThreshold = 0.4
// m_logFactor = 5.0

// for rand tracts:
// m_maxNbDist = 1.0
// m_tractThreshold = 0.0
// m_logFactor = 0.0


class randCnbTreeBuilder
{
public:

    explicit randCnbTreeBuilder( std::string roiFilename, bool verbose = true ):
        m_roiLoaded( false ), m_treeReady( false ), m_logfile( 0 ), m_maxNbDist( 1.0 ), m_numComps( 0 ), m_verbose( verbose )
    {
        readRoi( roiFilename );
    }

    ~randCnbTreeBuilder()
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


    size_t roiSize() const
    {
        return m_roi.size();
    }

    bool readRoi( const std::string roiFilename );

    void buildCentroiRand( const unsigned int nbLevel, const float memory, bool keepDiscarded, TC_GROWTYPE growType, size_t baseSize = 0);

    void writeTree() const;

private:
    bool m_verbose;
    bool m_roiLoaded;
    bool m_treeReady;
    WHtree m_tree;
    std::vector< WHcoord > m_roi;
    std::vector< compactTract > m_leafTracts;

    size_t m_numComps;

    WHcoord m_datasetSize;
    HC_GRID m_datasetGrid;
    std::ofstream *m_logfile;
    dist_t m_maxNbDist;
    std::string m_inputFolder;
    std::string m_outputFolder;

    protoNode* getProtoNode( const nodeID_t &thisNode, std::vector< protoNode > &protoLeaves,
                    std::vector< protoNode > &protoNodes ) const;
    WHnode* fetchNode( const nodeID_t &thisNode, std::vector< WHnode > &leaves, std::vector< WHnode > &nodes ) const;

    void loadTracts();

    std::list< WHcoord > initialize( const unsigned int &nbLevel, std::vector< protoNode > &protoLeaves );

    bool scanNbs( const size_t &currentSeedID, const compactTract &currentTract, std::map< size_t, dist_t > *nbLeaves,
                    std::vector< WHcoord > &nbCoords, std::vector< protoNode > &protoLeaves, std::map< WHcoord, size_t > &roimap);

    void writeBases( const std::vector< size_t > &baseNodes, const std::string filename ) const;


};

#endif  // RANDCNBTREEBUILDER_H
