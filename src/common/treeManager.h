#ifndef TREEMANAGER_H
#define TREEMANAGER_H

// parallel execution
#include <omp.h>

// std library
#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <ctime>

// boost library
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

// hClustering
#include "compactTract.h"
#include "WHcoord.h"
#include "distBlock.h"
#include "WHtree.h"
#include "listedCache.hpp"
#include "vistaManager.h"
#include "WHtreePartition.h"



class treeManager
{
public:
    explicit treeManager( WHtree* const tree, bool verbose = true) :
        m_tree( *tree ), m_logfile( 0 ), m_tractThreshold( 0.4 ), m_logFactor( 5.0 ), m_verbose( verbose )
    {
    }

    ~treeManager()
    {
    }

    bool setSingleTractFolder( std::string singleTractFolder )
    {
        m_singleTractFolder = singleTractFolder;
    }

    bool setMeanTractFolder( std::string meanTractFolder )
    {
        m_meanTractFolder = meanTractFolder;
    }

    bool setMaskFilename( std::string maskFilename )
    {
        m_maskFilename = maskFilename;
    }

    bool setFullTractFolder( std::string fullTractFolder )
    {
        m_fullTractFolder = fullTractFolder;
    }

    bool setDistMatrixFolder( std::string distMatrixFolder )
    {
        m_distMatrixFolder = distMatrixFolder;
    }

    bool setOutputFolder( std::string outputFolder )
    {
        m_outputFolder = outputFolder;
    }

    void log( std::ofstream* const logfile )
    {
        m_logfile = logfile;
    }

    void writeTree() const;
    void writeDebugTree() const;
    compactTract getMeanTract( const size_t inNode ) const;
    void writeMeanTracts( std::vector< size_t > inNodes ) const;
    void writeMeanTracts( std::vector< nodeID_t > inNodes ) const;
    void writeAllNodeTracts( const float memory ) const;

    void writeFullTract( const nodeID_t inNode, bool uFloat = true, bool doZip = false );
    void writeFullTract( std::vector< size_t > inNodes, bool uFloat = true, bool doZip = false );
    void writeFullTract( std::vector< nodeID_t > inNodes, bool uFloat = true, bool doZip = false );

    void writeFullBaseNodeTracts( bool uFloat = true, bool doZip = true, bool onlyMasks = false );

    float doCpcc();
    void flipX();
    void recaptureLeaves();

    WHtree &m_tree;

private:
    bool m_verbose;
    std::string m_treeFilename;
    std::string m_singleTractFolder;
    std::string m_meanTractFolder;
    std::string m_maskFilename;
    std::string m_fullTractFolder;
    std::string m_distMatrixFolder;
    std::string m_outputFolder;

    std::ofstream* m_logfile;

    float m_tractThreshold;
    float m_logFactor;

    void writeNodeTracts( std::vector< size_t >* const nodeVector, listedCache< compactTract >* const cache, size_t* const tractProg,
                    time_t* lastTime, const time_t &startTime ) const;

    void logWriteTract( const size_t &nodeID, compactTract &tract, const std::string &meanTractFolder,
                    volatile size_t* const threadcount ) const;

    void threadCout( const std::string &coutString, volatile size_t* const threadcounter );
};

#endif  // TREEMANAGER_H
