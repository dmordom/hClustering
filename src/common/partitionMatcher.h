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

class partitionMatcher
{
public:
    partitionMatcher( WHtree* const tree1, WHtree* const tree2, std::string matchFilename, bool verbose = true);

    ~partitionMatcher();

     void log( std::ofstream* const logfile )
    {
        m_logfile = logfile;
    }


    std::string reportBaseNodes() const;

    void findMatchingPartitions( const float lambda, const size_t levleDepth = 0);


    // returns if colors of tree1 were changed
    bool matchColors();


private:
    bool m_verbose;
    WHtree &m_tree1;
    WHtree &m_tree2;
    std::ofstream *m_logfile;
    std::vector< size_t > m_baseNodes1;
    std::vector< size_t > m_baseNodes2;
    std::vector< size_t > m_matchedBases1;
    std::vector< size_t > m_matchedBases2;
    std::vector< std::vector< size_t > > m_matchedBases4Node1;
    std::vector< std::vector< size_t > > m_matchedBases4Node2;

    std::vector< size_t > m_fullCorrespondence;

    void testBaseNodes();

    void loadCorrespondence( const std::string &matchFilename );

    size_t findRelativeBasenodeID( size_t absoluteID, const std::vector< size_t > &baseNodes ) const;

    std::vector< std::vector< bool > > getSignatureMatrix( const std::vector< size_t> &partition, const bool forTree1 ) const;

    float evalPartitionMatch( const float lambda, size_t size1, std::vector< std::vector< bool > > &signature1, size_t size2,  std::vector< std::vector< bool > > &signature2 ) const;


    size_t getClusterMatch( const size_t cluster1, const size_t cluster2 ) const;
    std::pair< std::vector< std::vector< size_t > >, std::vector< std::vector< size_t > > > getClusterMatchMatrix( const std::vector< size_t > &partition1, const std::vector< size_t > &partition2 ) const;
    std::pair< std::vector< size_t >, std::vector< size_t > > getClusterMatchTable( const std::vector< std::vector< size_t > > &matchMatrix ) const;
    double doPartMatch( const std::vector< size_t > &partition1, const std::vector< size_t > &partition2,
                       std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint1,
                       std::pair< std::vector< size_t >, std::vector< size_t > > *tablePoint2) const;

    size_t assignDepth( const size_t partSize );




    /// NEEDS IMPROVEMENT
    WHcoord shiftColor( const WHcoord &color, const size_t shiftIndex ) const;

    /// DEPRECATED
    std::vector< size_t> getRelMatchedBases( const size_t cluster, const bool forTree1 ) const;


    void searchPartition( size_t partSize, const std::vector< size_t > &thisPart, std::vector< std::vector< size_t > > &partVector);



};

#endif  // TREECOMPARER_H
