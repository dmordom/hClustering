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
#include <errno.h>

// boost library
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

// hClustering
#include "compactTract.h"
#include "compactTractChar.h"
#include "WHcoord.h"
#include "roiLoader.h"
#include "WHtree.h"
#include "protoNode.h"
#include "listedCache.hpp"
#include "fileManagerFactory.h"

#define DEBUG false

/**
 * defines the type of homogeneous merging restriction applied to the initial phase of the algorithm
 */
typedef enum
{
    TC_GROWOFF,
    TC_GROWNUM,
    TC_GROWSIZE
} TC_GROWTYPE;


/**
 * this class implements the main functionalities required for building and saving a centroid-neighborhood hierarchcial tree from tractography data:
 * reads a seed voxel coordinates list, builds a centroid hierarchical tree from tractography data and writes output files
 */
class CnbTreeBuilder
{
public:
    /**
     * Constructor
     * \param roiFilename file containing the list of seed voxels coordinates (from which the tractograms were computed)
     * \param thresholdRatio the relative number of streamlines that have to pass through a certain voxel in a tractogram to be allowed to contribute towards tract similarity, will determine threshold value (threshold is applied to filter out tractography noise)
     * \param maxNbDist the maximum distance(dissimilarity) value that a certain tractogram can have to its most similar neighbor in order not to be disregarded as an outlier
     * \param verbose the verbose output flag
     * \param noLog if true no logarithmic normalization will be used
     */
    explicit CnbTreeBuilder( const std::string& roiFilename, const float  thresholdRatio, const float  maxNbDist, bool verbose, bool noLog );

    //! Destructor
    ~CnbTreeBuilder() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    /**
     * sets the input folder
     * \param inputFolder path to the folder where tractogram files are located
     */
    inline void setInputFolder( const std::string& inputFolder ) { m_inputFolder = inputFolder; }

    /**
     * sets the output folder
     * \param outputFolder folder where to write the created tree files
     */
    inline void setOutputFolder( const std::string& outputFolder ) { m_outputFolder = outputFolder; }

    /**
     * sets the temporal folder to store mean tracts during tree building
     * \param tempFolder folder where to temporally write down mean tracts during tree building
     */
    inline void setTempFolder( const std::string& tempFolder ) { m_tempFolder = tempFolder; }

    /**
     * sets (or resets) the debug output flag in order to write additional result files with detailed information meant for debug purposes
     * \param debug the true/false flag to set the m_debug member to
     */
    inline void setDebugOutput( bool debug = true ) { m_debug = debug; }

    /**
     * sets (or resets) the verbose output flag in order to wrire progress information on the standard output
     * \param verbose the true/false flag to set the m_verbose member to
     */
    inline void setVerbose( bool verbose = true ) { m_verbose = verbose; }

    /**
     * queries whether the roi file was loaded and therefore the class is ready for tree building
     * \return the ready flag, if true roi file has been successfully loaded
     */
    inline bool ready() const { return m_roiLoaded; }

    /**
     * queriesfor the size of the currently loaded roi (the number of seed voxels from which the tree will be built)
     * \return the size of the vector containing the seed voxel coordinates
     */
    inline size_t roiSize() const { return m_roi.size(); }


    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * core function of this class, includes all the necessary tree building steps
     * \param nbLevel the seed voxel neighborhood level restriction to implement, valid values are: 6, 18, 26 (recommended), 32, 92 and 124
     * \param memory the amount of RAM memory to be used for tractogram cache to speed up building the tree, in GBs (recommended value around 50% of total machine memory)
     * \param growType defines the type of limit for the initial homogenoeus merging stage. TC_GROWNUM = # of clusters, TC_GROWSIZE = size of clusters, TC_GROWOFF = skip this stage.
     * \param baseSize cluster size/number limit to stop the initial homogeneous merging stage. If set to 1 or 0 (default) this stage will be skipped
     * \param keepDiscarded flag to keep track of discarded voxels in a special field in the tree file, default value true
     */
    void buildCentroid( const unsigned int nbLevel, const float memory, const TC_GROWTYPE growType,
                        const size_t baseSize = 0, const bool keepDiscarded = true );


private:
    // === PRIVATE DATA MEMBERS ===

    float           m_tractThreshold;    //!< The threshold to be temporarily applied to the tracts before computing dissimilarity measures, in logarithmic space. Calculated in class constructor
    float           m_logFactor;         //!< The Logarithmic factor to be applied when converting the tractogram data ftom log to natural units. Calculated in class constructor
    dist_t          m_maxNbDist;         //!< The maximum distance(dissimilarity) value that a certain tractogram can have to its most similar neighbor in order not to be disregarded as an outlier. Taken from input parameter
    std::string     m_inputFolder;       //!< The folder path that contains the seed voxel tractograms
    std::string     m_outputFolder;      //!< The folder path where to write the output files
    std::string     m_tempFolder;        //!< The folder path where to temporarily store the mean tractograms during the tree building process
    std::ofstream*  m_logfile;           //!< A pointer to the output log file stream

    WHtree          m_tree;              //!< The class that will hold the built tree
    WHcoord         m_datasetSize;       //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file
    HC_GRID         m_datasetGrid;       //!< Containe the type of coordinate frame (vista, nifti) the input data was stored in
    size_t          m_numStreamlines;    //!< number of stramlines generated form each seed voxel, needed to properly do the transformation between natural units and logarithmic units in the tractograms

    bool            m_niftiMode;         //!< The Nifti flag. If true, files are coordinates are in nifti reference frame, if False, in vista reference frame
    bool            m_roiLoaded;         //!< The ready flag. If true, the seed voxel list was successfully read and the class is ready to build the tree
    bool            m_treeReady;         //!< The tree building success flag. If true, the hierarchical centroid algorithm was successful and the class is ready to write the output data
    bool            m_debug;             //!< The debug output flag. If true, additional detailed outputs meant for debug will be written.
    bool            m_verbose;           //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    std::vector< WHcoord >  m_roi;       //!< A vector where the seed voxel coordinates are stored
    std::vector<size_t> m_trackids;      //!< Stores the ids of the seed tracts correesponding to each leaf
    std::vector< double >   m_leafNorms; //!< A vector storing the comoputed norms for each seed voxel tractograms
    std::vector< double >   m_nodeNorms; //!< A vector storing the computed norms for the mean tractograms (corresponding to tree nodes)

             size_t m_numComps;          //!< A variable to store the total number of tractogram dissimilarity comparisons done while building the tree, for post-analysis and optimizing purposes
    volatile size_t m_ncHits;            //!< A variable to store the total number of successful node tractogram hits in the cache list, for post-analysis and optimizing purposes
    volatile size_t m_ncMiss;            //!< A variable to store the total number of misses when searching for node tractogram in the cache list, for post-analysis and optimizing purposes
    volatile size_t m_lcHits;            //!< A variable to store the total number of successful leaf tractogram hits in the cache list, for post-analysis and optimizing purposes
    volatile size_t m_lcMiss;            //!< A variable to store the total number of misses when searching for leaf tractogram in the cache list, for post-analysis and optimizing purposes


    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Fetches a proto-node or proto-leaf from the appropiate vector in pointer form provided its ID. This proto-node/-leaf can be modified.
     * A proto-node/-leaf contains the information of the dissimilarity of a node/leaf to its neighbors during tree construction, before it has been insterted definitively in the hierarchy
     * \param thisNode the full-ID of the desired proto-node/-leaf
     * \param protoLeavesPointer a pointer to the vector containing the proto-leaves
     * \param protoNodesPointer a pointer to the vector containing the proto-nodes
     * \return a pointer to the desired protonode within the corresponding vector
     */
    protoNode* fetchProtoNode( const nodeID_t& thisNode, std::vector< protoNode >* protoLeavesPointer, std::vector< protoNode >* protoNodesPointer ) const;

    /**
     * Fetches a node or leaf from the appropiate vector in pointer form provided its ID. This node/-leaf can be modified.
     * \param thisNode the full-ID of the desired node/-leaf
     * \param leavesPointer a pointer to the vector containing the leaves
     * \param nodesPointer a pointer to the the vector containing the nodes
     * \return a pointer to the desired node within the corresponding vector
     */
    WHnode* fetchNode( const nodeID_t& thisNode, std::vector< WHnode >* leavesPointer, std::vector< WHnode >* nodesPointer ) const;

    /**
     * Computes the norms of all the seed voxel tractograms and stores them in the m_leafNorms vector
     */
    void computeNorms();

    /**
     * Finds out the neighbroghood relationships between seed voxels and calculates the tractogram dissimilarity between all neighbors, data is saved into the protoLeaves vector
     * Seed voxel tracts with dissimilarity to its most similar neighbor greater than m_maxNbDist are discarded
     * \param nbLevel the neighborhood level to be considered
     * \param cacheSize the maximum number of leaf tracts that can be saved into RAM cache
     * \param protoLeavesPointer a pointer to the vector of proto-leaves where the neghborhood information and distance to neighbors will be stored
     * \return a list containing the coordinates of the voxels that were discarded during the initialization process
     */
    std::list< WHcoord > initialize( const unsigned int nbLevel, const size_t cacheSize, std::vector< protoNode >* protoLeavesPointer );

    /**
     * Calculates the distance values between a given seed voxel tract and its seed voxel neighbors
     * \param currentSeedID the ID of the leaf we want to compute the neighbor distances of
     * \param currentTract the tractogram information corresponding to that seed voxel
     * \param protoLeaves the vector of proto-leaves where the neghborhood information and distance to neighbors is stored
     * \param nbIDs a vector with the IDs of the neighbors to the current seed voxel
     * \param nbLeavesPointer a pointer to a map structure where the distance to the neghbors will be stored with the neighbor ID as key
     * \param cachePointer a pointer to the cache holding the data of previously loaded tractograms
     * \return a bit indicating if the seed voxel is to be discarded due to the high dissimilarity to its neighbors (true) or accepted as valid (false)
     */
    bool scanNbs( const size_t currentSeedID, const compactTractChar* const currentTractPointer, const std::vector< protoNode >& protoLeaves,
                  const std::vector< size_t >& nbIDs, std::map< size_t, dist_t >* nbLeavesPointer, listedCache< compactTractChar >* cachePointer );

    /**
     * Fetches a node tractogram from cache if present. Otherwise, loads the tractogram from file into a tractogram class,
     * applies log transform, thresholds its values, adds the pre-computed norm value to the class and stores it in cache
     * \param nodeID the ID of the the corresponding node
     * \param nodeMngrPointer a pointer to the file manager that handles reading node tracts from file
     * \param nodesCachePointer a pointer to the node tractogram cache
     * \return a pointer to the compactTract object stored in chache with the loaded tractogram data
     */
    compactTract* loadNodeTract( const size_t nodeID, const fileManager* const nodeMngrPointer,
                                 listedCache< compactTract >* nodesCachePointer );

    /**
     * Fetches a leaf tractogram from cache if present. Otherwise, loads the tractogram from file into a tractogram class,
     * thresholds its values, adds the pre-computed norm value to the class and stores it in cache
     * \param leafID the ID of the the corresponding leaf
     * \param leafMngrPointer a pointer to the file manager that handles reading leaf leaf tracts from file
     * \param leavesCachePointer a pointer to the leaf tractogram cache
     * \return a pointer to the compactTractChar object stored in chache with the loaded tractogram data
     */
    compactTractChar* loadLeafTract( const size_t leafID, const fileManager* const leafMngrPointer,
                                     listedCache< compactTractChar >* leavesCachePointer );

    /**
     * Fetches a leaf or node tractogram from cache or file calling to either loadNodeTract or loadLeafTract members
     * \param nodeID the full-ID of the the corresponding leaf or node
     * \param leafMngrPointer a pointer to the file manager that handles reading leaf leaf tracts from file
     * \param nodeMngrPointer a pointer to the file manager that handles reading node tracts from file
     * \param leavesCachePointer a pointer to the leaf tractogram cache
     * \param nodesCachePointer a pointer to the node tractogram cache
     * \return a pointer to void with the address of the compactTractChar or compactTractChar object stored in chache with the loaded tractogram data
     */
    void* loadTract( const nodeID_t nodeID, const fileManager* const leafMngrPointer, const fileManager* const nodeMngrPointer,
                     listedCache< compactTractChar >* leavesCachePointer, listedCache< compactTract >* nodesCachePointer );


    /**
     * Deletes a tract file that was temporarily stored in the hard disk
     * \param nodeID ID of the node corresponding to the tract that is to be deleted
     * \param fileMngrPointer a pointer to the output node tractogram file manager
     * \param threadcounterPointer a pointer to the volatile variable keeping track of the number of active boost threads writing to disk
     */
    void deleteTract( const size_t& nodeID, const fileManager* const fileMngrPointer, volatile size_t* threadcounterPointer ) const;


    /**
     * Writes the data files from the computed trees to the output folder
     */
    void writeTree() const;
};

#endif  // CNBTREEBUILDER_H
