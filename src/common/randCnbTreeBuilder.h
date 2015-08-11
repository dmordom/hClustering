//---------------------------------------------------------------------------
//
// Project: hCLustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
//
// hClustering is free software: you can redistribute it and/or modify
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
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hCLustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------


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
#include "roiLoader.h"
#include "WHtree.h"
#include "protoNode.h"
#include "fileManager.h"
#include "cnbTreeBuilder.h"


#define DEBUG false


/**
 * This class simpler version of cbnTreeBuilder to build a centroid method hierarchical tree from small synthetic tractograms yoielding a random dissimilarity matrix
 * (no thresholding and no discarding outliers)
 * reads a seed voxel coordinates list, builds a centroid hierarchical tree from tractography data and writes output files
 */
class randCnbTreeBuilder
{
public:
    /**
     * Constructor
     * \param roiFilename file containing the list of seed voxels coordinates (for realistic neighborhood information)
     * \param verbose the verbose output flag
     */
    explicit randCnbTreeBuilder( std::string roiFilename, bool verbose = true );

    //! Destructor
    ~randCnbTreeBuilder() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets the output file stream for the program log file
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
    void buildRandCentroid( const unsigned int nbLevel, const float memory, const TC_GROWTYPE growType,
                        const size_t baseSize = 0, const bool keepDiscarded = true );


private:
    // === PRIVATE DATA MEMBERS ===

    dist_t          m_maxNbDist;         //!< The maximum distance(dissimilarity) value that a certain tractogram can have to its most similar neighbor in order not to be disregarded as an outlier. Taken from input parameter
    std::string     m_inputFolder;       //!< The folder path that contains the seed voxel tractograms
    std::string     m_outputFolder;      //!< The folder path where to write the output files
    std::ofstream*  m_logfile;           //!< A pointer to the output log file stream

    WHtree          m_tree;              //!< The class that will hold the built tree
    WHcoord         m_datasetSize;       //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file
    HC_GRID         m_datasetGrid;       //!< Containes the type of coordinate frame (vista, nifti) the input data was stored in
    bool            m_niftiMode;         //!< The Nifti flag. If true, files are coordinates are in vista reference frame, if False, in vista reference frame
    bool            m_roiLoaded;         //!< The ready flag. If true, the seed voxel list was successfully read and the class is ready to build the tree
    bool            m_treeReady;         //!< The tree building success flag. If true, the hierarchical centroid algorithm was successful and the class is ready to write the output data
    bool            m_debug;             //!< The debug output flag. If true, additional detailed outputs meant for debug will be written.
    bool            m_verbose;           //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    std::vector< WHcoord >  m_roi;       //!< A vector where the seed voxel coordinates are stored
    std::vector<size_t> m_trackids;      //!< Stores the ids of the seed tracts correesponding to each leaf

    size_t m_numComps;                   //!< A variable to store the total number of tractogram dissimilarity comparisons done while building the tree, for post-analysis and optimizing purposes
    std::vector< compactTract > m_leafTracts;  //!< A vector of all the tractogram objects with loaded data (only feasible if tracts are very small, this program is intended to be used with tractograms with less than 100 datapoints each)


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
     * Loads all tracts into a data member vector
     */
    void loadTracts();

    /**
     * Finds out the neighbroghood relationships between seed voxels and calculates the tractogram dissimilarity between all neighbors, data is saved into the protoLeaves vector
     * \param nbLevel the neighborhood level to be considered
     * \param protoLeavesPointer a pointer to the vector of proto-leaves where the neghborhood information and distance to neighbors will be stored
     * \return a list containing the coordinates of the voxels that were discarded during the initialization process (should be empty)
     */
    std::list< WHcoord > initialize( const unsigned int nbLevel, std::vector< protoNode >* protoLeavesPointer );

    /**
     * Calculates the distance values between a given seed voxel tract and its seed voxel neighbors
     * \param currentSeedID the ID of the leaf we want to compute the neighbor distances of
     * \param currentTract the tractogram information corresponding to that seed voxel
     * \param protoLeaves the vector of proto-leaves where the neghborhood information and distance to neighbors is stored
     * \param nbIDs a vector with the IDs of the neighbors to the current seed voxel
     * \param nbLeavesPointer a pointer to a map structure where the distance to the neghbors will be stored with the neighbor ID as key
     * \return a bit indicating if the seed voxel is to be discarded due to the high dissimilarity to its neighbors (true) or accepted as valid (false)
     */
    bool scanNbs( const size_t currentSeedID, const compactTract* const currentTractPointer, const std::vector< protoNode >& protoLeaves,
                  const std::vector< size_t >& nbIDs, std::map< size_t, dist_t >* nbLeavesPointer );

    /**
     * Writes a file with the base nodes (meta-leaves) obtained at the end of the homgeneous merging initial stage
     * \param baseNodes a vector with the base node IDs
     * \param filename the path to the file where the information will be written to
     */
    void writeBases( const std::vector< size_t >& baseNodes, const std::string& filename ) const;

    /**
     * Writes the data files from the computed trees to the output folder
     */
    void writeTree() const;

};

#endif  // RANDCNBTREEBUILDER_H
