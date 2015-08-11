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
// hCLustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------



#ifndef GRAPHTREEBUILDER_H
#define GRAPHTREEBUILDER_H

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

// hClustering
#include "WHcoord.h"
#include "distBlock.h"
#include "roiLoader.h"
#include "WHtree.h"
#include "fileManagerFactory.h"

/**
 * defines the type of graph linkage algorithm that  will be applied
 */
typedef enum
{
    TG_SINGLE,
    TG_COMPLETE,
    TG_AVERAGE,
    TG_WEIGHTED,
    TG_WARD
} TG_GRAPHTYPE;

/**
 * this class implements the main functionalities required for building and saving a graph-method-based hierarchcial tree from a precomputed distance matrix (using distBlocks program):
 * reads a seed voxel coordinates list, builds a graph hierarchical tree from distance datablocks and writes output files
 */
class graphTreeBuilder
{
public:
    /**
     * Constructor
     * \param roiFilename file containing the list of seed voxels coordinates (from which the tractograms were computed)
     * \param verbose the verbose output flag
     */
    explicit graphTreeBuilder( std::string roiFilename, bool verbose = true );

    //! Destructor
    ~graphTreeBuilder() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    /**
     * sets the input folder
     * \param inputFolder path to the folder where distance block files are located
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
     * \param graphMethod the graph linkage method to be used for tree building, possible options are: single, complete, average, weighted, ward
     */
    void buildGraph( const TG_GRAPHTYPE graphMethod );

private:
    // === PRIVATE DATA MEMBERS ===

    std::string     m_inputFolder;       //!< The folder path that contains the seed voxel tractograms
    std::string     m_outputFolder;      //!< The folder path where to write the output files
    std::ofstream*  m_logfile;           //!< A pointer to the output log file stream

    WHtree          m_tree;              //!< The class that will hold the built tree
    WHcoord         m_datasetSize;       //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file
    HC_GRID         m_datasetGrid;       //!< Containe the type of coordinate frame (vista, nifti) the input data was stored in
    size_t          m_numStreamlines;    //!< number of stramlines generated form each seed voxel

    bool            m_niftiMode;         //!< The Nifti flag. If true, files are coordinates are in vista reference frame, if False, in vista reference frame
    bool            m_roiLoaded;         //!< The ready flag. If true, the seed voxel list was successfully read and the class is ready to build the tree
    bool            m_treeReady;         //!< The tree building success flag. If true, the hierarchical centroid algorithm was successful and the class is ready to write the output data
    bool            m_debug;             //!< The debug output flag. If true, additional detailed outputs meant for debug will be written.
    bool            m_verbose;           //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    std::vector< WHcoord >  m_roi;       //!< A vector where the seed voxel coordinates are stored
    std::vector<size_t> m_trackids;      //!< Stores the ids of the seed tracts correesponding to each leaf


    /**
     * Fetches a node or leaf from the appropiate vector in pointer form provided its ID. This node/-leaf can be modified.
     * \param thisNode the full-ID of the desired node/-leaf
     * \param leavesPointer a pointer to the vector containing the leaves
     * \param nodesPointer a pointer to the the vector containing the nodes
     * \return a pointer to the desired node within the corresponding vector
     */
    WHnode* fetchNode( const nodeID_t thisNode, std::vector< WHnode >* leavesPointer, std::vector< WHnode >* nodesPointer ) const;

    /**
     * reconstructs the whole dist matrix into RAM memory from all the available distance Block.
     * WARNING: can be extremely memory intensive if distance matrix is big.
     * \param distMatrix a pointer to the matrix that will be filled with the loaded distance data
     */
    void loadDistMatrix( std::vector< std::vector< float > >* const distMatrix ) const;

    /**
     * Calculates a the distance between a newly merged cluster to another cluster from the distances of the pre-merge clusters to that cluster
     * WARNING: can be extremely memory intensive if distance matrix is big.
     * \param distance1 distance value from the first pre-merge cluster to a third cluster
     * \param distance2 distance value from the second pre-merge cluster to a third cluster
     * \param size1 size of first pre-merge cluster
     * \param size2 size of second pre-merge cluster
     * \param graphMethod the linkage method algorithm to be used to calculate the distance of the merged cluster to the third cluster
     * \return the distance value between the newly merged cluster to the third cluster
     */
    dist_t newGraphDist( const dist_t distance1, const dist_t distance2, const size_t size1, const size_t size2, const TG_GRAPHTYPE graphMethod ) const;

    /**
     * Writes the data files from the computed trees to the output folder
     */
    void writeTree() const;
};

#endif // GRAPHTREEBUILDER_H
