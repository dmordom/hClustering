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
#include "fileManagerFactory.h"
#include "WHtreePartition.h"

/**
 * this class implements tree related operations dealing with tractograms and distance blocks, such as cpcc computations
 */
class treeManager
{
public:
    /**
     * Constructor
     * \param tree a pointer to the tree object to work with
     * \param verbose the verbose output flag
     */
    explicit treeManager( WHtree* const tree, bool verbose );

    //! Destructor
    ~treeManager() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    /**
     * sets the single voxel (leaf) tract folder for reading / writing
     * \param leavesTractFolder folder where to read or write compact leaf tracts
     */
    inline void setSingleTractFolder( std::string leavesTractFolder ) { m_leavesTractFolder = leavesTractFolder; }

    /**
     * sets the compact mean tract folder for reading / writing
     * \param meanTractFolder folder where to read or write compact node tracts
     */
    inline void setMeanTractFolder( std::string meanTractFolder ) { m_meanTractFolder = meanTractFolder; }

    /**
     * sets the filename path where to find the white matter tract used when tracking
     * \param maskFilename the white matter mask filepath
     */
    inline void setMaskFilename( std::string maskFilename ) { m_maskFilename = maskFilename; }

    /**
     * sets the full tract folder
     * \param fullTractFolder folder where to write the full 3D image tractograms
     */
    inline void setFullTractFolder( std::string fullTractFolder ) { m_fullTractFolder = fullTractFolder; }

    /**
     * sets the distance matrix folder
     * \param distMatrixFolder folder where to find the distance matrix blocks
     */
    inline void setDistMatrixFolder( std::string distMatrixFolder ) { m_distMatrixFolder = distMatrixFolder; }

    /**
     * sets the output folder
     * \param outputFolder folder where to write the created files
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
     * sets the the zip flag in order to zip the files upon writing
     */
    inline void storeZipped() { m_zipFlag = true; }

    /**
     * resets the the zip flag in order not to zip the files upon writing
     */
    inline void storeUnzipped() { m_zipFlag = false; }

    /**
     * sets the the float flag in order to write tracts/images with float data type
     */
    inline void writeInFloat() { m_floatFlag = true; }

    /**
     * resets the the float flag in order to write tracts/images with unsigned char data type
     */
    inline void writeInChar() { m_floatFlag = false; }


    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Writes the data files from the computed trees to the output folder
     */
    void writeTree() const;

    /**
     * Writes the debug data files from the computed trees to the output folder
     */
    void writeDebugTree() const;

    /**
     * computes an average tractogram corresponding to a node in the tree
     * \param inNode the ID of the node to compute the tractogram from
     * \return the mean tractogram object
     */
    compactTract getMeanTract( const size_t inNode ) const;

    /**
     * computes and writes to file average tractogram corresponding to a set of nodes in the tree
     * \param inNodes the vector of node IDs to compute the tractograms from
     */
    void writeMeanTracts( std::vector< size_t > inNodes ) const;
    /**
     * \overload
     */
    void writeMeanTracts( std::vector< nodeID_t > inNodes ) const;

    /**
     * computes and writes to file the average tractograms corresponding to all nodes in the tree
     * \param memory the amount of available RAM memory in GBs
     */
    void writeAllNodeTracts( const float memory ) const;

    /**
     * computes and writes to file the full 3D image tractograms corresponding to a node in the tree
     * \param inNode the ID of the node to compute the tractogram from
     */
    void writeFullTract( const nodeID_t inNode );
    /**
     * computes and writes to file the full 3D image tractograms corresponding to a set of nodes in the tree
     * \param inNodes the vector of node IDs to compute the tractograms from
     */
    void writeFullTract( std::vector< size_t > inNodes );
    /**
     * \overload
     */
    void writeFullTract( std::vector< nodeID_t > inNodes );

    /**
     * computes and writes to file 3D masks of a node cluster cointained seed voxels mask
     * \param inNode the ID of the node to compute the cluster mask from
     */
    void writeClusterMasks( const nodeID_t inNode );
    /**
     * computes and writes to file 3D masks of a set of nodes' cluster cointained seed voxels masks
     * \param inNodes the IDs of the nodes to compute the cluster masks from
     */
    void writeClusterMasks( std::vector< size_t > inNodes );
    /**
     * \overload
     */
    void writeClusterMasks( std::vector< nodeID_t > inNodes );

    /**
     * computes and writes to file the full 3D image tractograms and cluster seed voxel masks corresponding to the base nodes in the tree (meta-leaves)
     * \param onlyMasks if set only the cluster seed voxel masks, and not the tractograms, will be written
     */
    void writeFullBaseNodeTracts( bool onlyMasks = false );

    /**
     * computes the cophenetic correlation coefficient of the currently loaded tree and adds it to the corresponding tree data member
     * \return the cpcc value of the loaded tree
     */
    float doCpcc();

    /**
     * flips the seed voxel coordinates of the currently loaded tree in the X axis, used to compare  right hemisphere trees with left hemisphere trees
     */
    void flipX();

    // === PUBLIC DATA MEMBERS ===

    WHtree &m_tree;     //!< a reference to the tree object to perform the operations on

private:
    // === PRIVATE DATA MEMBERS ===

    std::string     m_leavesTractFolder; //!< The folder to read / write leaf (single) tracts
    std::string     m_meanTractFolder;   //!< The folder to read / write node (mean) tracts
    std::string     m_maskFilename;      //!< The filepath to the white matter mask used when performing tractography
    std::string     m_fullTractFolder;   //!< The folder to write full 3D image tracts
    std::string     m_distMatrixFolder;  //!< The folder where to find the distance matrix blocks
    std::string     m_outputFolder;      //!< The folder to write output files

    std::ofstream*  m_logfile;           //!< A pointer to the output log file stream
    bool            m_niftiMode;         //!< The Nifti flag. If true, files are coordinates are in vista reference frame, if False, in vista reference frame
    bool            m_debug;             //!< The debug output flag. If true, additional detailed outputs meant for debug will be written.
    bool            m_verbose;           //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.

    bool            m_zipFlag;           //!< flag indicating if stored tractograms will be zipped or not
    bool            m_floatFlag;         //!< flag indicating tractograms will be written in char or float format

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * computes and write to file the mean tracts corresponding to a set of nodes in the tree
     * \param nodeVector a pointer to the vector containing the node IDs to compute tracts from
     * \param cache  a pointer to the cache memory where preloaded tracts are stored
     * \param tractProg a pointer to a variable to keep tract of the number of threads currently running
     * \param lastTime a pointer to a time object indicating the last time the output info was updated
     * \param startTime a time object indicating the time where the routine started
     */
    void writeNodeTracts( std::vector< size_t >* const nodeVector, listedCache< compactTract >* const cache, size_t* const tractProg,
                    time_t* lastTime, const time_t& startTime ) const;

    /**
     * log-transforms and writes a node tract to file (executed on a separate thread)
     * \param nodeID the ID of the node corresponding to the tract to be written
     * \param tract the tract to be written
     * \param meanTractFolder the folder where to write the tract
     * \param threadcount a pointer to a variable keeping count of the number of threads already running
     */
    void logWriteTract( const size_t &nodeID, compactTract* tract, const std::string& meanTractFolder,
                    volatile size_t* const threadcount ) const;

};

#endif  // TREEMANAGER_H
