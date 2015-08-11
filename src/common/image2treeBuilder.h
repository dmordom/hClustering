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


#ifndef IMAGE2TREEBUILDER_H
#define IMAGE2TREEBUILDER_H

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


// hClustering
#include "WHcoord.h"
#include "WFileParser.h"
#include "WHtree.h"
#include "fileManagerFactory.h"

/**
 * This class implements methods in order to accept a 3D partition label image and use this information to assign each base-node (meta-leaf) of a hierarchical tree
 * to one of the partition clusters, then create a new tree with the same base nodes as the orginal but with only one partition in the hierarchical structure:
 * the most similar to the one defined in the 3D partition label image. This single-partition tree can then be used to perform tree-comparison statistics between the
 * original tree and the single-partition tree.
 */
class image2treeBuilder
{
public:
    /**
     * Constructor
     * \param imageFilename file containing the 3D image with the partition labels
     * \param treeFilename file containing the tree with the base-nodes that must be matched to the partition
     * \param verbose the verbose output flag
     * \param baseFilename the file containing the base node information, leave empty if bases should simply be extracted from tree structure
     */
    explicit image2treeBuilder( std::string imageFilename, std::string treeFilename, bool verbose, std::string baseFilename = "" );

    //! Destructor
    ~image2treeBuilder() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets output file stream for the program log file
     * \param logfile a pointer to the output log file stream
     */
    inline void log( std::ofstream* const logfile ) { m_logfile = logfile; }

    /**
     * checks if all necessary data has been loaded
     * \return input ready flag, if true all is ready to generate the partition tree
     */
    inline bool inReady() const { return ( m_imageLoaded && m_treeLoaded && m_basesLoaded ); }

    /**
     * returns true if the new tree has been sucessfully generated
     * \return output ready flag, if true all partition tree has been sucessfully generated
     */
    inline bool outReady() const { return m_partTreeReady; }

    /**
     * sets the output folder
     * \param outputFolder folder where to write the created tree files
     */
    inline void setOutputFolder( std::string outputFolder ) { m_outputFolder = outputFolder; }

    /**
     * sets (or resets) the verbose output flag in order to wrire progress information on the standard output
     * \param verbose the true/false flag to set the m_verbose member to
     */
    inline void setVerbose( bool verbose = true ) { m_verbose = verbose; }

    /**
     * queriesfor the size of the currently loaded roi (the number of seed voxels from which the tree will be built)
     * \return the size of the vector containing the seed voxel coordinates
     */
    inline size_t roiSize() const { return m_tree.getNumLeaves(); }

    /**
     * carries out the base-node matching and partition tree generation process
     */
    void importImagePart();


private:
    // === PRIVATE DATA MEMBERS ===

    std::string     m_outputFolder;      //!< The folder path where to write the output files
    std::ofstream*  m_logfile;           //!< A pointer to the output log file stream

    bool m_imageLoaded;     //!< image loaded flag, if true the 3D partition image was successfully loaded.
    bool m_treeLoaded;      //!< tree loaded flag, if true the template tree was successfully loaded.
    bool m_basesLoaded;     //!< bases loaded flag, if true the base-nodes file was successfully generated, or bases were computed from tree.
    bool m_partTreeReady;   //!< output ready flag, if true the partition treewas successfully generated.
    bool m_verbose;         //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    bool m_niftiMode;       //!< The Nifti flag. If true, files are coordinates are in nifti reference frame, if False, in vista reference frame

    WHtree m_tree;                  //!< The class that will hold firstly the loaded tree, then the newly generated tree
    std::vector< WHcoord >  m_roi;  //!< A vector where the seed voxel coordinates are stored
    WHcoord m_datasetSize;          //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file
    HC_GRID m_datasetGrid;          //!< Containe the type of coordinate frame (vista, nifti) the input data was stored in

    std::vector< size_t > m_baseVector; //!< A vector with the base-node identifiers
    std::vector<std::vector<std::vector<size_t> > > m_partImage;    //!< A 3D matrix with the 3D partition labels loaded from the image

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Loads the 3D partition label image from file and stores the values in the corresponding data member matrix
     * \param imageFilename the path to the image file
     * \return the success bit, if true image was successfully loaded
     */
    bool loadImage( const std::string imageFilename );

    /**
     * Loads the base-node information from file
     * \param baseFilename the path to the bases file
     * \return the success bit, if true bases were successfully loaded
     */
    bool loadBases( const std::string baseFilename );

    /**
     * Writes the data files from the generated tree to the output folder
     */
    void writeTree() const;

};

#endif  // IMAGE2TREEBUILDER_H
