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

#ifndef FILEMANAGER_H
#define FILEMANAGER_H

// std library
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <errno.h>

// boost library
#include <boost/format.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/detail/singleton.hpp>
#include <boost/filesystem.hpp>

// hClustering
#include "compactTract.h"
#include "WHcoord.h"

#define NODE_COMPACT_FNAME        "compact_%06d"        // %06d will be replaced with the seed ID
#define NODE_COMPACT_PREFIX_SIZE  8                     // string size of "compact_" is 8
#define NODE_FULL_FNAME           "fulltract_%06d"      // %06d will be replaced with the seed ID
#define NODE_FULL_PREFIX_SIZE     10                    // string size of "fulltract_" is 10
#define NODE_CLUSTER_FNAME        "cluster_%06d"        // %06d will be replaced with the seed ID
#define NODE_CLUSTER_PREFIX_SIZE  8                     // string size of "cluster_" is 8
#define DISTBLOCK_FNAME           "dist_block_%03d_%03d"

#define NIFTI_LEAF_COMPACT_FNAME  "probtract_%d"        // %d will be replaced with the seed ID
#define NIFTI_LEAF_FULL_FNAME     "probtract_full_%06d"   // %06d will be replaced with the seed ID

#define VISTA_LEAF_COMPACT_FNAME  "connect_%s"          // %s will be replaced with the seed coordinate string: XXX_YYY_ZZZ
#define VISTA_LEAF_FULL_FNAME     "connect_full_%s"     // %s will be replaced with the seed coordinate string: XXX_YYY_ZZZ

#define NIFTI_EXT ".nii"
#define VISTA_EXT ".v"
#define COMPACT_EXT ".cmpct"

#define COMPACT_FLOAT 32
#define COMPACT_UINT8 8


/**
 * defines a data type (unsigned char, float, etc.)
 */
typedef enum {
    VTBit,
    VTUINT8,
    VTFloat32,
    VTError
} ValueType;

typedef enum {
    ETFull,
    ETCompact
} TractExtType;

/**
 * this class handles reading and writing of images, matrices and tractograms
 */
class fileManager
{
public:
    /**
     * Constructor
     * \param ioFolderInit input-output folder for the reading and writign of files
     */
    explicit fileManager( std::string ioFolderInit = "" ): m_ioFolder( ioFolderInit ), m_zipFlag( false ), m_floatFlag( true ), m_logFlag( true ), m_thresFlag( false ) {}

    //! Destructor
    ~fileManager() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets the input/output folder
     * \param ioFolder input-output folder for the file manager
     */
    inline void setFolder( const std::string& ioFolder ) { m_ioFolder = ioFolder; }

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

    /**
     * sets the the log flag in order to consider read tracts to be in logarithmic units
     */
    inline void readAsLog() { m_logFlag = true; }

    /**
     * resets the the log flag in order to consider read tracts to be in natural units
     */
    inline void readAsNat() { m_logFlag = false; }

    /**
     * sets the the threshold flag in order to consider read tracts to be thresholded
     */
    inline void readAsThres() { m_thresFlag = true; }

    /**
     * resets the the threshold flag in order to consider read tracts to be non-thresholded
     */
    inline void readAsUnThres() { m_thresFlag = false; }

    /**
     * returns size of compact tract according to white matter mask (needs mask to be pre-loaded)
     * \return tractogram size (number of voxels in the white matter)
     */
    inline size_t getTractSize() const { return m_flipVector.size(); }

    /**
     * returns the pre-loaded white matter mask image)
     * \return white matter mask
     */
    inline std::vector< std::vector< std::vector< bool > > > getMaskMatrix() const { return m_maskMatrix; }


    // MEMBER FUNCTIONS


    /**
     * returns the filename of the compact node tractogram defined by the ID
     * \param tractNode the node ID
     * \return a string indicating the filename corresponding to that compact tractogram
     */
    std::string getNodeTractFilename( const size_t tractNode ) const;

    /**
     * returns the filename of the full node tractogram defined by the ID
     * \param tractNode the node ID
     * \return a string indicating the filename corresponding to that full tractogram
     */
    std::string getFullNodeTractFilename( const size_t tractNode ) const;

    /**
     * returns the filename of the distance matrix block defined by the block location IDs inside the matrix
     * \param blockID1 the block row ID
     * \param blockID2 the block column ID
     * \return a string indicating the filename corresponding to that block
     */
    std::string getBlockFilename( const unsigned int blockID1, const unsigned int blockID2 ) const;

    /**
     * returns the filename of the cluster mask (an image with the voxels of contained leafs set to 1) by the block node ID
     * \param node the cluster node ID
     * \return a string indicating the filename corresponding to that cluster mask
     */
    std::string getClusterMaskFilename( const size_t node ) const;

    /**
     * extracts the numeric ID of a node from a full node tractogram filename
     * \param filename the tractogram filename
     * \return a size_t value indicating the ID of the node corresponding to that tract filename
     */
    size_t readFullNodeTractNumber(const std::string filename) const;

    /**
     * extracts the numeric ID of a node from a cluster mask filename
     * \param filename the cluster mask filename
     * \return a size_t value indicating the ID of the node corresponding to that cluster mask filename
     */
    size_t readClusterMaskNumber(const std::string filename) const;

    /**
     * rearranges a compact tractogram to correspond with a 3D image flip in the x-axis
     * \param tractogramPointer a pointer to the tract to be flipped (data will be changed)
     */
    void flipXtract(compactTract* tractogramPointer) const;

    /**
     * calculates the mean coordinate (in euclidean space) of the pre-loaded image (must be a binary mask to work properly)
     * used to calculate a cluster euclidean center from a mask of its leaf voxels
     * \return WHcoord object with the coordinates of the calculated center
     */
    WHcoord meanCoordFromMask() const;

    /**
     * blows a compact tractogram vector tract into a full 3D image matrix
     * \param compact the compact tractogram vector to be blown
     * \param fullTract a pointer to the 3D image matrix to hold the blown up tractogram
     */
    void compact2full(const std::vector< float>& compact, std::vector< std::vector< std::vector <float> > >* fullTract ) const;

    /**
     * compacts a full 3D image tractogram into a compact tractogram vector
     * \param fullTract the full 3D image tractogram vector to be compacted
     * \param compact a pointer to the vector to hold the compacted tractogram
     */
    void full2compact(const std::vector< std::vector< std::vector <float> > >& fullTract, std::vector<float>* compact ) const;

    /**
     * extracts node compact tractogram data from file and returns it in compactTract object
     * \param tractNode the tractogram node ID
     * \param tractPointer a pointer to the tractogram object to hold the read tract
     */
    void readNodeTract (const size_t tractNode, compactTract* tractPointer) const;

    /**
     * \overload
     * \param tractLeaf the tractogram leaf ID
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates
     * \param tractPointer a pointer to the tractogram object to hold the read tract
     */
    void readLeafTract (const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, compactTract* tractPointer) const;
    /**
     * \overload
     */
    void readLeafTract (const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, compactTractChar* tractPointer) const;

    /**
     * reads a tractogram file and returns it in compactTract or vector object properly normalized according to datatype
     * \param tractFilename the filepath where the tractogram file is located
     * \param tractPointer a pointer to the tractogram object or vector to hold the read tract
     */
    void readTract(const std::string& tractFilename, compactTract* tractPointer) const;
    /**
     * \overload
     */
    void readTract(const std::string& tractFilename, compactTractChar* tractPointer) const;
    /**
     * \overload
     */
    void readTract(const std::string& tractFilename, std::vector<float>* tractPointer) const;
    /**
     * \overload
     */
    void readTract(const std::string& tractFilename, std::vector<unsigned char>* tractPointer) const;

    /**
     * reads a distance matrix block and returns it in  matrix (vector of vectors) form
     * \param blockID a pair with the row and column block IDs in the matrix
     * \param dBlockPointer a pointer to the matrix to hold the read block
     */
    void readDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> >* dBlockPointer ) const;
    /**
     * \overload
     */
    void readDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> >* dBlockPointer ) const;

    /**
     * reads a 3D image tractogram file and returns it in compactTract object
     * \param fullTractFilename the filepath where the tractogram file is located
     * \param tractPointer a pointer to the tractogram object to hold the read tract
     */
    void readFullTract( const std::string fullTractFilename, compactTract* tractPointer ) const;

    /**
     * loads a 3D image usually to be used as white matter mask for tractogram conversion
     * if workingin nifti mode it also keeps the image header to use when writing further files
     * \param maskFilename the filepath where the 3D file is located
     */
    void loadMaskImage( std::string maskFilename);

    /**
     * writes a node tractogram to file
     * \param tractNode ID of the node corresponding to the tract to be written
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeNodeTract (const size_t tractNode, const compactTract& tractogram) const;

    /**
     * writes a leaf tractogram to file
     * \param leaf seed voxel coordinate corresponding to the tract to be written
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeLeafTract (const WHcoord& leaf, const compactTract& tractogram) const;
    /**
     * \overload
     */
    void writeLeafTract (const WHcoord& leaf,const compactTractChar& tractogram) const;
    /**
     * \overload
     * \param leaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const compactTract& tractogram) const;
    /**
     * \overload
     */
    void writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const compactTractChar& tractogram) const;
    /**
     * \overload
     * \param leaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTract& tractogram) const;
    /**
     * \overload
     */
    void writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTractChar& tractogram) const;

    /**
     * writes a node tractogram to file
     * \param tractFilename output filepath where the tragtogram will be stored
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeTract (const std::string& tractFilename, const compactTract& tractogram) const;
    /**
     * \overload
     */
    void writeTract (const std::string& tractFilename, std::vector<float> tractogram) const;

    /**
     * writes a distance matrix block matrix to file
     * \param blockID a pair with the row and column block ID in the matrix
     * \param dBlock the matrix holding the block data to be written
     */
    void writeDistBlock (const std::pair<unsigned int,unsigned int> blockID, const std::vector<std::vector<float> >& dBlock) const;
    /**
     * \overload
     * \param blockID1 the row block ID in the matrix
     * \param blockID2 the column block ID in the matrix
     * \param dBlock the matrix holding the block data to be written
     */
    void writeDistBlock (const unsigned int blockID1, const unsigned int blockID2, const std::vector<std::vector<float> >& dBlock) const;

    /**
     * writes a a binary 3D image to file
     * \param maskFilename output filepath where the mask will be stored
     * \param maskMatrix 3D matrix holding the image data to be written
     */
    void writeMask( const std::string& maskFilename, const std::vector<std::vector<std::vector<float> > >& maskMatrix) const;
    /**
     * \overload
     */
    void writeMask( const std::string& maskFilename, const std::vector<std::vector<std::vector<bool> > >& maskMatrix) const;

    /**
     * convert node compact tract into full image tractogram and write to file
     * \param tractNode ID of the node corresponding to the tractogram to be written
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeFullNodeTract (const size_t tractNode, const compactTract& tractogram) const;

    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeFullLeafTract (const size_t& tractLeaf, const std::vector< size_t >& indexVector, const compactTract& tractogram) const;
    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeFullLeafTract (const size_t& tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTract& tractogram) const;

    /**
     * convert compact tract into full image tractogram and write to file
     * \param tractFilename file path where to write the output tractogram file
     * \param tractogram tractogram object holding the tractogram data to be written
     */
    void writeFullTract (const std::string& tractFilename, const compactTract& tractogram) const;

    /**
     * deletes the file corresponding to the indicated node tractogram
     * \param tractNode ID of the node corresponding to the tractogram to be deleted
     * \returns return integer indicating success or failure of the executed system command
     */
    int deleteTractFile (const size_t tractNode) const;


    // VIRTUAL MEMBER FUNCTIONS

    /**
     * returns a string with the filename extension to be used in reading/writing files
     * \return the filename extension string
     */
    virtual std::string getFileExtension( TractExtType extType ) const = 0;

    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates, when tractogram naming is defined by seed ID rather than coordinate the vector should be empty (default value)
     */
    virtual std::string getLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector = std::vector< WHcoord >() ) const = 0;

    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates, when tractogram naming is defined by seed ID rather than coordinate the vector should be empty (default value)
     */
    virtual std::string getFullLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector = std::vector< WHcoord >() ) const = 0;


    /**
     * reads a 1D vector from file
     * \param vectorFilename file path to the 1D vector image file
     * \param vectorPointer a pointer to the standard vector to hold the read data
     * \return the original data type of the read values
     */
    virtual ValueType readVector(const std::string& vectorFilename, std::vector<float>* vectorPointer) const = 0;

    /**
     * reads a 2D matrix from file
     * \param matrixFilename file path to the 2D matrix image file
     * \param matrixPointer a pointer to the matrix structure to hold the read data
     * \return the original data type of the read values
     */
    virtual ValueType readMatrix( const std::string& matrixFilename, std::vector<std::vector<float> >* matrixPointer) const = 0;

    /**
     * reads a 3D image from file
     * \param imageFilename file path to the 3D image file
     * \param imagePointer a pointer to the image structure to hold the read data
     * \return the original data type of the read values
     */
    virtual ValueType readImage( const std::string& imageFilename, std::vector<std::vector<std::vector<float> > >* imagePointer) const = 0;

    /**
     * loads an image header and displayes its data in the standard output if so specified (only implemented in nifti mode)
     * \param filename file path to the image file containing the header
     * \param display if set the header info will be displayed in standard output, default value false
     */
    virtual void loadHeader( const std::string& filename, bool display = false ) = 0;

    /**
     * writes a 1D vector to file
     * \param vectorFilename file path to the 1D vector image to store the data in
     * \param dataValueType the data type to write the data in
     * \param vector the 1D vector holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    virtual void writeVector( const std::string& vectorFilename, const ValueType dataValueType, const std::vector<float>& vector, bool doZip = false ) const = 0;

    /**
     * writes a 2D matrix to file
     * \param matrixFilename file path to the 2D matrix image to store the data in
     * \param dataValueType the data type to write the data in
     * \param matrix the 2D matrix holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    virtual void writeMatrix( const std::string& matrixFilename, const ValueType dataValueType, const std::vector<std::vector<float> >& matrix, bool doZip = false ) const = 0;

    /**
     * writes a 3D image to file
     * \param imageFilename file path to the 3D image to store the data in
     * \param dataValueType the data type to write the data in
     * \param image the 3D matrix holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    virtual void writeImage( const std::string& imageFilename, const ValueType dataValueType, const std::vector<std::vector<std::vector<float> > >& image, bool doZip = false ) const = 0;


protected:
    // === PROTECTED DATA MEMBERS ===

    std::string m_ioFolder;     //!< string with the tractogram folder
    bool        m_zipFlag;      //!< flag indicating if stored tractograms will be zipped or not
    bool        m_floatFlag;    //!< flag indicating tractograms will be written in char or float format
    bool        m_logFlag;      //!< flag indicating whether read tracts will be marked as being in log units or natural
    bool        m_thresFlag;    //!< flag indicating whether read tracts will be marked as being thresholded or not
    std::vector<std::vector<std::vector<bool> > > m_maskMatrix;    //!< image mask indicating if the voxels in the image are part of the compressed tractograms
    std::vector<size_t> m_flipVector;    //!< vector with the positions of x-flipped elements

private:
};

#endif // FILEMANAGER_H
