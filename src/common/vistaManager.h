#ifndef VISTAMANAGER_H
#define VISTAMANAGER_H

// std library
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>

// boost library
#include <boost/format.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/detail/singleton.hpp>
#include <boost/filesystem.hpp>


// vista
#include <viaio/Vlib.h>
#include <viaio/VImage.h>

// hClustering
#include "compactTract.h"
#include "WHcoord.h"



class vistaManager
{
public:

    // constructor
     vistaManager (std::string tractFolderInit): m_tractFolder(tractFolderInit), m_zipFlag(false), m_floatFlag(true), m_logFlag(true), m_thresFlag(false) {}
    ~vistaManager () {}

    // set zip flag
    void storeZipped() {m_zipFlag=true;}

    // reset zip flag
    void storeUnzipped() {m_zipFlag=false;}

    // tracts will be written in float
    void writeInFloat() {m_floatFlag=true;}

    // tracts will be written in char
    void writeInChar() {m_floatFlag=false;}

    // tracts read will be marked as logarithmic units
    void readAsLog() {m_logFlag=true;}

    // tracts read will be marked as natural units
    void readAsNat() {m_logFlag=false;}

    // tracts read will be marked as thresholded
    void readAsThres() {m_thresFlag=true;}

    // tracts read will be marked as non-thresholded
    void readAsUnThres() {m_thresFlag=false;}

    // returns size of compact tract according to mask
    size_t getTractSize() const { return m_flipVector.size(); }

    // returns the mask matrix
    std::vector<std::vector<std::vector<bool> > > getMaskMatrix() const { return m_maskMatrix; }


    // "getTractFilename()": returns tract filename for compact tracts
    void getTractFilename(const size_t tractNode, std::string &filename) const;
    void getTractFilename(const WHcoord tractCoord, std::string &filename) const;

    // "getTractFilename()": returns tract filename    for full tracts
    void getFullTractFilename(const size_t tractNode, std::string &filename) const;
    void getFullTractFilename(const WHcoord tractCoord, std::string &filename) const;

    void getClusterMaskFilename(const size_t node, std::string* filename) const;
    void getClusterMaskFilenameZip(const size_t node, std::string* filename) const;

    size_t readClusterMaskNumber (const std::string filename) const;
    size_t readFullTractNumber (const std::string filename) const;





    // "readLeafTract()": extracts compact tractogram data from a vista file and returns it in compact Tract form
    void readLeafTract (const WHcoord tractCoord, compactTract &tractogram) const;
    void readLeafTract (const WHcoord tractCoord, compactTractChar &tractogram) const;


    // "readNodeTract()": extracts compact tractogram data from a vista file and returns it in compact Tract form
    void readNodeTract (const size_t tractNode, compactTract &tractogram) const;


    void readStringTract (const std::string &tractFilename, compactTract &tractogram) const;


    // "writeLeafTract()": write leaf tractogram data into vista format file
    void writeLeafTract (const WHcoord leaf, const compactTract &tractogram) const;
    void writeLeafTract (const WHcoord leaf, const compactTractChar &tractogram) const;


    // "writeNodeTract()": write node tractogram data into vista format file
    void writeNodeTract (const size_t tractNode, const compactTract &tractogram) const;

    void writeTract (const std::string tractFilename, const compactTract &tractogram) const;


    // "deleteTractFile()": delete tractogram file
    int deleteTractFile (const size_t tractNode) const;

    // "loadMask()": loads the (white matter) mask
    void loadMask( std::string maskFilename);

    // "writeMask()": writes a binary mask
    void writeMask( const std::string maskFilename, const std::vector<std::vector<std::vector<bool> > > &maskMatrix) const;

    // "readMatrix()": reads a 2D matrix
    void readMatrix( std::string matrixFilename, std::vector<std::vector<float> > *matrixPointer) const;

    // "writeMatrix()": writes a 2D matrix
    void writeMatrix( const std::string matrixFilename, const std::vector<std::vector<float> > &matrix ) const;

    // "readFullTract()": reads a full image tract and converts it to a compact tract
    void readFullTract( const std::string fullTractFilename, compactTract* tractPointer ) const;

    // "storeFulltract()": write full image tractogram into vista format file and compress
    void storeFullTract (const size_t tractNode, const compactTract &tractogram) const;
    void storeFullTract (const WHcoord tractCoord, const compactTract &tractogram) const;
    void storeFullTract (const std::string &tractFilename, const compactTract &tractogram) const;


    // "getBlockFilename()": returns tract filename
    void getBlockFilename(const unsigned int blockID1, const unsigned int blockID2, std::string &filename) const;

    // "readDistBlock()": extracts a distance block from a vista file and returns it in vector of vectors form
    void readDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> > &dBlock) const;
    void readDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> > &dBlock) const;

    // "writeDistBlock()": writes a distance block matrix to a vista file
    void writeDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> > &dBlock);
    void writeDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> > &dBlock);

    // "flipXtract()": flips the tractogram in the X axis
    void flipXtract(compactTract &tractogram);

    // "meanCoordFromMask()": claculatzes the mean coordinate of the loaded mask
    WHcoord meanCoordFromMask() const;







protected:
private:

    // string with the tractogram folder
    std::string m_tractFolder;

    // flag indicating if stored tractograms will be zipped or not
    bool m_zipFlag;

    // flag indicating tractograms will be written in char or float format
    bool m_floatFlag;

    // flag indicating whether read tracts will be marked as being in log units or natural
    bool m_logFlag;

    // flag indicating whether read tracts will be marked as being thresholded or not
    bool m_thresFlag;

    // image mask indicating if the voxels in the image are part of the compressed tractograms
    std::vector<std::vector<std::vector<bool> > > m_maskMatrix;

    // vector with the positions of x-flipped elements
    std::vector<size_t> m_flipVector;

    // "get tract information and put it in a vector
    void getTract(const std::string &tractFilename, std::vector<float> &tractogram) const;
    void getTract(const std::string &tractFilename, std::vector<unsigned char> &tractogram) const;



    // "WriteVImage()": write Vista image
    int WriteVImage( const VString Name, VImage &Image ) const;

    // "ReadImage()": read Vista image file
    int ReadImage( const VString Name, VImage &Image) const;


};

#endif // VISTAMANAGER_H
