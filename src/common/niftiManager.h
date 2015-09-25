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

#ifndef NIFTIMANAGER_H
#define NIFTIMANAGER_H

// std library
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <climits>
#include <fstream>
#include <stdint.h>

// boost library
#include <boost/format.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/detail/singleton.hpp>
#include <boost/filesystem.hpp>

// hClustering
#include "fileManager.h"

// libnifti
#include <nifti/nifti1_io.h>
#include <nifti/nifti1.h>


const size_t FLOAT_SIZE = sizeof(float);
const size_t UINT32_SIZE = sizeof(uint32_t);
const size_t UINT8_SIZE = sizeof(uint8_t);


/**
 * this class handles reading and writing of images, matrices and tractograms in nifti format
 */
class niftiManager
        : public fileManager
{
public:
    /**
     * Constructor
     * \param ioFolder input-output folder for the reading and writing of files
     */
    explicit niftiManager (std::string ioFolder="")
    {
        m_ioFolder = ioFolder;
        m_zipFlag = false;
        m_floatFlag = true;
        m_logFlag = true;
        m_thresFlag = false;
        m_header.sizeof_hdr = 0;
    }

    //! Destructor
    ~niftiManager () {}


    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * returns a string with the filename extension to be used in reading/writing files
     * \return the filename extension string
     */
    std::string getFileExtension( TractExtType extType ) const;

    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates, , when tractogram naming is defined by seed ID rather than coordinate the vector should be empty (default value)
     */
    std::string getLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector = std::vector< WHcoord >() ) const;

    /**
     * \overload
     * \param tractLeaf the ID of the leaf corresponding to the desired tract
     * \param indexVector vector holding the tract index corresponding to each leaf ID
     * \param coordVector vector holding the seed voxel coordinates, , when tractogram naming is defined by seed ID rather than coordinate the vector should be empty (default value)
     */
    std::string getFullLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector = std::vector< WHcoord >() ) const;


    /**
     * reads a 1D vector from file
     * \param vectorFilename file path to the 1D vector image file
     * \param vectorPointer a pointer to the standard vector to hold the read data
     * \return the original data type of the read values
     */
    ValueType readVector(const std::string& vectorFilename, std::vector<float>* vectorPointer ) const;

    /**
     * reads a 2D matrix from file
     * \param matrixFilename file path to the 2D matrix image file
     * \param matrixPointer a pointer to the matrix structure to hold the read data
     * \return the original data type of the read values
     */
    ValueType readMatrix( const std::string& matrixFilename, std::vector<std::vector<float> >* matrixPointer ) const;

    /**
     * reads a 3D image from file
     * \param imageFilename file path to the 3D image file
     * \param imagePointer a pointer to the image structure to hold the read data
     * \return the original data type of the read values
     */
    ValueType readImage( const std::string& imageFilename, std::vector<std::vector<std::vector<float> > >* imagePointer ) const;


    /**
     * loads an image header and displayes its data in the standard output if so specified
     * the header information is saved into a data member
     * \param filename file path to the image file containing the header
     * \param display if set the header info will be displayed in standard output, default value false
     */
    void loadHeader( const std::string& filename, bool display = false );


    /**
     * writes a 1D vector to file
     * \param vectorFilename file path to the 1D vector image to store the data in
     * \param dataValueType the data type to write the data in
     * \param vector the 1D vector holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    void writeVector( const std::string& vectorFilename, const ValueType dataValueType, const std::vector<float>& vector, bool doZip = false ) const;

    /**
     * writes a 2D matrix to file
     * \param matrixFilename file path to the 2D matrix image to store the data in
     * \param dataValueType the data type to write the data in
     * \param matrix the 2D matrix holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    void writeMatrix( const std::string& matrixFilename, const ValueType dataValueType, const std::vector<std::vector<float> >& matrix, bool doZip = false ) const;

    /**
     * writes a 3D image to file
     * \param imageFilename file path to the 3D image to store the data in
     * \param dataValueType the data type to write the data in
     * \param image the 3D matrix holding the data to be written
     * \param doZip if set written files will be gzipped, default value false
     */
    void writeImage( const std::string& imageFilename, const ValueType dataValueType, const std::vector<std::vector<std::vector<float> > >& image, bool doZip = false ) const;


protected:
private:

    // === PRIVATE DATA MEMBERS ===

    nifti_1_header m_header; //!< nifti header structure to copy when writing full images


    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * display pre-loaded header information into standard output
     */
    void displayHeader() const;

    /**
     * read nifti header information from file
     * \param imageFilenameRef file containing the header to be read
     * \return nifti_1_header structure containing the header information
     */
    nifti_1_header readHeader( const std::string& imageFilenameRef ) const;

    /**
     * generates from scratch a standard header to use when writing a nifti image
     * \param dimx x-axis dimension of the image to be written
     * \param dimy y-axis dimension of the image to be written
     * \param dimz z-axis dimension of the image to be written
     * \param dataValueType data type that the values to be written will have in file
     * \param imageHeader pointer to the header structure to hold the generated header information
     */
    void generateHeader( const size_t dimx, const size_t dimy, const size_t dimz, ValueType dataValueType, nifti_1_header* imageHeader ) const;
};

#endif // NIFTIMANAGER_H
