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

#ifndef ROILOADER_H
#define ROILOADER_H

// std library
#include <vector>
#include <string>

// hClustering
#include "WHcoord.h"
#include "WFileParser.h"
#include "WStringUtils.h"


/**
 * this class implements a roi loading routine (reading dataset size, track ids and seed coordinates)
 * that will be shared by many classes and programs
 */
class RoiLoader
{
public:
    /**
      * Constructor
      * \param niftiMode if true working grid is nifti, otherwise, vista
      * \param autoFitGrid if true read grid will be automatically reconverted to workign grid if necessary
      */
    explicit RoiLoader( bool niftiMode, bool autoFitGrid = true ): m_niftiMode( niftiMode ), m_autoFitGrid( autoFitGrid ) {}

    //! Destructor
    ~RoiLoader() {}

    /**
     * parses the seed voxel file and saves the coordinates in a vector
     * \param roiFilename the file path to the seed voxel coordinate file
     * \param datasetGridPtr a pointer to save the type of coordinate frame (vista, nifti) the input data was stored in
     * \param datasetSizePtr a pointer to save the size in voxels of the dataset where the seed voxel coordinates correspond.
     * \param numStreamlinesPtr a pointer to save the number of streamlines generated per seed voxel.
     * \param coordinatesPtr a pointer to save the vector where the seed voxel coordinates are stored
     * \param trackidsPtr a pointer to save the vector with the ids of the seed tracts correesponding to each leaf
     * \return the success flag, if true roi file has been successfully loaded
     */
    bool readRoi( const std::string& roiFilename, HC_GRID* datasetGridPtr, WHcoord* datasetSizePtr, size_t* numStreamlinesPtr, std::vector< WHcoord >*  coordinatesPtr, std::vector<size_t>* trackidsPtr );

private:
    bool m_niftiMode;   //!< The Nifti flag. If true, files are coordinates are in nifti reference frame, if False, in vista reference frame
    bool m_autoFitGrid; //!< The autofit flag. If true, grid will  be automatically converted to that defined by m_niftiMode if necessary
};

#endif // ROILOADER_H
