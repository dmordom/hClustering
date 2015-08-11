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

#ifndef SEEDMATCHER_H
#define SEEDMATCHER_H


// std library
#include <vector>
#include <list>
#include <string>
#include <map>
#include <fstream>
#include <sstream>

// hClustering
#include "WHcoord.h"
#include "WFileParser.h"
#include "compactTract.h"
#include "fileManagerFactory.h"
#include "roiLoader.h"



/**
 * this class handles interpolated projection of seed voxel tractograms onto the vertices of a surface
 */
class surfProjecter
{
public:
    /**
     * Constructor
     * \param roiFfilename the file with the coordinates and ids of the tractograms to be projected
     * \param surfFfilename the file with the coordinates of the surface vertices
     * \param verbose the file verbose output flag
     */
    surfProjecter( std::string roiFfilename, std::string surfFfilename, bool verbose );

    //! Destructor
    ~surfProjecter() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * returns the ready flag
     * \return true if both seed roi and surface vertices were successfully loaded
     */
    inline bool ready() const { return m_loaded; }

    /**
     * sets the averaging kernel interpolation mode
     * \param kernelRadius the radius around the surf coordinate to average from
     */
    inline void kernelMean( float kernelRadius ) { m_gaussKernel = false; m_kernelRadius = kernelRadius; }

    /**
     * sets the gaussian kernel interpolation mode
     * \param kernelRadius the radius around the surf coordinate to consider for interpolation
     * \param fwhm the full-width half-maximum of the desired gaussian curve
     */
    inline void kernelGauss( float kernelRadius, float fwhm ) { m_gaussKernel = true;
                                                                m_kernelRadius = kernelRadius;
                                                                m_sigma = fwhm / 2.35482; }

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Computes the projection with simple nearest neighbor interpolation
     */
    void matchCoordsNearestNb();

    /**
     * Computes the projection with the current active kernel mode
     */
    void matchCoordsKernel();

    /**
     * Writes to file the computed mean tracts corresponding to each surface voxel
     * \param leafTractFolder the folder with the leaf tractograms
     * \param outputTractFolder the folder where to write the output tracts
     * \param useFloat the float flag. If set tracts will be written in float32, otherwise in uint8
     * \param doZip the zip flag. If set tracts will be zipped after being written.
     */
    void writeMeanTracts( std::string leafTractFolder, std::string outputTractFolder, bool useFloat, bool doZip);



private:
    // === PRIVATE DATA MEMBERS ===

    bool m_loaded;                      //!< The ready to compute flag
    bool m_verbose;                     //!< The verbose output flag
    bool m_gaussKernel;                 //!< The gauss kernel flag. If set gaussian kernel will be used, otherwise averaging kernel will.
    bool            m_niftiMode;        //!< The Nifti flag. If true, files are coordinates are in nifti reference frame, if False, in vista reference frame.
    size_t          m_numStreamlines;   //!< Number of stramlines generated form each seed voxel, needed to properly do the transformation between natural units and logarithmic units in the tractograms.
    float           m_logFactor;        //!< The Logarithmic factor to be applied when converting the tractogram data ftom log to natural units. Calculated in class constructor.
    float m_kernelRadius;               //!< The radius of the interpolation kernel (in voxel distance units).
    float m_sigma;                      //!< The standard deviation of the gaussian function.
    HC_GRID         m_roiDatasetGrid;   //!< Contains the type of coordinate frame (vista, nifti) the input data was stored in.
    HC_GRID         m_surfDatasetGrid;  //!< Contains the type of coordinate frame (vista, nifti, RAS) the surface vertex data was stored in.
    WHcoord         m_roiDatasetSize;   //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file.
    WHcoord         m_surfDatasetSize;  //!< Contains the size in voxels of the dataset where the surface coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file.
    std::string m_leafTractFolder;      //!< The folder containing the leaf tractograms.
    std::vector<size_t> m_trackids;     //!< Stores the ids of the seed tracts correesponding to each leaf.
    std::vector< WHcoord > m_roiCoords; //!< A vector where the seed voxel coordinates are stored.
    std::vector< WHcoord > m_surfCoords;//!< A vector where the surface vertex coordinates are stored.
    std::vector< std::vector< WHcoord > > m_coordMatch; //!< The roi coordinates matched to every surface vertex.
    std::vector< std::vector< std::pair< float, size_t > > > m_matchDists;  //!<  The IDs matched to every surface vector and the physical distances to them (for kernel interpolation purposes).

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * Computes defined by the active kernel and the provided IDs and distances
     */
    compactTract getMeanTract( std::vector< std::pair< float, size_t > > inIDs ) const;

    /**
     * Matches a surface vertes with the neares roi voxel
     */
    std::pair< float, size_t > matchSurfNearest( WHcoord surfCoord );






};

#endif  // SEEDMATCHER_H
