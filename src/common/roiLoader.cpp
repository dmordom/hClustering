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
// hClustering is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
//
// For more reference on the underlying algorithm and research they have been used for refer to:
// - Moreno-Dominguez, D., Anwander, A., & Knösche, T. R. (2014).
//   A hierarchical method for whole-brain connectivity-based parcellation.
//   Human Brain Mapping, 35(10), 5000-5025. doi: http://dx.doi.org/10.1002/hbm.22528
// - Moreno-Dominguez, D. (2014).
//   Whole-brain cortical parcellation: A hierarchical method based on dMRI tractography.
//   PhD Thesis, Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig.
//   ISBN 978-3-941504-45-5
//// (at your option) any later version.
// http://creativecommons.org/licenses/by-nc/3.0
//
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------
#include "roiLoader.h"


bool RoiLoader::readRoi(const std::string &roiFilename, HC_GRID *datasetGridPtr, WHcoord *datasetSizePtr, size_t* numStreamlinesPtr, std::vector< WHcoord > *coordinatesPtr, std::vector< size_t > *trackidsPtr )
{
    HC_GRID& datasetGrid( *datasetGridPtr );
    WHcoord& datasetSize( *datasetSizePtr );
    size_t& numStreamlines( *numStreamlinesPtr );
    std::vector< WHcoord >& coordinates( *coordinatesPtr);
    std::vector< size_t >& trackids( *trackidsPtr );

    coordinates.clear();
    trackids.clear();
    bool surf( false );

    WFileParser parser( roiFilename );
    if( !parser.readFile() )
    {
        std::cerr << "ERROR @ roiLoader::readRoi(): Parser error" << std::endl;
        return false;
    }
    std::vector< std::string > lines = parser.getRawLines();
    if( lines.size() == 0 )
    {
        std::cerr << "ERROR @ roiLoader::readRoi(): File is empty" << std::endl;
        return false;
    }

    {
        std::vector< std::vector< std::string > > datasetStrings = parser.getLinesForTagSeparated( "imagesize" );
        if( datasetStrings.size() == 0 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): Dataset size was not found in tree file" << std::endl;
            return false;
        }
        if( datasetStrings.size() > 1 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): Dataset attribute had multiple lines" << std::endl;
            return false;
        }
        WHcoord thisDatasetSize( string_utils::fromString< coord_t >( datasetStrings[0][0] ), string_utils::fromString< coord_t >(
                        datasetStrings[0][1] ), string_utils::fromString< coord_t >( datasetStrings[0][2] ) );

        std::string gridString( datasetStrings[0][3] );
        if( gridString == getGridString( HC_VISTA ) )
        {
            datasetGrid = HC_VISTA;
        }
        else if( gridString == getGridString( HC_NIFTI ) )
        {
            datasetGrid = HC_NIFTI;
        }
        else if( gridString == getGridString( HC_SURF ) )
        {
            datasetGrid = HC_SURF;
            surf = true;
        }
        else
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): Dataset grid type string \"" << gridString
                            << "\" could not be identified" << std::endl;
            return false;
        }
        datasetSize = thisDatasetSize;
    }

    {
        std::vector< std::vector< std::string > > streamNumberStrings = parser.getLinesForTagSeparated( "streams" );
        if( streamNumberStrings.size() == 0 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): tracking streams number was not found in roi file" << std::endl;
            return false;
        }
        if( streamNumberStrings.size() > 1 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): tracking streams number attribute has multiple lines" << std::endl;
            return false;
        }
        if( streamNumberStrings[0].size() > 1 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): tracking streams number attribute has multiple elements" << std::endl;
            return false;
        }
        numStreamlines = string_utils::fromString< size_t >(streamNumberStrings[0][0]);
    }


    {
        std::vector< std::vector< std::string > > coordStrings = parser.getLinesForTagSeparated( "roi" );
        if( coordStrings.size() == 0 )
        {
            std::cerr << "ERROR @ roiLoader::readRoi(): no roi coordinates in roi file (lacking #roi tag?)" << std::endl;
            return false;
        }
        coordinates.reserve( coordStrings.size() );
        for( size_t i = 0; i < coordStrings.size(); ++i )
        {
            WHcoord tempCoord( string_utils::fromString< coord_t >( coordStrings[i][0] ), string_utils::fromString< coord_t >(
                            coordStrings[i][1] ), string_utils::fromString< coord_t >( coordStrings[i][2] ) );
            coordinates.push_back( tempCoord );
        }
        if( m_autoFitGrid )
        {
            if( datasetGrid == HC_VISTA )
            {
                if( m_niftiMode )
                {
                    for( size_t i = 0; i<coordinates.size(); ++i )
                    {
                        coordinates[i] = coordinates[i].vista2nifti( datasetSize );
                    }
                    datasetGrid = HC_NIFTI;
                }
            }
            else if( datasetGrid == HC_NIFTI )
            {
                if( !m_niftiMode )
                {
                    for( size_t i = 0; i<coordinates.size(); ++i )
                    {
                        coordinates[i] = coordinates[i].nifti2vista( datasetSize );
                    }
                    datasetGrid = HC_VISTA;
                }
            }
            else
            {
                if( m_niftiMode )
                {
                    for( size_t i = 0; i<coordinates.size(); ++i )
                    {
                        coordinates[i] = coordinates[i].surf2nifti( datasetSize );
                    }
                    datasetGrid = HC_NIFTI;
                }
                else
                {
                    for( size_t i = 0; i<coordinates.size(); ++i )
                    {
                        coordinates[i] = coordinates[i].surf2vista( datasetSize );
                    }
                    datasetGrid = HC_VISTA;
                }
            }
        }
        else
        {
            if( datasetGrid == HC_VISTA )
            {
                if( m_niftiMode )
                {
                    std::cerr << "ERROR @ roiLoader::readRoi(): " << getGridString( datasetGrid ) << " format indicated by roi file does not coincide with active fileManager format (nifti) " << std::endl;
                    return false;
                }
            }
            else if( datasetGrid == HC_NIFTI )
            {
                if( !m_niftiMode )
                {
                    std::cerr << "ERROR @ roiLoader::readRoi(): " << getGridString( datasetGrid ) << " format indicated by roi file does not coincide with active fileManager format (vista) " << std::endl;
                    return false;
                }
            }
            else
            {
                std::cerr << "ERROR @ roiLoader::readRoi(): format indicated by roi is neither nifti nor vista and autofit was not set" << std::endl;
                return false;
            }
        }
    }

    {
        std::vector< std::vector< std::string > > indexStrings = parser.getLinesForTagSeparated( "trackindex" );
        if( indexStrings.empty() )
        {
            if( datasetGrid == HC_NIFTI && !surf )
            {
                std::cerr << "ERROR @ roiLoader::readRoi(): no tract ids in roi file, necessary to work on nifti mode" << std::endl;
                return false;
            }
            else
            {
                trackids.reserve( coordinates.size() );
                for( size_t i = 0; i < coordinates.size(); ++i )
                {
                    trackids.push_back( i );
                }
            }
        }
        else
        {
            trackids.reserve( indexStrings.size() );
            for( size_t i = 0; i < indexStrings.size(); ++i )
            {
                size_t tempIndex( string_utils::fromString< size_t >( indexStrings[i][0] ) );
                trackids.push_back( tempIndex );
            }
        }
    }

    if(  !trackids.empty() && ( coordinates.size() != trackids.size() ) )
    {
        std::cerr << "ERROR @ roiLoader::readRoi(): coordinate list and track id list have different lengths " << std::endl;
        return false;
    }
    return true;


}
