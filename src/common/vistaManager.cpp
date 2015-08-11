//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
// This file is also part of OpenWalnut ( http://www.openwalnut.org ).
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

#include "vistaManager.h"
#include <errno.h>

// PUBLIC MEMBERS

/*
// vistaManager::getTractFilename(): returns tract filename
std::string vistaManager::getLeafTractFilename (const WHcoord& tractCoord ) const
{
    std::string filename = str(boost::format( VISTA_LEAF_COMPACT_FNAME ) % tractCoord.getNameString() );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension();
    return filepath;
}// end vistaManager::getTractFilename() -----------------------------------------------------------------
*/
std::string vistaManager::getLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector ) const
{
    if( coordVector.empty() )
    {
        VError( "getLeafTractFilename(): coordVector is empty" );
    }
    if( tractLeaf >= coordVector.size() )
    {
        VError( "getLeafTractFilename(): leaf ID provided is higher than coordinate vector length" );
    }
    WHcoord tractCoord( coordVector[tractLeaf] );   
    std::string filename = str(boost::format( VISTA_LEAF_COMPACT_FNAME ) % tractCoord.getNameString() );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETCompact );
    return filepath;
}// end vistaManager::getTractFilename() -----------------------------------------------------------------

/*
// vistaManager::getFullTractFilename(): returns tract filename
std::string vistaManager::getFullLeafTractFilename(const WHcoord& tractCoord) const
{
    std::string filename = str(boost::format( VISTA_LEAF_FULL_FNAME ) % tractCoord.getNameString() );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension();
    return filepath;
}// end vistaManager::getFullTractFilename() -----------------------------------------------------------------
*/
std::string vistaManager::getFullLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector ) const
{
    if( coordVector.empty() )
    {
        VError( "getLeafTractFilename(): coordVector is empty" );
    }
    if( tractLeaf >= coordVector.size() )
    {
        VError( "getLeafTractFilename(): leaf ID provided is higher than coordinate vector length" );
    }
    WHcoord tractCoord( coordVector[tractLeaf] );
    std::string filename = str(boost::format( VISTA_LEAF_FULL_FNAME ) % tractCoord.getNameString() );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETFull );
    return filepath;
}// end vistaManager::getFullTractFilename() -----------------------------------------------------------------


// "readVector()": reads a 1D vector
ValueType vistaManager::readVector(const std::string& vectorFilenameRef, std::vector<float>* vectorPointer) const {


    std::vector<float>& vector = *vectorPointer;
    std::string vectorFilename( vectorFilenameRef );
    std::string vectorFilenameZipped( vectorFilenameRef );

    std::string extension( boost::filesystem::path(vectorFilename).extension().string() );


    if ( extension == ".gz" )
    {
        vectorFilename.resize(vectorFilename.size()-3);
        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + vectorFilenameZipped + " > " + vectorFilename );
        int success(system(sysCommand.c_str()));
    }
    else if ( extension == getFileExtension( ETCompact ) )
    {
    }
    else
    {
        std::cerr << "File \"" << vectorFilenameZipped << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return VTError;
    }

    // read vector image
    VImage vectorImage;
    if ( ! readVista( (char *)vectorFilename.c_str(), &vectorImage ) )
        VError( "readVector(): Failed to read vector image" );

    VRepnKind vectorRepnType(VPixelRepn(vectorImage));
    ValueType vectorValueType;
    if(vectorRepnType == VUByteRepn )
    {
        vectorValueType = VTUINT8;
    }
    else if(vectorRepnType == VFloatRepn )
    {
        vectorValueType = VTFloat32;
    }
    else if(vectorRepnType == VBitRepn )
    {
        vectorValueType = VTBit;
    }
    else
    {
        VError( "readVector(): vector representation type not recognized (neither VFloat nor VUByte)" );
    }

    const size_t Bands   = VImageNBands( vectorImage );
    const size_t Rows    = VImageNRows( vectorImage );

    if (Bands!=1 || Rows!= 1)
        VError("ERROR @ vistaManager::readVector(): vector image must have 1 row and 1 band only");

    const size_t nElements = VImageNColumns(vectorImage);
    vector.reserve(nElements); // Reserve memory in vector to increase speed

    for (size_t i=0 ; i<nElements ; ++i)
    {
        if (vectorRepnType==VFloatRepn)
        {
            vector.push_back( VPixel(vectorImage,0,0,i,VFloat) );
        }
        else if (vectorRepnType==VUByteRepn)
        {
            vector.push_back( (float)VPixel(vectorImage,0,0,i,VUByte) );
        }
        else if (vectorRepnType==VBitRepn)
        {
            vector.push_back( VPixel(vectorImage,0,0,i,VBit) );
        }
        else
            VError( "readVector(): vector representation type not recognized (neither VFloat nor VUByte nor VBit)" );
        }

    // clean up
    VDestroyImage( vectorImage );

    if (vector.empty())
        std::cerr<<" WARNING: vector is empty"<<std::endl;

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + vectorFilename );
        int success(system(sysCommand.c_str()));
    }
    return vectorValueType;
}// end vistaManager::readVector() -----------------------------------------------------------------

// "readMatrix()": reads a 2D matrix
ValueType vistaManager::readMatrix(const  std::string& matrixFilenameRef, std::vector<std::vector<float> >* matrixPointer ) const
{
    std::vector<std::vector<float> >& matrix = *matrixPointer;
    std::string matrixFilename( matrixFilenameRef );
    std::string matrixFilenameZipped( matrixFilenameRef );

    std::string extension( boost::filesystem::path(matrixFilename).extension().string() );


    if ( extension == ".gz" )
    {
        matrixFilename.resize(matrixFilename.size()-3);
        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + matrixFilenameZipped + " > " + matrixFilename );
        int success(system(sysCommand.c_str()));
    }
    else if ( extension == getFileExtension( ETFull ) )
    {
    }
    else
    {
        std::cerr << "File \"" << matrixFilenameZipped << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return VTError;
    }

    // read matrix image
    VImage matrixImage;
    if ( ! readVista( (char *)matrixFilename.c_str(), &matrixImage ) )
        VError( "loadMatrix(): Failed to read matrix image" );

    VRepnKind maskRepnType(VPixelRepn(matrixImage));
    ValueType maskValueType;
    if(maskRepnType == VBitRepn )
    {
        maskValueType = VTBit;
    }
    else if(maskRepnType == VUByteRepn )
    {
        maskValueType = VTUINT8;
    }
    else if(maskRepnType == VFloatRepn )
    {
        maskValueType = VTFloat32;
    }
    else
    {
        VError( "loadMatrix(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }

    const size_t Bands   = VImageNBands( matrixImage );
    if (Bands != 1)
        VError( "loadMatrix(): matrix image has multiple bands" );

    const size_t dimy = VImageNRows( matrixImage );
    const size_t dimx = VImageNColumns( matrixImage );
    matrix.clear();
    {
        std::vector<float> yvect(dimy,false);
        matrix.resize(dimx,yvect);
    }

    for (size_t i=0 ; i<dimy ; ++i)
    {
        for (size_t j=0 ; j<dimx ; ++j)
        {
            if (maskRepnType==VBitRepn)
                matrix[j][i] = VPixel(matrixImage,0,i,j,VBit);
            else if (maskRepnType==VFloatRepn)
                matrix[j][i] = VPixel(matrixImage,0,i,j,VFloat);
            else if (maskRepnType==VUByteRepn)
                matrix[j][i] = (float)VPixel(matrixImage,0,i,j,VUByte);
            else
                VError( "loadMatrix(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );
        }
    }

    // clean up
    VDestroyImage( matrixImage );

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + matrixFilename );
        int success(system(sysCommand.c_str()));
    }
    return maskValueType;
}// end vistaManager::readMatrix() -----------------------------------------------------------------

// "readImage()": reads a 3D image
ValueType vistaManager::readImage( const std::string& imageFilenameRef, std::vector<std::vector<std::vector<float> > >* imagePointer ) const
{
    std::vector<std::vector<std::vector<float> > >& imageRef( *imagePointer );
    std::string imageFilename( imageFilenameRef );
    std::string imageFilenameZipped( imageFilenameRef );

    std::string extension( boost::filesystem::path(imageFilename).extension().string() );


    if ( extension == ".gz" )
    {
        imageFilename.resize(imageFilename.size()-3);

        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + imageFilenameZipped + " > " + imageFilename );
        int success(system(sysCommand.c_str()));
    }
    else if ( extension == getFileExtension( ETFull ) )
    {
    }
    else
    {
        std::cerr << "File \"" << imageFilenameZipped << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return VTError;
    }

    // read mask image
    VImage thisImage;
    if ( ! readVista( (char *)imageFilename.c_str(), &thisImage ) )
        VError( "readImage(): Failed to read image" );

    VRepnKind ImageRepnType(VPixelRepn(thisImage));
    ValueType imageValueType;
    if(ImageRepnType == VBitRepn )
    {
        imageValueType = VTBit;
    }
    else if(ImageRepnType == VUByteRepn )
    {
        imageValueType = VTUINT8;
    }
    else if(ImageRepnType == VFloatRepn )
    {
        imageValueType = VTFloat32;
    }
    else
    {
        VError( "readImage(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }

    const size_t dimz = VImageNBands( thisImage );  // bands are the z coordinate
    const size_t dimy = VImageNRows( thisImage );  // rows are the y coordinate
    const size_t dimx = VImageNColumns( thisImage ); // columns are the x coordinate
    imageRef.clear();
    {
        std::vector<float> zvect(dimz,0);
        std::vector<std::vector<float> >yzmatrix(dimy,zvect);
        imageRef.resize(dimx,yzmatrix);
    }
    for (size_t i=0 ; i<dimz ; ++i)
    {
        for (size_t j=0 ; j<dimy ; ++j)
        {
            for (size_t k=0 ; k<dimx ; ++k)
            {

                if (ImageRepnType==VBitRepn)
                    imageRef[k][j][i] = VPixel(thisImage,i,j,k,VBit);
                else if (ImageRepnType==VFloatRepn)
                    imageRef[k][j][i] = VPixel(thisImage,i,j,k,VFloat);
                else
                    imageRef[k][j][i] = (float)VPixel(thisImage,i,j,k,VUByte);
            }
        }
    }

    // clean up
    VDestroyImage( thisImage );

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + imageFilename );
        int success(system(sysCommand.c_str()));
    }
    return imageValueType;
}// end vistaManager::readImage() -----------------------------------------------------------------


// "writeVector()": writes a 1D vector
void vistaManager::writeVector( const std::string& vectorFilename, const ValueType dataValueType, const std::vector<float>& vector, bool doZip ) const
{
    if ( vector.empty() )
    {
        std::cerr<< "ERROR @ vistaManager::writeVector(): vector is empty, it has not been stored" <<std::endl;
        return;
    }

    const size_t Columns = vector.size();

    VRepnKind vectorRepnType;
    if( dataValueType == VTBit )
    {
        vectorRepnType = VBitRepn;
    }
    else if(dataValueType == VTUINT8 )
    {
        vectorRepnType = VUByteRepn ;
    }
    else if(dataValueType == VTFloat32 )
    {
        vectorRepnType = VFloatRepn ;
    }
    else
    {
        VError( "writeVector(): vector representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }

    VImage vistaImage;
    vistaImage = VCreateImage(1,1,Columns,vectorRepnType);


    for( size_t i=0 ; i<Columns ; ++i )
    {
        if( dataValueType == VTBit )
        {
            VPixel(vistaImage,0,0,i,VBit) = ( vector[i] != 0 );
        }
        else if(dataValueType == VTUINT8 )
        {
            VPixel(vistaImage,0,0,i,VUByte) = (unsigned char)vector[i];
        }
        else if(dataValueType == VTFloat32 )
        {
            VPixel(vistaImage,0,0,i,VFloat) = vector[i];
        }
        else
        {
            VError( "writeVector(): vector representation type not recognized (neither VFloat nor VUByte nor VBit)" );
        }
    }

    // write vector to file
    if ( ! writeVista( (char *)vectorFilename.c_str() , vistaImage ) )
        VError( "writeVector(): Failed to open output vector file " );

    // clean up
    VDestroyImage( vistaImage );

    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + vectorFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end vistaManager::writeVector() -----------------------------------------------------------------


// "writeMatrix()": writes a 2D matrix
void vistaManager::writeMatrix( const std::string& matrixFilename, const ValueType dataValueType, const std::vector<std::vector<float> >& matrix, bool doZip ) const
{
    if (matrix.empty())
    {
        std::cerr<< "ERROR @ vistaManager::writeMatrix(): matrix is empty, it has not been stored" <<std::endl;
        return;
    }

    const size_t dimx = matrix.size();
    const size_t dimy = matrix[0].size();

    VRepnKind matrixRepnType;
    if( dataValueType == VTBit )
    {
        matrixRepnType = VBitRepn;
    }
    else if(dataValueType == VTUINT8 )
    {
        matrixRepnType = VUByteRepn ;
    }
    else if(dataValueType == VTFloat32 )
    {
        matrixRepnType = VFloatRepn ;
    }
    else
    {
        VError( "writeMatrix(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }

    VImage vistaImage;
    vistaImage = VCreateImage(1,dimy,dimx,matrixRepnType);

    for( size_t i=0 ; i<dimy ; ++i )
    {
        for( size_t j=0 ; j<dimx ; ++j )
        {
            if( dataValueType == VTBit )
            {
                VPixel(vistaImage,0,i,j,VBit) = ( matrix[j][i] != 0 );
            }
            else if(dataValueType == VTUINT8 )
            {
                VPixel(vistaImage,0,i,j,VUByte) = (unsigned char)matrix[j][i];
            }
            else if(dataValueType == VTFloat32 )
            {
                VPixel(vistaImage,0,i,j,VFloat) = matrix[j][i];
            }
            else
            {
                VError( "writeMatrix(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );
            }
        }
    }

    // write matrix to file
    if ( ! writeVista( (char *)matrixFilename.c_str() , vistaImage ) )
        VError( "writeMatrix(): Failed to open output matrix file " );

    // clean up
    VDestroyImage( vistaImage );

    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + matrixFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end vistaManager::writeMatrix() -----------------------------------------------------------------


// "writeImage()": writes a 3D image
void vistaManager::writeImage( const std::string& imageFilename, const ValueType dataValueType, const std::vector<std::vector<std::vector<float> > >& image, bool doZip ) const
{
    if (image.empty())
    {
        std::cerr<< "ERROR @ vistaManager::writeImage(): image matrix is empty, image has not been written" <<std::endl;
        return;
    }

    const size_t dimx = image.size();
    const size_t dimy = image[0].size();
    const size_t dimz = image[0][0].size();

    VRepnKind ImageRepnType;
    if( dataValueType == VTBit )
    {
        ImageRepnType = VBitRepn;
    }
    else if(dataValueType == VTUINT8 )
    {
        ImageRepnType = VUByteRepn ;
    }
    else if(dataValueType == VTFloat32 )
    {
        ImageRepnType = VFloatRepn ;
    }
    else
    {
        VError( "writeImage(): image representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }

    VImage vistaImage;
    vistaImage = VCreateImage(dimz,dimy,dimx,ImageRepnType);
    VFillImage(vistaImage,VAllBands,0);

    for( size_t i=0 ; i<dimz ; ++i )
    {
        for( size_t j=0 ; j<dimy ; ++j )
        {
            for( size_t k=0 ; k<dimx ; ++k )
            {
                float datapoint = image[k][j][i];
                if (datapoint)
                {
                    if( dataValueType == VTBit )
                    {
                        VPixel(vistaImage,i,j,k,VBit) = 1;
                    }
                    else if(dataValueType == VTUINT8 )
                    {
                        VPixel(vistaImage,i,j,k,VUByte) = (unsigned char)datapoint;
                    }
                    else if(dataValueType == VTFloat32 )
                    {
                        VPixel(vistaImage,i,j,k,VFloat) = datapoint;
                    }
                    else
                    {
                        VError( "writeImage(): image representation type not recognized (neither VFloat nor VUByte nor VBit)" );
                    }
                }
            }
        }
    }

    // write image to file
    if ( ! writeVista( (char *)imageFilename.c_str() , vistaImage ) )
        VError( "writeImage(): Failed to open output image file " );

    // clean up
    VDestroyImage( vistaImage );

    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + imageFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end vistaManager::writeImage() -----------------------------------------------------------------



// PROTECTED MEMBERS



// "readVista()": read Vista image file
int vistaManager::readVista( const VString& Name, VImage* ImagePointer ) const
{

    VImage& Image = *ImagePointer;

   VAttrList     list;   // attribute list
   VAttrListPosn pos;    // position in list


   // initialize results
   Image = NULL;


   boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
   ioMutex.lock();
   {
       // read file
       FILE* file;   // input file
       file = VOpenInputFile (Name, TRUE);
       list = VReadFile (file, NULL);
       if (!file) {
           fclose (file);
           VError ("ReadImage(): Failed to open input file '%s'", Name);
       }
       fclose (file);
   }
   ioMutex.unlock();

   if (!list)
   {
      VError ("readVista(): Failed to read input file '%s'", Name);
   }

   // extract image
   for (VFirstAttr (list, &pos); VAttrExists (&pos); VNextAttr (&pos))
   {
      if (VGetAttrRepn (&pos) != VImageRepn) continue;
      if (Image)
      {
         VDestroyImage (Image);
         Image = NULL;
         break;
      }
      VGetAttrValue (&pos, NULL, VImageRepn, &Image);
      VSetAttrValue (&pos, NULL, VImageRepn, NULL);
   }
   if (! (Image) &&  VAttrExists (&pos))
   {
       VError ("readVista(): Input file '%s' contains multiple images",  Name);
   }
   if (! (Image) && !VAttrExists (&pos))
   {
       VError ("readVista(): Input file '%s' does not contain an image", Name);
   }

   // clean-up
   VDestroyAttrList (list);

    return 1;

} // end "readVista()" -----------------------------------------------------------------



// "writeVista()": write Vista image
int vistaManager::writeVista( const VString& Name, const VImage& Image )  const {

    VAttrListPosn pos;       /* position in list */
    VBoolean      success;   /* success flag     */
    VAttrList list;

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;

    ioMutex.lock();
    {
        FILE*  file;      /* output file      */
        list = VCreateAttrList();
        VAppendAttr( list, "image", NULL, VImageRepn, Image );
        /* open file */
        file = fopen (Name, "w");
        if (!file) {
            std::string errorstring( "writeVista(): Failed to open output vista file. Code: "+ boost::lexical_cast<std::string>(errno) + " " + strerror(errno) );
            VError ( errorstring.c_str() );
        }
        /* write file */
        success = VWriteFile (file, list);
        fclose (file);
    }
    ioMutex.unlock();

   if (!success)
   { VError ("writeVista(): Failed to write output file '%s'", Name);
   }

   /* remove images */
   for (VFirstAttr (list, &pos); VAttrExists (&pos); VNextAttr (&pos))
      if (VGetAttrRepn (&pos) == VImageRepn)
         VSetAttrValue (&pos, NULL, VImageRepn, NULL);

   /* clean-up*/
   VDestroyAttrList (list);

    return 1;

} // end "writeVista()" -----------------------------------------------------------------

