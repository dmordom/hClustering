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

#include "niftiManager.h"
#include <nifti/nifti1_io.h>


// PUBLIC MEMBERS


std::string niftiManager::getFileExtension( TractExtType extType ) const
{
    if( extType == ETFull )
    {
        return NIFTI_EXT;
    }
    else if( extType == ETCompact )
    {
        return COMPACT_EXT;
    }
    else
    {
        throw std::runtime_error("ERROR @ niftiManager::getFileExtension(): extension type not recognized");
    }
}

/*
// niftiManager::getTractFilename(): returns tract filename
std::string niftiManager::getLeafTractFilename (const WHcoord& tractCoord ) const
{
    throw std::runtime_error( "niftiManager::getLeafTractFilename(): getLeafTractFilename(const WHcoord& tractCoord ) not implemented in nifti mode" );
    return "";
}// end niftiManager::getTractFilename() -----------------------------------------------------------------
*/
std::string niftiManager::getLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector ) const
{
    if( tractLeaf >= indexVector.size() )
    {
        throw std::runtime_error( "getLeafTractFilename(): leaf ID is higher than index vector" );
        return "";
    }
    std::string filename = str(boost::format( NIFTI_LEAF_COMPACT_FNAME ) % indexVector[tractLeaf] );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETCompact );
    return filepath;
}// end niftiManager::getTractFilename() -----------------------------------------------------------------

/*
// niftiManager::getFullTractFilename(): returns tract filename
std::string niftiManager::getFullLeafTractFilename(const WHcoord& tractCoord) const
{
    throw std::runtime_error( "niftiManager::getFullLeafTractFilename(): getFullLeafTractFilename(const WHcoord& tractCoord ) not implemented in nifti mode" );
    return "";
}// end niftiManager::getFullTractFilename() -----------------------------------------------------------------
*/

std::string niftiManager::getFullLeafTractFilename( const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector  ) const
{
    if( tractLeaf >= indexVector.size() )
    {
        throw std::runtime_error( "getFullLeafTractFilename(): leaf ID is higher than index vector" );
        return "";
    }
    std::string filename = str(boost::format( NIFTI_LEAF_FULL_FNAME ) % indexVector[tractLeaf] );
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETFull );
    return filepath;
}// end niftiManager::getFullTractFilename() -----------------------------------------------------------------



// niftiManager::readVector():reads 1D nifti vector
ValueType niftiManager::readVector(const std::string& vectorFilenameRef, std::vector<float>* vectorPointer) const
{
    std::vector<float>& vector = *vectorPointer;
    std::string vectorFilename( vectorFilenameRef );
    std::string vectorFilenameZipped( vectorFilenameRef );
    std::string extension;

    std::string zipExtension( boost::filesystem::path(vectorFilename).extension().string() );


    if ( zipExtension == ".gz" )
    {
        vectorFilename.resize(vectorFilename.size()-3);
        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + vectorFilenameZipped + " > " + vectorFilename );
        int success(system(sysCommand.c_str()));
        extension = ( boost::filesystem::path(vectorFilename).extension().string() );

    }
    else if ( zipExtension == getFileExtension( ETFull ) || zipExtension == getFileExtension( ETCompact ) )
    {
        extension = zipExtension;
    }
    else
    {
        std::cerr << "File \"" << vectorFilenameZipped << "\" has no recognized extension (\"" << zipExtension << "\") stopping." << std::endl;
        return VTError;
    }


    ValueType vectorValueType;

    if( extension == getFileExtension( ETFull ) )
    {
        nifti_image*  niftiImage = NULL;

        boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
        ioMutex.lock();
        {
            niftiImage = nifti_image_read(vectorFilename.c_str(), 1);
        }
        ioMutex.unlock();

        if( !niftiImage )
        {
            std::cerr<< "ERROR @ niftiManager::readVector(): there was an error calling nifti_image_read() on vector "<< vectorFilename <<std::endl;
            return VTError;
        }

        if( niftiImage->ndim != 1 )
        {
            std::cerr << "ERROR @ niftiManager::readVector(): nifti file has more than 1 dimension (" << niftiImage->ndim <<"): "<< niftiImage->nx <<" "<< niftiImage->ny <<" "<< niftiImage->nz <<" "<< niftiImage->nt <<" "<< niftiImage->nu <<" "<< niftiImage->nv <<std::endl;
            return VTError;
        }
        const size_t repType( niftiImage->datatype );
        const size_t repByteSize( niftiImage->nbyper );
        const size_t dataSize( niftiImage->nvox );


        if( repType == DT_FLOAT32 )
        {
            vectorValueType = VTFloat32;
        }
        else if ( repType == DT_UINT8 )
        {
            vectorValueType = VTUINT8;
        }
        else
        {
            std::cerr<< "ERROR @ niftiManager::readVector(): vector representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
            return VTError;
        }

        const size_t dimx( niftiImage->nx );
        if( dimx != dataSize )
        {
            std::cerr<< "ERROR @ niftiManager::readVector(): vector dimension does not match data size" << std::endl;
            return VTError;
        }

        vector.clear();
        vector.resize(dataSize,0);

        void* niftiData(niftiImage->data);

        for (int i=0 ; i<dataSize ; ++i) {

            size_t byteOffset( i * repByteSize );
            if( byteOffset >= dataSize )
            {
                std::cerr<< "ERROR @ niftiManager::readVector():: pointer offset was too high when loading nifti data" << std::endl;
                return VTError;
            }
            void* dataPos = static_cast<unsigned char*>(niftiData) + byteOffset;

            if (repType == DT_UINT8)
            {
                unsigned char datapoint = *(static_cast<unsigned char*>(dataPos));
                vector[i] = (float)datapoint;
            }
            else if (repType == DT_FLOAT32)
            {
                float datapoint = *(static_cast<float*>(dataPos));
                vector[i] = datapoint;
            }
            else
            {
                std::cerr<< "ERROR @ niftiManager::readVector(): vector representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
                return VTError;
            }
        }
        // clean up
        nifti_image_free( niftiImage );
    }
    else if ( extension == getFileExtension( ETCompact ) )
    {

        std::ifstream inFile;
        inFile.open( vectorFilename.c_str(), std::ios::in | std::ios::binary | std::ios::ate );
        if( !inFile.is_open() )
        {
            std::stringstream errormessage;
            errormessage << "ERROR @ niftiManager::readVector(): there was an error opening the input image file "<< vectorFilename <<std::endl;
            throw std::runtime_error(  errormessage.str() );
        }
        size_t dataSize = inFile.tellg();
        void* data = malloc( dataSize );
        inFile.seekg( 0, std::ios::beg );

        boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
        ioMutex.lock();
        {
            inFile.read( static_cast<char*>(data), dataSize );
        }
        ioMutex.unlock();

        inFile.close();

        const size_t headerSize = UINT32_SIZE + UINT32_SIZE;
        void* dataPos = static_cast<unsigned char*>(data);
        const size_t repType = *(static_cast<uint32_t*>(dataPos));
        dataPos = static_cast<unsigned char*>(data) + UINT32_SIZE;
        const size_t dimx = *(static_cast<uint32_t*>(dataPos));
        size_t repByteSize( 0 );

        if( repType == (FLOAT_SIZE * CHAR_BIT) )
        {
            vectorValueType = VTFloat32;
            repByteSize = FLOAT_SIZE;

        }
        else if ( repType == (UINT8_SIZE * CHAR_BIT) )
        {
            vectorValueType = VTUINT8;
            repByteSize = UINT8_SIZE;
        }
        else
        {
            std::cerr<< "ERROR @ niftiManager::readVector(): vector representation type not recognized: " << repType << " (not FLOAT32 nor UNIT8)" << std::endl;
            return VTError;
        }

        const size_t dataSizeTest = headerSize + ( dimx * repByteSize );
        if( dataSize != dataSizeTest )
        {
            std::stringstream errormessage;
            errormessage << "ERROR @ niftiManager::readVector(): file data as read: "<< dataSize << " is different as file data according to header: " << dataSizeTest <<std::endl;
            errormessage<<" headerSize: "<< headerSize <<" dimx: "<< dimx <<" repSize: "<<repByteSize << std::endl;
            throw std::runtime_error(  errormessage.str() );
        }

        vector.clear();
        vector.resize(dimx,0);

        for( size_t i = 0; i < dimx; ++i )
        {
            size_t byteOffset( headerSize + ( i * repByteSize ) );
            if( byteOffset >= dataSize )
            {
                std::cerr<< "ERROR @ niftiManager::readVector():: pointer offset was too high when loading nifti data" << std::endl;
                return VTError;
            }
            dataPos = static_cast<unsigned char*>(data) + byteOffset;

            if ( vectorValueType == VTUINT8 )
            {
                uint8_t datapoint = *(static_cast<uint8_t*>(dataPos));
                vector[i] = (float)datapoint;
            }
            else if ( vectorValueType == VTFloat32 )
            {
                float datapoint = *(static_cast<float*>(dataPos));
                vector[i] = datapoint;
            }
            else
            {
                std::cerr<< "ERROR @ niftiManager::readVector(): vector representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
                return VTError;
            }
        }
        // clean up
        free( data );

    }
    else
    {
        std::cerr << "File \"" << vectorFilename << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return VTError;
    }

    if ( zipExtension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + vectorFilename );
        int success(system(sysCommand.c_str()));
    }
    return vectorValueType;
}// end niftiManager::readVector() ----------------------------------------------------------------


// niftiManager::readMatrix():reads 2D nifti matrix
ValueType niftiManager::readMatrix( const std::string& matrixFilenameRef, std::vector<std::vector<float> >* matrixPointer ) const
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

    nifti_image*  niftiImage = NULL;

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        niftiImage = nifti_image_read(matrixFilename.c_str(), 1);
    }
    ioMutex.unlock();

    if( !niftiImage )
    {
        std::cerr<< "ERROR @ niftiManager::readMatrix(): there was an error calling nifti_image_read() on matrix "<< matrixFilename <<std::endl;
        return VTError;
    }
    if( niftiImage->ndim > 2 )
    {
        std::cerr << "ERROR @ niftiManager::readMatrix(): nifti file has more than 2 dimensions (" << niftiImage->ndim <<"): "<< niftiImage->nx <<" "<< niftiImage->ny <<" "<< niftiImage->nz <<" "<< niftiImage->nt <<" "<< niftiImage->nu <<" "<< niftiImage->nv <<std::endl;
        return VTError;
    }
    const size_t repType( niftiImage->datatype );
    const size_t repByteSize( niftiImage->nbyper );
    const size_t dataSize( niftiImage->nvox );

    ValueType matrixValueType;

    if( repType == DT_FLOAT32 )
    {
        matrixValueType = VTFloat32;
    }
    else if ( repType == DT_UINT8 )
    {
        matrixValueType = VTUINT8;
    }
    else
    {
        std::cerr<< "ERROR @ niftiManager::readMatrix(): matrix representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
        return VTError;
    }

    const size_t dimx( niftiImage->nx );
    const size_t dimy( niftiImage->ny );

    matrix.clear();
    {
        std::vector<float> yvect(dimy,0);
        matrix.resize(dimx,yvect);
    }
    void* niftiData(niftiImage->data);

    for (int i=0 ; i<dimy ; ++i)
    {
        for (int j=0 ; j<dimx ; ++j)
        {
            size_t voxOffset( (i*dimx) + j );
            if( voxOffset >= dataSize )
            {
                std::cerr<< "ERROR @ niftiManager::readDistBlock():: pointer offset was too high when loading nifti data" << std::endl;
                return VTError;
            }
            size_t byteOffset( voxOffset * repByteSize );
            void* dataPos = static_cast<unsigned char*>(niftiData) + byteOffset;

            if( repType == DT_FLOAT32 )
            {
                float datapoint = *(static_cast<float*>(dataPos));
                matrix[j][i] = datapoint;
            }
            else if ( repType == DT_UINT8 )
            {
                unsigned char datapoint = *(static_cast<unsigned char*>(dataPos));
                matrix[j][i] = (float)datapoint;
            }
            else
            {
                std::cerr<< "ERROR @ niftiManager::readMatrix(): matrix representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
                return VTError;
            }
        }
    }
    // clean up
    nifti_image_free( niftiImage );

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + matrixFilename );
        int success(system(sysCommand.c_str()));
    }
    return matrixValueType;
}// end niftiManager::readMatrix() -----------------------------------------------------------------


// "readImage()": reads a 3D image
ValueType niftiManager::readImage( const std::string& imageFilenameRef, std::vector<std::vector<std::vector<float> > >* imagePointer ) const
{
    std::vector<std::vector<std::vector<float> > >& image( *imagePointer );
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

    nifti_image*  niftiImage = NULL;

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        niftiImage = nifti_image_read(imageFilename.c_str(), 1);
    }
    ioMutex.unlock();

    if( !niftiImage )
    {
        std::cerr<< "ERROR @ niftiManager::readImage(): there was an error calling nifti_image_read() on image file "<< imageFilename <<std::endl;
        return VTError;
    }

    if( niftiImage->ndim > 3 )
    {
        std::cerr << "ERROR @ niftiManager::readImage(): nifti file has more than 3 dimensions (" << niftiImage->ndim <<"): "<< niftiImage->nx <<" "<< niftiImage->ny <<" "<< niftiImage->nz <<" "<< niftiImage->nt <<" "<< niftiImage->nu <<" "<< niftiImage->nv <<std::endl;
        return VTError;
    }
    const size_t repType( niftiImage->datatype );
    const size_t repByteSize( niftiImage->nbyper );
    const size_t dataSize( niftiImage->nvox );

    ValueType imageValueType;

    if(repType == DT_UINT8 )
    {
        imageValueType = VTUINT8;
    }
    else if(repType == DT_FLOAT32 )
    {
        imageValueType = VTFloat32;
    }
    else
    {
        std::cerr<< "ERROR @ niftiManager::readImage(): image representation type not recognized (neither UINT8 nor FLOAT32)" << std::endl;
        return VTError;
    }

    const size_t dimx( niftiImage->nx );
    const size_t dimy( niftiImage->ny );
    const size_t dimz( niftiImage->nz);

    image.clear();
    {
        std::vector<float> zvector(dimz,0);
        std::vector<std::vector<float> >yzmatrix(dimy,zvector);
        image.resize(dimx,yzmatrix);
    }
    void* niftiData(niftiImage->data);

    for (int i=0 ; i<dimz ; ++i)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            for (int k=0 ; k<dimx ; ++k)
            {
                size_t voxOffset( (i*dimy*dimx) + (j*dimx) + k );
                if( voxOffset >= dataSize )
                {
                    std::cerr<< "ERROR @ niftiManager::readImage():: pointer offset was too high when loading nifti data" << std::endl;
                    return VTError;
                }
                size_t byteOffset( voxOffset * repByteSize );
                void* dataPos = static_cast<unsigned char*>(niftiData) + byteOffset;

                if (repType == DT_UINT8)
                {
                    unsigned char datapoint = *(static_cast<unsigned char*>(dataPos));
                    image[k][j][i] = (float)datapoint;
                }
                else if (repType == DT_FLOAT32)
                {
                    float datapoint = *(static_cast<float*>(dataPos));
                    image[k][j][i] = datapoint;
                }
                else
                {
                    std::cerr<< "ERROR @ niftiManager::readImage(): image representation type not recognized (neither UINT8 nor FLOAT32)" << std::endl;
                    return VTError;
                }
            }
        }
    }
    // clean up
    nifti_image_free( niftiImage );

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + imageFilename );
        int success(system(sysCommand.c_str()));
    }
    return imageValueType;
}// end niftiManager::readImage() -----------------------------------------------------------------


// "loadHeader()": load header info into member variable
void niftiManager::loadHeader( const std::string& imageFilename, bool display )
{
    m_header = readHeader( imageFilename );
    if( display )
    {
        displayHeader();
    }
    return;
}// end niftiManager::loadHeader() -----------------------------------------------------------------




// "writeVector()": writes a 1D vector in binary format
void niftiManager::writeVector( const std::string& vectorFilename, const ValueType dataValueType, const std::vector<float>& vector, bool doZip ) const
{
    if (vector.empty())
    {
        std::cerr<< "ERROR @ niftiManager::writeVector(): vector is empty, image has not been written" <<std::endl;
        return;
    }
    const size_t dimx = vector.size();
    size_t repByteSize( 0 );
    size_t repType( 0 );


    if( dataValueType == VTFloat32 )
    {
        repType = FLOAT_SIZE * CHAR_BIT;
        repByteSize = FLOAT_SIZE;
    }
    else if ( dataValueType == VTUINT8 )
    {
        repType = UINT8_SIZE * CHAR_BIT;
        repByteSize = UINT8_SIZE ;
    }
    else
    {
        std::cerr<< "ERROR @ niftiManager::writeVector(): vector representation type not recognized (not FLOAT32 nor UNIT8)" << std::endl;
        return;
    }

    const size_t headerSize = UINT32_SIZE + UINT32_SIZE;
    const size_t dataSize = headerSize + (dimx*repByteSize);
    void* data = malloc (dataSize );

    //write header
    void* dataPos = static_cast<unsigned char*>(data);
    *(static_cast<uint32_t*>(dataPos)) = repType;
    dataPos = static_cast<unsigned char*>(data) + UINT32_SIZE;
    *(static_cast<uint32_t*>(dataPos)) = dimx;

    //write data
    for (int i=0 ; i<dimx ; ++i)
    {
        size_t dataOffset = headerSize + ( i * repByteSize );
        if( dataOffset >= dataSize )
        {
            throw std::runtime_error( "ERROR @ niftiManager::writeVector():: pointer offset was too high when loading nifti data");
        }

        dataPos = static_cast<unsigned char*>(data) + dataOffset;

        if (dataValueType == VTUINT8)
        {
            unsigned char datapoint = (unsigned char)vector[i];
            *(static_cast<uint8_t*>(dataPos)) = datapoint;
        }
        else if (dataValueType == VTFloat32)
        {
            float datapoint = vector[i];
            *(static_cast<float*>(dataPos)) = datapoint;
        }
        else
        {
            throw std::runtime_error("ERROR @ niftiManager::writeVector(): image representation type not recognized (neither UINT8 nor FLOAT32)");
        }
    }



    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        std::ofstream outFile;
        outFile.open( vectorFilename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary );

        if( !outFile.is_open() )
        {
            std::stringstream errormessage;
            errormessage << "ERROR @ niftiManager::writeVector(): there was an error opening the output image file "<< vectorFilename <<std::endl;
            throw std::runtime_error(  errormessage.str() );
        }
        outFile.write( static_cast<const char*>(data), dataSize );
        outFile.close();
    }
    ioMutex.unlock();


    // clean up
    free( data );


    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + vectorFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end niftiManager::writeVector() -----------------------------------------------------------------




// "writeMatrix()": writes a 2D matrix
void niftiManager::writeMatrix( const std::string& matrixFilename, const ValueType dataValueType, const std::vector<std::vector<float> >& matrix, bool doZip ) const
{
    if (matrix.empty())
    {
        std::cerr<< "ERROR @ niftiManager::writeMatrix(): matrix is empty, image has not been written" <<std::endl;
        return;
    }
    const size_t dimx = matrix.size();
    const size_t dimy = matrix[0].size();

    nifti_1_header matrixHeader;
    matrixHeader.sizeof_hdr = 0;
    generateHeader( dimx, dimy, 1, dataValueType, &matrixHeader);

    nifti_image* niftiImage = nifti_convert_nhdr2nim( matrixHeader, NULL);

    const size_t repByteSize( niftiImage->nbyper );
    const size_t dataSize( niftiImage->nvox );
    const size_t repType( niftiImage->datatype );

    if ( dataSize != (dimx*dimy) )
    {
        throw std::runtime_error( "niftiManager::writeMatrix(): matrix header nvox and matrix size do not match" );
    }
    if ( repType != matrixHeader.datatype )
    {
        throw std::runtime_error( "niftiManager::writeMatrix(): created header datatype and input argument dont match" );
    }

    if(niftiImage->data != NULL)
    {
        throw std::runtime_error( "niftiManager::writeVector(): nifti_convert_nhdr2nim had already allocated memory" );
    }
    niftiImage->data = malloc (dataSize * repByteSize);

    for (int i=0 ; i<dimy ; ++i)
    {
        for (int j=0 ; j<dimx ; ++j)
        {
            size_t voxOffset( (i*dimx) + j );
            if( voxOffset >= dataSize )
            {
                throw std::runtime_error( "ERROR @ niftiManager::writeMatrix():: pointer offset was too high when loading nifti data");
            }
            size_t byteOffset( voxOffset * repByteSize );
            void* dataPos = static_cast<unsigned char*>(niftiImage->data) + byteOffset;

            if (repType == DT_UINT8)
            {
                unsigned char datapoint = (unsigned char)matrix[j][i];
                *(static_cast<unsigned char*>(dataPos)) = datapoint;
            }
            else if (repType == DT_FLOAT32)
            {
                float datapoint = matrix[j][i];
                *(static_cast<float*>(dataPos)) = datapoint;
            }
            else
            {
                throw std::runtime_error("ERROR @ niftiManager::writeMatrix(): image representation type not recognized (neither UINT8 nor FLOAT32)");
            }
        }
    }

        if( nifti_set_filenames(niftiImage, matrixFilename.c_str(), 0, 1) )
    {
        std::stringstream errormessage;
        errormessage << "ERROR @ niftiManager::writeMatrix(): there was an error calling nifti_set_filenames() on output image file "<< matrixFilename <<std::endl;
        throw std::runtime_error(  errormessage.str() );
    }

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        nifti_image_write( niftiImage );
    }
    ioMutex.unlock();

    // clean up
    nifti_image_free( niftiImage );

    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + matrixFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end niftiManager::writeMatrix() -----------------------------------------------------------------



// "writeImage()": writes a 3D image
void niftiManager::writeImage( const std::string& imageFilename, const ValueType dataValueType, const std::vector<std::vector<std::vector<float> > >& image, bool doZip ) const
{
    if (image.empty())
    {
        std::cerr<< "ERROR @ niftiManager::writeImage(): image matrix is empty, image has not been written" <<std::endl;
        return;
    }
    const size_t dimx = image.size();
    const size_t dimy = image[0].size();
    const size_t dimz = image[0][0].size();

    nifti_1_header imageHeader;
    imageHeader.sizeof_hdr = 0;
    if(m_header.sizeof_hdr != 0)
    {
        imageHeader = m_header;
    }
    generateHeader( dimx, dimy, dimz, dataValueType, &imageHeader);

    nifti_image* niftiImage = nifti_convert_nhdr2nim( imageHeader, NULL);

    const size_t repByteSize( niftiImage->nbyper );
    const size_t dataSize( niftiImage->nvox );
    const size_t repType( niftiImage->datatype );

    if ( dataSize != (dimx*dimy*dimz) )
    {
        throw std::runtime_error( "niftiManager::writeImage(): image header nvox and matrix size do not match" );
    }
    if ( repType != imageHeader.datatype )
    {
        throw std::runtime_error( "niftiManager::writeImage(): created header datatype and input argument dont match" );
    }

    if(niftiImage->data != NULL)
    {
        throw std::runtime_error( "niftiManager::writeVector(): nifti_convert_nhdr2nim had already allocated memory" );
    }
    niftiImage->data = malloc (dataSize * repByteSize);

    for (int i=0 ; i<dimz ; ++i)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            for (int k=0 ; k<dimx ; ++k)
            {
                size_t voxOffset( (i*dimy*dimx) + (j*dimx) + k );
                if( voxOffset >= dataSize )
                {
                    throw std::runtime_error( "ERROR @ niftiManager::writeImage():: pointer offset was too high when loading nifti data");
                }
                size_t byteOffset( voxOffset * repByteSize );
                void* dataPos = static_cast<unsigned char*>(niftiImage->data) + byteOffset;
                if (repType == DT_UINT8)
                {
                    unsigned char datapoint = (unsigned char)image[k][j][i];
                    *(static_cast<unsigned char*>(dataPos)) = datapoint;
                }
                else if (repType == DT_FLOAT32)
                {
                    float datapoint = image[k][j][i];
                    *(static_cast<float*>(dataPos)) = datapoint;
                }
                else
                {
                     throw std::runtime_error("ERROR @ niftiManager::writeImage(): image representation type not recognized (neither UINT8 nor FLOAT32)");
                }
            }
        }
    }


    if( nifti_set_filenames(niftiImage, imageFilename.c_str(), 0, 1) )
    {
        std::stringstream errormessage;
        errormessage << "ERROR @ niftiManager::writeImage(): there was an error calling nifti_set_filenames() on output image file "<< imageFilename <<std::endl;
        throw std::runtime_error(  errormessage.str() );
    }

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        nifti_image_write( niftiImage );
    }
    ioMutex.unlock();

    // clean up
    nifti_image_free( niftiImage );

    if (doZip)
    { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + imageFilename);
        int success(system(sysCommand.c_str()));
    }
    return;
}// end niftiManager::writeImage() -----------------------------------------------------------------




//////////////// PRIVATE MEMBER FUNCTIONS ///////////////////



// "displayHeader()": display headerinfo on standard output
void niftiManager::displayHeader() const
{
    if(m_header.sizeof_hdr == 0)
    {
        std::cerr<< "WARNING @ niftiManager::displayHeader():no saved image header found. cannot write info... " <<std::endl;
        return;
    }
    char*  infoChar;
    disp_nifti_1_header( infoChar, &m_header);
//    std::string infoString(infoChar);
//    std::cout << std::endl << std::endl << " DONE. WRITING INFOSTRING TO STD::COUT..." << std::endl << std::endl;
//    std::cout << infoString << std::endl;
    return;
}// end niftiManager::displayHeader() -----------------------------------------------------------------


// "readHeader()": read nifti header information from file
nifti_1_header niftiManager::readHeader( const std::string& imageFilenameRef ) const
{
    std::string imageFilenameZipped( imageFilenameRef );
    std::string imageFilename( imageFilenameRef );

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
        exit(1);
    }

    nifti_image*  niftiImage = NULL;

    boost::mutex& ioMutex(boost::detail::thread::singleton<boost::mutex>::instance()) ;
    ioMutex.lock();
    {
        niftiImage = nifti_image_read(imageFilename.c_str(), 0);  // load_data = 0, only header and not image data is read
    }
    ioMutex.unlock();

    if( !niftiImage )
    {
        std::stringstream errormessage;
        errormessage << "ERROR @ niftiManager::readHeader(): there was an error calling nifti_image_read() on image file "<< imageFilename <<std::endl;
        throw std::runtime_error(  errormessage.str() );
    }

    if( niftiImage->ndim > 3 )
    {
        std::stringstream errormessage;
        errormessage << "ERROR @ niftiManager::readHeader(): nifti file has more than 3 dimensions (" << niftiImage->ndim <<"): "<< niftiImage->nx <<" "<< niftiImage->ny <<" "<< niftiImage->nz <<" "<< niftiImage->nt <<" "<< niftiImage->nu <<" "<< niftiImage->nv <<std::endl;
        throw std::runtime_error(  errormessage.str() );
    }

    nifti_1_header header = nifti_convert_nim2nhdr( niftiImage );
    // clean up
    nifti_image_free( niftiImage );

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + imageFilename );
        int success(system(sysCommand.c_str()));
    }
    return header;
}// end niftiManager::readHeader() -----------------------------------------------------------------


// "niftiManager::generateHeader()": generates a header for writing a nifti image
void niftiManager::generateHeader( const size_t dimx, const size_t dimy, const size_t dimz, const ValueType dataValueType, nifti_1_header* imageHeaderPointer ) const
{
    if( dimx > SHRT_MAX || dimy > SHRT_MAX || dimz > SHRT_MAX )
    {
        std::stringstream errormessage;
        errormessage << "ERROR @ niftiManager::generateHeader(): dimensions are too long for a nifti image (" << dimx <<" "<< dimy <<" "<< dimz <<") MAX value for any dimension is: " << SHRT_MAX <<std::endl;
        throw std::runtime_error(  errormessage.str() );
    }

    nifti_1_header& imageHeader = *( imageHeaderPointer );
    if( dimy == 1 || dimz ==1 )
    {
        imageHeader.sizeof_hdr = 0;
    }

    if(imageHeader.sizeof_hdr != 0)
    {
        if(imageHeader.dim[0] > 3 )
        {
            std::stringstream errormessage;
            errormessage << "ERROR @ niftiManager::generateHeader(): nifti stored header more than 3 dimensions (" << imageHeader.dim[0] <<"): "<< imageHeader.dim[1] <<" "<< imageHeader.dim[2] <<" "<< imageHeader.dim[3] <<" "<< imageHeader.dim[4] <<" "<< imageHeader.dim[5] <<" "<<imageHeader.dim[6] <<std::endl;
            throw std::runtime_error(  errormessage.str() );
        }
        if( ( imageHeader.dim[3] != dimz ) || ( imageHeader.dim[2] != dimy ) || ( imageHeader.dim[1] != dimx ) )
        {
            std::stringstream errormessage;
            errormessage << "ERROR @ niftiManager::generateHeader(): image matrix and stored header dimensions dont match (" << dimx <<" "<< dimy <<" "<< dimz <<") ("<< imageHeader.dim[1] <<" "<< imageHeader.dim[2] <<" "<< imageHeader.dim[3] <<")" <<std::endl;
            throw std::runtime_error(  errormessage.str() );
        }
    }
    else
    {
        if( dimx > 1 && dimy > 1 && dimz > 1)
        {
            std::cerr<< "====="<<std::endl;
            std::cerr<< "WARNING @ niftiManager::generateHeader():writing a 3D image but no saved 3D image header found."<< std::endl;
            std::cerr<< "Creating a new header from scratch. Output image might have wrong orientation when read by another program ... " <<std::endl;
            std::cerr<< "====="<<std::endl;
        }
        imageHeader.sizeof_hdr = 348;
        for( size_t i=0; i<10; ++i)
        {
            imageHeader.data_type[i] = 0;
        }
        for( size_t i=0; i<18; ++i)
        {
            imageHeader.db_name[i] = 0;
        }
        imageHeader.extents = 0;
        imageHeader.session_error = 0;
        imageHeader.regular = 0;
        imageHeader.dim_info = NIFTI_SLICE_UNKNOWN;
        imageHeader.dim[0] = 3;
        if( dimz == 1 )
        {
            imageHeader.dim[0] = 2;
        }
        if( dimy == 1 )
        {
            imageHeader.dim[0] = 1;
        }
        imageHeader.dim[1] = dimx;
        imageHeader.dim[2] = dimy;
        imageHeader.dim[3] = dimz;
        for(size_t i=4; i<8; ++i)
        {
            imageHeader.dim[i] = 1;
        }
        imageHeader.intent_p1 = 0;
        imageHeader.intent_p2 = 0;
        imageHeader.intent_p3 = 0;
        imageHeader.intent_code = 0;
        imageHeader.slice_start = 0;
        imageHeader.pixdim[0] = 1;
        imageHeader.pixdim[1] = 1;
        imageHeader.pixdim[2] = 0;
        imageHeader.pixdim[3] = 0;
        if( dimy > 1 )
        {
            imageHeader.pixdim[2] = 1;
        }
        if( dimz > 1 )
        {
            imageHeader.pixdim[3] = 1;
        }
        for(size_t i=4; i<8; ++i)
        {
            imageHeader.pixdim[i] = 0;
        }
        imageHeader.vox_offset = 0;
        imageHeader.scl_slope = 0;
        imageHeader.scl_inter = 0;
        imageHeader.xyzt_units = NIFTI_UNITS_MM;
        imageHeader.cal_max = 0;
        imageHeader.cal_min = 0;
        imageHeader.slice_duration = 1;
        imageHeader.toffset = 0;
        imageHeader.glmax = 0;
        imageHeader.glmin = 0;
        imageHeader.descrip[0] = '\0';
        imageHeader.aux_file[0] = '\0';
        imageHeader.qform_code = NIFTI_XFORM_UNKNOWN;
        imageHeader.sform_code = 0;
        imageHeader.quatern_b = 0;
        imageHeader.quatern_c = 0;
        imageHeader.quatern_d = 0;
        imageHeader.qoffset_x = 0;
        imageHeader.qoffset_y = 0;
        imageHeader.qoffset_z = 0;
        for( size_t i=0; i<4; ++i)
        {
            imageHeader.srow_x[i] = 0;
            imageHeader.srow_y[i] = 0;
            imageHeader.srow_z[i] = 0;
        }
        imageHeader.intent_name[0] = '\0';
        imageHeader.magic[0] = '\0';
    }


    if( dataValueType == VTBit )
    {
        imageHeader.datatype = DT_UINT8;
        imageHeader.bitpix= 8;
    }
    else if(dataValueType == VTUINT8 )
    {
        imageHeader.datatype = DT_UINT8;
        imageHeader.bitpix= 8;
    }
    else if(dataValueType == VTFloat32 )
    {
        imageHeader.datatype = DT_FLOAT32;
        imageHeader.bitpix = 32;
    }
    else
    {
        throw std::runtime_error( "niftiManager::generateHeader(): image representation type not recognized (neither VFloat nor VUByte nor VBit)" );
    }
    return;
}// end niftiManager::generateHeader() -----------------------------------------------------------------
