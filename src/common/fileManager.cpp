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


#include "fileManager.h"

#include <errno.h>

// PUBLIC MEMBERS

// fileManager::getTractFilename(): returns tract filename
std::string fileManager::getNodeTractFilename(const size_t tractNode ) const
{
    std::string filename = str(boost::format( NODE_COMPACT_FNAME ) % tractNode);
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETCompact );
    return filepath;
}// end fileManager::getTractFilename() -----------------------------------------------------------------

std::string fileManager::getFullNodeTractFilename(const size_t tractNode ) const
{
    std::string filename = str(boost::format( NODE_FULL_FNAME ) % tractNode);
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETFull );
    return filepath;
}// end fileManager::getFullTractFilename() -----------------------------------------------------------------

std::string fileManager::getClusterMaskFilename(const size_t node ) const
{
    std::string filename = str(boost::format( NODE_CLUSTER_FNAME ) % node);
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETFull);
    return filepath;
}// end fileManager::getClusterMaskFilename() -----------------------------------------------------------------

std::string fileManager::getBlockFilename(const unsigned int blockID1, const unsigned int blockID2) const
{
    std::string filename = str(boost::format( DISTBLOCK_FNAME ) % blockID1 % blockID2);
    std::string filepath = m_ioFolder + "/" + filename + getFileExtension( ETFull);
    return filepath;
}// end fileManager::getBlockFilename() -----------------------------------------------------------------


size_t fileManager::readClusterMaskNumber (const std::string filename) const
{
    std::string fileNameStem( boost::filesystem::path(filename).stem().string() );

    std::string clusterNumberString( fileNameStem.begin() + NODE_CLUSTER_PREFIX_SIZE, fileNameStem.end() );
    size_t clusterNumber( boost::lexical_cast<size_t>( clusterNumberString ));
    return clusterNumber;
}// end fileManager::readClusterMaskNumber() -----------------------------------------------------------------

size_t fileManager::readFullNodeTractNumber (const std::string filename) const
{
    std::string fileNameStem( boost::filesystem::path(filename).stem().string() );

    std::string clusterNumberString( fileNameStem.begin() + NODE_FULL_PREFIX_SIZE, fileNameStem.end() );
    size_t clusterNumber( boost::lexical_cast<size_t>( clusterNumberString ));
    return clusterNumber;
}// end fileManager::readClusterMaskFilename() -----------------------------------------------------------------

/*
// "fileManager::readLeafTract()": extracts compact tractogram data from a file and returns it in compact Tract form
void fileManager::readLeafTract (const WHcoord tractCoord, compactTract* tractPointer ) const
{
    compactTract& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename(tractCoord);
    readTract( tractFilename, &( tractogram.m_tract ) );
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}
void fileManager::readLeafTract (const WHcoord tractCoord, compactTractChar* tractPointer) const
{
    compactTractChar& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename(tractCoord);
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readLeafTract() -----------------------------------------------------

void fileManager::readLeafTract (const size_t tractLeaf, compactTract* tractPointer) const
{
    compactTract& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename(tractLeaf);
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;}//end fileManager::readLeafTract() -----------------------------------------------------

void fileManager::readLeafTract (const size_t tractLeaf, compactTractChar* tractPointer) const
{
    compactTractChar& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename(tractLeaf);
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readLeafTract() -----------------------------------------------------
*/

void fileManager::readLeafTract (const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, compactTract* tractPointer) const
{
    compactTract& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename(tractLeaf, indexVector, coordVector );
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readLeafTract() -----------------------------------------------------
void fileManager::readLeafTract (const size_t tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, compactTractChar* tractPointer) const
{
    compactTractChar& tractogram = *tractPointer;
    std::string tractFilename = getLeafTractFilename( tractLeaf, indexVector,  coordVector );
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readLeafTract() -----------------------------------------------------

// "fileManager::readNodeTract()": extracts compact tractogram data from a file and returns it in compact Tract form
void fileManager::readNodeTract (const size_t tractNode, compactTract* tractPointer) const
{
    compactTract& tractogram = *tractPointer;
    std::string tractFilename = getNodeTractFilename(tractNode);
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readNodeTract() -----------------------------------------------------


void fileManager::readTract (const std::string& tractFilename, compactTract* tractPointer) const
{
    compactTract& tractogram = *tractPointer;
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readNodeTract() -----------------------------------------------------

void fileManager::readTract (const std::string& tractFilename, compactTractChar* tractPointer) const
{
    compactTractChar& tractogram = *tractPointer;
    readTract(tractFilename, &( tractogram.m_tract ) );
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end fileManager::readNodeTract() -----------------------------------------------------


// "fileManager::writeNodeTract()": write tractogram data into file
void fileManager::writeNodeTract (const size_t tractNode, const compactTract& tractogram) const
{
    std::string tractFilename = getNodeTractFilename(tractNode);
    writeTract(tractFilename, tractogram);
    return;
}//end fileManager::writeNodeTract() -----------------------------------------------------


// "deleteTractFile()": delete tractogram file
int fileManager::deleteTractFile ( const size_t tractNode ) const
{

    std::string tractFilename = getNodeTractFilename(tractNode);

    // Variables for calling system commands
    std::string delString("rm -f ");
    std::string sysCommand(delString +tractFilename);

    if(m_zipFlag)
    {
        sysCommand = sysCommand +".gz";
    }

    // remove file
    return system(sysCommand.c_str());

}// end fileManager::deleteTractFile(" -----------------------------------------------------------------




// "full2compact()": compacts a 3D image tract into vector form
void fileManager::full2compact(const std::vector< std::vector< std::vector <float> > >& fullTract, std::vector<float>* compactPointer ) const
{
    if (m_maskMatrix.empty())
    {
        throw std::runtime_error("ERROR @ fileManager::full2compact(): Mask hast not been loaded,  tract has not been processed");
    }
    const size_t dimx = fullTract.size();
    const size_t dimy = fullTract[0].size();
    const size_t dimz = fullTract[0][0].size();
    if ( dimx != m_maskMatrix.size() || dimy != m_maskMatrix[0].size() || dimz!= m_maskMatrix[0][0].size() )
    {
        throw std::runtime_error("ERROR @ fileManager::full2compact(): mask and full tract dimensions do not match");
    }
    std::vector<float>& compact( *compactPointer );
    compact.clear();
    compact.reserve( getTractSize() );
    size_t inBounds( 0 ), outOfBounds( 0 );
    for( int i=0 ; i<dimz ; ++i )
    {
        for( int j=0 ; j<dimy ; ++j )
        {
            for( int k=0 ; k<dimx ; ++k )
            {
                float datapoint( fullTract[k][j][i] );
                if( m_maskMatrix[k][j][i] )
                {
                    compact.push_back( datapoint );
                    if( datapoint )
                    {
                        ++inBounds;
                    }
                }
                else
                {
                    if( datapoint )
                    {
                        ++outOfBounds;
                    }
                }
            }
        }
    }

    if( outOfBounds > 0 )
    {
        std::cerr << "WARNING @ fileManager::full2compact(): full tract had " << outOfBounds << " non-zero voxels outside of mask (and " << inBounds << " non-zero voxels within the mask)" << std::endl;
    }
    if( getTractSize() != compact.size() )
    {
        std::cerr << "Tract lenght: " << compact.size() << ". Mask length: " << getTractSize() << std::endl;
        throw std::runtime_error( "fileManager::full2compact(): Mask and tractogram sizes do not match" );
    }
    return;
}// end fileManager::full2compact() -----------------------------------------------------------------




// "compact2full()": blows a compact vector tract into a full 3D image
void fileManager::compact2full(const std::vector<float>& compact, std::vector< std::vector< std::vector <float> > >* fullTractPointer ) const
{
    if (m_maskMatrix.empty())
    {
        throw std::runtime_error("ERROR @ fileManager::storeFullTract(): Mask hast not been loaded, tract has not been processed");
    }
    const size_t dimx   = m_maskMatrix.size();
    const size_t dimy    = m_maskMatrix[0].size();
    const size_t dimz = m_maskMatrix[0][0].size();

    std::vector<std::vector<std::vector<float> > >& fullTract( *fullTractPointer );
    fullTract.resize( dimx, std::vector<std::vector<float> >( dimy, ( std::vector<float>( dimz, 0 ) ) ) );

    std::vector<float>::const_iterator tractIter(compact.begin());
    for( int i=0 ; i<dimz ; ++i )
    {
        for( int j=0 ; j<dimy ; ++j )
        {
            for( int k=0 ; k<dimx ; ++k )
            {
                if (m_maskMatrix[k][j][i])
                {
                    if (tractIter == compact.end())
                    {
                        std::cerr << "Tract lengtt: " << compact.size() << std::endl;
                        throw std::runtime_error( "fileManager::compact2full(): Mask and tractogram sizes do not match, tractogram is shorter" );
                    }
                    fullTract[k][j][i] = *tractIter++;
                }
            }
        }
    }
    if (tractIter != compact.end())
    {
        std::cerr << "Tract lenght: " << compact.size() << ". Mask length: " << getTractSize() << std::endl;
        throw std::runtime_error( "fileManager::compact2full(): Mask and tractogram sizes do not match, tractogram is longer" );
    }
    return;
}// end fileManager::compact2full() -----------------------------------------------------------------



// "writeFullTract()": convert compact tract into full image tractogram and write to file
void fileManager::writeFullNodeTract (const size_t tractNode, const compactTract& tractogram) const
{
    std::string tractFilename = getFullNodeTractFilename(tractNode);
    writeFullTract(tractFilename,tractogram);
    return;
}// end fileManager::writeFullNodeTract() -----------------------------------------------------------------

/*
void fileManager::writeFullLeafTract (const WHcoord& tractCoord, const compactTract& tractogram) const
{
    std::string tractFilename = getFullLeafTractFilename(tractCoord);
    writeFullTract(tractFilename,tractogram);
    return;
}// end fileManager::writeFullLeafTract() -----------------------------------------------------------------
*/
void fileManager::writeFullLeafTract (const size_t& tractLeaf, const std::vector< size_t >& indexVector, const compactTract& tractogram) const
{
    std::string tractFilename = getFullLeafTractFilename(tractLeaf, indexVector);
    writeFullTract(tractFilename,tractogram);
    return;
}// end fileManager::writeFullLeafTract() -----------------------------------------------------------------
void fileManager::writeFullLeafTract (const size_t& tractLeaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTract& tractogram) const
{
    std::string tractFilename = getFullLeafTractFilename(tractLeaf, indexVector, coordVector);
    writeFullTract(tractFilename,tractogram);
    return;
}// end fileManager::writeFullLeafTract() -----------------------------------------------------------------

void fileManager::writeFullTract (const std::string& tractFilename, const compactTract& tractogram) const
{
    std::vector< std::vector< std::vector <float> > > fullTract;
    compact2full( tractogram.m_tract, &fullTract );

    ValueType tractValueType;
    if(m_floatFlag)
    {
        tractValueType = VTFloat32;
    }
    else
    {
        tractValueType = VTUINT8;
        const size_t dimx   = fullTract.size();
        const size_t dimy    = fullTract[0].size();
        const size_t dimz = fullTract[0][0].size();
        for( int i=0 ; i<dimx ; ++i )
        {
            for( int j=0 ; j<dimy ; ++j )
            {
                for( int k=0 ; k<dimz ; ++k )
                {
                    if (fullTract[i][j][k] != 0)
                    {
                        fullTract[i][j][k] = fullTract[i][j][k] * 255.;
                    }
                }
            }
        }
    }
    writeImage( tractFilename, tractValueType, fullTract, m_zipFlag );
    return;
}// end fileManager::writeFullTract() -----------------------------------------------------------------


void fileManager::readFullTract( const std::string fullTractFilename, compactTract* tractPointer ) const
{
    compactTract& tractogram = *tractPointer;
    std::vector<std::vector<std::vector<float> > > fullTractMatrix;
    ValueType fulltractValueType = readImage(fullTractFilename,&fullTractMatrix);

    if ( fulltractValueType == VTError )
    {
        throw std::runtime_error( "ERROR @ fileManager::readFullTract(): Failed to read full tract image" );
    }
    if ( fulltractValueType != VTFloat32 && fulltractValueType != VTUINT8 )
    {
        throw std::runtime_error("ERROR @ fileManager::readFullTract(): tract representation type not recognized (neither Float32 nor UINT8)");
    }

    full2compact( fullTractMatrix, &(tractogram.m_tract) );
    tractogram.m_inLogUnits = m_logFlag;
    tractogram.m_norm = 0;
    tractogram.m_normReady = false;
    tractogram.m_thresholded = m_thresFlag;

    if ( fulltractValueType == VTUINT8 )
    {
        for( int i = 0; i < tractogram.m_tract.size(); ++i )
        {
            tractogram.m_tract[i] = tractogram.m_tract[i]/255.;
        }
    }
    return;
}// end fileManager::readFullTract() -----------------------------------------------------------------














void fileManager::flipXtract(compactTract* tractogramPointer) const
{
    compactTract& tractogram = *tractogramPointer;

    if (m_flipVector.empty()) {
        std::cerr<< "ERROR @ fileManager::flipXtract(): Mask hast not been loaded, tractogram has not been flipped"<<std::endl;
        return;
    }
    if (m_flipVector.size()!=tractogram.size()) {
        std::cerr<< "ERROR @ fileManager::flipXtract(): flip vector and tractogram sizes dont match"<<std::endl;
        return;
    }

    compactTract flippedTract(tractogram);
    flippedTract.m_tract.assign(m_flipVector.size(),0);
    for (size_t i=0; i<m_flipVector.size(); ++i)
        flippedTract.m_tract[i]=tractogram.m_tract[m_flipVector[i]];

    tractogram=flippedTract;

    return;

}// end fileManager::flipXtract() -----------------------------------------------------------------

WHcoord fileManager::meanCoordFromMask() const
{
    if (m_maskMatrix.empty()) {
        std::cerr<< "ERROR @ fileManager::storeFullTract(): Mask hast not been loaded, returning 0 coordinate"<<std::endl;
        return WHcoord();
    }
    size_t sumX( 0 ), sumY( 0 ), sumZ( 0 ), sumElements( 0 );
    for( int i=0 ; i<m_maskMatrix.size() ; ++i )
    {
        for( int j=0 ; j< m_maskMatrix[i].size() ; ++j )
        {
            for( int k=0 ; k<m_maskMatrix[i][j].size() ; ++k )
            {
                if (m_maskMatrix[i][j][k])
                {
                    sumX += i;
                    sumY += j;
                    sumZ += k;
                    ++sumElements;
                }
            }
        }
    }
    size_t meanX( sumX/sumElements );
    size_t meanY( sumY/sumElements );
    size_t meanZ( sumZ/sumElements );
    WHcoord meanCoord(meanX,meanY,meanZ);

    return meanCoord;
} // end "meanCoordFromMask()" -----------------------------------------------------------------








/*
void fileManager::writeLeafTract (const WHcoord& leaf, const compactTractChar& tractogram) const
{
    compactTract floatTract(tractogram);
    writeLeafTract( leaf, floatTract);
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------
void fileManager::writeLeafTract (const WHcoord& leaf, const compactTract& tractogram) const
{
    std::string tractFilename = getLeafTractFilename(leaf);
    writeTract( tractFilename, tractogram.m_tract );
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------
*/
void fileManager::writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const compactTract& tractogram) const
{
    std::string tractFilename = getLeafTractFilename( leaf, indexVector );
    writeTract( tractFilename, tractogram.m_tract );
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------
void fileManager::writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const compactTractChar& tractogram) const
{
    compactTract floatTract(tractogram);
    std::string tractFilename = getLeafTractFilename( leaf, indexVector );
    writeTract( tractFilename, floatTract.m_tract );
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------
void fileManager::writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTract& tractogram) const
{
    std::string tractFilename = getLeafTractFilename( leaf, indexVector, coordVector );
    writeTract( tractFilename, tractogram.m_tract );
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------
void fileManager::writeLeafTract (const size_t leaf, const std::vector< size_t >& indexVector, const std::vector< WHcoord >& coordVector, const compactTractChar& tractogram) const
{
    compactTract floatTract(tractogram);
    std::string tractFilename = getLeafTractFilename( leaf, indexVector, coordVector );
    writeTract( tractFilename, floatTract.m_tract );
    return;
}//end fileManager::writeLeafTract() -----------------------------------------------------

void fileManager::writeTract (const std::string& tractFilename, const compactTract& tractogram) const
{
    writeTract( tractFilename, tractogram.m_tract );
    return;
}//end fileManager::writeTract() -----------------------------------------------------
void fileManager::writeTract (const std::string& tractFilename, std::vector<float> tractogram ) const
{
    ValueType tractValueType;
    if(m_floatFlag)
    {
        tractValueType = VTFloat32;
    }
    else
    {
        tractValueType = VTUINT8;
        for( size_t i=0; i<tractogram.size(); ++i)
        {
            tractogram[i] = tractogram[i] * 255.;
        }
    }
    writeVector( tractFilename, tractValueType, tractogram, m_zipFlag );
    return;
}//end fileManager::writeTract() -----------------------------------------------------




void fileManager::writeDistBlock (const std::pair<unsigned int,unsigned int> blockID, const std::vector<std::vector<float> >& dBlock) const {
    writeDistBlock(blockID.first, blockID.second, dBlock);
    return;
}// end fileManager::writeDistBlock() -----------------------------------------------------------------
void fileManager::writeDistBlock (const unsigned int blockID1, const unsigned int blockID2, const std::vector<std::vector<float> >& dBlock) const {

    if( dBlock.empty()){
        std::cerr<< "ERROR @ fileManager::writeDistBlock(): block is empty, no file was written"<<std::endl;
        return;
    }
    std::string blockFilename = getBlockFilename( blockID1, blockID2 );
    writeMatrix( blockFilename, VTFloat32, dBlock, m_zipFlag );
    return;
}// end fileManager::writeDistBlock() -----------------------------------------------------------------


void fileManager::writeMask( const std::string& maskFilename, const std::vector<std::vector<std::vector<float> > >& maskMatrix) const
{
    writeImage( maskFilename, VTBit, maskMatrix, m_zipFlag );
    return;
}// end fileManager::writeMask() -----------------------------------------------------------------
void fileManager::writeMask( const std::string& maskFilename, const std::vector<std::vector<std::vector<bool> > >& maskMatrix) const
{
    const size_t dimx = maskMatrix.size();
    const size_t dimy = maskMatrix[0].size();
    const size_t dimz = maskMatrix[0][0].size();
    std::vector<std::vector<std::vector<float> > > maskMatrixFloat( dimx, std::vector<std::vector<float> >( dimy, std::vector<float>( dimz, 0 ) ) );

    for (int i=0 ; i<dimx ; ++i)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            for (int k=0 ; k<dimz ; ++k)
            {
                maskMatrixFloat[i][j][k] = maskMatrix[i][j][k];
            }
        }
    }
    writeMask( maskFilename, maskMatrixFloat );
    return;
}// end fileManager::writeMask() -----------------------------------------------------------------


// fileManager::loadMaskImage(): loads the full tractogram mask in memory
void fileManager::loadMaskImage( std::string maskFilename) {

    std::vector<std::vector<std::vector<float> > > maskMatrix;

    ValueType maskValueType = readImage( maskFilename, &maskMatrix );
    loadHeader(maskFilename);

    const size_t dimx = maskMatrix.size();
    const size_t dimy = maskMatrix[0].size();
    const size_t dimz = maskMatrix[0][0].size();

    m_maskMatrix.clear();

    {
        std::vector<bool> zvect(dimz,false);
        std::vector<std::vector<bool> >yzmatrix(dimy,zvect);
        m_maskMatrix.resize(dimx,yzmatrix);
    }

    size_t maskSum=0;
    for (int i=0 ; i<dimx ; ++i)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            for (int k=0 ; k<dimz ; ++k)
            {
                m_maskMatrix[i][j][k] = maskMatrix[i][j][k];
                if(m_maskMatrix[i][j][k])
                    ++maskSum;
            }
        }
    }

    m_flipVector.clear();
    m_flipVector.resize(maskSum,0);
    for (size_t i=0; i<maskSum; ++i)
    {
        m_flipVector[i]=i;
    }
    size_t indexer(0);
    for (int k=0 ; k<dimz ; ++k)
    {
        for (int j=0 ; j<dimy ; ++j)
        {
            size_t valuesPerRow(0);
            for (int i=0 ; i<dimx ; ++i)
            {
                if(m_maskMatrix[i][j][k])
                {
                    ++valuesPerRow;
                }
            }
            std::reverse(m_flipVector.begin()+indexer,m_flipVector.begin()+indexer+valuesPerRow);
            indexer+=valuesPerRow;
        }
    }
    return;
}// end fileManager::loadMaskImage() -----------------------------------------------------------------

void fileManager::readDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> >* dBlockPointer ) const
{
    readDistBlock(blockID.first,blockID.second,dBlockPointer);
    return;
}// end fileManager::readDistBlock() -----------------------------------------------------------------


void fileManager::readDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> >* dBlockPointer ) const
{
    std::string blockFilename = getBlockFilename(blockID1, blockID2);

    ValueType blockRepnType = readMatrix(blockFilename,dBlockPointer);

    if (blockRepnType==VTError)
    {
        throw std::runtime_error(  "ERROR @ fileManager::readDistBlock(): there was an error when reading the distance block"  );
    }

    if (blockRepnType!=VTFloat32)
    {
        throw std::runtime_error(  "ERROR @ fileManager::readDistBlock(): there was an error when reading the distance block"  );
    }

    return;
}// end fileManager::readDistBlock() -----------------------------------------------------------------

void fileManager::readTract(const std::string& tractFilename, std::vector<float>* tractPointer) const
{
    std::vector<float>& tract = *tractPointer;
    std::vector<float> vector;
    ValueType tractValueType = readVector(tractFilename, &vector);

    if (tractValueType==VTError)
    {
        throw std::runtime_error(  "ERROR @ fileManager::readTract(): there was an error when reading the tractogram"  );
    }

    tract.clear();
    tract.reserve( vector.size() );
    if(tractValueType == VTUINT8 )
    {
        for( size_t i=0; i < vector.size(); ++i )
        {
            tract.push_back( vector[i] / 255. );
        }
    }
    else if(tractValueType == VTFloat32 )
    {
        tract = vector;
    }
    else
    {
        throw std::runtime_error( "readTract(): tract representation type not recognized (neither VFloat nor VUByte)" );
    }
    return;
}// end fileManager::readTract() -----------------------------------------------------------------

void fileManager::readTract(const std::string& tractFilename, std::vector<unsigned char>* tractPointer) const
{
    std::vector<unsigned char>& tract = *tractPointer;
    std::vector<float> vector;
    ValueType tractValueType = readVector(tractFilename, &vector);

    if (tractValueType==VTError)
    {
        throw std::runtime_error(  "ERROR @ fileManager::readTract(): there was an error when reading the tractogram"  );
    }

    tract.clear();
    tract.reserve( vector.size() );
    if(tractValueType == VTUINT8 )
    {
        for( size_t i=0; i < vector.size(); ++i )
        {
            tract.push_back( (unsigned char)vector[i] );
        }
    }
    else if(tractValueType == VTFloat32 )
    {
        for( size_t i=0; i < vector.size(); ++i )
        {
            tract.push_back( (unsigned char)(vector[i] * 255.) );
        }
    }
    else
    {
        throw std::runtime_error( "readTract(): tract representation type not recognized (neither VFloat nor VUByte)" );
    }
    return;
}// end fileManager::readTract() -----------------------------------------------------------------
