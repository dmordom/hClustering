#include "vistaManager.h"
#include <errno.h>

#define NEWBOOST true

// PUBLIC MEMBERS

// vistaManager::getTractFilename(): returns tract filename
void vistaManager::getTractFilename(const WHcoord tractCoord, std::string &filename) const {
    filename = (m_tractFolder + "/connect_" + tractCoord.getNameString() + ".v");
    return;
}// end vistaManager::getTractFilename() -----------------------------------------------------------------

// vistaManager::getTractFilename(): returns tract filename
void vistaManager::getTractFilename(const size_t tractNode, std::string &filename) const {
    filename = str(boost::format("%s/compact_%06d.v") % m_tractFolder % tractNode);
    return;
}// end vistaManager::getTractFilename() -----------------------------------------------------------------
// vistaManager::getFullTractFilename(): returns tract filename
void vistaManager::getFullTractFilename(const WHcoord tractCoord, std::string &filename) const {
    filename = (m_tractFolder + "/fulltract_" + tractCoord.getNameString() + ".v");
    return;
}// end vistaManager::getFullTractFilename() -----------------------------------------------------------------

void vistaManager::getFullTractFilename(const size_t tractNode, std::string &filename) const {
    filename = str(boost::format("%s/fulltract_%06d.v") % m_tractFolder % tractNode);
    return;
}// end vistaManager::getFullTractFilename() -----------------------------------------------------------------

void vistaManager::getClusterMaskFilename(const size_t node, std::string* filename) const
{
    *filename = str(boost::format("%s/cluster_%06d.v") % m_tractFolder % node);
    return;
}// end vistaManager::getClusterMaskFilename() -----------------------------------------------------------------

void vistaManager::getClusterMaskFilenameZip(const size_t node, std::string* filename) const
{
    getClusterMaskFilename(node, filename);
    *filename += ".gz";
    return;
}// end vistaManager::getClusterMaskFilename() -----------------------------------------------------------------

size_t vistaManager::readClusterMaskNumber (const std::string filename) const
{
#if NEWBOOST
    std::string fileNameStem( boost::filesystem::path(filename).stem().string() );
#else
    std::string fileNameStem( boost::filesystem::path(filename).stem() );
#endif
    std::string clusterNumberString( fileNameStem.begin()+8, fileNameStem.end() ); // "cluster_" has 8 charachters
    size_t clusterNumber( boost::lexical_cast<size_t>( clusterNumberString ));
    return clusterNumber;
}// end vistaManager::readClusterMaskNumber() -----------------------------------------------------------------

size_t vistaManager::readFullTractNumber (const std::string filename) const
{
#if NEWBOOST
    std::string fileNameStem( boost::filesystem::path(filename).stem().string() );
#else
    std::string fileNameStem( boost::filesystem::path(filename).stem() );
#endif
    std::string clusterNumberString( fileNameStem.begin()+10, fileNameStem.end() ); // "fulltract_" has 10 charachters
    size_t clusterNumber( boost::lexical_cast<size_t>( clusterNumberString ));
    return clusterNumber;
}// end vistaManager::readClusterMaskFilename() -----------------------------------------------------------------


// "vistaManager::readLeafTract()": extracts compact tractogram data from a vista file and returns it in compact Tract form
void vistaManager::readLeafTract (const WHcoord tractCoord, compactTract &tractogram) const {
    std::string tractFilename;
    getTractFilename(tractCoord, tractFilename);
    getTract(tractFilename,tractogram.m_tract);
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}
void vistaManager::readLeafTract (const WHcoord tractCoord, compactTractChar &tractogram) const {
    std::string tractFilename;
    getTractFilename(tractCoord, tractFilename);
    getTract(tractFilename,tractogram.m_tract);
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end vistaManager::readLeafTract() -----------------------------------------------------

// "vistaManager::readNodeTract()": extracts compact tractogram data from a vista file and returns it in compact Tract form
void vistaManager::readNodeTract (const size_t tractNode, compactTract &tractogram) const
{
    std::string tractFilename;
    getTractFilename(tractNode, tractFilename);
    getTract(tractFilename,tractogram.m_tract);
    tractogram.m_inLogUnits=m_logFlag;
    tractogram.m_thresholded=m_thresFlag;
    return;
}//end vistaManager::readNodeTract() -----------------------------------------------------


void vistaManager::readStringTract (const std::string &tractFilename, compactTract &tractogram) const
{
getTract(tractFilename,tractogram.m_tract);
tractogram.m_inLogUnits=m_logFlag;
tractogram.m_thresholded=m_thresFlag;
return;
}//end vistaManager::readNodeTract() -----------------------------------------------------


void vistaManager::writeLeafTract (const WHcoord leaf, const compactTract &tractogram) const {
    std::string tractFilename;
    getTractFilename(leaf,tractFilename);
    const int columns(tractogram.m_tract.size());

    VImage vtract;

    if(m_floatFlag) {
        vtract = VCreateImage(1,1,columns,VFloatRepn);
        for (size_t i=0; i<columns; ++i) {
            VPixel(vtract,0,0,i,VFloat) = tractogram.m_tract[i];
        }
    } else {
        vtract = VCreateImage(1,1,columns,VUByteRepn);
        for (size_t i=0; i<columns; ++i) {
            VPixel(vtract,0,0,i,VUByte) = tractogram.m_tract[i]*255;
        }
    }


    // write tractogram to file
    if ( ! WriteVImage( (char *)tractFilename.c_str() , vtract ) )
        VError( "write_Vtract(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( vtract );

    return;

}
void vistaManager::writeLeafTract (const WHcoord leaf, const compactTractChar &tractogram) const {
    std::string tractFilename;
    getTractFilename(leaf,tractFilename);
    const int columns(tractogram.m_tract.size());

    VImage vtract;

    vtract = VCreateImage(1,1,columns,VUByteRepn);
    for (size_t i=0; i<columns; ++i) {
        VPixel(vtract,0,0,i,VUByte) = tractogram.m_tract[i];
    }

    // write tractogram to file
    if ( ! WriteVImage( (char *)tractFilename.c_str() , vtract ) )
        VError( "write_Vtract(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( vtract );

    return;

}//end vistaManager::writeLeafTract() -----------------------------------------------------


// "vistaManager::writeNodeTract()": write tractogram data into vista format file
void vistaManager::writeNodeTract (const size_t tractNode, const compactTract &tractogram) const {

    std::string tractFilename;
    getTractFilename(tractNode,tractFilename);

    writeTract(tractFilename, tractogram);

    return;

}//end vistaManager::writeNodeTract() -----------------------------------------------------

void vistaManager::writeTract (const std::string tractFilename, const compactTract &tractogram) const
{
        const int columns(tractogram.m_tract.size());

    VImage vtract;

    if(m_floatFlag) {
        vtract = VCreateImage(1,1,columns,VFloatRepn);
        for (size_t i=0; i<columns; ++i) {
            VPixel(vtract,0,0,i,VFloat) = tractogram.m_tract[i];
        }
    } else {
        vtract = VCreateImage(1,1,columns,VUByteRepn);
        for (size_t i=0; i<columns; ++i) {
            VPixel(vtract,0,0,i,VUByte) = tractogram.m_tract[i]*255;
        }
    }


    // write tractogram to file
    if ( ! WriteVImage( (char *)tractFilename.c_str() , vtract ) )
        VError( "write_Vtract(): Failed to open output tractogram file " );

    if (m_zipFlag) { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + tractFilename);
        int success(system(sysCommand.c_str()));
    }

    // clean up
    VDestroyImage( vtract );

    return;

}//end vistaManager::writeTract() -----------------------------------------------------


// "deleteTractFile()": delete tractogram file
int vistaManager::deleteTractFile ( const size_t tractNode ) const {

    std::string tractFilename;
    getTractFilename(tractNode,tractFilename);

    // Variables for calling system commands
    std::string delString("rm -f ");
    std::string sysCommand(delString +tractFilename);

    // remove file
    return system(sysCommand.c_str());

}// end vistaManager::deleteTractFile(" -----------------------------------------------------------------




// vistaManager::loadMask(): loads the full tractogram mask in memory
void vistaManager::loadMask( std::string maskFilename) {

    std::string maskFilenameZipped( maskFilename );

#if NEWBOOST
    std::string extension( boost::filesystem::path(maskFilename).extension().string() );
#else
    std::string extension( boost::filesystem::path(maskFilename).extension() );
#endif


    if ( extension == ".gz" )
    {
        maskFilename.resize(maskFilename.size()-3);

        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + maskFilenameZipped + " > " + maskFilename );
        int success(system(sysCommand.c_str()));
    }
    else if ( extension == ".v" )
    {
    }
    else
    {
        std::cerr << "File \"" << maskFilenameZipped << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return;
    }


    // read mask image
    VImage mask_image;


    if ( ! ReadImage( (char *)maskFilename.c_str(), mask_image ) )
        VError( "store_Vtract_Image(): Failed to read mask image" );

    VRepnKind maskRepnType(VPixelRepn(mask_image));
    if (maskRepnType!=VBitRepn && maskRepnType!=VFloatRepn && maskRepnType!=VUByteRepn)
        VError( "storeFullVtract(): Mask representation type not recognized (neither VFloat nor VUByte nor VBit)" );

    const size_t Bands   = VImageNBands( mask_image );
    const size_t Rows    = VImageNRows( mask_image );
    const size_t Columns = VImageNColumns( mask_image );


    m_maskMatrix.clear();

    {
        std::vector<bool> column(Columns,false);
        std::vector<std::vector<bool> >slice(Rows,column);
        m_maskMatrix.resize(Bands,slice);
    }
    size_t maskSum=0;
    for (int i=0 ; i<Bands ; ++i) {
        for (int j=0 ; j<Rows ; ++j) {
            for (int k=0 ; k<Columns ; ++k) {

                if (maskRepnType==VBitRepn)
                    m_maskMatrix[i][j][k] = VPixel(mask_image,i,j,k,VBit);
                else if (maskRepnType==VFloatRepn)
                    m_maskMatrix[i][j][k] = VPixel(mask_image,i,j,k,VFloat);
                else
                    m_maskMatrix[i][j][k] = VPixel(mask_image,i,j,k,VUByte);

                if(m_maskMatrix[i][j][k])
                    ++maskSum;
            }
        }
    }

    // clean up
    VDestroyImage( mask_image );

    m_flipVector.clear();
    m_flipVector.resize(maskSum,0);
    for (size_t i=0; i<maskSum; ++i)
        m_flipVector[i]=i;
    size_t indexer(0);
    for (int i=0 ; i<Bands ; ++i) {
        for (int j=0 ; j<Rows ; ++j) {
            size_t valuesPerRow(0);
            for (int k=0 ; k<Columns ; ++k) {
                if(m_maskMatrix[i][j][k])
                    ++valuesPerRow;
            }
            std::reverse(m_flipVector.begin()+indexer,m_flipVector.begin()+indexer+valuesPerRow);
            indexer+=valuesPerRow;
        }
    }

    if ( extension == ".gz" )
    {
        // Variables for calling system commands
        std::string sysCommand( "rm -f " + maskFilename );
        int success(system(sysCommand.c_str()));
    }


}// end vistaManager::loadMask() -----------------------------------------------------------------


void vistaManager::writeMask( const std::string maskFilename, const std::vector<std::vector<std::vector<bool> > > &maskMatrix ) const
{
    if (maskMatrix.empty()) {
        std::cerr<< "ERROR @ vistaManager::writeMask(): Mask matrix is empty, mask has not been stored" <<std::endl;
        return;
    }


    const size_t Bands   = maskMatrix.size();
    const size_t Rows    = maskMatrix[0].size();
    const size_t Columns = maskMatrix[0][0].size();

    VImage maskImage;

    maskImage = VCreateImage(Bands,Rows,Columns,VBitRepn);

    VFillImage(maskImage,VAllBands,0);

    for( int i=0 ; i<Bands ; ++i )
    {
        for( int j=0 ; j<Rows ; ++j )
        {
            for( int k=0 ; k<Columns ; ++k )
            {
                if (maskMatrix[i][j][k])
                {
                    VPixel(maskImage,i,j,k,VBit) = 1;
                }
            }
        }
    }

    // write mask to file
    if ( ! WriteVImage( (char *)maskFilename.c_str() , maskImage ) )
        VError( "writeMask(): Failed to open output mask file " );

    // clean up
    VDestroyImage( maskImage );

    if (m_zipFlag) { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + maskFilename);
        int success(system(sysCommand.c_str()));
    }

    return;
}// end vistaManager::writeMask() -----------------------------------------------------------------

void vistaManager::readMatrix( std::string matrixFilename, std::vector<std::vector<float> > *matrixPointer ) const
{
    std::vector<std::vector<float> > &matrix = *matrixPointer;

    std::string matrixFilenameZipped( matrixFilename );

#if NEWBOOST
    std::string extension( boost::filesystem::path(matrixFilename).extension().string() );
#else
    std::string extension( boost::filesystem::path(matrixFilename).extension() );
#endif

    if ( extension == ".gz" )
    {
        matrixFilename.resize(matrixFilename.size()-3);

        // Variables for calling system commands
        std::string sysCommand( "gzip -dcf " + matrixFilenameZipped + " > " + matrixFilename );
        int success(system(sysCommand.c_str()));
    }
    else if ( extension == ".v" )
    {
    }
    else
    {
        std::cerr << "File \"" << matrixFilenameZipped << "\" has no recognized extension (\"" << extension << "\") stopping." << std::endl;
        return;
    }


    // read matrix image
    VImage matrixImage;


    if ( ! ReadImage( (char *)matrixFilename.c_str(), matrixImage ) )
        VError( "loadMatrix(): Failed to read matrix image" );

    VRepnKind maskRepnType(VPixelRepn(matrixImage));
    if (maskRepnType!=VBitRepn && maskRepnType!=VFloatRepn && maskRepnType!=VUByteRepn)
        VError( "loadMatrix(): matrix representation type not recognized (neither VFloat nor VUByte nor VBit)" );

    const size_t Bands   = VImageNBands( matrixImage );

    if (Bands != 1)
        VError( "loadMatrix(): matrix image has multiple bands" );


    const size_t Rows    = VImageNRows( matrixImage );
    const size_t Columns = VImageNColumns( matrixImage );


    matrix.clear();

    {
        std::vector<float> column(Columns,false);
        matrix.resize(Rows,column);
    }

    for (int i=0 ; i<Rows ; ++i) {
        for (int j=0 ; j<Columns ; ++j) {
            if (maskRepnType==VBitRepn)
                matrix[i][j] = VPixel(matrixImage,0,i,j,VBit);
            else if (maskRepnType==VFloatRepn)
                matrix[i][j] = VPixel(matrixImage,0,i,j,VFloat);
            else
                matrix[i][j] = VPixel(matrixImage,0,i,j,VUByte);
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
    return;

}// end vistaManager::readMatrix() -----------------------------------------------------------------


void vistaManager::writeMatrix( const std::string matrixFilename, const std::vector<std::vector<float> > &matrix ) const
{
    if (matrix.empty()) {
        std::cerr<< "ERROR @ vistaManager::writeMatrix(): matrix is empty, it has not been stored" <<std::endl;
        return;
    }


    const size_t Rows    = matrix.size();
    const size_t Columns = matrix[0].size();

    VImage matrixImage;

    matrixImage = VCreateImage(1,Rows,Columns,VFloatRepn);



    for( int i=0 ; i<Rows ; ++i )
    {
        for( int j=0 ; j<Columns ; ++j )
        {
            VPixel(matrixImage,0,i,j,VFloat) = matrix[i][j];
        }
    }

    // write matrix to file
    if ( ! WriteVImage( (char *)matrixFilename.c_str() , matrixImage ) )
        VError( "writeMatrix(): Failed to open output matrix file " );

    // clean up
    VDestroyImage( matrixImage );

    if (m_zipFlag) { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + matrixFilename);
        int success(system(sysCommand.c_str()));
    }

    return;
}// end vistaManager::writeMatrix() -----------------------------------------------------------------


void vistaManager::readFullTract( const std::string fullTractFilename, compactTract* tractPointer ) const
{
    compactTract &tractogram = *tractPointer;
    if (m_maskMatrix.empty()) {
        std::cerr<< "ERROR @ vistaManager::readFullTract(): Mask hast not been loaded, tractogram has not been read"<<std::endl;
        return;
    }

    VImage fullTractImage;

    if ( ! ReadImage( (char *)fullTractFilename.c_str(), fullTractImage ) )
        VError( "readFullTract(): Failed to read full tract image" );

    VRepnKind tractRepnType(VPixelRepn(fullTractImage));
    if ( tractRepnType!=VFloatRepn && tractRepnType!=VUByteRepn)
        VError( "storeFullVtract(): tract representation type not recognized (neither VFloat nor VUByte)" );

    const size_t Bands   = VImageNBands( fullTractImage );
    const size_t Rows    = VImageNRows( fullTractImage );
    const size_t Columns = VImageNColumns( fullTractImage );

    if ( Bands != m_maskMatrix.size() || Rows != m_maskMatrix[0].size() || Columns!= m_maskMatrix[0][0].size() )
        VError( "storeFullVtract(): mask and full tract dimensions do not match" );

    tractogram.m_inLogUnits = m_logFlag;
    tractogram.m_norm = 0;
    tractogram.m_normReady = false;
    tractogram.m_thresholded = m_thresFlag;
    tractogram.m_tract.empty();
    tractogram.m_tract.reserve( getTractSize() );

    for( int i = 0; i < Bands; ++i )
    {
        for( int j = 0; j < Rows; ++j )
        {
            for( int k = 0; k < Columns; ++k )
            {
                if ( m_maskMatrix[i][j][k] )
                {
                    if ( tractRepnType == VFloatRepn )
                    {
                        tractogram.m_tract.push_back( VPixel( fullTractImage, i, j, k, VFloat ) );
                    }
                    else
                    {
                        tractogram.m_tract.push_back( VPixel( fullTractImage, i, j, k, VUByte ) );
                    }
                }
            }
        }
    }

    // clean up
    VDestroyImage( fullTractImage );

    return;
}// end vistaManager::readFullTract() -----------------------------------------------------------------



// "storeFulltract()": write full image tractogram into vista format file and compress
void vistaManager::storeFullTract (const size_t tractNode, const compactTract &tractogram) const {
    std::string tractFilename;
    getFullTractFilename(tractNode,tractFilename);
    storeFullTract(tractFilename,tractogram);
    return;
}
void vistaManager::storeFullTract (const WHcoord tractCoord, const compactTract &tractogram) const {
    std::string tractFilename;
    getFullTractFilename(tractCoord,tractFilename);
    storeFullTract(tractFilename,tractogram);
    return;
}
void vistaManager::storeFullTract (const std::string &tractFilename, const compactTract &tractogram) const {

    if (m_maskMatrix.empty()) {
        std::cerr<< "ERROR @ vistaManager::storeFullTract(): Mask hast not been loaded, file has not been stored"<<std::endl;
        return;
    }

    const size_t Bands   = m_maskMatrix.size();
    const size_t Rows    = m_maskMatrix[0].size();
    const size_t Columns = m_maskMatrix[0][0].size();

    VImage fullTractImage;

    if(m_floatFlag)
        fullTractImage = VCreateImage(Bands,Rows,Columns,VFloatRepn);
    else
        fullTractImage = VCreateImage(Bands,Rows,Columns,VUByteRepn);

    VFillImage(fullTractImage,VAllBands,0);

    std::vector<float>::const_iterator tractIter(tractogram.m_tract.begin());


    for( int i=0 ; i<Bands ; ++i )
    {
        for( int j=0 ; j<Rows ; ++j )
        {
            for( int k=0 ; k<Columns ; ++k )
            {
                if (m_maskMatrix[i][j][k])
                {
                    if (tractIter == tractogram.m_tract.end())
                        VError( "store_Vtract_Image(): Mask and tractogram sizes do not match" );

                    if (m_floatFlag)
                    {
                        VPixel(fullTractImage,i,j,k,VFloat) = *tractIter++;
                    }
                    else
                    {
                        VPixel(fullTractImage,i,j,k,VUByte) = (*tractIter++)*255;
                    }
                }
            }
        }
    }


    if (tractIter != tractogram.m_tract.end())
        VError( "store_Vtract_Image(): Mask and tractogram sizes do not match" );

    // write tractogram to file
    if ( ! WriteVImage( (char *)tractFilename.c_str() , fullTractImage ) )
        VError( "store_Vtract_Image(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( fullTractImage );

    if (m_zipFlag) { // zip file
        // Variables for calling system commands
        std::string gzipString("gzip -f ");
        std::string sysCommand(gzipString + tractFilename);
        int success(system(sysCommand.c_str()));
    }

    return;

}// end vistaManager::storeFulltract() -----------------------------------------------------------------

void vistaManager::getBlockFilename(const unsigned int blockID1, const unsigned int blockID2, std::string &filename) const {
    filename = str(boost::format("%s/dist_block_%03d_%03d.v") % m_tractFolder % blockID1 % blockID2);
    return;
}// end vistaManager::getBlockFilename() -----------------------------------------------------------------

void vistaManager::readDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> > &dBlock) const {
    readDistBlock(blockID.first,blockID.second,dBlock);
    return;
}// end vistaManager::readDistBlock() -----------------------------------------------------------------


void vistaManager::readDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> > &dBlock) const {
    std::string blockFilename;
    getBlockFilename(blockID1, blockID2, blockFilename);

    VImage blockImage;
    dBlock.clear();

    if ( ! ReadImage( (char *)blockFilename.c_str(), blockImage ) )
        VError( "ERROR @ vistaManager::readDistBlock(): Failed to read distance block image" );

    VRepnKind blockRepnType(VPixelRepn(blockImage));

    if (blockRepnType!=VFloatRepn)
        VError("ERROR @ vistaManager::readDistBlock(): distance block image must be of type float");
    if (VImageNBands(blockImage)!=1)
        VError("ERROR @ vistaManager::readDistBlock(): distance block image must have 1 band only");
    size_t nColumns = VImageNColumns(blockImage);
    size_t nRows = VImageNRows(blockImage);

    {
        std::vector<float> tempRow(nColumns,0);
        dBlock.resize(nRows,tempRow);
    }

    for (size_t i=0; i<nRows; ++i) {
        for (size_t j=0; j<nColumns; ++j) {
            dBlock[i][j] = VPixel(blockImage,0,i,j,VFloat);
        }
    }

    // clean up
    VDestroyImage(blockImage);

    return;
}// end vistaManager::readDistBlock() -----------------------------------------------------------------

void vistaManager::writeDistBlock (const std::pair<unsigned int,unsigned int> blockID, std::vector<std::vector<float> > &dBlock) {
    writeDistBlock(blockID.first, blockID.second, dBlock);
    return;
}// end vistaManager::writeDistBlock() -----------------------------------------------------------------

void vistaManager::writeDistBlock (const unsigned int blockID1, const unsigned int blockID2, std::vector<std::vector<float> > &dBlock) {

    if( dBlock.empty()){
        std::cerr<< "ERROR @ vistaManager::writeDistBlock(): block is empty, no file was written"<<std::endl;
        return;
    }

    std::string blockFilename;
    getBlockFilename(blockID1, blockID2, blockFilename);

    size_t rows(dBlock.size());
    size_t columns(dBlock[0].size());

    VImage blockImage = VCreateImage(1,rows,columns,VFloatRepn);
    VFillImage(blockImage,VAllBands,0);

    for (size_t i=0; i<rows; ++i) {
        for (size_t j=0; j<dBlock[i].size(); ++j) {
            VPixel(blockImage,0,i,j,VFloat) = dBlock[i][j];
        }
    }

    // write block to file
    if ( ! WriteVImage( (char *)blockFilename.c_str() , blockImage ) )
        VError( "writeDistBlock(): Failed to open output tractogram file " );

    // clean up
    VDestroyImage( blockImage );

    return;

}// end vistaManager::writeDistBlock() -----------------------------------------------------------------


void vistaManager::flipXtract(compactTract &tractogram) {

    if (m_flipVector.empty()) {
        std::cerr<< "ERROR @ vistaManager::flipXtract(): Mask hast not been loaded, tractogram has not been flipped"<<std::endl;
        return;
    }
    if (m_flipVector.size()!=tractogram.size()) {
        std::cerr<< "ERROR @ vistaManager::flipXtract(): flip vector and tractogram sizes dont match"<<std::endl;
        return;
    }

    compactTract flippedTract(tractogram);
    flippedTract.m_tract.assign(m_flipVector.size(),0);
    for (size_t i=0; i<m_flipVector.size(); ++i)
        flippedTract.m_tract[i]=tractogram.m_tract[m_flipVector[i]];

    tractogram=flippedTract;

    return;

}// end vistaManager::flipXtract() -----------------------------------------------------------------

WHcoord vistaManager::meanCoordFromMask() const
{
    if (m_maskMatrix.empty()) {
        std::cerr<< "ERROR @ vistaManager::storeFullTract(): Mask hast not been loaded, returning 0 coordinate"<<std::endl;
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
                    sumX += k;
                    sumY += j;
                    sumZ += i;
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



// PROTECTED MEMBERS


// vistaManager::getTract(): extracts compact tractogram data from a vista file and returns it in compact Tract form
void vistaManager::getTract(const std::string &tractFilename, std::vector<float> &tractogram) const {

    VImage tractImage;
    tractogram.clear();

    if ( ! ReadImage( (char *)tractFilename.c_str(), tractImage ) )
        VError( "ERROR @ vistaManager::getTract(): Failed to read mask image" );

    VRepnKind tractRepnType(VPixelRepn(tractImage));

    if (tractRepnType!=VFloatRepn && tractRepnType!=VUByteRepn)
        VError("ERROR @ vistaManager::getTract(): tractogram image must be of type char or float");

    if (VImageNBands(tractImage)!=1 && VImageNRows(tractImage)!= 1)
        VError("ERROR @ vistaManager::getTract(): tractogram image must have 1 row and 1 band only");

    size_t nElements = VImageNColumns(tractImage);
    tractogram.reserve(nElements); // Reserve memory in vector to increase speed

    // Copy tractogram data
    if (tractRepnType==VUByteRepn) {
        //it is a seed tractogram in char format
        for (size_t i=0; i<nElements; ++i)
            tractogram.push_back(VPixel(tractImage,0,0,i,VUByte)/255.);
    } else {
        //it is a node tractogram in float format
        for (int i=0; i<nElements; ++i)
            tractogram.push_back(VPixel(tractImage,0,0,i,VFloat));
    }

    // clean up
    VDestroyImage(tractImage);

    return;

}
void vistaManager::getTract(const std::string &tractFilename, std::vector<unsigned char> &tractogram) const {

    VImage tractImage;
    tractogram.clear();

    if ( ! ReadImage( (char *)tractFilename.c_str(), tractImage ) )
        VError( "ERROR @ vistaManager::getTract(): Failed to read mask image" );

    VRepnKind tractRepnType(VPixelRepn(tractImage));

    if (tractRepnType!=VUByteRepn)
        VError("ERROR @ vistaManager::getTract(): tractogram image must be of type char");

    if (VImageNBands(tractImage)!=1 && VImageNRows(tractImage)!= 1)
        VError("ERROR @ vistaManager::getTract(): tractogram image must have 1 row and 1 band only");

    size_t nElements = VImageNColumns(tractImage);
    tractogram.reserve(nElements); // Reserve memory in vector to increase speed

    // Copy tractogram data
    for (size_t i=0; i<nElements; ++i)
    {
        tractogram.push_back(VPixel(tractImage,0,0,i,VUByte));
    }

    // clean up
    VDestroyImage(tractImage);

    return;

}// end vistaManager::getTract() -----------------------------------------------------------------


// "WriteVImage()": write Vista image
int vistaManager::WriteVImage( const VString Name, VImage &Image )  const {

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
            std::string errorstring( "WriteVImage(): Failed to open output vista file. Code: "+ boost::lexical_cast<std::string>(errno) + " " + strerror(errno) );
            VError ( errorstring.c_str() );
        }
        /* write file */
        success = VWriteFile (file, list);
        fclose (file);
    }
    ioMutex.unlock();

   if (!success)
   { VError ("WriteVImage(): Failed to write output file '%s'", Name);
   }

   /* remove images */
   for (VFirstAttr (list, &pos); VAttrExists (&pos); VNextAttr (&pos))
      if (VGetAttrRepn (&pos) == VImageRepn)
         VSetAttrValue (&pos, NULL, VImageRepn, NULL);

   /* clean-up*/
   VDestroyAttrList (list);

    return 1;

} // end "WriteVImage()" -----------------------------------------------------------------



// "ReadImage()": read Vista image file
int vistaManager::ReadImage( const VString Name, VImage &Image ) const {


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
      VError ("ReadImage(): Failed to read input file '%s'", Name);
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
       VError ("ReadImage(): Input file '%s' contains multiple images",  Name);
   }
   if (! (Image) && !VAttrExists (&pos))
   {
       VError ("ReadImage(): Input file '%s' does not contain an image", Name);
   }

   // clean-up
   VDestroyAttrList (list);

    return 1;

} // end "ReadImage()" -----------------------------------------------------------------

