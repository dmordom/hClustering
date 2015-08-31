#include "distMatComputer.h"


distMatComputer::distMatComputer( const std::string& roiFilename, const float thresholdRatio, const bool verbose, const bool noLog )
{
    m_verbose = verbose;
    m_veryVerbose = false;
    m_ready2go = false;
    m_blockSize = 0;
    m_blocksPerRow = 0;
    m_subBlockSize = 0;
    m_subBlocksPerBlock = 0;
    m_zipFlag = true;


    fileManagerFactory fileMF;
    m_niftiMode = fileMF.isNifti();
    RoiLoader roiLoader( m_niftiMode, true );
    m_roiLoaded = roiLoader.readRoi( roiFilename, &m_datasetGrid, &m_datasetSize, &m_numStreamlines, &m_coordinates, &m_trackids );
    std::cout << "Roi loaded, " << m_coordinates.size() << " seed voxels" << std::endl;


    if( noLog )
    {
        m_logFactor = 0;
    }
    else
    {
        if( m_numStreamlines == 0 )
        {
            std::cerr << "WARNING @ distMatComputer::distMatComputer(): provided a number of stramlines per voxel of 0, interpreting it as requesting no logarithmic normalization of tracts (input tracts must also be in natural units)" << std::endl;
            m_logFactor = 0;
        }
        else
        {
            m_logFactor = log10( m_numStreamlines );
        }
    }

    if( thresholdRatio <= 0 || thresholdRatio >= 1 )
    {
        if( thresholdRatio != 0 )
        {
            std::cerr << "WARNING @ distMatComputer::distMatComputer(): threshold ratio provided (" << thresholdRatio << ") is out of bounds [0,1), using a value of 0.0 (no thresholding)" << std::endl;
        }
        m_tractThreshold = 0;
    }
    else if ( m_logFactor == 0 )
    {
        m_tractThreshold = thresholdRatio; // if using natural units the normalized tract threshold is the same as the threshold ratio
    }
    else
    {
        m_tractThreshold = log10( m_numStreamlines * thresholdRatio ) / m_logFactor; // if using log units the threshold must be calculated
    }
    if( verbose )
    {
        std::cout << "Final normalized threshold: " << m_tractThreshold << std::endl;
    }

    m_startingBlock = std::make_pair< size_t, size_t >( 0, 0 );
    m_finishBlock = m_startingBlock;

}

void distMatComputer::setBlockSize( float memory, size_t blockSize )
{
    if ( !m_roiLoaded )
    {
        std::cerr << "ERROR: Roi was not loaded." << std::endl;
        return;
    }
    if ( m_inputFolder.empty() )
    {
        std::cerr << "ERROR: Input folder was not set." << std::endl;
        return;
    }
    if ( m_outputFolder.empty() )
    {
        std::cerr << "ERROR: Output folder was not set." << std::endl;
        return;
    }
    // computed this way to avoid overflow
    size_t byteMemory =  1024 * memory;
    byteMemory =  1024 * byteMemory;
    byteMemory =  1024 * byteMemory;

    size_t maxDistBlockSize(sqrt((byteMemory)/(sizeof(dist_t)*2.))); // max block is half the available memory
    maxDistBlockSize = ( ( size_t ) ( ( maxDistBlockSize ) / MIN_BLOCK_SIZE ) ) * MIN_BLOCK_SIZE; // adapt max samples to be a multiple of the sample unit

    if( blockSize > 0 && blockSize < MIN_BLOCK_SIZE )
    {
        if( m_verbose )
        {
            std::cerr<<"WARNING: indicated block size ( " << blockSize << " ) is smaller than minimum block size, setting it to the minimum of: " << MIN_BLOCK_SIZE <<" elements. "<<std::endl;
        }
        blockSize = MIN_BLOCK_SIZE;
    }

    if( blockSize == 0 || blockSize > maxDistBlockSize )
    {
        if( m_verbose )
        {
            std::cerr<<"WARNING: indicated block size ( " << blockSize << " ) is 0 or bigger than available memory, setting block size to maximum size: " << maxDistBlockSize <<" elements. "<<std::endl;
            blockSize = maxDistBlockSize;
        }
    }

    if( blockSize > m_coordinates.size() )
    {
        if( m_verbose )
        {
            std::cout<<"block size is bigger than seed set. Matrix will have a single block of size " << m_coordinates.size() << "x" << m_coordinates.size() << std::endl;
        }
        m_blockSize = m_coordinates.size();
        m_blocksPerRow = 1;
    }
    else
    {
        m_blockSize = blockSize;
        m_blocksPerRow = ( m_coordinates.size( ) / m_blockSize );
        if ( m_coordinates.size() % m_blockSize )
        {
            ++m_blocksPerRow;
        }
        if( m_verbose )
        {
            std::cout<< m_blocksPerRow <<"x"<< m_blocksPerRow <<" blocks of size "<< m_blockSize <<"x"<< m_blockSize <<std::endl;
        }
    }

    size_t distBlockBytes( m_blockSize * m_blockSize * sizeof( dist_t ) );

    // remaining available memory
    size_t remainingMemoryBytes = byteMemory - distBlockBytes;

    // get the size of a single tract
    size_t tractBytes;

    compactTractChar testTract;
    fileManagerFactory testFMF( m_inputFolder );
    fileManager& testFM( testFMF.getFM() );
    testFM.readLeafTract(0, m_trackids, m_coordinates, &testTract);
    m_trackSize = testTract.size();
    tractBytes = m_trackSize*sizeof( testTract.tract()[0] );

    if( m_verbose )
    {
        std::cout <<"Tractogram size: "<< m_trackSize <<" elements ("<< tractBytes / (1024.*1024.) << " MBytes)"<< std::endl;
    }
    size_t maxSubBlockSize( remainingMemoryBytes / ( 2 * tractBytes ) );

    if ( maxSubBlockSize < MIN_SUB_BLOCK_SIZE )
    {
        maxSubBlockSize = MIN_SUB_BLOCK_SIZE;
        if( m_verbose )
        {
            std::cout <<"Memory restrictions are too strict, not enough memory for minimum sub_block, nevertheless setting to predefined minimum sub-block size: "<< maxSubBlockSize << std::endl;
        }
    }
    if ( m_blockSize > maxSubBlockSize )
    {
        m_subBlocksPerBlock = 1;
        m_subBlockSize = m_blockSize;
        while( m_subBlockSize > maxSubBlockSize )
        {
            ++m_subBlocksPerBlock;
            if ( ( m_blockSize % m_subBlocksPerBlock ) == 0 )
            {
                m_subBlockSize = m_blockSize / m_subBlocksPerBlock;
            }
        }
    }
    else
    {
        m_subBlockSize = m_blockSize;
        m_subBlocksPerBlock = 1;
    }

    size_t tractBlockBytes( tractBytes * m_subBlockSize );
    size_t memoryUsageBytes( distBlockBytes + ( 2 * tractBlockBytes ) + sizeof(*this) + ( sizeof( double ) * m_coordinates.size() ) );

    if( m_verbose )
    {
        std::cout <<"Using "<< m_subBlocksPerBlock <<"x"<< m_subBlocksPerBlock <<" tractogram sub-blocks of "<< m_subBlockSize <<" tracts for each distance block."<< std::endl;
        std::cout << "Expected memory usage: " << ( memoryUsageBytes / (1024 * 1024 ) / 1024.) << " GBytes (" << memoryUsageBytes / (1024 * 1024) << " MBytes). [ distBlock: " << distBlockBytes/(1024 * 1024) << " MB. tract subBlocks: 2 x " << tractBlockBytes/(1024 * 1024) << " MB. ]" << std::endl;
    }

    if( m_finishBlock.first == 0 && m_finishBlock.second == 0 )
    {
        m_finishBlock.first = m_blocksPerRow - 1;
        m_finishBlock.second = m_blocksPerRow - 1;
    }
    else if( m_finishBlock.first > 0 || m_finishBlock.second > 0 )
    {
        setFinishBlock( m_finishBlock.first, m_finishBlock.second );
    }

    if( m_startingBlock.first > 0 || m_startingBlock.second > 0 )
    {
        setStartingBlock( m_startingBlock.first, m_startingBlock.second );
    }

    m_ready2go = true;
    return;
}

void distMatComputer::setStartingBlock( size_t start_row, size_t start_column )
{
    if( start_column < start_row )
    {
        size_t newColTemp( start_row );
        start_row = start_column;
        start_column = newColTemp;
    }

    if( m_ready2go )
    {
        if( start_row >= m_blocksPerRow )
        {
            std::cerr << "WARNING @ distMatComputer::setStartingBlock: starting block indices are higher than number of matrix blocks, setting to 0;" << std::endl;
            m_startingBlock.first = 0;
            m_startingBlock.second = 0;
        }
        else
        {
            if( start_row > m_finishBlock.first || ( start_row == m_finishBlock.first && start_column >= m_finishBlock.second ) )
            {
                std::cerr << "WARNING @ distMatComputer::setStartingBlock: starting block indices are higher than preset finishing block indices, setting to 0;" << std::endl;
                m_startingBlock.first = 0;
                m_startingBlock.second = 0;
            }
            m_startingBlock.first = start_row;
            m_startingBlock.second = start_column;
        }
    }
    else
    {
        m_startingBlock.first = start_row;
        m_startingBlock.second = start_column;
    }
}

void distMatComputer::setFinishBlock( size_t finish_row, size_t finish_column )
{
    if( finish_column < finish_row )
    {
        size_t newColTemp( finish_row );
        finish_row = finish_column;
        finish_column = newColTemp;
    }

    if( m_ready2go )
    {
        if( finish_row >= m_blocksPerRow )
        {
            std::cerr << "WARNING @ distMatComputer::setFinishBlock: finish block indices are higher than number of matrix blocks, setting to final block indices;" << std::endl;
            m_finishBlock.first = m_blocksPerRow - 1;
            m_finishBlock.second = m_blocksPerRow - 1;
        }
        else
        {
            if( finish_row < m_startingBlock.first || ( finish_row == m_startingBlock.first && finish_column <= m_startingBlock.second ) )
            {
                std::cerr << "WARNING @ distMatComputer::setFinishBlock: finish block indices are higher than preset starting block indices, setting to max;" << std::endl;
                m_finishBlock.first = m_blocksPerRow - 1;
                m_finishBlock.second = m_blocksPerRow - 1;
            }
            m_finishBlock.first = finish_row;
            m_finishBlock.second = finish_column;
        }
    }
    else
    {
        m_finishBlock.first = finish_row;
        m_finishBlock.second = finish_column;
    }
}



void distMatComputer::doDistBlocks()
{
    //variables to keep tract of min and max values
    dist_t maxValue(-1), minValue(2);

    // variables for progress printout
    time_t distmatStartTime(time(NULL)), currentTime(time(NULL));
    size_t blockProgress(0), progressInt(0), elapsedTime(0), expectedRemain(0), totalBlocks(0);
    float progress(0);

    // write index file
    writeIndex();

    // compute tractogram norms
    computeNorms();

    // just to test the number of blocks that will be computed
    for (size_t row = m_startingBlock.first ; row <= m_finishBlock.first ; ++row)
    {
        for (size_t column = row ; column < m_blocksPerRow ; ++column)
        {
            if( row == m_startingBlock.first && column < m_startingBlock.second )
            {
                continue;
            }
            if( row == m_finishBlock.first && column > m_finishBlock.second )
            {
                continue;
            }
            ++totalBlocks;
        }
    }

    // loop through rows of the distance matrix (each element will be a separate distance block)
    for (size_t row = m_startingBlock.first ; row <= m_finishBlock.first ; ++row)
    {
        // loop through columns of the distance matrix (each element will be a separate distance block)
        for (size_t column = row ; column < m_blocksPerRow ; ++column)
        {
            if( row == m_startingBlock.first && column < m_startingBlock.second )
            {
                continue;
            }

            if( row == m_finishBlock.first && column > m_finishBlock.second )
            {
                 continue;
            }

            if( m_verbose )
            {
                std::cout << "Computing Block " << row << "-" << column << ". " << std::flush;
            }

            std::pair< dist_t, dist_t > blockMinMax = computeDistBlock( row, column );

            if( blockMinMax.first < minValue )
            {
                minValue = blockMinMax.first;
            }
            if( blockMinMax.second > maxValue )
            {
                maxValue = blockMinMax.second;
            }
            //update progress

            ++blockProgress;


            if( m_verbose )
            {
                currentTime = (time(NULL));

                // printout progress
                progress = ( blockProgress * 100.) / totalBlocks;
                progressInt = progress;

                currentTime = ( time( NULL ) );
                elapsedTime = ( difftime( currentTime, distmatStartTime ) );
                expectedRemain = ( elapsedTime * ( ( 100.-progress ) / progress ));

                std::stringstream message;
                message << "\rCompleted block " << row << "-" << column << ". " << progressInt <<  "% completed (" << blockProgress << " of " << totalBlocks << "). ";
                message << "Elapsed: " << elapsedTime/3600 <<"h "<<  (elapsedTime%3600)/60 <<"' "<< ((elapsedTime%3600)%60) <<"\". ";
                message << "Remaining: "<< expectedRemain/3600 <<"h "<<  (expectedRemain%3600)/60 <<"' "<< ((expectedRemain%3600)%60) <<"\"           ";
                std::cout << message.str() <<std::endl;

            }
        }
    }

    currentTime = ( time( NULL ) );
    elapsedTime = ( difftime( currentTime, distmatStartTime ) );
    std::cout << "100% of blocks completed (" << totalBlocks << " of matrix total " << (m_blocksPerRow*(m_blocksPerRow+1))/2. << "). ";
    std::cout << "Elapsed time: " << elapsedTime/3600 <<"h "<<  (elapsedTime%3600)/60 <<"' "<< ((elapsedTime%3600)%60) <<"\". " << std::endl;
    std::cout<<"Total MAX value: "<< maxValue <<". Total min value: "<< minValue <<std::endl;

}// end "doDistBlocks()" -----------------------------------------------------------------


void distMatComputer::computeNorms()
{
    // loop  through all the seed voxels and compute tractogram norms
    std::cout << "Precomputing tractogram norms" << std::endl;
    time_t loopStart( time( NULL ) ), lastTime( time( NULL ) );
    m_leafNorms.assign( m_coordinates.size(), 0 );
    size_t progCount( 0 );


    fileManagerFactory fileSingleMF(m_inputFolder);
    fileManager& fileSingle(fileSingleMF.getFM());
    fileSingle.readAsUnThres();
    fileSingle.readAsLog();

#pragma omp parallel for schedule( static )
    for( size_t i = 0; i < m_coordinates.size(); ++i )
    {
        compactTractChar thisTract;

        fileSingle.readLeafTract( i, m_trackids, m_coordinates, &thisTract );
        thisTract.threshold( m_tractThreshold );
        m_leafNorms[i] = thisTract.computeNorm();

#pragma omp atomic
        ++progCount;

#pragma omp single nowait // only one thread executes output
        if( m_verbose )
        {
            time_t currentTime( time( NULL ) );
            if( currentTime - lastTime > 1 )
            {
                lastTime = currentTime;
                size_t currentCount( progCount );
                float progress( currentCount * 100. / m_coordinates.size() );
                size_t intProgress( progress );
                size_t elapsedTime( difftime( currentTime, loopStart ) );
                std::stringstream message;
                message << "\r" << intProgress << " % of norms computed (" << currentCount << " tracts). ";
                if( progress > 0 )
                {
                    size_t expectedRemain( elapsedTime * ( ( 100. - progress ) / progress ) );
                    message << "Expected remaining time: ";
                    message << expectedRemain / 3600 << "h ";
                    message << ( expectedRemain % 3600 ) / 60 << "' ";
                    message << ( expectedRemain % 3600 ) % 60 << "\". ";
                }
                message << "Elapsed time: ";
                message << elapsedTime / 3600 << "h " << ( elapsedTime % 3600 ) / 60 << "' ";
                message << ( elapsedTime % 3600 ) % 60 << "\". ";
                std::cout << message.str() <<std::flush;
            }
        } // end verbose
    } // end  for


    int timeTaken = difftime( time( NULL ), loopStart );
    if( m_verbose )
    {
        std::cout << "\r" << std::flush << "100 % of norms computed. Time taken: " << timeTaken / 3600 << "h " << ( timeTaken
                        % 3600 ) / 60 << "' " << ( ( timeTaken % 3600 ) % 60 ) << "\"    " << std::endl;
    }

    if( m_verbose )
    {
        for( size_t i = 0; i < m_leafNorms.size(); ++i )
        {
            if( m_leafNorms[i] == 0)
            {
                std::cerr << "track " << i << " has norm 0" <<std::endl;
            }
        }
    }
} // end distMatComputer::computeNorms() -------------------------------------------------------------------------------------


std::pair< dist_t, dist_t > distMatComputer::computeDistBlock( size_t row, size_t column )
{
    //variables to keep tract of min and max values
    dist_t maxValue(-1), minValue(2);

    // get dimensions of the block and ID vectors for the row and column seeds
    size_t firstRowSeedID( row * m_blockSize ), postlastRowSeedID( ( row + 1 ) * m_blockSize );
    if ( postlastRowSeedID > m_coordinates.size() )
    {
        postlastRowSeedID = m_coordinates.size();
    }
    size_t thisBlockRowSize( postlastRowSeedID - firstRowSeedID );
    std::vector<size_t> blockRowIDs;
    blockRowIDs.reserve( thisBlockRowSize );
    for( size_t i = firstRowSeedID; i < postlastRowSeedID; ++i )
    {
        blockRowIDs.push_back( i );
    }

    size_t firstColumnSeedID( column * m_blockSize ), postlastColumnSeedID( ( column + 1 ) * m_blockSize );
    if ( postlastColumnSeedID > m_coordinates.size() )
    {
        postlastColumnSeedID = m_coordinates.size();
    }
    size_t thisBlockColumnSize( postlastColumnSeedID - firstColumnSeedID );
    std::vector<size_t> blockColumnIDs;
    blockColumnIDs.reserve( thisBlockColumnSize );
    for( size_t i = firstColumnSeedID; i < postlastColumnSeedID; ++i )
    {
        blockColumnIDs.push_back( i );
    }

    // reserve memory and initialize the distance block
    std::vector< std::vector< dist_t > > distBlockValues;
    {
        std::vector< dist_t > blockRowTemp( thisBlockColumnSize, 0 );
        std::vector< std::vector< dist_t > > blockTemp( thisBlockRowSize, blockRowTemp );
        distBlockValues.swap( blockTemp );
    }

    // loop over sub block rows
    for( size_t subRow = 0; subRow < m_subBlocksPerBlock; ++subRow )
    {
        size_t firstSubRowSeedPos( subRow * m_subBlockSize ), postlastSubRowSeedPos( ( subRow + 1 ) * m_subBlockSize );
        if ( postlastSubRowSeedPos > thisBlockRowSize )
        {
            postlastSubRowSeedPos = thisBlockRowSize;
        }
        size_t subBlockRowSize( postlastSubRowSeedPos - firstSubRowSeedPos );
        std::vector< double > subRowNorms( m_leafNorms.begin()+firstRowSeedID+firstSubRowSeedPos,  m_leafNorms.begin()+firstRowSeedID+postlastSubRowSeedPos );


        if( m_subBlocksPerBlock == 1 )
        {
            if( m_verbose )
            {
                std::cout << " Loading row tracts..." << std::flush;
            }
        }
        else if( m_veryVerbose )
        {
            std::cout << std::endl << "\t Loading sub-row " << subRow << " tracts..." << std::flush;
        }

        //load row tracts
        unsigned char* rowTracts;
        loadTractSet( blockRowIDs[firstSubRowSeedPos], blockRowIDs[postlastSubRowSeedPos-1], rowTracts, false );
        if( m_verbose )
        {
            std::cout << "Done. " << std::flush;

        }

        // loop over sub block columns
        for( size_t subColumn = 0; subColumn < m_subBlocksPerBlock; ++subColumn )
        {

            size_t firstSubColumnSeedPos( subColumn * m_subBlockSize ), postlastSubColumnSeedPos( ( subColumn + 1 ) * m_subBlockSize );
            if ( postlastSubColumnSeedPos > thisBlockColumnSize )
            {
                postlastSubColumnSeedPos = thisBlockColumnSize;
            }
            size_t subBlockColumnSize( postlastSubColumnSeedPos - firstSubColumnSeedPos );
            std::vector< double > subColumnNorms( m_leafNorms.begin()+firstColumnSeedID+firstSubColumnSeedPos,  m_leafNorms.begin()+firstColumnSeedID+postlastSubColumnSeedPos );


            if( m_subBlocksPerBlock == 1 )
            {
                if( m_verbose )
                {
                    std::cout << "Loading column tracts..." << std::flush;
                }
            }
            else if( m_veryVerbose )
            {
                std::cout << std::endl << "\t Loading sub-column " << subColumn << " tracts..." << std::flush;
            }

            unsigned char* colTracts;
            if( row == column )
            {
                // if its a block in the diagonal do not compute sub-blocks under the diagonal
                if( subRow > subColumn )
                {
                    continue;
                }
                else if( subRow == subColumn )
                {
                    // if its a sub-block in the diagonal sets are equal
                    transposeSet( subBlockColumnSize, rowTracts, colTracts );
                }
            }
            else
            {
                // load tractograms in transposed position
                loadTractSet( blockColumnIDs[firstSubColumnSeedPos], blockColumnIDs[postlastSubColumnSeedPos-1], colTracts, true );
            }
            if( m_verbose )
            {
                std::cout << "Done. " << std::flush;

            }
            if( m_subBlocksPerBlock == 1 )
            {
                if( m_verbose )
                {
                    std::cout << "Computing distances..." << std::flush;
                }
            }
            else if( m_veryVerbose )
            {
                std::cout << "Computing distaces..." << std::flush;
            }

            computeDistances( subRowNorms, rowTracts, subColumnNorms, colTracts, &distBlockValues, firstSubRowSeedPos, firstSubColumnSeedPos );
            if( m_verbose )
            {
                std::cout << "Done. " << std::flush;


            }
            /////////////////////////// progress report

            // cleanup
            delete [] colTracts;
        }
        // cleanup
        delete [] rowTracts;
    }

    /////////////////////////// progress report

    // get min and max values
    for( size_t i = 0; i < distBlockValues.size(); ++i )
    {
        size_t jInit( 0 );
        if( row == column )
        {
            jInit = i + 1;
        }

        for( size_t j = jInit; j < distBlockValues[i].size(); ++j )
        {
            dist_t value(distBlockValues[i][j]);
            if( minValue > value )
            {
                minValue = value;
            }
            if( maxValue < value )
            {
                maxValue = value;
            }
        }
    }

    // write to file
    fileManagerFactory blockFMF( m_outputFolder );
    fileManager& blockFM(blockFMF.getFM() );
    blockFM.writeInFloat();
    if( m_zipFlag )
    {
        blockFM.storeZipped();
    }
    else
    {
        blockFM.storeUnzipped();
    }
    blockFM.writeDistBlock( row, column, distBlockValues );

    return std::make_pair< dist_t, dist_t >( minValue, maxValue );
}// end "computeDistBlock()" -----------------------------------------------------------------

void distMatComputer::writeIndex()
{
    std::vector< std::pair< size_t, size_t > > roiBlockIndex;
    roiBlockIndex.reserve( m_coordinates.size() );

    for (size_t row=0 ; row < m_blocksPerRow ; ++row)
    {
        size_t firstRowSeedID( row * m_blockSize ), postlastRowSeedID( ( row + 1 ) * m_blockSize );
        if ( postlastRowSeedID > m_coordinates.size() )
        {
            postlastRowSeedID = m_coordinates.size();
        }
        size_t blockSize( postlastRowSeedID - firstRowSeedID );

        for (size_t i=0; i < blockSize; ++i)
        {
            roiBlockIndex.push_back(std::make_pair( row, i ) );
        }
    }

    std::string indexFilename = m_outputFolder + "/" + MATRIX_INDEX_FILENAME;
    std::ofstream outFileStream( indexFilename.c_str() );
    if( !outFileStream )
    {
        std::cerr << "ERROR: unable to open output index file: \"" << indexFilename << "\"" << std::endl;
        throw std::runtime_error("ERROR: unable to open output index file");
    }

    if( m_verbose )
    {
        std::cout<< "Writing distance matrix index file in \""<< indexFilename <<"\""<< std::endl;
    }

    outFileStream << "#distindex" << std::endl;

    for( size_t i = 0; i < m_coordinates.size(); ++i )
    {
        outFileStream << m_coordinates[i];
        outFileStream << boost::format( " b %03d i %04d" ) % roiBlockIndex[i].first  % roiBlockIndex[i].second;
        outFileStream << std::endl;
    }

    outFileStream << "#enddistindex" << std::endl;

}// end "writeIndex()" -----------------------------------------------------------------

void distMatComputer::computeDistances(
                       std::vector< double >& rowNorms, const unsigned char* rowTractSet,
                       std::vector< double >& columnNorms, const unsigned char* columnTractSet,
                       std::vector< std::vector< dist_t > >* distBlockPointer,
                       size_t blockRowOffset, size_t blockColumnOffset )
{

    std::vector< std::vector< dist_t > >& distBlockValues( *distBlockPointer );
    const double normalizer( 255. * 255. );
    size_t rowBlockSize( rowNorms.size() ), columnBlockSize( columnNorms.size() );

    #pragma omp parallel for // loop through tractograms (use parallel threads)
    for(size_t i=0 ; i<rowBlockSize ; ++i)
    {
        // if a sqrsum is 0 then distance is 1
        if (rowNorms[i]!=0.)
        {
            double* dotprodArray = new double[ columnBlockSize ];

            //initialize
            for(size_t init=0; init < columnBlockSize; ++init)
            {
                    dotprodArray[init]=0;
            }

            size_t rowTractArrayOffset( i * m_trackSize );
            double value1( 0 );
            size_t k(0), j(0);
            const unsigned char* p_value2;
            double* p_result;

            for(k=0 ; k<m_trackSize ; ++k)
            {

                p_result = dotprodArray;

                value1 = rowTractSet[rowTractArrayOffset+k]/normalizer;
                if (value1)
                {
                    p_value2 = columnTractSet + ( k * columnBlockSize );

                    for (j=0 ; j<columnBlockSize ;++j)
                    {

                        *(p_result++) += *(p_value2++) * value1;
                    }
                }
            }
            size_t jInit( 0 );

            for (j=0 ; j<columnBlockSize ;++j)
            {

                double vprod(0);

                if (columnNorms[j]!=0.)
                {
                    vprod = dotprodArray[j]/(rowNorms[i] * columnNorms[j]);
                }

                //insert value in distance matrix
                distBlockValues[blockRowOffset+i][blockColumnOffset+j] = 1 - (vprod);
            }


            delete [] dotprodArray;
        }
        else
        {
            for (int j=0 ; j<columnBlockSize ;++j)
            {
                distBlockValues[blockRowOffset+i][blockColumnOffset+j] = 1;
            }
        }
    }

}// end "computeDistances()" -----------------------------------------------------------------

void distMatComputer::loadTractSet( size_t firstID, size_t lastID, unsigned char* tractSet, bool transposed )
{
    size_t setOffset(0);
    size_t setElements( m_trackSize * ( lastID + 1 - firstID ) );
    tractSet = new unsigned char [ setElements ];
    for( size_t i = 0; i <  setElements; ++i )
    {
        tractSet[i]=0;
    }

    for(size_t i = firstID; i <= lastID; ++i )
    {
        fileManagerFactory tractFMF( m_inputFolder );
        fileManager& tractFM( tractFMF.getFM() );
        tractFM.readAsLog();
        tractFM.readAsUnThres();
        compactTractChar tract;
        tractFM.readLeafTract( i, m_trackids, m_coordinates, &tract );
        tract.threshold( m_tractThreshold );
        std::vector< unsigned char > tractValues( tract.tract() );
        if( transposed )
        {
            for( size_t tractPos = 0; tractPos < m_trackSize; ++tractPos )
            {
                tractSet[ ( tractPos * m_trackSize ) + setOffset ] = tractValues[ tractPos ];
            }
        }
        else
        {
            for( size_t tractPos = 0; tractPos < m_trackSize; ++tractPos )
            {
                tractSet[ tractPos + ( setOffset * m_trackSize) ] = tractValues[ tractPos ];
            }
        }
        ++setOffset;
    }
    return;
}

void distMatComputer::transposeSet( size_t setSize, unsigned char* originalSet, unsigned char* transposedSet )
{
    transposedSet = new unsigned char [ m_trackSize * setSize ];
    for( size_t i = 0; i <  m_trackSize * setSize; ++i )
    {
        transposedSet[i] = 0;
    }

    for(size_t setOffset = 0; setOffset < setSize; ++setOffset )
    {
        for( size_t tractPos = 0; tractPos < m_trackSize; ++tractPos )
        {
            transposedSet[ ( tractPos * m_trackSize ) + setOffset ] = originalSet[ tractPos + (setOffset * m_trackSize) ];
        }
    }
    return;
}



/*
void distMatComputer::doDistBlocks()
{


    //variables to keep tract of min and max values
    float max_value(-1), min_value(2);


    // variables for progress printout
    time_t dist_matrix_start(time(NULL)), dist_matrix_current_time(time(NULL));
    size_t block_prog(0), progress_int(0), block_prog_time(0), expected_remain(0);
    float progress(0);

    // loop through rows of the distance matrix (each element will be a separate distance block)
    for (size_t row = m_startingBlock.first ; row <= m_finishBlock.first ; ++row)
    {
        // get subset of row seeds
        size_t first_row_seed( row * m_blockSize ), postlast_row_seed( ( row + 1 ) * m_blockSize );
        if ( postlast_row_seed > m_coordinates.size() )
        {
            postlast_row_seed = m_coordinates.size();
        }

        if ( m_blockSize == m_subBlockSize )
        {
            time_t block_start(time(NULL));
            std::cout <<"Loading row block "<<row<<"..."<< std::flush;


            //allocate memory for tractogram block, square sum and create cropped seed vector
            unsigned char *row_tract_block = new unsigned char [m_trackSize*seed_set_row.size()];
            double *row_sqrsum = new double[seed_set_row.size()];

            //load data into block
            load_tract_block(tract_dir, seed_set_row, seed_set_row.size(), m_trackSize, row_tract_block);

            //compute square sum
            do_sqrsum(row_sqrsum, row_tract_block, seed_set_row.size(), m_trackSize);


            // loop through columns of the distance matrix (each element will be a separate distance block)
            for (int col=row ; col < m_blocksPerRow ; ++col)
            {

                std::cout <<std::endl <<"Computing block "<<row<<"-"<<col<<" in one go..."<< std::flush;

                // reserve space for the distance block matrix
                std::vector<std::vector<float> > dist_block;

                if (m_coordinates.size()%m_blockSize)
                { // not all blocks wil be square
                    if ( (row<m_blocksPerRow-1)&&(col<m_blocksPerRow-1) )
                    { //its one of the full square blocks
                        std::vector<float> test_vect(m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                    else if ( (row==m_blocksPerRow-1)&&(col==m_blocksPerRow-1) )
                    { // its the square block in the lower right corner
                        std::vector<float> test_vect(m_coordinates.size()%m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_coordinates.size()%m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                    else
                    { // its one of the rectangular blocks in the final column
                        std::vector<float> test_vect(m_coordinates.size()%m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                }
                else
                { // all blocks will be square
                    std::vector<float> test_vect(m_blockSize,0);
                    std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                    dist_block.swap(test_block);
                }


                if (row==col)
                { // block is on the diagonal

                    std::cout <<"transposing..."<< std::flush;


                    // get transposed version of the block
                    unsigned char *row_tract_block_tp = new unsigned char [m_trackSize*seed_set_row.size()];

                    for(size_t i=0; i<seed_set_row.size();++i)
                    {
                        for(size_t j=0; j<m_trackSize;++j)
                        {
                            row_tract_block_tp[(j*seed_set_row.size())+i]=row_tract_block[(i*m_trackSize)+j];
                        }
                    }


                    // Correlate row block with itself
                    compute_dist_block(dist_block, row_tract_block, row_tract_block_tp, row_sqrsum, row_sqrsum, 0, 0, seed_set_row.size(), seed_set_row.size(), m_trackSize, threads);

                    delete [] row_tract_block_tp;

                    for(size_t i=0; i<dist_block.size();++i)
                    {
                        for(size_t j=0; j<=i;++j)
                        {
                            dist_block[i][j]=0;
                        }
                    }



                }
                else
                {


                    block_start = (time(NULL));


                    // get subsets of seeds
                    size_t first_col_seed(col*m_blockSize), postlast_col_seed((col+1)*m_blockSize);
                    if (postlast_col_seed>m_coordinates.size())
                    {
                        postlast_col_seed = m_coordinates.size();
                    }

                    std::vector<coordinate> seed_set_col(&m_coordinates[first_col_seed], &m_coordinates[postlast_col_seed]);


                    //allocate memory for tractogram block, square sum and create cropped seed vector
                    unsigned char *col_tract_block = new unsigned char [m_trackSize*seed_set_col.size()];
                    double *col_sqrsum = new double[seed_set_col.size()];


                    //load data into block (dta will be loaded transposed, one tractogram per column, to speed up processing of distances later on)
                    load_tract_block_transposed(tract_dir, seed_set_col, seed_set_col.size(), m_trackSize, col_tract_block);


                    //compute square sum
                    do_sqrsum_transposed(col_sqrsum, col_tract_block, seed_set_col.size(), m_trackSize, threads);

                    // Correlate row block with with col block
                    compute_dist_block(dist_block, row_tract_block, col_tract_block, row_sqrsum, col_sqrsum, 0, 0, seed_set_row.size(), seed_set_col.size(), m_trackSize, threads);


                    // clear memory of local col block
                    delete [] col_tract_block;
                    delete [] col_sqrsum;

                }

                //update progress
                ++block_prog;

                //get minimum and maximum values for this block and update total maximum and minimum
                float block_max(-1), block_min(2);
                size_t col_begin(0);
                for (int i=0; i<dist_block.size(); ++i)
                {
                    if (row==col)
                    {
                        col_begin=i+1;
                    }
                    for (int j=col_begin; j<dist_block[i].size(); ++j)
                    {
                        if (dist_block[i][j]<0)
                        {
                            if (dist_block[i][j]<-0.0001)
                            {
                                std::cerr<< "WARNING: negative value: " << dist_block[i][j] << std::endl;
                            }
                            dist_block[i][j]=0;
                        }
                        if (dist_block[i][j]>block_max)
                        {
                            block_max=dist_block[i][j];
                        }
                        if (dist_block[i][j]<block_min)
                        {
                            block_min=dist_block[i][j];
                        }
                    }
                }
                if (block_max>max_value)
                {
                    max_value=block_max;
                }
                if(block_min<min_value)
                {
                    min_value=block_min;
                }

                dist_matrix_current_time = (time(NULL));
                time_t blockTime = ( difftime(dist_matrix_current_time,block_start) );

                // printout progress
                progress = (block_prog*100.)/(m_blocksPerRow*(m_blocksPerRow+1)/2);
                progress_int = progress;
                std::cout<<"Saving to disk..."<< std::flush;

                // write down block
                std::string test_name("/dist_block_000_000.v");
                char name[test_name.size()];
                sprintf(name, "/dist_block_%03d_%03d.v",row,col);
                std::string dist_block_filename(out_dir + name);
                write_dist_block (dist_block_filename, dist_block);

                // printout progress
                dist_matrix_current_time = (time(NULL));
                block_prog_time = ( difftime(dist_matrix_current_time,dist_matrix_start) );
                expected_remain = (block_prog_time*((100.-progress)/progress));

                std::cout << std::endl <<"Block "<<row<<"-"<<col<<" computed. in: ";
                std::cout<<(blockTime)/60 <<"' "<< ((blockTime)%60) <<"\" ";

                std::cout <<"Max: "<< block_max <<". min: "<< block_min <<". Total "<<progress_int<<"% done."<<std::flush;
                if (!((row==m_blocksPerRow-1)&&(col==m_blocksPerRow-1)))
                {
                    std::cout<<" Expected total time left: "<< expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\"                         ";
                }
                std::cout<<"                                                  "<<std::endl;




            }


            // clear memory of local row block
            delete [] row_tract_block;
            delete [] row_sqrsum;

        }
        else
        {


            // loop through columns of the distance matrix (each element will be a separate distance block)
            for (int col=row ; col < m_blocksPerRow ; ++col)
            {

#if 0
                //if( row <3 || (row == 3 && col<12) )
                if( row != col )
                {
                    continue;
                }
#endif




                std::cout <<"\r"<< std::flush <<"Computing block "<<row<<"-"<<col<<"..."<< std::flush;

                // reserve space for the distance block matrix
                std::vector<std::vector<float> > dist_block;

                if (m_coordinates.size()%m_blockSize)
                { // not all blocks wil be square
                    if ( (row<m_blocksPerRow-1)&&(col<m_blocksPerRow-1) ) { //its one of the full square blocks
                        std::vector<float> test_vect(m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                    else if
                            ( (row==m_blocksPerRow-1)&&(col==m_blocksPerRow-1) ) { // its the square block in the lower right corner
                        std::vector<float> test_vect(m_coordinates.size()%m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_coordinates.size()%m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                    else
                    { // its one of the rectangular blocks in the final column
                        std::vector<float> test_vect(m_coordinates.size()%m_blockSize,0);
                        std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                        dist_block.swap(test_block);
                    }
                }
                else
                { // all blocks will be square
                    std::vector<float> test_vect(m_blockSize,0);
                    std::vector<std::vector<float> > test_block(m_blockSize,test_vect);
                    dist_block.swap(test_block);
                }




                if (row==col)
                { // block is on the diagonal




                    std::vector<coordinate> empty_coord;
                    std::vector<std::vector<float> > empty_tracts;


                    fill_dist_block(tract_dir, seed_set_row, empty_coord, dist_block, m_subBlockSize, m_trackSize, threads, verbose, veryvb);

                    //update progress
                    ++block_prog;


                }
                else
                { // block is over the diagonal

                    // get subsets of seeds
                    size_t first_col_seed(col*m_blockSize), postlast_col_seed((col+1)*m_blockSize);
                    if (postlast_col_seed>m_coordinates.size())
                        postlast_col_seed = m_coordinates.size();

                    std::vector<coordinate> seed_set_col(&m_coordinates[first_col_seed], &m_coordinates[postlast_col_seed]);

                    std::vector< std::vector<float> > rand_tracts_col;
                    if (rand_mode)
                        rand_tracts_col.assign(&rand_tracts[first_col_seed],&rand_tracts[postlast_col_seed]);

                    //fill distance block matrix

                     fill_dist_block(tract_dir, seed_set_row, seed_set_col, dist_block, m_subBlockSize,m_trackSize, threads, verbose, veryvb);

                    //update progress
                    ++block_prog;
                    ++block_prog;

                }






                //get minimum and maximum values for this block and update total maximum and minimum
                float block_max(-1), block_min(2);
                size_t col_begin(0);
                for (int i=0; i<dist_block.size(); ++i) {
                    if (row==col)
                        col_begin=i+1;
                    for (int j=col_begin; j<dist_block[i].size(); ++j) {
                        if (dist_block[i][j]<0)
                        {
                            if (dist_block[i][j]<-0.0001)
                            {
                                std::cerr<< "WARNING: negative value: " << dist_block[i][j] << std::endl;
                            }
                            dist_block[i][j]=0;
                        }
                        if (dist_block[i][j]>block_max)
                            block_max=dist_block[i][j];
                        if (dist_block[i][j]<block_min)
                            block_min=dist_block[i][j];
                    }
                }
                if (block_max>max_value)
                    max_value=block_max;
                if(block_min<min_value)
                    min_value=block_min;

                // printout progress
                progress = (block_prog*100.)/(m_blocksPerRow*(m_blocksPerRow));
                progress_int = progress;
                std::cout <<"\r"<< std::flush <<"Block "<<row<<"-"<<col<<" computed. Saving to disk..."<< std::flush;

                // write down block
                std::string test_name("/dist_block_000_000.v");
                char name[test_name.size()];
                sprintf(name, "/dist_block_%03d_%03d.v",row,col);
                std::string dist_block_filename(out_dir + name);
                write_dist_block (dist_block_filename, dist_block);

                // printout progress
                dist_matrix_current_time = (time(NULL));
                block_prog_time = ( difftime(dist_matrix_current_time,dist_matrix_start) );
                expected_remain = (block_prog_time*((100.-progress)/progress));

                std::cout <<"\r"<< std::flush <<"Block "<<row<<"-"<<col<<" computed. Max: "<< block_max <<". min: "<< block_min <<". Total "<<progress_int<<"% done."<<std::flush;
                if (!((row==m_blocksPerRow-1)&&(col==m_blocksPerRow-1)))
                    std::cout<<" Expected total time left: "<< expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\"                         ";
                std::cout<<"                                                  "<<std::endl;
            } // end for cols
        } // end if else
    }// end for rows
    std::cout<<"Total MAX value: "<< max_value <<". Total min value: "<< min_value <<std::endl;

}// end "do_dist_blocks()" -----------------------------------------------------------------




// "fill_dist_block()": compute the distances from a distance block, computation will be divided into different sub-blocks due to memory issues
void distMatComputer::fill_dist_block(std::string &tract_dir,
                     std::vector<coordinate> &row_dist_seed_set,
                     std::vector<coordinate> &col_dist_seed_set,
                     std::vector<std::vector<float> > &dist_block,
                     const size_t &tract_block_size,
                     const size_t &tract_length,
                     const size_t &threads,
                     const bool &verbose,
                     const bool &veryvb) {


    std::cerr <<" 0 % of current block completed."<< std::flush;
    // progress printing variables
    time_t block_start(time(NULL)), block_current_time(time(NULL)), block_last_time(time(NULL));
    size_t block_prog(0), progress_int(0), block_prog_time(0), update_time(0), expected_remain(0);
    float progress(0);

    // get number of blocks to be used (can be lower if its the last row)
    if (row_dist_seed_set.empty())
        throw std::runtime_error ("ERROR [fill_dist_block()]: seed set is empty");

    size_t num_row_tblocks(row_dist_seed_set.size()/tract_block_size);
    if (row_dist_seed_set.size()%tract_block_size)
        ++num_row_tblocks;

    size_t num_col_tblocks(col_dist_seed_set.size()/tract_block_size);
    if (col_dist_seed_set.size()%tract_block_size)
        ++num_col_tblocks;
    if (col_dist_seed_set.empty())
        num_col_tblocks = 0;

    // looping through row subsets
    for (int row=0; row<num_row_tblocks ; ++row) {

        // get block size (can be lower if its the last row)
        size_t row_tblock_size(tract_block_size);
        if ((row==num_row_tblocks-1)&&(row_dist_seed_set.size()%tract_block_size))
            row_tblock_size = (row_dist_seed_set.size()%tract_block_size);
        size_t row_distb_offset (row*tract_block_size);

        //allocate memory for tractogram block, square sum and create cropped seed vector
        unsigned char *row_tract_block = new unsigned char [tract_length*row_tblock_size];
        double *row_sqrsum = new double[row_tblock_size];
        std::vector<coordinate> row_tblock_seed_set(&row_dist_seed_set[row_distb_offset],&row_dist_seed_set[row_distb_offset+row_tblock_size]);

        //load data into block
        load_tract_block(tract_dir, row_tblock_seed_set, row_tblock_size, tract_length, row_tract_block);

        //compute square sum
        do_sqrsum(row_sqrsum, row_tract_block, row_tblock_size, tract_length, threads);


        if (col_dist_seed_set.empty()) {   // distance block is part of the diagonal, only compute upper diagonal values

            // looping through column subsets
            for (int col=row; col<num_row_tblocks ; ++col) {

                if (col==row) { // sub-block is also on the diagonal

                    // get transposed version of the block
                    unsigned char *row_tract_block_tp = new unsigned char [tract_length*row_tblock_size];

                    for(size_t i=0; i<row_tblock_size;++i) {
                        for(size_t j=0; j<tract_length;++j) {
                            row_tract_block_tp[(j*row_tblock_size)+i]=row_tract_block[(i*tract_length)+j];
                        }
                    }

                    // Correlate row block with itself
                    compute_dist_block(dist_block, row_tract_block, row_tract_block_tp, row_sqrsum, row_sqrsum, row_distb_offset, row_distb_offset, row_tblock_size, row_tblock_size, tract_length, threads);

                    delete [] row_tract_block_tp;

                    // uodate progress
                    block_prog += 1;
                }


                else { // sub-block is not on the diagonal


                    // get block size (can be lower if its the last column)
                    size_t col_tblock_size(tract_block_size);
                    if ((col==num_row_tblocks-1)&&(row_dist_seed_set.size()%tract_block_size))
                        col_tblock_size = (row_dist_seed_set.size()%tract_block_size);
                    size_t col_distb_offset (col*tract_block_size);

                    //allocate memory for tractogram block, square sum and create cropped seed vector
                    unsigned char *col_tract_block = new unsigned char [tract_length*col_tblock_size];
                    double *col_sqrsum = new double[col_tblock_size];
                    std::vector<coordinate> col_tblock_seed_set(&row_dist_seed_set[col_distb_offset],&row_dist_seed_set[col_distb_offset+col_tblock_size]);

                    //load data into block (dta will be loaded transposed, one tractogram per column, to speed up processing of distances later on)
                    load_tract_block_transposed(tract_dir, col_tblock_seed_set, col_tblock_size, tract_length, col_tract_block);

                    //compute square sum
                    do_sqrsum_transposed(col_sqrsum, col_tract_block, col_tblock_size, tract_length, threads);

                    // Correlate row block with with col block
                    compute_dist_block(dist_block, row_tract_block, col_tract_block, row_sqrsum, col_sqrsum, row_distb_offset, col_distb_offset, row_tblock_size, col_tblock_size, tract_length, threads);

                    // uodate progress
                    block_prog += 2;

                    // clear memory of local col block
                    delete [] col_tract_block;
                    delete [] col_sqrsum;

                }

                // write progress
                progress = (block_prog*100.)/(num_row_tblocks*num_row_tblocks);
                progress_int = progress;
                block_last_time = block_current_time;
                block_current_time = (time(NULL));
                update_time = difftime(block_current_time,block_last_time);
                block_prog_time = ( difftime(block_current_time,block_start) );
                expected_remain = (block_prog_time*((100.-progress)/progress));
                std::cerr <<"\r"<< std::flush << progress_int <<" % of current block completed. Expected remaining time (for this block): ";
                std::cerr<<expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\" (update time: "<< update_time*2/60 <<"')    "<< std::flush;


            }

        } else { // block is above the diagonal, compute all distances


            for (int col=0; col<num_col_tblocks ; ++col) {


                // get block size (can be lower if its the last column)
                size_t col_tblock_size(tract_block_size);
                if ((col==num_col_tblocks-1)&&(col_dist_seed_set.size()%tract_block_size))
                    col_tblock_size = (col_dist_seed_set.size()%tract_block_size);
                size_t col_distb_offset (col*tract_block_size);

                //allocate memory for tractogram block, square sum and create cropped seed vector
                unsigned char *col_tract_block = new unsigned char [tract_length*col_tblock_size];
                double *col_sqrsum = new double[col_tblock_size];
                std::vector<coordinate> col_tblock_seed_set(&col_dist_seed_set[col_distb_offset],&col_dist_seed_set[col_distb_offset+col_tblock_size]);

                //load data into block
                load_tract_block_transposed(tract_dir, col_tblock_seed_set, col_tblock_size, tract_length, col_tract_block);

                //compute square sum
                do_sqrsum_transposed(col_sqrsum, col_tract_block, col_tblock_size, tract_length, threads);

                // Correlate row block with with col block
                compute_dist_block(dist_block, row_tract_block, col_tract_block, row_sqrsum, col_sqrsum, row_distb_offset, col_distb_offset, row_tblock_size, col_tblock_size, tract_length, threads);

                // uodate progress
                block_prog += 1;

                // clear memory of local col block
                delete [] col_tract_block;
                delete [] col_sqrsum;



                // write progress
                progress = (block_prog*100.)/(num_row_tblocks*num_col_tblocks);
                progress_int = progress;
                block_last_time = block_current_time;
                block_current_time = (time(NULL));
                update_time = difftime(block_current_time,block_last_time);
                block_prog_time = ( difftime(block_current_time,block_start) );
                expected_remain = (block_prog_time*((100.-progress)/progress));
                std::cerr <<"\r"<< std::flush << progress_int <<" % of current block completed. Expected remaining time (for this block): ";
                std::cerr<<expected_remain/3600 <<"h "<<  (expected_remain%3600)/60 <<"' "<< ((expected_remain%3600)%60) <<"\" (update time: "<< update_time*2/60 <<"')    "<< std::flush;



            }

        }

        // clear memory of local row block
            delete [] row_tract_block;
            delete [] row_sqrsum;


    }
    block_current_time = (time(NULL));
    block_prog_time = ( difftime(block_current_time,block_start) );
    std::cerr <<"\r"<< std::flush <<"100 % of block completed. Matrix size: "<<dist_block.size()<<"x"<< dist_block[dist_block.size()-1].size()<< std::flush;
    std::cerr<<". Time taken: "<<block_prog_time/3600 <<"h "<<  (block_prog_time%3600)/60 <<"' "<< ((block_prog_time%3600)%60) <<"\"         "<< std::flush;
    if (verbose)
        std::cout<<std::endl;

    // Print out to check
    if (veryvb) {
        std::cout <<"distance matrix:"<< std::endl;
        for (std::vector<std::vector<dist_t> >::const_iterator iter_row(dist_block.begin()); iter_row!=dist_block.end(); ++iter_row) {

            for (std::vector<dist_t>::const_iterator iter_col(iter_row->begin()); iter_col!=iter_row->end(); ++iter_col) {
                printf("%1.3f ",*iter_col);
            }
            std::cout<<std::endl;
        }
    }

    return;

}// end "fill_dist_block()" -----------------------------------------------------------------





// "compute_dist_block()":
void distMatComputer::compute_dist_block(std::vector<std::vector<float> > &dist_block,
                        unsigned char* row_tract_block,
                        unsigned char* col_tract_block,
                        double* row_sqrsum,
                        double* col_sqrsum,
                        const size_t  row_distb_offset,
                        const size_t  col_distb_offset,
                        const size_t  row_tblock_size,
                        const size_t  col_tblock_size,
                        const size_t &tract_length,
                        const size_t &threads){

    std::cout <<" Computing distances..."<< std::flush;

    const double normalizer(65025);


    omp_set_num_threads(threads);
    #pragma omp parallel for // loop through tractograms (use parallel threads)
    for(int i=0 ; i<row_tblock_size ; ++i) {

        if (row_sqrsum[i]!=0.){// if a sqrsum is 0 then distance is 1


            double* cov_array = new double[col_tblock_size];

            //initialize
            for(size_t init=0; init < col_tblock_size; ++init)
                    cov_array[init]=0;

            size_t tb1_offset(i*tract_length);

            double value1;
            size_t k(0), j(0);
            unsigned char* p_value2;
            double* p_result;

            for(k=0 ; k<tract_length ; ++k){

                p_result = cov_array;

                value1 = row_tract_block[tb1_offset+k]/normalizer;
                if (value1) {
                    p_value2 = col_tract_block + (k*col_tblock_size);

                    for (j=0 ; j<col_tblock_size ;++j) {

                        *(p_result++) += *(p_value2++) * value1;
                    }
                }
            }

            for (int j=0 ; j<col_tblock_size ;++j) {

                double vprod(0);

                if (col_sqrsum[j]!=0.)
                    vprod = cov_array[j]/sqrt( (row_sqrsum[i]) *(col_sqrsum[j]));

                //insert value in distance matrix
                dist_block[row_distb_offset+i][col_distb_offset+j] = 1 - (vprod); // 65025 = 255 * 255
            }

            delete [] cov_array;
        } else {
            for (int j=0 ; j<col_tblock_size ;++j)
                dist_block[row_distb_offset+i][col_distb_offset+j] = 1;
        }
    }
    std::cout <<"Done."<< std::flush;


}// end "compute_dist_block()" Two blocks -----------------------------------------------------------------

// fill_rand_block(): load distance block with random data (upper triangle)
void distMatComputer::fill_rand_block(std::vector<std::vector<float> > rand_tracts_row,
                     std::vector<std::vector<float> > rand_tracts_col,
                     std::vector<std::vector<float> > &dist_block,
                     const size_t rand_tract_length,
                     const size_t &threads,
                     const bool &verbose,
                     const bool &veryvb) {


    std::cerr <<" 0 % of current block completed."<< std::flush;
    // progress printing variables
    time_t block_start(time(NULL)), block_current_time(time(NULL)), block_last_time(time(NULL));
    size_t block_prog(0), progress_int(0), block_prog_time(0), update_time(0), expected_remain(0);
    float progress(0);


    if (rand_tracts_col.empty()) {
        omp_set_num_threads(threads);
        #pragma omp parallel for // loop through tractograms (use parallel threads)
        for (int i=0; i< rand_tracts_row.size(); ++i) {
            for(int j=i+1; j<rand_tracts_row.size(); j++) {
                dist_block[i][j] = tract_distance(rand_tracts_row[i],rand_tracts_row[j]);
            }
        }


    } else {

        omp_set_num_threads(threads);
        #pragma omp parallel for // loop through tractograms (use parallel threads)
        for (int i=0; i< rand_tracts_row.size(); ++i) {
            for(int j=0; j<rand_tracts_col.size(); j++) {
                dist_block[i][j] = tract_distance(rand_tracts_row[i],rand_tracts_col[j]);
            }
        }





    }



    block_current_time = (time(NULL));
    block_prog_time = ( difftime(block_current_time,block_start) );
    std::cerr <<"\r"<< std::flush <<"100 % of block completed. Matrix size: "<<dist_block.size()<<"x"<< dist_block[dist_block.size()-1].size()<< std::flush;
    std::cerr<<". Time taken: "<<block_prog_time/3600 <<"h "<<  (block_prog_time%3600)/60 <<"' "<< ((block_prog_time%3600)%60) <<"\"         "<< std::flush;
    if (verbose)
        std::cout<<std::endl;

    // Print out to check
    if (veryvb) {
        std::cout <<"distance matrix:"<< std::endl;
        for (std::vector<std::vector<dist_t> >::const_iterator iter_row(dist_block.begin()); iter_row!=dist_block.end(); ++iter_row) {

            for (std::vector<dist_t>::const_iterator iter_col(iter_row->begin()); iter_col!=iter_row->end(); ++iter_col) {
                printf("%1.3f ",*iter_col);
            }
            std::cout<<std::endl;
        }
    }

    return;
}// end "fill_rand_block()" -----------------------------------------------------------------

*/
