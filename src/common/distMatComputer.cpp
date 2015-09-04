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
    size_t blockProgress(0), progressInt(0), elapsedTime(0), expectedRemain(0), totalBlocks(0), totalSubBlocks( 0 );
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
            if( row == column )
            {
                totalSubBlocks += ( m_subBlocksPerBlock * ( m_subBlocksPerBlock + 1) / 2 );
            }
            else
            {
                totalSubBlocks += m_subBlocksPerBlock * m_subBlocksPerBlock;
            }
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


            std::pair< dist_t, dist_t > blockMinMax = computeDistBlock( row, column );

            if( blockMinMax.first < minValue )
            {
                minValue = blockMinMax.first;
            }
            if( blockMinMax.second > maxValue )
            {
                maxValue = blockMinMax.second;
            }
            //update progress (number of subblocks computed)
            if( row == column )
            {
                blockProgress += ( m_subBlocksPerBlock * ( m_subBlocksPerBlock + 1) / 2 );
            }
            else
            {
                blockProgress += m_subBlocksPerBlock * m_subBlocksPerBlock;
            }


            if( m_verbose )
            {
                currentTime = (time(NULL));

                // printout progress
                progress = ( blockProgress * 100.) / totalSubBlocks;
                progressInt = progress;

                currentTime = ( time( NULL ) );
                elapsedTime = ( difftime( currentTime, distmatStartTime ) );
                expectedRemain = ( elapsedTime * ( ( 100.-progress ) / progress ));

                std::stringstream message;
                message << "\rCompleted block " << row << "-" << column << ". " << progressInt <<  "% completed (" << blockProgress << " of " << totalSubBlocks << " sub-blocks). ";
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
    // variables for progress printout
    size_t blockProgress(0), progressInt(0), totalSubBlocks( 0 );
    float progress(0);

    if( row == column )
    {
        totalSubBlocks = m_subBlocksPerBlock * ( m_subBlocksPerBlock +1 ) / 2;
    }
    else
    {
        totalSubBlocks = m_subBlocksPerBlock * m_subBlocksPerBlock;
    }


    std::string blockMessage("\r\t Block " + string_utils::toString( row ) + "-" + string_utils::toString( column ) + ". ");
    if( m_verbose )
    {
        std::cout << blockMessage << std::flush;
    }
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


        if( m_verbose )
        {
            if ( m_subBlocksPerBlock > 1 )
            {
                progress = blockProgress * 100. / totalSubBlocks;
                progressInt = progress;
                std::cout << blockMessage << progressInt << "% complete. ";
                std::cout << "Loading sub-row " << subRow << " tracts..." << std::flush;
            }
            else
            {
                std::cout << " Loading row tracts..." << std::flush;
            }
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


            if( m_verbose )
            {
                if ( m_subBlocksPerBlock > 1 )
                {
                    progress = blockProgress * 100. / totalSubBlocks;
                    progressInt = progress;
                    std::cout << blockMessage << progressInt << "% complete. Sub-block " << subRow << "-" << subColumn <<". ";
                }
                std::cout << " Loading column tracts..." << std::flush;
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
                if ( m_subBlocksPerBlock > 1 )
                {
                    progress = blockProgress * 100. / totalSubBlocks;
                    progressInt = progress;
                    std::cout << blockMessage << progressInt << "% complete. Sub-block " << subRow << "-" << subColumn <<". ";
                }
                std::cout << " Computing distances..." << std::flush;
            }


            computeDistances( subRowNorms, rowTracts, subColumnNorms, colTracts, &distBlockValues, firstSubRowSeedPos, firstSubColumnSeedPos );
            if( m_verbose )
            {
                std::cout << "Done. " << std::flush;


            }
            ++blockProgress;

            // cleanup
            delete [] colTracts;
        }
        // cleanup
        delete [] rowTracts;
    }


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

    if( m_verbose )
    {
        std::cout << blockMessage <<" 100% completed. Writing block to file..." << std::flush;
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
    if( m_verbose )
    {
        std::cout << "Done." << std::flush;
    }

    if( m_veryVerbose )
    {
        std::cout << blockMessage <<" 100% completed. min: " << minValue << ". Max: " << maxValue << std::endl;
    }

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


