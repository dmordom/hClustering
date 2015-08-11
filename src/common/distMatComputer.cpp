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
    m_ready2go = true;
    return;
}



/*
// "do_dist_blocks()": compute and write down the distance blocks
void distMatComputer::doDistBlocks()
{


    //variables to keep tract of min and max values
    float max_value(-1), min_value(2);


    // variables for progress printout
    time_t dist_matrix_start(time(NULL)), dist_matrix_current_time(time(NULL));
    size_t block_prog(0), progress_int(0), block_prog_time(0), expected_remain(0);
    float progress(0);

    // loop through rows of the distance matrix (each element will be a separate distance block)
    for (int row=0 ; row < m_blocksPerRow ; ++row)
    {

#if 0
        if( row < 3)
        {
           continue;
        }
#endif

        // get subset of row seeds
        size_t first_row_seed(row*m_blockSize), postlast_row_seed((row+1)*m_blockSize);
        if (postlast_row_seed>m_coordinates.size())
        {
            postlast_row_seed = m_coordinates.size();
        }
        std::vector<WHcoord> seed_set_row(&m_coordinates[first_row_seed], &m_coordinates[postlast_row_seed]);


        if (m_blockSize==m_subBlockSize)
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

#if 0
                if((row == 3 && ( col>row && col<15)))
                {
                   continue;
                }
#endif


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






// "load_tract_block()": load single tractograms into a block
void distMatComputer::load_tract_block(std::string tract_dir,
                      std::vector<coordinate> &seed_set,
                      const size_t  tract_block_size,
                      const size_t &tract_length,
                      unsigned char* tract_block) {

    std::cout <<" Loading tract set..."<< std::flush;

    for(int i=0; i < tract_block_size; ++i) {

        coordinate current_coord(seed_set[i]);
        get_Vtract_th(get_Vname(tract_dir,current_coord),&tract_block[i*tract_length],tract_length, threshold);

    }
    std::cout <<"Done."<< std::flush;
    return;
}// end "load_tract_block()" -----------------------------------------------------------------

// "load_tract_block_transposed()": load single transposed tractograms into a block
void distMatComputer::load_tract_block_transposed(std::string tract_dir,
                                 std::vector<coordinate> seed_set,
                                 const size_t  tract_block_size,
                                 const size_t &tract_length,
                                 unsigned char* tract_block) {

    std::cout <<" Loading tp tract set..."<< std::flush;

    unsigned char * current_tract = new unsigned char [tract_length];


    for(int i=0; i < tract_block_size; ++i) {

        coordinate current_coord(seed_set[i]);
        get_Vtract_th(get_Vname(tract_dir,current_coord),current_tract,tract_length, threshold);
        for (int j=0; j < tract_length; ++j )
            tract_block[(j*tract_block_size)+i] = current_tract[j];
    }
    delete current_tract;

    std::cout <<"Done."<< std::flush;
    return;
}// end "load_tract_block_transposed()" -----------------------------------------------------------------


// "do_sqrsum()": precoumpute squared sums of tract blocks
void distMatComputer::do_sqrsum(double* sqrsum_array,
               unsigned char* tract_block,
               const size_t  tract_block_size,
               const size_t &tract_length,
               const size_t &threads) {

    std::cout <<" Computing square sum..."<< std::flush;

    //initialize
    for(int i=0; i < tract_block_size; ++i)
        sqrsum_array[i]=0;

    // compute
    omp_set_num_threads(threads);
    #pragma omp parallel for // loop through tractograms (use parallel threads)
    for(int i=0; i < tract_block_size; ++i) {
        size_t start(i*tract_length);
        size_t end(start+tract_length);
        double value;
        for(size_t j=start; j < end; ++j) {
            value = tract_block[j]/255.;
            sqrsum_array[i] += value*value;
        }
    }

    std::cout <<"Done."<< std::flush;
    return;
}// end "do_sqrsum()" -----------------------------------------------------------------


// "do_sqrsum_transposed()": precoumpute squared sums of transposed tract blocks
void distMatComputer::do_sqrsum_transposed(double* sqrsum_array,
                          unsigned char* tract_block,
                          const size_t  tract_block_size,
                          const size_t &tract_length,
                          const size_t &threads) {

    std::cout <<" Computing tp square sum..."<< std::flush;

    //initialize
    for(int i=0; i < tract_block_size; ++i)
        sqrsum_array[i]=0;

    // compute
   // omp_set_num_threads(threads);
   // #pragma omp parallel for // loop through tractograms (use parallel threads)
    for(int i=0; i < tract_length; ++i) {
        size_t offset(i*tract_block_size);
        double value;
        for(int j=0; j < tract_block_size; ++j) {
            value = tract_block[offset+j]/255.;
            sqrsum_array[j] += value * value;
        }
    }
    std::cout <<"Done."<< std::flush;
    return;
}// end "do_sqr_sum()" -----------------------------------------------------------------





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
