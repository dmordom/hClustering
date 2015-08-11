#ifndef DISTMATCOMPUTER_H
#define DISTMATCOMPUTER_H

// std library
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <algorithm>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>



// parallel execution
#include <omp.h>

// classes
#include "WHnode.h"
#include "WHcoord.h"
#include "fileManagerFactory.h"
#include "compactTract.h"
#include "roiLoader.h"


#define MIN_BLOCK_SIZE 500
#define MIN_SUB_BLOCK_SIZE 10

/**
 * This class computes a distance matrix from compact probabilistic tracts
 */
class distMatComputer
{
public:
    /**
     * Constructor
     * \param roiFilename file containing the list of seed voxels coordinates (from which the tractograms were computed)
     * \param thresholdRatio the relative number of streamlines that have to pass through a certain voxel in a tractogram to be allowed to contribute towards tract similarity, will determine threshold value (threshold is applied to filter out tractography noise)
     * \param verbose the verbose output flag
     * \param noLog if true no logarithmic normalization will be used
     */
    explicit distMatComputer( const std::string& roiFilename, const float thresholdRatio, const bool verbose, const bool noLog );

    //! Destructor
    ~distMatComputer() {}


    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets the input folder
     * \param inputFolder path to the folder where tractogram files are located
     */
    inline void setInputFolder( const std::string& inputFolder ) { m_inputFolder = inputFolder; }

    /**
     * sets the output folder
     * \param outputFolder folder where to write the created distance matrix files
     */
    inline void setOutputFolder( const std::string& outputFolder ) { m_outputFolder = outputFolder; }

    /**
     * sets (or resets) the verbose output flag in order to wrire progress information on the standard output
     * \param verbose the true/false flag to set the m_verbose member to
     */
    inline void setVerbose( bool verbose = true ) { m_verbose = verbose; }

    /**
     * sets (or resets) the very verbose output flag in order to wrire extra progress information on the standard output
     * \param veryVerbose the true/false flag to set the m_veryVerbose member to
     */
    inline void setVeryVerbose( bool veryVerbose = true ) { m_veryVerbose = veryVerbose; }

    /**
     * queries whether the roi file was loaded and therefore the class is ready for distance matrix computation
     * \return the ready flag, if true class
     */
    inline bool ready() const { return m_ready2go; }

    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * Calculates memory usage and size of the tractogram sub blocks to be loaded, and sub blocks to divide the matrix in
     * \param memory the available RAM memory in GBytes
     * \param blockSize the desired size for the distBlocks, if 0 maximum size for the available memory will be used
     */
    void setBlockSize( float memory, size_t blockSize = 0 );

    /**
     * Computes and writes the distance matrix
     */
    void doDistBlocks();




private:
    // === PRIVATE DATA MEMBERS ===

    bool            m_verbose;           //!< The verbose output flag. If true, additional and progress information will be shown through the standard output on execution.
    bool            m_niftiMode;         //!< The Nifti flag. If true, files are coordinates are in nifti reference frame, if False, in vista reference frame
    bool            m_roiLoaded;         //!< The roi loaded flag. If true, the seed voxel list was successfully read
    bool            m_veryVerbose;       //!< The very verbose output flag. If true, comprehensive progress information will be shown through the standard output on execution.
    bool            m_ready2go;          //!< The ready flag. If true, the class is ready to compute the distance matrix
    float           m_tractThreshold;    //!< The threshold to be applied to the tracts before computing dissimilarity measures, in normalized space. Calculated in class constructor
    float           m_logFactor;         //!< The Logarithmic factor to be applied when converting the tractogram data ftom log to natural units. Calculated in class constructor
    std::string     m_inputFolder;       //!< The folder path that contains the seed voxel tractograms
    std::string     m_outputFolder;      //!< The folder path where to write the output files

    HC_GRID         m_datasetGrid;       //!< Containe the type of coordinate frame (vista, nifti) the input data was stored in
    WHcoord         m_datasetSize;       //!< Contains the size in voxels of the dataset where the seed voxel coordinates correspond. Necessary for proper coordinate fam conversion. Taken from file
    size_t          m_numStreamlines;    //!< number of stramlines generated form each seed voxel, needed to properly do the transformation between natural units and logarithmic units in the tractograms

    std::vector< WHcoord >  m_coordinates; //!< A vector where the seed voxel coordinates are stored
    std::vector<size_t> m_trackids;        //!< Stores the ids of the seed tracts correesponding to each leaf
    std::vector< double >   m_leafNorms;   //!< A vector storing the comoputed norms for each seed voxel tractograms

    size_t m_blockSize;             //!< The number of elements in a row of each distance block
    size_t m_blocksPerRow;          //!< The number of distance blocks in a block-row of the distance matrix
    size_t m_subBlockSize;          //!< The number of tractograms that will be loaded at the same time
    size_t m_subBlocksPerBlock;     //!< The number of tract-subblocks in each distance block
    size_t m_trackSize;             //!< The number of points in a compact tractogram

    // === PRIVATE MEMBER FUNCTIONS ===


    /*

    // "fill_dist_block()": compute the distances from a distance block, computation will be divided into different sub-blocks due to memory issues
    void fill_dist_block(std::string &tract_dir,
                         std::vector<WHcoord> &seed_set_row,
                         std::vector<WHcoord> &seed_set_col,
                         std::vector<std::vector<float> > &dist_block,
                         const size_t &tract_block_size,
                         const size_t &tract_length,
                         const size_t &threads,
                         const bool &verbose,
                         const bool &veryvb);

    // "load_tract_block()": load single tractograms into a block
    void load_tract_block(std::string tract_dir,
                          std::vector<WHcoord> &seed_set,
                          const size_t  tract_block_size,
                          const size_t &tract_length,
                          unsigned char* tract_block);

    // "load_tract_block_transposed()": load single transposed tractograms into a block
    void load_tract_block_transposed(std::string tract_dir,
                                     std::vector<WHcoord> seed_set,
                                     const size_t  tract_block_size,
                                     const size_t &tract_length,
                                     unsigned char* tract_block);

    // "do_sqrsum()": precoumpute squared sums of tract blocks
    void do_sqrsum(double* sqrsum_array,
                   unsigned char* tract_block,
                   const size_t  tract_block_size,
                   const size_t &tract_length,
                   const size_t &threads);

    // "do_sqrsum_transposed()": precoumpute squared sums of transposed tract blocks
    void do_sqrsum_transposed(double* sqrsum_array,
                              unsigned char* tract_block,
                              const size_t  tract_block_size,
                              const size_t &tract_length,
                              const size_t &threads);


    // "compute_dist_block()":
    void compute_dist_block(std::vector<std::vector<float> > &dist_block,
                            unsigned char* row_tract_block,
                            unsigned char* col_tract_block,
                            double* row_sqrsum,
                            double* col_sqrsum,
                            const size_t  row_distb_offset,
                            const size_t  col_distb_offset,
                            const size_t  row_tblock_size,
                            const size_t  col_tblock_size,
                            const size_t &tract_length,
                            const size_t &threads);


    // fill_rand_block(): load distance block with random data (upper triangle)
    void fill_rand_block(std::vector<std::vector<float> > rand_tracts_row,
                         std::vector<std::vector<float> > rand_tracts_col,
                         std::vector<std::vector<float> > &dist_block,
                         const size_t rand_tract_length,
                         const size_t &threads,
                         const bool &verbose,
                         const bool &veryvb);

                         */

};


#endif // DISTMATCOMPUTER_H