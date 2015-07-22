#ifndef DO_DIST_BLOCKS_H
#define DO_DIST_BLOCKS_H

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

//functions
#include "do_Vista.h"
#include "get_names.h"
#include "tract_dist.h"

extern float threshold;

// "do_dist_blocks()": compute and write down the distance blocks
void do_dist_blocks(std::string tract_dir,
                    std::string out_dir,
                    const size_t roiNum,
                    const size_t bsize,
                    const size_t tract_block_size,
                    const size_t tract_length,
                    const size_t threads,
                    const bool rand_mode,
                    const bool verbose,
                    const bool veryvb);

// "fill_dist_block()": compute the distances from a distance block, computation will be divided into different sub-blocks due to memory issues
void fill_dist_block(std::string &tract_dir,
                     std::vector<size_t> &seed_set_row,
                     std::vector<size_t> &seed_set_col,
                     std::vector<std::vector<float> > &dist_block,
                     const size_t &tract_block_size,
                     const size_t &tract_length,
                     const size_t &threads,
                     const bool &verbose,
                     const bool &veryvb);

// "load_tract_block()": load single tractograms into a block
void load_tract_block(std::string tract_dir,
                      std::vector<size_t> &seed_set,
                      const size_t  tract_block_size,
                      const size_t &tract_length,
                      unsigned char* tract_block);

// "load_tract_block_transposed()": load single transposed tractograms into a block
void load_tract_block_transposed(std::string tract_dir,
                                 std::vector<size_t> seed_set,
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


#endif // DO_DIST_BLOCKS_H
