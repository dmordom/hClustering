//---------------------------------------------------------------------------
//
// Project: hCLustering
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
// hCLustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------


#ifndef DISTBLOCK_H
#define DISTBLOCK_H

// std library
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>

// boost library
#include <boost/lexical_cast.hpp>

// hClustering
#include "WHcoord.h"
#include "WFileParser.h"
#include "fileManagerFactory.h"


#define MATRIX_INDEX_FILENAME "roi_index.txt"

/**
 * this class manages the reading and writing of blocks belonging to a dissimilarity matrix
 */
class distBlock
{
public:
    /**
     * Constructor
     * \param distBlockFolderInit folder where the distance matrix blocks are located / will be written to
     */
    explicit distBlock( const std::string& distBlockFolderInit );

    //! Destructor
    ~distBlock() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * returns size of of the loaded distance block in number of rows/columns
     * \return block size (number of rows/columns)
     */
    inline size_t size() const { return m_block.size(); }

    /**
     * returns the size of the complete distance matrix the loaded block is part of in number of rows/columns
     * \return full matrix size (number of rows/columns)
     */
    inline size_t matrixSize() const { return m_fullIndex.size(); }

    /**
     * returns the index ready flag, if true the matrix index file has been successfully loaded
     * \return index ready flag
     */
    inline bool indexReady() const { return m_indexReady; }

    /**
     * returns the block ready flag, if true the distance block has been successfully loaded
     * \return block ready flag
     */
    inline bool blockReady() const { return m_blockReady; }

    /**
     * returns the loaded block ID (a pair of integers indicating position within the full matrix)
     * \return block ID integer pair
     */
    inline std::pair< unsigned int, unsigned int > blockID() const { return m_blockID; }

    /**
     * returns the maximum block ID a block can have in the current distance matrix
     * \return maximum block ID of current matrix
     */
    inline unsigned int topBlock() const { return m_maxBlockID; }

    /**
     * returns the number of blocks contained in the current distance matrix
     * \return total number of existing blocks
     */
    inline unsigned int numBlocks() const { return ( ( ( m_maxBlockID + 1 ) * ( m_maxBlockID + 2 ) ) / 2 ); }

    // MEMBER FUNCTIONS


    /**
     * loads the block index that links a position in the block to a specific seed voxel coordinate
     * \return success bit, if true index was loaded successfully
     */
    bool readIndex();

    /**
     * loads the the distance block dissimilarity values and stores them in the corresponding data member
     * \param blockID integer pair with the block ID that is to be loaded
     */
    void loadBlock( std::pair< unsigned int, unsigned int > blockID );
    /**
     * \overload
     * \param blockID1 row position integer of the block ID
     * \param blockID2 column position integer of the block ID
     */
    void loadBlock( unsigned int blockID1, unsigned int blockID2 );
    /**
     * \overload
     * loads the block containing the distance value between two tracts defined by their seed coordinates
     * assumes that tractogram IDs follows same ordering as per seed voxel coordinates (z->y->x)
     * \param coord1 seed coordinate of first tract
     * \param coord2 seed coordinate of second tract
     */
    void loadBlock( const WHcoord& coord1, const WHcoord& coord2 );

    /**
     * fetches the distance value between two seed voxel tracts
     * \param coord1 seed coordinate of first tract
     * \param coord2 seed coordinate of second tract
     * \return distance value
     */
    float getDistance( const WHcoord& coord1, const WHcoord& coord2 );

    /**
     * returns the ranges of seed coordinates with contained within the currently loaded block
     * assumes that tractogram IDs follows same ordering as per seed voxel coordinates (z->y->x)

     * \return a pair of coordinate ranges (for rows and columns) each consisting of a pair of first and last coordinate of range
     */
    std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > getBlockRange();

    /**
     * writes the loaded block into the folder specified at object creation
     */
    void writeBlock();

    friend class fileManager;
    friend class vistaManager;
    friend class niftiManager;
    friend class cpcc;

protected:
private:
    // === PRIVATE DATA MEMBERS ===

    bool m_indexReady; //!< index successfully loaded flag
    bool m_blockReady; //!< block successfully loaded flag
    unsigned int m_maxBlockID;      //!< maximum block ID of the distance matrix
    std::string m_distBlockFolder;  //!< folder containing the distance matrix blocks (or where these blocks will be written to)
    std::map< WHcoord, std::pair< unsigned int, unsigned int > > m_fullIndex; //!< index map linking each seed voxel coordinate to a block and position within the block
    std::pair< unsigned int, unsigned int > m_blockID;  //!< block ID within the full matrix
    std::map< WHcoord, unsigned int > m_blockIndex1;    //!< index map linking each seed voxel coordinate a specific row position in the current block
    std::map< WHcoord, unsigned int > m_blockIndex2;    //!< index map linking each seed voxel coordinate a specific column position in the current block
    std::vector< std::vector< float > > m_block;        //!< block pairwise distance information matrix

    // === PRIVATE MEMBER FUNCTIONS ===

    /**
     * retrieve the path for the index file corresponding to the current distance matrix
     * \return the index file path
     */
    std::string getIndexFilename();

    /**
     * calculates based on the index information the block ID where the distance information of two pairs of seed coordinates is located
     * \param coord1 seed coordinate of first tract
     * \param coord2 seed coordinate of second tract
     * \return the block ID pair pointing to the block with the required information
     */
    std::pair< unsigned int, unsigned int > whichBlock( const WHcoord& coord1, const WHcoord& coord2 );
};

#endif  // DISTBLOCK_H
