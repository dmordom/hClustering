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
#include "vistaManager.h"

class distBlock
{
public:
    explicit distBlock( std::string distBlockFolderInit ) :
        m_indexReady( false ), m_blockReady( false ), m_maxBlockID( 0 ), m_distBlockFolder( distBlockFolderInit ), m_blockID(
                        std::make_pair( 0, 0 ) )
    {
        readIndex();
    }
    ~distBlock()
    {
    }

    size_t size() const
    {
        return m_block.size();
    }
    size_t matrixSize() const
    {
        return m_fullIndex.size();
    }

    bool indexReady() const
    {
        return m_indexReady;
    }

    bool blockReady() const
    {
        return m_blockReady;
    }

    std::pair< unsigned int, unsigned int > blockID() const
    {
        return m_blockID;
    }

    unsigned int topBlock() const
    {
        return m_maxBlockID;
    }

    unsigned int numBlocks() const
    {
        return ( ( ( m_maxBlockID + 1 ) * ( m_maxBlockID + 2 ) ) / 2 );
    }

    bool readIndex();

    void loadBlock( std::pair< unsigned int, unsigned int > blockID );
    void loadBlock( unsigned int blockID1, unsigned int blockID2 );
    void loadBlock( WHcoord coord1, WHcoord coord2 );

    float getDistance( WHcoord coord1, WHcoord coord2 );

    std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > getBlockRange();

    void writeBlock();

    friend class vistaManager;
    friend class cpcc;

protected:
private:
    bool m_indexReady;
    bool m_blockReady;
    unsigned int m_maxBlockID;
    std::string m_distBlockFolder;
    std::map< WHcoord, std::pair< unsigned int, unsigned int > > m_fullIndex;
    std::pair< unsigned int, unsigned int > m_blockID;
    std::map< WHcoord, unsigned int > m_blockIndex1;
    std::map< WHcoord, unsigned int > m_blockIndex2;
    std::vector< std::vector< float > > m_block;

    void getIndexFilename( std::string* indexFilename );

    std::pair< unsigned int, unsigned int > whichBlock( WHcoord coord1, WHcoord coord2 );
};

#endif  // DISTBLOCK_H
