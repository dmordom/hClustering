// std library
#include <vector>
#include <utility>
#include <map>
#include <string>

#include "distBlock.h"

// PUBLIC FUNCTIONS

bool distBlock::readIndex()
{

    m_fullIndex.clear();
    m_blockIndex1.clear();
    m_blockIndex1.clear();
    m_indexReady = false;
    m_blockReady = false;
    m_maxBlockID = 0;

    std::string indexFilename;
    getIndexFilename( &indexFilename );

    WFileParser parser( indexFilename );
    if( !parser.readFile() )
    {
        std::cerr << "ERROR @ distBlock::readIndex(): Parser error" << std::endl;
        return false;
    }

    std::vector< std::string > lines = parser.getRawLines();
    if( lines.size() == 0 )
    {
        std::cerr << "ERROR @ distBlock::readIndex(): Index file is empty" << std::endl;
        return false;
    }

    std::vector< std::vector< std::string > > indexStrings = parser.getLinesForTagSeparated( "distindex" );
    for( size_t i = 0; i < indexStrings.size(); ++i )
    {
        if( indexStrings[i].size() != 7 || indexStrings[i][3] != "b" || indexStrings[i][5] != "i" )
        {
            std::cerr << "ERROR @ distBlock::readIndex(): format of index file is not correct" << std::endl;
            return false;
        }
        WHcoord tempCoord( boost::lexical_cast< coord_t >( indexStrings[i][0] ), boost::lexical_cast< coord_t >(
                        indexStrings[i][1] ), boost::lexical_cast< coord_t >( indexStrings[i][2] ) );
        std::pair< unsigned int, unsigned int > tempPair( boost::lexical_cast< unsigned int >( indexStrings[i][4] ),
                        boost::lexical_cast< unsigned int >( indexStrings[i][6] ) );
        m_fullIndex[tempCoord] = tempPair;
        if( tempPair.first > m_maxBlockID )
            m_maxBlockID = tempPair.first;
    }

    m_indexReady = true;
    return true;
} // end "readIndex()" -----------------------------------------------------------------

void distBlock::loadBlock( std::pair< unsigned int, unsigned int > blockID )
{
    loadBlock( blockID.first, blockID.second );
    return;
} // end "loadblock()" -----------------------------------------------------------------

void distBlock::loadBlock( unsigned int blockID1, unsigned int blockID2 )
{
    if( !m_indexReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::loadblock(): Index was not previously loaded" );
    }
    else if( ( blockID1 > m_maxBlockID ) || ( blockID2 > m_maxBlockID ) )
    {
        throw std::runtime_error( "ERROR @ distBlock::loadblock(): blockID is out of bounds" );
    }
    else
    {
        if( blockID1 > blockID2 )
        {
            unsigned int tempID( blockID2 );
            blockID2 = blockID1;
            blockID1 = tempID;
        }

        // load distanceblock
        vistaManager vMngr( m_distBlockFolder );
        vMngr.readDistBlock( blockID1, blockID2, m_block );
        m_blockID = std::make_pair( blockID1, blockID2 );
        m_blockIndex1.clear();
        m_blockIndex2.clear();

        // fill block indices
        for( std::map< WHcoord, std::pair< unsigned int, unsigned int > >::const_iterator iter( m_fullIndex.begin() ); iter
                        != m_fullIndex.end(); ++iter )
        {
            if( iter->second.first == blockID1 )
                m_blockIndex1[iter->first] = iter->second.second;
            if( iter->second.first == blockID2 )
                m_blockIndex2[iter->first] = iter->second.second;
        }

        //set flags
        m_blockReady = true;
    }
    return;
} // end "loadblock()" -----------------------------------------------------------------

void distBlock::loadBlock( WHcoord coord1, WHcoord coord2 )
{
    std::pair< unsigned int, unsigned int > blockID = whichBlock( coord1, coord2 );
    loadBlock( blockID );
    return;
} // end "loadblock()" -----------------------------------------------------------------

float distBlock::getDistance( WHcoord coord1, WHcoord coord2 )
{
    unsigned int index1( 0 ), index2( 0 );

    if( !m_indexReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): Index was not previously loaded" );
    }
    else if( !m_blockReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): no block is loaded" );
    }
    else
    {
        std::map< WHcoord, unsigned int >::iterator findIter( m_blockIndex1.find( coord1 ) );
        if( findIter == m_blockIndex1.end() )
        {
            throw std::runtime_error( "ERROR @ distBlock::getDistance(): first coordinate is out of bounds" );
        }
        else
        {
            index1 = findIter->second;
        }

        findIter = ( m_blockIndex2.find( coord2 ) );
        if( findIter == m_blockIndex2.end() )
        {
            throw std::runtime_error( "ERROR @ distBlock::getDistance(): second coordinate is out of bounds" );
        }
        else
        {
            index2 = findIter->second;
        }

        return m_block[index1][index2];
    }
} // end "getDistance()" -----------------------------------------------------------------

std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > distBlock::getBlockRange()
{
    if( !m_indexReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): Index was not previously loaded" );
    }
    else if( !m_blockReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): no block is loaded" );
    }
    else
    {
        std::pair< std::pair< WHcoord, WHcoord >, std::pair< WHcoord, WHcoord > > outRange;
        outRange.first.first = ( m_blockIndex1.begin() )->first;
        outRange.first.second = ( --m_blockIndex1.end() )->first;
        outRange.second.first = ( m_blockIndex2.begin() )->first;
        outRange.second.second = ( --m_blockIndex2.end() )->first;
        return outRange;
    }
} // end "getBlockRange()" -----------------------------------------------------------------


void distBlock::writeBlock()
{
    if( !m_indexReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): Index was not previously loaded" );
    }
    else if( !m_blockReady )
    {
        throw std::runtime_error( "ERROR @ distBlock::getDistance(): no block is loaded" );
    }
    else
    {
        vistaManager vMngr( m_distBlockFolder );
        vMngr.writeDistBlock( m_blockID, m_block );
    }
    return;
} // end "writeBlock()" -----------------------------------------------------------------


// PRIVATE FUNCTIONS

void distBlock::getIndexFilename( std::string* indexFilename )
{
    *indexFilename = m_distBlockFolder + "/roi_index.txt";
    return;
} // end "getIndexFilename()" -----------------------------------------------------------------

std::pair< unsigned int, unsigned int > distBlock::whichBlock( WHcoord coord1, WHcoord coord2 )
{
    unsigned int rowBlockID, colBlockID;
    std::map< WHcoord, std::pair< unsigned int, unsigned int > >::iterator findIter( m_fullIndex.find( coord1 ) );
    if( findIter == m_fullIndex.end() )
    {
        throw std::runtime_error( "ERROR @ distBlock::whichBlock(): first coordinate is out of bounds" );
    }
    else
        rowBlockID = findIter->second.first;

    findIter = ( m_fullIndex.find( coord2 ) );
    if( findIter == m_fullIndex.end() )
    {
        throw std::runtime_error( "ERROR @ distBlock::whichBlock(): second coordinate is out of bounds" );
    }
    else
    {
        colBlockID = findIter->second.first;
    }

    if( colBlockID < rowBlockID )
    {
        unsigned int tempID( rowBlockID );
        rowBlockID = colBlockID;
        colBlockID = tempID;
    }

    return std::make_pair( rowBlockID, colBlockID );
} // end "whichBlock()" -----------------------------------------------------------------


