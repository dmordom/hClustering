#ifndef LISTEDCACHE_H
#define LISTEDCACHE_H

// std library
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <stdexcept>

template< typename T > class listedCache
{
    typedef typename std::list< T >::iterator valueIter_t;
    typedef typename std::list< size_t >::iterator fifoIter_t;

public:
    listedCache( size_t listSize, size_t sizeLimitInit = 0 ) :
        m_sizeLimit( sizeLimitInit )
    {
        m_tracker.resize( listSize, std::make_pair( false, std::make_pair( m_values.begin(), m_fifo.begin() ) ) );
    }
    ~listedCache()
    {
    }

    // returns number of elements currently stored in cache
    size_t size() const
    {
        return m_values.size();
    }

    // returns the maximum number of objects to be stored in cache after a call to cleanup()
    size_t limit() const
    {
        return m_sizeLimit;
    }

    // sets the maximum number of objects to be stored in cache after a call to cleanup()
    void setLimit( size_t sizeLimit )
    {
        m_sizeLimit = sizeLimit;
    }

    // returns the index of the oldest accessed element
    size_t oldest()
    {
        // return index to oldest object
        return m_fifo.front();
    }

    // returns a pointer to the element associated with index
    T* get( const size_t index )
    {
        if( index >= m_tracker.size() )
            throw std::runtime_error( "ERROR @ listedCache::get(): index is out of bounds" );
        // if object is not in the cache, return null pointer
        if( !m_tracker[index].first )
            return 0;
        //move index to back of fifo (it was used)
        m_fifo.splice( m_fifo.end(), m_fifo, m_tracker[index].second.second );
        m_tracker[index].second.second = ( --m_fifo.end() );
        // return pointer to the object
        return &( *( m_tracker[index].second.first ) );
    }

    // returns a pointer to the element associated with index but doesnt update its udage
    T* getNoUpdate( const size_t index ) const
    {
        if( index >= m_tracker.size() )
            throw std::runtime_error( "ERROR @ listedCache::get(): index is out of bounds" );
        // if object is not in the cache, return null pointer
        if( !m_tracker[index].first )
            return 0;
        // return pointer to the object
        return &( *( m_tracker[index].second.first ) );
    }

    // insterts a new object and returns a pointer to the stored value
    T* insert( const size_t index, const T &value )
    {
        if( index >= m_tracker.size() )
            throw std::runtime_error( "ERROR @ listedCache::insert(): index is out of bounds" );
        if( has( index ) )
        {
            std::cerr << "WARNING @ listedCache::insert(): element is already loaded in cache, new element not inserted"
                            << std::endl;
            return 0;
        }
        // insert object
        m_tracker[index].second.second = m_fifo.insert( m_fifo.end(), index );
        m_tracker[index].second.first = m_values.insert( m_values.end(), value );
        m_tracker[index].first = true;
        return &( *( m_tracker[index].second.first ) );
    }

    // erases the element associated with index
    void erase( const size_t index )
    {
        if( index >= m_tracker.size() )
            throw std::runtime_error( "ERROR @ listedCache::insert(): index is out of bounds" );
        if( !has( index ) )
            return;
        // eliminate object from fifo and cache
        m_fifo.erase( m_tracker[index].second.second );
        m_values.erase( m_tracker[index].second.first );
        m_tracker[index].first = false;
        return;
    }

    // erases all the oldest elements while the size is over the specified limit
    void cleanup()
    {
        if( m_sizeLimit == 0 && m_values.size() > 0 )
        {
            clear();
            return;
        }
        while( m_values.size() > m_sizeLimit )
        {
            // eliminate object from fifo and cache
            size_t index( m_fifo.front() );
            m_fifo.pop_front();
            m_values.erase( m_tracker[index].second.first );
            m_tracker[index].first = false;
        }
        return;
    }

    // erases all the elements and resets the tracker vector
    void clear()
    {
        m_fifo.clear();
        m_values.clear();
        m_tracker.assign( m_tracker.size(), std::make_pair( false, std::make_pair( m_values.begin(), m_fifo.begin() ) ) );
        return;
    }

    // erases all the elements and frees the memory of the tracker vector, this cache object may not be used any more
    void shutdown()
    {
        m_fifo.clear();
        m_values.clear();
        {
            std::vector< std::pair< bool, std::pair< valueIter_t, fifoIter_t > > > emptytracker;
            m_tracker.swap( emptytracker );
        }
        return;
    }

private:
    size_t m_sizeLimit;
    std::list< size_t > m_fifo;
    std::list< T > m_values;
    std::vector< std::pair< bool, std::pair< valueIter_t, fifoIter_t > > > m_tracker;

    bool has( const size_t index )
    {
        return ( m_tracker[index].first );
    }

};

#endif // LISTEDCACHE_H
