//---------------------------------------------------------------------------
//
// Project: hClustering
//
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// d.mor.dom@gmail.com
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno//
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
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------

#ifndef LISTEDCACHE_H
#define LISTEDCACHE_H

// std library
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <stdexcept>

/**
 * This class implements a time and memory efficient cache list template where once the max size is reached, the entry that was least recently accessed is eliminated when a new item is added.
 * the data is saved in a list for fast adding/deleting, but each item's memory address is saved in a vector (where item ID = position in the vector) fast keyword search and memory access
 * a separate fifo list keeps track of the last usage order of each item. (items memory adress in the fifo list is also kpt in the lookup vector mentioned before)
 */
template< typename T > class listedCache
{
    typedef typename std::list< T >::iterator valueIter_t;
    typedef typename std::list< size_t >::iterator fifoIter_t;

public:
    /**
     * Constructor
     * \param listSize total number of different entries that may at some point be introduced in the cache list (as in number of different possible IDs)
     * \param sizeLimitInit maximum size (as in maximum number of entries) that the cache list is allowed to have
     */
    listedCache( const size_t listSize, const size_t sizeLimitInit = 0 ) : m_sizeLimit( sizeLimitInit )
    {
        m_tracker.resize( listSize, std::make_pair( false, std::make_pair( m_values.begin(), m_fifo.begin() ) ) );
    }

     //! Destructor
    ~listedCache() {}

    // === IN-LINE MEMBER FUNCTIONS ===

    /**
     * sets the maximum number of objects to be stored in cache [after a call to cleanup()]
     * \param sizeLimit new max cache size
     */
    inline void setLimit( size_t sizeLimit ) { m_sizeLimit = sizeLimit; }

    /**
     * returns number of elements currently stored in cache (can be temporarily greater than the size limit)
     * \return current cache size
     */
    inline size_t size() const { return m_values.size(); }

    /**
     * returns the maximum number of objects to be stored in cache [after a call to cleanup()]
     * \return max cache size
     */
    inline size_t limit() const { return m_sizeLimit; }

    /**
     * returns the index of the oldest accessed (least recently accessed) element
     * \return least recently accessed element ID
     */
    inline size_t oldest() { return m_fifo.front(); }

    /**
     * checks whether a specific element is contained in the cache list
     * \param index an index identifying the desired entry (usually an ID number)
     * \return true if the element is stored in the list, false otherwise
     */
    inline bool has( const size_t index ) { return ( m_tracker[index].first ); }


    // === PUBLIC MEMBER FUNCTIONS ===

    /**
     * fetches the object identified by the index parameter and updates its relative usage status in the fifo
     * \param index an index identifying the desired entry (usually an ID number)
     * \return a pointer to the element associated with index
     */
    T* const get( const size_t index )
    {
        T* const objecTpointer = getNoUpdate( index );
        if( objecTpointer != NULL )
        {
            //move index to back of fifo (it was used)
            m_fifo.splice( m_fifo.end(), m_fifo, m_tracker[index].second.second );
            m_tracker[index].second.second = ( --m_fifo.end() );
        }
        return objecTpointer;
    }

    /**
     * fetches the object identified by the index parameter but does not update its relative usage status
     * \param index an index identifying the desired entry (usually an ID number)
     * \return a pointer to the element associated with index
     */
    T* const getNoUpdate( const size_t index ) const
    {
        if( index >= m_tracker.size() )
            throw std::runtime_error( "ERROR @ listedCache::get(): index is out of bounds" );
        // if object is not in the cache, return null pointer
        if( !m_tracker[index].first )
        {
            return NULL;
        }
        else
        {
            // return pointer to the object
            return &( *( m_tracker[index].second.first ) );
        }
    }

    /**
     * inserts a new object and returns a pointer to the stored value
     * \param index an index identifying the new entry (usually an ID number)
     * \param value the data of the new entry
     * \return a pointer to the new element inside the cache
     */
    T* const insert( const size_t index, const T& value )
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

    /**
     * erases the element entry associated with index
     * \param index an index identifying the entry to be deleted (usually an ID number)
     */
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

    /**
     * if the cache list size is over the specified limit, iteratively erases the least recently accessed elements until the list size matches the limit
     */
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

    /**
     * erases all the elements in the list and resets the tracker vector
     */
    void clear()
    {
        m_fifo.clear();
        m_values.clear();
        m_tracker.assign( m_tracker.size(), std::make_pair( false, std::make_pair( m_values.begin(), m_fifo.begin() ) ) );
        return;
    }

    /**
     * erases all the elements and frees the memory of the tracker vector, this cache object may not be used any more
     */
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
    // === PRIVATE DATA MEMBERS ===

    size_t              m_sizeLimit; //!< max cache size
    std::list< size_t > m_fifo;      //!< fifo list for relative usage information
    std::list< T >      m_values;    //!< list containing the element data values
    std::vector< std::pair< bool, std::pair< valueIter_t, fifoIter_t > > > m_tracker;  //!< lookup vector storing for each ID the information on  whether that element is stored in the list, the address of the data value and the address of the related fifo entry
};

#endif // LISTEDCACHE_H
