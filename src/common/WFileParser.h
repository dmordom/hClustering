//---------------------------------------------------------------------------
//
// Project: OpenWalnut ( http://www.openwalnut.org )
//
// Copyright 2009 OpenWalnut Community, BSV@Uni-Leipzig and CNCF@MPI-CBS
// For more information see http://www.openwalnut.org/copying
//
// This file is part of OpenWalnut.
//
// OpenWalnut is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// OpenWalnut is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenWalnut. If not, see <http://www.gnu.org/licenses/>.
//
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//
// Project: hClustering
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
// hClustering is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
//---------------------------------------------------------------------------

#ifndef WFILEPARSER_H
#define WFILEPARSER_H

// std library
#include <fstream>
#include <string>
#include <vector>

// boost library
// Use filesystem version 2 for compatibility with newer boost versions.
#if 0
#ifndef BOOST_FILESYSTEM_VERSION
    #define BOOST_FILESYSTEM_VERSION 2
#endif
#endif

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>


/**
 * this class implements text file loading and several convinience methods to access the recovered text
 */
class WFileParser
{
public:
    /**
     * constructor
     * \param fileName
     */
    explicit WFileParser( const std::string fileName );

    /**
     * destructor
     */
    ~WFileParser();

    /**
     * helper function to read a text file
     * \return string containing the file
     */
    bool readFile();

    /**
     * getter
     * \return the content of the loaded as vector of strings for each line
     */
    std::vector<std::string>getRawLines();

    /**
     * getter
     * \param tag tag marking a certain type of content
     * \return vector of strings for each line for the tag
     */
    std::vector<std::string>getLinesForTag( std::string tag );

    /**
     * getter
     * \param tag tag marking a certain type of content
     * \return same as getLinesForTag but each line represented as vector of strings
     */
    std::vector< std::vector<std::string> >getLinesForTagSeparated( std::string tag );

protected:
private:
    std::string m_fileName; //!< the file name of the file to parse

    std::vector<std::string>m_rawLines; //!< vector of every line in the file

    std::string m_tagIndicator; //!< string marking a line as tag

    std::string m_endIndicator; //!< string marking the end of a tagged area

    std::string m_delimiter; //!< delimiter for entries in a line
};

inline std::vector<std::string> WFileParser::getRawLines()
{
    return m_rawLines;
}

#endif  // WFILEPARSER_H
