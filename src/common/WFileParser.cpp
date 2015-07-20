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
// This file is also part of the
// Whole-Brain Connectivity-Based Hierarchical Parcellation Project
// David Moreno-Dominguez
// moreno@cbs.mpg.de
// www.cbs.mpg.de/~moreno
//
//---------------------------------------------------------------------------

#include <string>
#include <vector>

#include "WFileParser.h"


WFileParser::WFileParser( const std::string fileName ) :
    m_fileName( fileName ),
    m_tagIndicator( "#" ),
    m_endIndicator( "end" ),
    m_delimiter( " " )
{
}

WFileParser::~WFileParser()
{
}

bool WFileParser::readFile()
{
    using namespace boost::filesystem; //NOLINT

    if( !exists( m_fileName ) )
    {
        //debugLog() << "file doesn't exist";
        return false;
    }

    std::ifstream ifs( m_fileName.c_str(), std::ifstream::in );

    std::string line;

    if( !ifs.is_open() )
    {
        //debugLog() << "file open failed";
        ifs.close();
        return false;
    }

    while( !ifs.eof() )
    {
        getline( ifs, line );

        m_rawLines.push_back( std::string( line ) );
    }

    ifs.close();

    return true;
}

std::vector<std::string>WFileParser::getLinesForTag( std::string tag )
{
    std::string startTag = m_tagIndicator + tag;
    std::string endTag = m_tagIndicator + m_endIndicator + tag;

    std::vector<std::string>returnVector;

    size_t i = 0;
    while( i < m_rawLines.size() )
    {
        if( m_rawLines[i] == startTag )
        {
            //debugLog() << "coordinates tag at line " << i;
            ++i;
            break;
        }
        else
        {
            ++i;
        }
    }
    while( i < m_rawLines.size() )
    {
        if( m_rawLines[i] == endTag )
        {
            //debugLog() << "endcoordinates tag at line " << i;
            ++i;
            break;
        }
        else
        {
            returnVector.push_back( m_rawLines[i] );
            ++i;
        }
    }
    return returnVector;
}

std::vector< std::vector<std::string> >WFileParser::getLinesForTagSeparated( std::string tag )
{
    std::string startTag = m_tagIndicator + tag;
    std::string endTag = m_tagIndicator + m_endIndicator + tag;

    std::vector<std::vector<std::string > >returnVector;

    size_t i = 0;
    while( i < m_rawLines.size() )
    {
        if( m_rawLines[i] == startTag )
        {
            //debugLog() << "coordinates tag at line " << i;
            ++i;
            break;
        }
        else
        {
            ++i;
        }
    }

    while( i < m_rawLines.size() )
    {
        if( m_rawLines[i] == endTag )
        {
            //debugLog() << "endcoordinates tag at line " << i;
            ++i;
            break;
        }
        else
        {
            std::vector<std::string>svec;
            boost::regex reg( m_delimiter );
            boost::sregex_token_iterator it( m_rawLines[i].begin(), m_rawLines[i].end(), reg, -1 );
            boost::sregex_token_iterator end;
            while( it != end )
            {
                svec.push_back( *it++ );
            }

            returnVector.push_back( svec );
            ++i;
        }
    }
    return returnVector;
}
