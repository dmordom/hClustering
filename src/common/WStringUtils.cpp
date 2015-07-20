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

#include <algorithm>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "WStringUtils.h"

const std::string string_utils::WHITESPACE( "\r\n\t " ); //!< These characters are regarded as whitespaces

std::string string_utils::rTrim( const std::string &source, const std::string& t )
{
    std::string str = source;
    str.erase( str.find_last_not_of( t ) + 1 );
    return str;
}

std::string string_utils::lTrim( const std::string& source, const std::string& t )
{
    std::string str = source;
    str.erase( 0 , source.find_first_not_of( t ) );
    return str;
}

std::string string_utils::trim( const std::string& source, const std::string& t )
{
    std::string str = source;
    return lTrim( rTrim( str , t) , t );
}

std::string string_utils::toUpper( const std::string& source )
{
    std::string str = source;
    std::transform( source.begin(), source.end(), str.begin(), ::toupper );
    return str;
}

std::string string_utils::toLower( const std::string& source )
{
    std::string str = source;
    std::transform( source.begin(), source.end(), str.begin(), ::tolower );
    return str;
}

std::vector< std::string > string_utils::tokenize( const std::string& source,
                                                   const std::string& delim,
                                                   bool compress )
{
    std::vector< std::string > result;
    namespace ba = boost::algorithm;
    ba::token_compress_mode_type compression = ba::token_compress_on;
    if( !compress )
    {
        compression = ba::token_compress_off;
    }
    ba::split( result, source, ba::is_any_of( delim ), compression );
    if( !result.empty() )
    {
        // NOTE: moved the back() command to another if as if compiled on Windows, OW crashes since the compiler does not stop evaluation of the
        // condition after the first statement evaluates to false.

        if( result.back() == "" )
        {
            result.pop_back();
        }
    }
    return result;
}
