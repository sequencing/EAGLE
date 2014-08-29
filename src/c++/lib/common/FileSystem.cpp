/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Time display.
 **
 ** \author Mauricio Varea
 **/

#include <fstream>
#include <algorithm>
#include <boost/format.hpp>

#include "common/Logger.hh"
#include "common/FileSystem.hh"

namespace eagle
{
namespace common
{

/*
 * \brief Non-recursive GLOB operation
 */
std::vector<boost::filesystem::path> Glob::glob( const boost::filesystem::path& dir ) const
{
    std::vector<boost::filesystem::path> files;
    if (boost::filesystem::is_directory(dir))
    {
        for( boost::filesystem::directory_iterator it(dir), end;
             it != end;
             ++it)
        {
            boost::smatch match;
            if (boost::regex_match( it->path().leaf().string(), match, pattern_ ))
            {
                EAGLE_DEBUG(4, "... " << it->path().leaf() );
                files.push_back( it->path() );
            }
        }
    } else {
        boost::smatch match;
        if (boost::regex_match( dir.leaf().string(), match, pattern_ ))
        {
            EAGLE_DEBUG(4, "... " << dir.leaf() );
            // ... return vector with only one element
            files.push_back( dir );
        }
    }
    if (files.empty())
    {
        EAGLE_WARNING( (boost::format("Regex \"%s\" did not match any files in %s") % pattern_ % dir ).str() );
    } else {
        std::stable_sort(files.begin(),files.end());
    }
    return files;
}



} // namespace common
} // namespace eagle
