/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description A preprocessor-based Logger + time display.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_COMMON_FILE_SYSTEM_HH
#define EAGLE_COMMON_FILE_SYSTEM_HH

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>


namespace eagle
{
namespace common
{

class Glob
{
public:
//    Glob( boost::regex pattern ) : pattern_(pattern) {}
    Glob( std::string pattern = ".*" ) : pattern_(pattern, boost::regex_constants::mod_x) {}

    std::vector<boost::filesystem::path> glob( const boost::filesystem::path& dir ) const;
private:
    boost::regex pattern_;
};


} // namespace common
} // namespace eagle

#endif // EAGLE_COMMON_FILE_SYSTEM_HH
