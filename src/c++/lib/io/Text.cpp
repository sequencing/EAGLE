/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component to read/write text files (specific protocols to be derived from here).
 **
 ** \author Mauricio Varea
 **/

#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "io/Text.hh"
#include "common/Exceptions.hh"
#include "common/Logger.hh"


namespace eagle
{
namespace io
{


DsvReader::DsvReader(const std::vector<boost::filesystem::path> &pathList)
: pathList_(pathList)
, thisPath_(pathList_.begin())
, lineCount_(0)
{
    if (pathList.size() > 0)
    {
        openFirstFile();
    }
}

DsvReader::DsvReader(const boost::filesystem::path &path)
: pathList_(std::vector<boost::filesystem::path>(1,path))
, thisPath_(pathList_.begin())
, lineCount_(0)
{
    openFirstFile();
}

void DsvReader::openFirstFile()
{
    open(thisPath_->string().c_str());
    if (!is_open())
    {
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, (boost::format("Failed to open file %s for reading") % *thisPath_).str()));
    }
}

bool DsvReader::openNextFile()
{
    close();
    if (pathList_.end() != thisPath_)
    {
        ++thisPath_;
    }
    if (pathList_.end() != thisPath_)
    {
        open(thisPath_->string().c_str());
        lineCount_ = 0;
    }
    return is_open();
}

std::string DsvReader::getNextLine(char comment)
{
    std::string line;
    do {
        std::getline(*this,line);
        while (fail())
        {
            if (openNextFile())
            {
                std::getline(*this,line);
            }
            else
            {
                EAGLE_WARNING_IF( pathList_.end() != thisPath_,
                                  (boost::format("Could not read beyond %s:%lu") % *(thisPath_ - 1) % lineCount_ ).str() );
                return "";
            }
        }
        ++lineCount_;
    } while (line.empty() || comment == line[0]);
    return line;
}


void DsvWriter::open(unsigned int i)
{
    close();
    if (i >= pathList_.size())
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException( (boost::format("Tried to open stream #%d, but only %d files were provided")
                                                           % (i + 1) % pathList_.size()).str()) );
    }
    thisPath_ = pathList_.begin() + i;
    if ( boost::filesystem::exists(*thisPath_) )
    {
        if (overwrite_)
        {
            EAGLE_WARNING( "Overwriting " << *thisPath_ << " due to the --force switch." );
        } else {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Cannot write to %s: File already exists!") % *thisPath_).str()));
        }
    }
    std::ofstream::open(thisPath_->string().c_str(),std::ios_base::out | std::ios_base::trunc);
    if (!std::ofstream::is_open())
    {
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, (boost::format("Failed to open file %s for writing") % *thisPath_).str()));
    }
}



} // namespace io
} // namespace eagle
