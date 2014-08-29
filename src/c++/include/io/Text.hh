/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Base for many text (ie., non-binary) files.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_IO_TEXT_HH
#define EAGLE_IO_TEXT_HH

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/utility.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "common/Logger.hh"


namespace eagle
{
namespace io
{

/*
 * \brief Generic delimiter-separated-values file
 */
class DsvReader: public std::ifstream, boost::noncopyable
{
public:
    DsvReader(const std::vector<boost::filesystem::path> &pathList);
    DsvReader(const boost::filesystem::path &path);
    DsvReader() : lineCount_(0) {}
    size_t size() const {return pathList_.size();}
    boost::filesystem::path file(unsigned int i) const {assert(pathList_.size() > i); return pathList_[i];}
    std::string getNextLine(char comment = '#');

    template<char C>  // C = delimiter
    bool getNextLineFields(std::vector<std::string> &fields, char comment = '#')
    {
        std::string line = getNextLine( comment );
        if (!line.empty())
        {
            boost::split(fields, line, boost::is_any_of( std::string(1,C) ));
            return true;
        }
        return false;
    }

protected:
    bool openNextFile();

    const std::vector<boost::filesystem::path> pathList_;
    std::vector<boost::filesystem::path>::const_iterator thisPath_;
    unsigned long lineCount_;

private:
    void openFirstFile();
};


/*
 * \brief Generic delimiter-separated-values file
 */
class DsvWriter: public std::ofstream, boost::noncopyable
{
public:
    DsvWriter(const std::vector<boost::filesystem::path> pathList, const bool overwrite = false)
    : pathList_(pathList), overwrite_(overwrite) {}
    DsvWriter(const boost::filesystem::path singlePath, const bool overwrite = false)
    : pathList_(1,singlePath), overwrite_(overwrite)   {}
    DsvWriter(const bool overwrite = false) : overwrite_(overwrite) {}
    void open(unsigned int i);
    boost::filesystem::path file(unsigned int i) const {assert(pathList_.size() > i); return pathList_[i];}
    std::vector<boost::filesystem::path>::const_iterator begin() const {return pathList_.begin();}
    std::vector<boost::filesystem::path>::const_iterator end() const {return pathList_.end();}
protected:
    const std::vector<boost::filesystem::path> pathList_;
    const bool overwrite_;
    std::vector<boost::filesystem::path>::const_iterator thisPath_;
};


} // namespace io
} // namespace eagle

#endif // EAGLE_IO_TEXT_HH
