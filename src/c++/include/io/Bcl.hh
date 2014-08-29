/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Reader/Writer component for BCL files.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_IO_BCL_HH
#define EAGLE_IO_BCL_HH


#include <string>
#include <vector>
#include "common/Exceptions.hh"


namespace eagle
{
namespace io
{


class BclTile
{
public:
    BclTile( const unsigned long long expectedReadCount, const unsigned int clusterLength, const std::string &filenameTemplate, const std::string &statsFilenameTemplate, const std::string &filterFilename, const std::string &clocsFilename, const std::string &controlFilename, const bool verbose=true );
    void addClusterToRandomLocation( const char *bufCluster, const bool isPassingFilter = true );
    void flushToDisk() const;

private:
    void writeBclFile( const unsigned int cycle ) const;
    void writeStatsFile( const unsigned int cycle ) const;
    void writeFilterFile() const;
    void writeClocsFile() const;
    void writeControlFile() const;

    unsigned long long expectedReadCount_;
    unsigned int clusterLength_;
    std::string filenameTemplate_;
    std::string statsFilenameTemplate_;
    std::string filterFilename_;
    std::string clocsFilename_;
    std::string controlFilename_;
    //    vector<bool> usedLocations;
    std::vector<std::vector<unsigned int> > stats_;
    char *ramTile_;
    std::vector<char> passFilter_;
};



/*

#include <fstream>
#include <memory>

// silence boost::sysem error code warnings
//#ifndef BOOST_SYSTEM_NO_DEPRECATED 
//#define BOOST_SYSTEM_NO_DEPRECATED
//#endif

#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

//#include "flowcell/TileMetadata.hh"


class BclReader : public boost::ptr_vector<std::ifstream>
{
public:
    //BclReader(const isaac::flowcell::Tile &tile);
    BclReader(const isaac::flowcell::TileMetadata &tile, unsigned int firstCycle, unsigned int lastCycle);
    void next(std::string &read);
    static unsigned int getClusterCount(const boost::filesystem::path &bclFilePath);
private:
    const std::vector<boost::filesystem::path> bclFilePathList_;
    static std::ifstream &open(std::ifstream &is, const boost::filesystem::path &bclFilePath);
    static unsigned int getClusterCount(std::ifstream &is, const boost::filesystem::path &bclFilePath);
};
*/

} // namespace io
} // namespace eagle

#endif // EAGLE_IO_BCL_HH
