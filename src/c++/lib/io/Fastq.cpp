/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Writer component for FASTQ files.
 **
 ** \author Lilian Janin
 **/

#include <iostream>
#include <fstream>
#include <cerrno>
#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include "io/Fastq.hh"

using namespace std;


namespace eagle
{
namespace io
{


FastqTile::FastqTile( const unsigned long long expectedReadCount, const unsigned int clusterLength, const string &filenameTemplate, const string &statsFilenameTemplate, const string &filterFilename, const string &clocsFilename, const string &controlFilename, const bool verbose )
        : expectedReadCount_     ( expectedReadCount )
        , clusterLength_         ( clusterLength )
        , filenameTemplate_      ( filenameTemplate )
        , statsFilenameTemplate_ ( statsFilenameTemplate )
        , filterFilename_        ( filterFilename )
        , clocsFilename_         ( clocsFilename )
        , controlFilename_       ( controlFilename )
        , stats_                 ( clusterLength_, std::vector<unsigned int>(4, 0) )
        , passFilter_            ( expectedReadCount_, '\0' )
    {
        if (verbose)
        {
            clog << (boost::format("Creating new Fastq tile as %s, expecting %d reads") % filenameTemplate_ % expectedReadCount_).str() << endl;
        }
        assert (expectedReadCount_ < static_cast<unsigned int>(0xFFFFFFFF) && "Tile too large: FASTQ filter files can only contain 2^32 entries per tile");
        //boost::filesystem::create_directories( path_ );

        ramTile_ = new char[clusterLength_*expectedReadCount_];
        assert (ramTile_ != 0);
    }

    void FastqTile::addClusterToRandomLocation( const char *bufCluster, const bool isPassingFilter )
    {
        static unsigned long long nextPos = 0; // TODO: make this random
        if (nextPos >= expectedReadCount_)
        {
            BOOST_THROW_EXCEPTION( eagle::common::OutOfLimitsException( "Trying to add a cluster to a full tile" ) );
        }
        for (unsigned int i=0; i<clusterLength_; ++i)
        {
            ramTile_[nextPos+expectedReadCount_*i] = bufCluster[i];
        }
        if (isPassingFilter) {
            passFilter_[nextPos] = '\1';
        }
        nextPos++;
    }

    void FastqTile::flushToDisk() const
    {
        clog << "Flushing tile to disk" << endl;
        for (unsigned int i=0; i<clusterLength_; ++i)
        {
            writeFastqFile(i);
            writeStatsFile(i);
        }
        writeFilterFile();
        writeClocsFile();
        writeControlFile();
    }

    void FastqTile::writeFastqFile( const unsigned int cycle ) const
    {
        string filename = (boost::format( filenameTemplate_) % (cycle+1)).str();
        ofstream os( filename.c_str() );
        if( !os.good() )
        {
            cerr << "Can't write to " << filename << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
        os.write( (const char*)&expectedReadCount_, 4);
        os.write( &ramTile_[expectedReadCount_*cycle], expectedReadCount_);
    }

    void FastqTile::writeStatsFile( const unsigned int cycle ) const
    {
        //        clog << "Writing stats file" << endl;
        string filename = (boost::format( statsFilenameTemplate_) % (cycle+1)).str();
        ofstream statsWriter( filename.c_str() );
        if( !statsWriter.good() )
        {
            cerr << "Can't write to " << filename << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
        double intensity = 25.0;
        int zero = 0;
        const char * tmp = reinterpret_cast<const char *>(&clusterLength_);
        statsWriter.write(tmp, sizeof(clusterLength_));
        tmp = reinterpret_cast<char *>(&intensity);

        for (unsigned i = 0; 9 > i; ++i)
            statsWriter.write(tmp, sizeof(intensity));

        for (unsigned i = 0; 4 > i; ++i)
        {
            tmp = reinterpret_cast<const char *>(&stats_[clusterLength_ - 1][i]);
            statsWriter.write(tmp, sizeof(stats_[clusterLength_ - 1][i]));
        }

        tmp = reinterpret_cast<char *>(&zero);
        //std::cerr << "0 stats, clusterLength_ " << clusterLength_ << std::endl;
        statsWriter.write(tmp, sizeof(zero));
        //std::cerr << "1 stats, clusterLength_ " << clusterLength_ << std::endl;

        for (unsigned i = 0; 4 > i; ++i)
        {
            tmp = reinterpret_cast<const char *>(&stats_[clusterLength_ - 1][i]);
            statsWriter.write(tmp, sizeof(stats_[clusterLength_ - 1][i]));
        }
    }

    void FastqTile::writeFilterFile() const
    {
        clog << "Writing filter file" << endl;

        ofstream os( filterFilename_.c_str() );
        if( !os.good() )
        {
            cerr << "Can't write to " << filterFilename_ << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
        unsigned int header1 = 0, header2 = 3;
        assert(sizeof(unsigned int) == 4);
        os.write( (const char*)&header1, 4);
        os.write( (const char*)&header2, 4);
        os.write( (const char*)&expectedReadCount_, 4);
        os.write( &passFilter_[0], expectedReadCount_);
    }

    void FastqTile::writeClocsFile() const
    {
        // This is the .clocs format
        // max 2048x20000 positions divided in bins of 25x25 positions => 82 bins wide x 800 bins high
        // each bin must contain no more than 255 clusters
        // what we store are the position offsets in each bin
        clog << "Writing clocs file" << endl;

        ofstream os( clocsFilename_.c_str() );
        if( !os.good() )
        {
            cerr << "Can't write to " << clocsFilename_ << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
        unsigned int leftToProcess = expectedReadCount_;
        assert(sizeof(unsigned int) == 4);
        const unsigned char header1 = 1;
        const unsigned int binCount = (expectedReadCount_ / 255) + ((expectedReadCount_ % 255) != 0);
        os.write( (const char*)&header1, 1);
        os.write( (const char*)&binCount, 4);
        for (unsigned int i=0; i<binCount; ++i)
        {
            const unsigned char clusterCount = min<unsigned int>( leftToProcess, 255 );
            os.write( (const char*)&clusterCount, 1);
            for (unsigned char j=0; j<clusterCount; ++j)
            {
                const unsigned char one = 1;
                const unsigned char zero = 0;
                os.write( (const char*)&one, 1);
                os.write( (const char*)&zero, 1);
                --leftToProcess;
            }
        }

        if( leftToProcess != 0 )
        {
            cerr << "leftToProcess=" << leftToProcess << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Error encoding clocs" ) );
        }

    }

    void FastqTile::writeControlFile() const
    {
        clog << "Writing control file" << endl;

        ofstream os( controlFilename_.c_str() );
        if( !os.good() )
        {
            cerr << "Can't write to " << controlFilename_ << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
        assert(sizeof(unsigned int) == 4);
        const unsigned int header1 = 0, header2 = 2;
        const unsigned int header3 = expectedReadCount_;
        os.write( (const char*)&header1, 4);
        os.write( (const char*)&header2, 4);
        os.write( (const char*)&header3, 4);
        for (unsigned int i=0; i<expectedReadCount_; ++i)
        {
            const unsigned int zero = 0;
            os.write( (const char*)&zero, 2);
        }
    }

} // namespace io
} // namespace eagle
