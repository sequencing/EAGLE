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


FastqTile::FastqTile( const unsigned long long expectedReadCount, const unsigned int clusterLength, const string &read1FastqFilename, const string &read2FastqFilename, const RunInfo &runInfo, const int lane, const unsigned int tileId, const bool verbose )
        : expectedReadCount_ ( expectedReadCount )
        , clusterLength_     ( clusterLength )
        , read1FastqFilename_( read1FastqFilename )
        , read2FastqFilename_( read2FastqFilename )
        , read1FastqFile_    ( read1FastqFilename.c_str() )
        , read2FastqFile_    ( read2FastqFilename.c_str() )
        , infoFile_          ( (read1FastqFilename_ + ".info").c_str() )
        , readNamePrefix_    ( (boost::format("@EAGLE:%s:%s:%d:%04g:") % runInfo.runNumber % runInfo.flowcell % lane % tileId).str() ) // @HSQ1004:200:C0M7LACXX:1:1101:1187:2070 1:N:0:1
        , totalReadCount_    ( 0 )
        , passedFilterReadCount_( 0 )
    {
        if (verbose)
        {
            clog << (boost::format("Creating new Fastq tile as (%s, %s), expecting %d reads") % read1FastqFilename % read2FastqFilename % expectedReadCount_).str() << endl;
        }

        if( !read1FastqFile_.good() )
        {
            cerr << "Can't write to " << read1FastqFilename << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }

        if( read2FastqFilename != "" && !read2FastqFile_.good() )
        {
            cerr << "Can't write to " << read2FastqFilename << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }

        if( !infoFile_.good() )
        {
            cerr << "Can't write to " << read1FastqFilename << ".info" << endl;
            BOOST_THROW_EXCEPTION( eagle::common::IoException( errno, "Cannot create file" ) );
        }
   }

    void FastqTile::addCluster( const string &read1Nucleotides, const string &read1Qualities, const string &read2Nucleotides, const string &read2Qualities, const bool isPassingFilter, const unsigned long coordX, const unsigned long coordY )
    {
        string readNamePart1 = readNamePrefix_ + (boost::format("%d:%d") % coordX % coordY).str();

        string read1Buf = readNamePart1 + " 1:" + (isPassingFilter?"N":"Y") + ":0:1\n" + read1Nucleotides + "\n+\n" + read1Qualities + "\n";
        read1FastqFile_.write(read1Buf.c_str(), read1Buf.size());

        string read2Buf = readNamePart1 + " 2:" + (isPassingFilter?"N":"Y") + ":0:1\n" + read2Nucleotides + "\n+\n" + read2Qualities + "\n";
        read2FastqFile_.write(read2Buf.c_str(), read2Buf.size());

        totalReadCount_++;
        if (isPassingFilter)
        {
            passedFilterReadCount_++;
        }
    }

    void FastqTile::finaliseAndWriteInfo()
    {
        read1FastqFile_.close();
        read2FastqFile_.close();

        string info = (boost::format("TotalReadsRaw\t%lld\nTotalReadsPF\t%lld") % totalReadCount_ % passedFilterReadCount_).str();
        infoFile_.write( info.c_str(), info.length() );

        infoFile_.close();
    }


} // namespace io
} // namespace eagle
