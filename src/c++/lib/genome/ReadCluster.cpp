/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "model/Nucleotides.hh"
#include "genome/SharedFastaReference.hh"
#include "genome/ReadCluster.hh"


using namespace std;


namespace eagle
{
namespace genome
{


ReadClusterSharedData::ReadClusterSharedData( const unsigned int clusterLength, const eagle::io::RunInfo &runInfo, const boost::filesystem::path& sampleGenomeDir, const vector<boost::filesystem::path>& qualityTableFiles, const boost::filesystem::path& mismatchTableFile, const boost::filesystem::path& homopolymerIndelTableFile, const boost::filesystem::path& motifQualityDropTableFile, const boost::filesystem::path& qqTableFile, const unsigned int userRandomSeed, const std::vector< std::string >& errorModelOptions )
    : clusterLength_  ( clusterLength )
    , runInfo_        ( runInfo )
    , errorModel_     ( qualityTableFiles, mismatchTableFile, homopolymerIndelTableFile, motifQualityDropTableFile, qqTableFile, errorModelOptions )
    , userRandomSeed_ ( userRandomSeed )
{
    SharedFastaReference::init( sampleGenomeDir );

    FragmentStructureType2Generic fragmentStructure1;
    unsigned int readNum = 1;
    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo.reads)
    {
        if (rd.isIndex)
        {
            fragmentStructure1.addBarcode();
        }
        else
        {
            fragmentStructure1.addRead( readNum++ );
        }
    }
    multiplexedFragmentStructures.push_back( fragmentStructure1 );

    FragmentStructureType2GenericReverse fragmentStructure2;
    readNum = 1;
    BOOST_FOREACH(const eagle::io::ReadDescription &rd, runInfo.reads)
    {
        if (rd.isIndex)
        {
            fragmentStructure2.addBarcode();
        }
        else
        {
            fragmentStructure2.addRead( readNum++ );
        }
    }
    multiplexedFragmentStructures.push_back( fragmentStructure2 );

}


ReadClusterWithErrors::ReadClusterWithErrors(ReadClusterSharedData &sharedData, const EnrichedFragment& eFragment, boost::shared_ptr<boost::mt19937> randomGen)
    : sharedData_(sharedData)
    , randomGen_(randomGen)
    , eFragment_(eFragment)
    , lazyEvaluationDone_( false )
    , buf_(sharedData_.clusterLength_,' ')
                                            /*
    , cigar_()
    , read1PosInBuf_(0)
    , read1Length_(0)
    , read2PosInBuf_(0)
    , read2Length_(0)
                                            */
{
}

const char *ReadClusterWithErrors::getBclCluster( bool generateCigar, const bool dropLastBase )
{
    unsigned int posInCluster = 0;
    unsigned int readNum = 0;
    ClusterErrorModelContext clusterErrorModelContext;

    lazyEvaluationDone_ = true;

    BOOST_FOREACH(const eagle::io::ReadDescription &rd, sharedData_.runInfo_.reads)
    {
        clusterErrorModelContext.initialiseForNewRead();

        unsigned int lastCigarOp = 0;
        unsigned int lastCigarOpCount = 0;
        if (generateCigar)
        {
            cigar_.resize( readNum+1 );
            usedDnaLength_.resize( readNum+1 );
        }

//        unsigned int readLength = rd.lastCycle - rd.firstCycle + 1; //eFragment_.getReadLength( readNum );
        for (unsigned int cycle=rd.firstCycle, posToRead=0; cycle<=rd.lastCycle; ++cycle, ++posToRead)
        {
            char base = eFragment_.getBase( readNum, posToRead );
            unsigned int quality, randomErrorType;
            char bclBase;
            unsigned int newCigarOp = 0;
            sharedData_.errorModel_.getQualityAndRandomError( *randomGen_, cycle, base, quality, randomErrorType, bclBase, clusterErrorModelContext );

            switch (randomErrorType)
            {
            case ErrorModel::NoError:
            case ErrorModel::BaseSubstitution:
            case ErrorModel::BaseInsertion:
                {
                    char bclBaseAndQ = bclBase | (quality<<2);
                    assert( posInCluster < buf_.size() );
                    buf_[posInCluster++] = bclBaseAndQ;
                    //            ++stats[cycle][bclBase & 3];

                    if (randomErrorType == ErrorModel::BaseInsertion)
                    {
                        --posToRead; // next cycle will read from the same pos
                        newCigarOp = 1; // 1='I'
                    }
                    else
                    {
                        newCigarOp = 0; // 0='M'
                    }
                }
                break;
            case ErrorModel::BaseDeletion:
                {
                    --cycle; // repeat same cycle but will read next pos
                    newCigarOp = 2; // 2='D'
                }
                break;
            default:
                assert( false );
            }

            if (generateCigar && (!dropLastBase || cycle<rd.lastCycle))
            {
                if (newCigarOp == lastCigarOp)
                {
                    ++lastCigarOpCount;
                }
                else
                {
                    if (lastCigarOpCount > 0)
                    {
                        cigar_[readNum].push_back( lastCigarOpCount << 4 | lastCigarOp );
                    }
                    lastCigarOp = newCigarOp;
                    lastCigarOpCount = 1;
                }

                switch (newCigarOp)
                {
                case 0: // 'M'
                case 2: // 'D'
                    ++usedDnaLength_[readNum];
                    break;
                case 1: // 'I'
                    break;
                default:
                    assert( false && "shouldn't reach here");
                }
            }
        }

        // Pushes last CIGAR entry
        if (generateCigar)
        {
            cigar_[readNum].push_back( lastCigarOpCount << 4 | lastCigarOp );
        }

        readNum++;
    }
    assert( posInCluster == buf_.size() );
    return buf_.c_str();
}


const std::vector<unsigned int>& ReadClusterWithErrors::getCigar( unsigned int readNum, const bool dropLastBase )
{
    if (!lazyEvaluationDone_)
    {
        getBclCluster( true, dropLastBase );
    }
    assert (!cigar_.empty() && "CIGAR string didn't get generated as it should have");
    assert (!cigar_[readNum].empty() && "CIGAR string didn't get generated as it should have");
    return cigar_[readNum];
}

unsigned int ReadClusterWithErrors::getUsedDnaLength( unsigned int readNum, const bool dropLastBase )
{
    if (!lazyEvaluationDone_)
    {
        getBclCluster( true, dropLastBase );
    }
    assert (!usedDnaLength_.empty() && "DNA length didn't get calculated as it should have");
    return usedDnaLength_[readNum];
}

const string ReadClusterWithErrors::getNucleotideOrQualitySequenceForRead( unsigned int readNum, bool getNucleotides, bool revComp, const bool dropLastBase )
{
    eagle::model::IUPAC converter;

    if (!lazyEvaluationDone_)
    {
        getBclCluster( false, dropLastBase );
    }

    string result;
//    BOOST_FOREACH(const eagle::io::ReadDescription &rd, sharedData_.runInfo_.reads) //TODO: see if we can optimise this by using the fragmentStructures directly
    const eagle::io::ReadDescription &rd = sharedData_.runInfo_.reads[readNum];
    if (!rd.isIndex)
    {
        unsigned int readLength = rd.lastCycle - rd.firstCycle + 1;
        unsigned int lastCycle = rd.lastCycle;
        if (dropLastBase)
        {
            --lastCycle;
            --readLength;
        }
        unsigned int pos = revComp?readLength-1:0;
        result.resize( readLength );
        for (unsigned int cycle=rd.firstCycle; cycle<=lastCycle; ++cycle)
        {
            unsigned char bclBase = buf_[cycle-1];
//            clog << cycle << "," << (unsigned int)bclBase << endl;
            if (getNucleotides)
            {
                result[pos] = (bclBase>>2)?converter.normFromBcl( revComp?~bclBase:bclBase ):'N'; // Even though converter.normFromBcl handles quality==bclBase>>2==0, it doesn't handle bclBase==~0, hence the test
            }
            else
            {
                result[pos] = bclBase/4 + 33;
            }
//                    clog << result[pos] << endl;
            if (revComp) { --pos; } else { ++pos; }
        }
//                clog << result << endl;
    }
    return result.empty()?"*":result;
}


ReadClusterFactory::ReadClusterFactory( const eagle::io::RunInfo &runInfo, const boost::filesystem::path& sampleGenomeDir, const vector<boost::filesystem::path>& qualityTableFiles, const boost::filesystem::path& mismatchTableFile, const boost::filesystem::path& homopolymerIndelTableFile, const boost::filesystem::path& motifQualityDropTableFile, const boost::filesystem::path& qqTableFile, const unsigned int userRandomSeed, const std::vector< std::string >& errorModelOptions )
    : runInfo_(runInfo)
    , sharedData_( runInfo_.getClusterLength(), runInfo_, sampleGenomeDir, qualityTableFiles, mismatchTableFile, homopolymerIndelTableFile, motifQualityDropTableFile, qqTableFile, userRandomSeed, errorModelOptions )
{
}

ReadClusterWithErrors ReadClusterFactory::getReadClusterWithErrors( const eagle::model::Fragment &f ) {
    unsigned int seed = ((f.fragmentNum_+1) * sharedData_.userRandomSeed_) ^ 0x9e3779b9; // Avoid too many zeroes in the non-random initial state, as they delay real randomness
    boost::shared_ptr< boost::mt19937 > randomGen( new boost::mt19937(seed) );
    randomGen->discard( 10 ); // Generates several random values to make the state more random

    ReadClusterWithErrors cluster( sharedData_, EnrichedFragment( f, sharedData_.multiplexedFragmentStructures, randomGen->operator()()%2 ), randomGen );
    return cluster;
}


} // namespace genome
} // namespace eagle
