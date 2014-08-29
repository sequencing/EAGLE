/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <fstream>
#include <boost/iostreams/device/file.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "genome/SharedFastaReference.hh"
#include "genome/BamMetadata.hh"
#include "io/BamIndexer.hh"
#include "model/Nucleotides.hh"
#include "HomopolymerIndelVcfGenerator.hh"


using namespace std;
//#define DEBUG


namespace eagle
{
namespace main
{

HomopolymerIndelVcfGenerator::HomopolymerIndelVcfGenerator (const HomopolymerIndelVcfGeneratorOptions &options )
    : options_( options )
{
    // Parse table file
    io::DsvReader tsvReader( options.indelProbabilitiesFile );
    vector<string> tokens;
    while ( tsvReader.getNextLineFields<'\t'>(tokens) )
    {
        if ( tokens.size() == 0 ) { continue; }
        assert ( tokens.size() == 4 && "There should be 4 entries per line" );
        unsigned int homopolymerLength = boost::lexical_cast<unsigned int,std::string>( tokens[0] );
        unsigned int indelLength       = boost::lexical_cast<unsigned int,std::string>( tokens[2] );
        double proba          = boost::lexical_cast<double,std::string>( tokens[3] );
        if (tokens[1] == "ins")
        {
            if (homopolymerLength >= insertionProbabilities_.size())
            {
                insertionProbabilities_.resize( homopolymerLength+1 );
            }
            if (indelLength >= insertionProbabilities_[homopolymerLength].size())
            
            insertionProbabilities_[homopolymerLength].resize ( indelLength+1 );
            insertionProbabilities_[homopolymerLength][indelLength] = proba;
        }
        else if (tokens[1] == "del")
        {
            if (homopolymerLength >= deletionProbabilities_.size())
            {
                deletionProbabilities_.resize( homopolymerLength+1 );
            }
            if (indelLength >= deletionProbabilities_[homopolymerLength].size())
            
            deletionProbabilities_[homopolymerLength].resize ( indelLength+1 );
            deletionProbabilities_[homopolymerLength][indelLength] = proba;
        }
        else
        {
            EAGLE_ERROR("Error while reading table: column 2 should be ins or del");
        }
    }
}

void HomopolymerIndelVcfGenerator::run()
{
//        char base = genome::SharedFastaReference::get()->get( globalPos, i, overlapContigBoundary );
    genome::SharedFastaReference::init( options_.referenceGenome );
    const vector<string> contigNames = genome::SharedFastaReference::get()->allContigNames();

    unsigned long globalPos = 0;
    BOOST_FOREACH( const string& contigName, contigNames )
    {
        unsigned long contigLength = genome::SharedFastaReference::get()->getContigLength(contigName);
        char lastBase = 0;
        unsigned int homopolymerLength = 0;
        for (unsigned long i=1; i<=contigLength; ++i, ++globalPos)
        {
            bool overlapContigBoundary;
            char base = genome::SharedFastaReference::get()->get( globalPos, 0, overlapContigBoundary );
            base = toupper( base );
            if (base == lastBase && base != 'N')
            {
                ++homopolymerLength;
            }
            else
            {
                processHomopolymer( homopolymerLength, /*globalPos - homopolymerLength,*/ contigName, i - homopolymerLength, base );
                lastBase = base;
                homopolymerLength = 1;
            }
        }
    }
}

void HomopolymerIndelVcfGenerator::processHomopolymer( const unsigned int homopolymerLength, const string &contigName, const unsigned long startPos, const char base )
{
    if (homopolymerLength < 5) { return; }
//    cout << contigName << "\t" << startPos << "\t" << homopolymerLength << endl;
    static int variantNum = 0;

    double random = randomGen_() / (double) randomGen_.max();
    int insertionTableIndex = min<int>(homopolymerLength,insertionProbabilities_.size()-1);
    for (unsigned int i=1; i<insertionProbabilities_[insertionTableIndex].size(); ++i)
    {
        const double prob = insertionProbabilities_[insertionTableIndex][i];
        random -= prob;
        if (random < 0)
        {
//            cout << " Insertion of length " << i << endl;
            cout << contigName << "\t" << startPos << "\thomo" << contigName << "Indel" << ++variantNum << "\t" << base << "\t";
            for (unsigned int j=0; j<=i; ++j)
            {
                cout << base;
            }
            cout << "\t.\tPASS\tDB\tGT\t";

            double random2 = randomGen_() / (double) randomGen_.max();
            if (random2 < 0.45)
            {
                cout << "1/0";
            }
            else if (random2 < 0.9)
            {
                cout << "0/1";
            }
            else
            {
                cout << "1/1";
            }
            cout << endl;
            return;
        }
    }

    int deletionTableIndex = min<int>(homopolymerLength,deletionProbabilities_.size()-1);
    for (unsigned int i=1; i<deletionProbabilities_[deletionTableIndex].size(); ++i)
    {
        const double prob = deletionProbabilities_[deletionTableIndex][i];
        random -= prob;
        if (random < 0)
        {
//            cout << " Deletion of length " << i << endl;
            cout << contigName << "\t" << startPos << "\thomo" << contigName << "Indel" << ++variantNum << "\t";
            for (unsigned int j=0; j<=i; ++j)
            {
                cout << base;
            }
            cout << "\t" << base << "\t.\tPASS\tDB\tGT\t";

            double random2 = randomGen_() / (double) randomGen_.max();
            if (random2 < 0.45)
            {
                cout << "1/0";
            }
            else if (random2 < 0.9)
            {
                cout << "0/1";
            }
            else
            {
                cout << "1/1";
            }
            cout << endl;
            return;
        }
    }
}


} // namespace main
} // namespace eagle
