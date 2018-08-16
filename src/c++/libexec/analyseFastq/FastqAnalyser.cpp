/**
 ** Copyright (c) 2018 Illumina, Inc.
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
//#include "genome/FastqMetadata.hh"
//#include "io/FastqIndexer.hh"
#include "model/Nucleotides.hh"
#include "FastqAnalyser.hh"


using namespace std;
//#define DEBUG


namespace eagle
{
namespace main
{

FastqAnalyser::FastqAnalyser (const FastqAnalyserOptions &options )
    : options_( options )
{
    //genome::SharedFastaReference::init( options_.referenceGenome );
}

void FastqAnalyser::run()
{
    string line_readname, line_bases, line_sep, line_quals;
    ifstream inp( options_.fastqFile.string().c_str() );
    while (getline(inp, line_readname))
    {
        getline(inp, line_bases);
        getline(inp, line_sep);
        getline(inp, line_quals);
        processRead(line_readname, line_bases, line_quals);
    }
    endOfFastq();
}

void FastqAnalyser::processRead(string line_readname, string line_bases, string line_quals)
{
    unsigned int seqLength = line_quals.size();//alignment.lSeq;

    // Convert letter [!-...] to qscore value 0..50
    for (unsigned int i=0; i<seqLength; ++i) {
        line_quals[i] -= 33;
    }

    const char *quals = line_quals.c_str();//alignment.getQual();
    int qSum = 0;
    for (unsigned int i=0; i<seqLength; ++i) {
        int qual = quals[i];
        qSum += qual;
    }
    int averageQ = qSum / seqLength;

#define MIN_QSCORE 2
#define MAX_QSCORE 50

    unsigned int templateNum = max(0, min(MAX_QSCORE, averageQ));
    if (qualityTable_.size() <= templateNum)
    {
        qualityTable_.resize( templateNum+1 );
    }
    if (qualityTable_[templateNum].size() < seqLength)
    {
        qualityTable_[templateNum].resize( seqLength );
        for (unsigned int i=0; i<seqLength; ++i)
            qualityTable_[templateNum][i].resize( MAX_QSCORE+1 );
    }
    for (unsigned int i=0; i<seqLength; ++i) {
        int qual = quals[i];
        qual = max(MIN_QSCORE, min(MAX_QSCORE, qual));
        qualityTable_[templateNum][i][qual]++;
    }

    // Count quals
    for (unsigned int i=0; i<seqLength; ++i) {
        unsigned int qual = quals[i];
        if (qual >= qualCount_.size())
        {
            qualCount_.resize( qual+1 );
        }
        qualCount_[qual]++;
    }

}

void FastqAnalyser::endOfFastq()
{
    cout << " *** Quality table ***" << endl;
    cout << "Output to QualityTable.Rx.qtable2" << endl;
    ofstream ofs( "QualityTable.Rx.qtable2" );
    ofs << "#templateNum\tcycle\tQ1:#Q1\tQ2:#Q2\t..." << endl;

    ofs << "0\t0";
    for (unsigned int templateNum = 1; templateNum < qualityTable_.size(); ++templateNum)
    {
        unsigned long long count = 0;
        //for (unsigned int cycle = 0; cycle < qualityTable_[templateNum].size(); ++cycle)
        if (qualityTable_[templateNum].size() > 1)
        {
            unsigned int cycle = 1; // any cycle would do
            for (unsigned int qscore = 0; qscore < qualityTable_[templateNum][cycle].size(); ++qscore)
            {
                count += qualityTable_[templateNum][cycle][qscore];
            }
        }

        if (count)
        {
            ofs << '\t' << templateNum << ':' << count;
        }
    }
    ofs << endl;

    for (unsigned int templateNum = 1; templateNum < qualityTable_.size(); ++templateNum)
    {
        for (unsigned int cycle = 0; cycle < qualityTable_[templateNum].size(); ++cycle)
        {
            ofs << templateNum << '\t' << (cycle + 1);
            for (unsigned int qscore = 0; qscore < qualityTable_[templateNum][cycle].size(); ++qscore)
            {
                unsigned long long count = qualityTable_[templateNum][cycle][qscore];
                if (count)
                    ofs << '\t' << qscore << ':' << count;
            }
            ofs << endl;
        }
    }
    ofs << endl;

    // output qualcount
    ofstream ofs2( "QualityTable.Rx.counts" );
    cout << "# Count per qscore: See file QualityTable.Rx.counts" << endl;
    for (unsigned int i=0; i<qualCount_.size(); ++i)
    {
        if (qualCount_[i])
        ofs2 << i << '\t' << qualCount_[i] << endl;
    }
}

} // namespace main
} // namespace eagle
