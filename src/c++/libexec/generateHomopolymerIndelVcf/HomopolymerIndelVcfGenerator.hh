/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_HOMOPOLYMER_INDEL_VCF_GENERATOR_HH
#define EAGLE_MAIN_HOMOPOLYMER_INDEL_VCF_GENERATOR_HH

#include <vector>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "genome/Reference.hh"
#include "io/BamParserFilter.hh"
#include "HomopolymerIndelVcfGeneratorOptions.hh"


namespace eagle
{
namespace main
{


class HomopolymerIndelVcfGenerator
{
public:
    HomopolymerIndelVcfGenerator (const HomopolymerIndelVcfGeneratorOptions &options);
    void run();

private:
    const HomopolymerIndelVcfGeneratorOptions &options_;
    boost::mt19937 randomGen_;
    std::vector< std::vector< double > > insertionProbabilities_, deletionProbabilities_;

    void processHomopolymer( const unsigned int homopolymerLength, const std::string &contigName, const unsigned long startPos, const char base );
};


} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_HOMOPOLYMER_INDEL_VCF_GENERATOR_HH
