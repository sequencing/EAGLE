/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "genome/SharedFastaReference.hh"


namespace eagle
{
namespace genome
{

PreferredFastaReader* SharedFastaReference::fastaReference_ = 0;

std::vector< boost::filesystem::path >  SharedFastaReference::sampleGenomeDir_;
std::vector< PreferredFastaReader* >    SharedFastaReference::fastaReferenceArray_;
std::map< std::string, unsigned int >   SharedFastaReference::dictionary_;


} // namespace genome
} // namespace eagle
