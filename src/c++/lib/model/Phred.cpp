/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Pre-computed Phred quality scores
 **
 ** \author Lilian Janin
 **/

#include "model/Phred.hh"

namespace eagle
{
namespace model
{

std::vector< double > Phred::qualityToProbability_;


} // namespace model
} // namespace eagle
