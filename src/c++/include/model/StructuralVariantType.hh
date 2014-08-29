/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component that deals with variant events.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MODEL_STRUCTURAL_VARIANT_TYPE_HH
#define EAGLE_MODEL_STRUCTURAL_VARIANT_TYPE_HH

#include <string>
#include <bitset>

namespace eagle
{
namespace model
{
namespace variant
{


typedef std::bitset<7> Type;
const Type Undefined(0);
const Type SNP(0x1);
const Type INS(0x2);
const Type DEL(0x4);
const Type InDel(INS | DEL);  // 0x6
const Type DUP(0x8);
const Type INV(0x10);
const Type XOVER(0x20);
const Type Translocation(0x40);

// initialization via specialization
template <typename T>  T initialize(std::string chr, unsigned long pos, std::string ref, std::string alt, unsigned int altGtIndex=1);


} // namespace variant
} // namespace model
} // namespace eagle

#endif // EAGLE_GENOME_STRUCTURAL_VARIANT_TYPE_HH
