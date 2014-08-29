/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **/

#include <iostream>
#include <string>
#include <vector>
#include <boost/assign.hpp>

using namespace std;
using boost::assign::list_of;

#include "Helpers.hh"

#include "RegistryName.hh"
#include "testVcf.hh"

using eagle::io::VcfVariant;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestVcf, registryName("Vcf"));


void TestVcf::setUp()
{
}

void TestVcf::tearDown()
{
}


void TestVcf::testVcfParserSnp()
{
    VcfVariant variant;

    // Keep this as reference
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr", "pos", "id", "ref", "alt", "qual=.", "filter=PASS", "info=" ), eagle::common::CorruptedFileException );

    // Typical SNP
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C" ) );

    // Minimal SNP
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "1", "0", ".", "G", "T" ) );

    // Too minimal: missing fields
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "" , "0", ".", "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "" , ".", "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", "" , "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", ".", "" , "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", ".", "G", ""  ), eagle::common::CorruptedFileException );

    // Too minimal: missing fields filled with a space
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  " ", "0", ".", "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", " ", ".", "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", " ", "G", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", ".", " ", "T" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant(  "1", "0", ".", "G", " " ), eagle::common::CorruptedFileException );

    // Invalid position
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "-1", "VariantId1", "A", "C" ), eagle::common::CorruptedFileException );

    // Invalid combination of letters
    //    (don't we want to support this for centromeres? SAGE-27)
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "A" ), eagle::common::CorruptedFileException );

    // Should we throw on "strange" nucleotide letters?
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "B", "D" ) );

    // Should we throw on these special "N" SNPs?
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "N" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "N", "A" ) );

    // We should throw on non-existing nucleotide letters
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "X", "Z" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "@", "&" ), eagle::common::CorruptedFileException );

    // Do we support the '.' "any base" notation? (is it part of VCF standard?)
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", ".", "A" ) );

    // We allow lowercase bases
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "a", "c" ) );
}

void TestVcf::testVcfParserInsertions()
{
    VcfVariant variant;

    // Typical insertion
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "AACGT" ) );

    // This notation is NOT SUPPORTED YET
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<INS>", ".", "PASS", "?" ), eagle::common::CorruptedFileException );

    // This notation is NOT SUPPORTED YET (contigs from assembly files)
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "A[<ctg1>:1[", ".", "PASS", "?" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId2", "C", "]<ctg1>:1000]C", ".", "PASS", "?" ), eagle::common::CorruptedFileException );
}

void TestVcf::testVcfParserDeletions()
{
    VcfVariant variant;

    // Typical deletion
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "AACGT", "A" ) );

    // This notation is NOT SUPPORTED YET
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DEL>", ".", "PASS", "END=23456;SVLEN=11111" ), eagle::common::CorruptedFileException );

    // Check for invalid redundant information (END-START != SVLEN)
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DEL>", ".", "PASS", "END=23456;SVLEN=99999" ), eagle::common::CorruptedFileException );
}

void TestVcf::testVcfParserIndels()
{
    VcfVariant variant;

    // Typical indels
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "ACCCC", "AGGGG" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "AC", "AGGGGTTT" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "ACCCCTTT", "AG" ) );

    // Should we throw on an indel that should be described as a SNP? (indel leads to more conflicts)
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "AC", "AG" ), eagle::common::CorruptedFileException );
}

void TestVcf::testVcfParserTranslocations()
{
    VcfVariant variant;

    // Typical translocations
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "A[chr2:54321[" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "C", "C]chr2:54321]" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "G", "[chr2:54321[G" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "T", "]chr2:54321]T" ) );

    // Typical translocations but common base not matching
    //    (at the moment, this is treated as *Translocation with SNP*)
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C[chr2:54321[" ), eagle::common::CorruptedFileException );
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "C", "G]chr2:54321]" ), eagle::common::CorruptedFileException );
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "G", "[chr2:54321[T" ), eagle::common::CorruptedFileException );
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "T", "]chr2:54321]A" ), eagle::common::CorruptedFileException );

    // Translocation with insertions are NOT supported yet
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "ACGTCGT[chr2:54321[" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "]chr2:54321]CGTCGTA" ) );

    // We need to check that the REF base appears at the correct place in the ALT field (these tests produce a WARNING at the moment)
/*
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "CCGTCGA[chr2:54321[" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "]chr2:54321]AGTCGTC" ), eagle::common::CorruptedFileException );
*/
}

void TestVcf::testVcfParserInversions()
{
    VcfVariant variant;

    // Typical inversion is NOT supported
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<INV>" ), eagle::common::CorruptedFileException );

    // Translocation-like notation
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "A]chr1:12400]" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12346", "VariantId2", "C", "[chr1:12401[C" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12401", "VariantId3", "G", "[chr1:12346[G" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12400", "VariantId4", "T", "T]chr1:12345]" ) );
}

void TestVcf::testVcfParserDuplications()
{
    VcfVariant variant;

    // Typical duplication is NOT supported, and will probably never be as it's an imprecise variant
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DUP>" ), eagle::common::CorruptedFileException );

    // Tandem duplication is NOT supported yet, but should be
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DUP:TANDEM>", ".", "PASS", "END=23456;SVLEN=11111" ), eagle::common::CorruptedFileException );

    // (uncomment when supported) Tandem duplication is NOT supported yet - but when we start supporting it, we should allow only 1 field of the redundant END/SVLEN sub-fields, and check for incompatible values
/*
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DUP:TANDEM>", ".", "PASS", "END=23456" ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DUP:TANDEM>", ".", "PASS", "SVLEN=11111" ) );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DUP:TANDEM>", ".", "PASS", "END=23456;SVLEN=999" ), eagle::common::CorruptedFileException );
*/
}

void TestVcf::testVcfParserQualityField()
{
    VcfVariant variant;

    // Quality should be "." or a positive number
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", "." ) );
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", "12" ) );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", "-1" ), eagle::common::CorruptedFileException );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", "ABC" ), eagle::common::CorruptedFileException );
}

void TestVcf::testVcfParserFilterField()
{
    VcfVariant variant;

    // Filter should be "PASS"?
    CPPUNIT_ASSERT_NO_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", ".", "PASS" ) );
    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "C", ".", " ; " ), eagle::common::CorruptedFileException );
}

void TestVcf::testVcfParserInfoField()
{
    VcfVariant variant;

    // We do NOT want to support imprecise variants (there is still quite some work to do in parsing the INFO field, yet)
//    CPPUNIT_ASSERT_THROW( variant = VcfVariant( "chr1", "12345", "VariantId1", "A", "<DEL>", ".", "PASS", "IMPRECISE" ), eagle::common::CorruptedFileException );
}

