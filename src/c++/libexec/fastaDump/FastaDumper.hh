/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to dump a reference.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MAIN_FASTA_DUMPER_HH
#define EAGLE_MAIN_FASTA_DUMPER_HH

#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

#include "common/Logger.hh"
#include "model/Struct.hh"
#include "genome/Reference.hh"


// TODO: get this from  eagle::io::Fasta::ContigWidth  or something like that
#define CONTIG_WIDTH 70

//


namespace eagle
{
namespace main
{


class MetaLocus : public eagle::model::Locus
{
public:
    MetaLocus(unsigned long pos) : eagle::model::Locus(" ", pos), global(true) {}  // FASTA cannot define a contig name as ' ', so
                                                                                   //    we are using this to indicate it's a global position
    MetaLocus(std::string chr, unsigned long pos) : eagle::model::Locus(chr,pos), global(false) {}
    MetaLocus(std::string chr) : eagle::model::Locus(chr), global(false) {}
    MetaLocus() {}
    bool global;
};

class FastaDumper
{
public:
    FastaDumper(  // SAFE_MODE
        const std::vector<boost::filesystem::path> inputFiles,
        const std::string position,
        const unsigned long size
        )
    : reference_( inputFiles )
    , location_( seek(position) )
    , size_(size)
    {}
    FastaDumper(  // WHOLE_DIR
        const boost::filesystem::path inputDir,
        const std::string position,
        const unsigned long size
        )
    : reference_( inputDir )
    , location_( seek(position) )
    , size_(size)
    {}
    void run();
private:
    void display(eagle::model::Locus location, unsigned int nameWidth=20, unsigned int posWidth=8)
    {
        std::cout << std::setw( nameWidth ) << std::right << location.chr() << ":";
        std::cout << std::setw( posWidth )  << std::left  << location.pos() << "| ";
    }
    MetaLocus seek(std::string ml);

    eagle::genome::MultiFastaReference reference_;
    const MetaLocus location_;
    const unsigned long size_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_FASTA_DUMPER_HH
