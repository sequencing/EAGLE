/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'fastaDump'
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MAIN_FASTA_DUMPER_OPTIONS_HH
#define EAGLE_MAIN_FASTA_DUMPER_OPTIONS_HH

#include <string>
#include <vector>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class FastaDumperOptions : public eagle::common::Options
{
public:
    enum Modes
    {
         UNDEFINED,
         SAFE_MODE,
         WHOLE_DIR
    };
    FastaDumperOptions();
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       fastaDump <fasta1.fa> [<fasta2.fa> [... <fastaN.fa>]]  [options]\n")
                                          + std::string("Or:\n")
                                          + std::string("       fastaDump <fastaDir>  [options]");}
    void postProcess(boost::program_options::variables_map &vm);

public:
    std::vector< boost::filesystem::path > fastaFiles;
    std::string position;
    unsigned long size;

    Modes mode;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_FASTA_DUMPER_OPTIONS_HH
