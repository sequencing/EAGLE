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

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "common/Exceptions.hh"
#include "FastaDumperOptions.hh"

namespace eagle
{
namespace main
{


namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

FastaDumperOptions::FastaDumperOptions()
    : eagle::common::Options(false)
    , position("1")
    , size(0UL)
    , mode(UNDEFINED)
{
    namedOptions_.add_options()
        ("position,p",   bpo::value< std::string >(&position)->default_value(position),
                                     "Position to start dumping from. Can be either:\n"
                                     " a)        Number => a global position\n"
                                     " b) String:Number => contig name followed by\n"
                                     "                     a local position\n"
                                     " c) String        => just a contig name\n"
                                     "                     (starts from 1st position\n"
                                     "                     in that chromosome)\n")
        ("size,n",       bpo::value< unsigned long >(&size),
                                     "Amount of bases to dump\n(defaults to until-the-end behaviour)")
        ;

    unnamedOptions_.add_options()
        ("positional", bpo::value< std::vector< bfs::path > >(&fastaFiles), "list of files, or just 1 directory")
        ;
    positionalOptions_.add("positional",-1);
}

void FastaDumperOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    if (fastaFiles.empty())
    {
        throw bpo::validation_error(bpo::validation_error::at_least_one_value_required, "", "positional");
    }

    check.addPathOptions(fastaFiles,"positional");
    check.inputPathsExist();

    if ( 1 == fastaFiles.size() && bfs::is_directory(fastaFiles[0]) )
    {
        mode = WHOLE_DIR;
    } else {
        for (unsigned int i = 0; i < fastaFiles.size(); i++)
        {
            if (bfs::is_directory(fastaFiles[i]))
            {
                const boost::format message = boost::format("\n   *** FASTA file #%d has an invalid value: ***"
                                                            "\n   ***       It should point to a file, but a directory already exists with name %s ***\n")
                                                            % i % fastaFiles[i];
                BOOST_THROW_EXCEPTION(eagle::common::InvalidOptionException(message.str()));
            }
        }
        mode = SAFE_MODE;
    }

    check.inRange<unsigned long>(std::make_pair(size,"size"),1);  // 1 <= size < inf
                                                                  // (size == 0) not allowed, so internally used to represent 'until the end'
}

} //namespace main
} // namespace eagle
