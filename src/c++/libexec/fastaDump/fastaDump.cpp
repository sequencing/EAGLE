/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 **/


#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>


#include "FastaDumper.hh"
#include "FastaDumperOptions.hh"


static void fastaDumperLauncher(const eagle::main::FastaDumperOptions &options)
{
    if (eagle::main::FastaDumperOptions::WHOLE_DIR == options.mode)
    {
        std::clog << (boost::format("Looking for a reference genome in %s ...") % options.fastaFiles[0] ).str() << std::endl;
        eagle::main::FastaDumper fastaDumper(
            options.fastaFiles[0],
            options.position,
            options.size
            );
        fastaDumper.run();
    } else {
        eagle::main::FastaDumper fastaDumper(
            options.fastaFiles,
            options.position,
            options.size
            );
        fastaDumper.run();
    }
}

int main(int argc, char *argv[])
{
    eagle::common::run(fastaDumperLauncher, argc, argv);
}
