/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Time display.
 **
 ** \author Mauricio Varea
 **/

#include <fstream>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Logger.hh"

namespace eagle
{
namespace common
{

std::string displayTime(size_t time)
{   // time in microSeconds
    size_t mili   =    time / 1000;
    size_t hour   =    mili / (60 * 60 * 1000);
    size_t r_hour =    mili % (60 * 60 * 1000);
    size_t min    =  r_hour / (60 * 1000);
    size_t r_min  =  r_hour % (60 * 1000);
    size_t sec    =   r_min / (1000);
    return boost::lexical_cast<std::string>(mili) + "ms"
         + ((mili < 1000) ? "" : (boost::format(" (%uh:%um:%us)") % hour % min % sec).str());
}

std::string displayTime(size_t time, size_t &acc)
{
    acc += time;
    return displayTime(time);
}


} // namespace common
} // namespace eagle
