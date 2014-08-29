/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description A preprocessor-based Logger + time display.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_COMMON_LOGGER_HH
#define EAGLE_COMMON_LOGGER_HH

#include <iostream>
#include <string>

#include "config.h"
#include "common/Exceptions.hh"

#ifdef EAGLE_DEBUG_MODE
#define EAGLE_DEBUG(x,y) std::clog << "* Debug *: " << std::string((x),' ') << y << std::endl
#define EAGLE_DEBUG_IF(x,y,z) if((x)) {std::clog << "* Debug *: " << std::string((y),' ') << z << std::endl;}
#undef  EAGLE_SILENT_MODE
#else
#define EAGLE_DEBUG(x,y)
#define EAGLE_DEBUG_IF(x,y,z)
#endif

#ifndef EAGLE_SILENT_MODE
#define EAGLE_PRINT(x) std::clog << x << std::endl
#else
#define EAGLE_PRINT(x)
#endif

#define EAGLE_WARNING(x)                  std::cerr << "** Warning:" << __FILE__ << ":" << __LINE__ << ": **" << std::endl << "** Warning **: " << x << std::endl
#define EAGLE_WARNING_IF(x,y)      if(x) {std::cerr << "** Warning:" << __FILE__ << ":" << __LINE__ << ": **" << std::endl << "** Warning **: " << y << std::endl;}
#define EAGLE_WARNING_CONT(x)             std::cerr << "** Warning **: " << x << std::endl
#define EAGLE_WARNING_CONT_IF(x,y) if(x) {std::cerr << "** Warning **: " << y << std::endl;}

/* Always prefer BOOST_THROW_EXCEPTION() directly, as it allows you to throw the appropriate EagleException(). */
/* This is just a lazy default way to terminate the execution */
#define EAGLE_ERROR(M) BOOST_THROW_EXCEPTION(eagle::common::EagleException(0, std::string("*** ERROR *** :\n") + (M) ))


namespace eagle
{
namespace common
{

std::string displayTime(size_t time);
std::string displayTime(size_t time, size_t &acc);

} // namespace common
} // namespace eagle

#endif // EAGLE_COMMON_LOGGER_HH
