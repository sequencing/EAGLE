/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Implementation of the common exception mechanism.
 **
 ** \author Mauricio Varea
 **/

#include <cstring>
#include <cerrno>
#include <boost/date_time.hpp>

#include "common/Exceptions.hh"

namespace eagle
{
namespace common
{

ExceptionData::ExceptionData(int errorNumber, const std::string &message) : boost::exception(),
            errorNumber_(errorNumber), message_(message)
{
}

std::string ExceptionData::getContext() const
{
    const std::string now = boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());
    return now + ": " + std::string(errorNumber_?strerror(errorNumber_):"") + ": " + boost::diagnostic_information(*this);
}

IoException::IoException(int errorNumber, const std::string &message)
    : std::ios_base::failure(message)
    , ExceptionData(errorNumber, message)
{
}

ResourceException::ResourceException(int errorNumber, const std::string &message)
    : ExceptionData(errorNumber, message)
{
}

MemoryException::MemoryException(const std::string &message)
    : std::bad_alloc(),
      ExceptionData(ENOMEM, message)
{
}

CorruptedFileException::CorruptedFileException(const std::string &type, const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, "Corrupt '" + type + "' file: " + message)
{
}

UnsupportedVersionException::UnsupportedVersionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidParameterException::InvalidParameterException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidOptionException::InvalidOptionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PreConditionException::PreConditionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PostConditionException::PostConditionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

OutOfLimitsException::OutOfLimitsException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}


} // namespace common
} // namespace eagle
