/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Declaration of the common exception mechanism.
 **
 ** All exceptions must carry the same data (independently of the
 ** exception type) to homogenize the reporting and processing of
 ** errors.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_COMMON_EXCEPTIONS_HH
#define EAGLE_COMMON_EXCEPTIONS_HH

#include <string>
#include <stdexcept>
#include <ios>
#include <boost/cerrno.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>


namespace eagle
{
namespace common
{

/**
 ** \brief Virtual base class to all the exception classes in EAGLE.
 **
 ** Use BOOST_THROW_EXCEPTION to get the context info (file, function, line)
 ** at the throw site.
 **/
class ExceptionData : public boost::exception
{
public:
    ExceptionData(int errorNumber=0, const std::string &message="");
    ExceptionData(const ExceptionData &e) : boost::exception(e), errorNumber_(e.errorNumber_), message_(e.message_) {}
    virtual ~ExceptionData() throw ()
    {
    }
    int getErrorNumber() const
    {
        return errorNumber_;
    }
    std::string getMessage() const
    {
        return message_;
    }
    std::string getContext() const;
private:
    const int errorNumber_;
    const std::string message_;
    ExceptionData &operator=(const ExceptionData &);
};

class EagleException: public std::exception, public ExceptionData
{
public:
    EagleException(int errorNumber, const std::string &message) : ExceptionData(errorNumber, message) {}
    EagleException(const EagleException &e) : std::exception(e), ExceptionData(e) {}
private:
    EagleException &operator=(const EagleException &);
};

/**
 ** \brief Exception thrown when a file cannot be partially/totally accessed.
 **
 **/
class IoException: public std::ios_base::failure, public ExceptionData
{
public:
    IoException(int errorNumber, const std::string &message);
};

/**
 * \brief Exception thrown when there is insufficient resources to perform an operation. For example
 *        if the adjusting the soft ulimit fails due to a set hard limit
 */
class ResourceException: public std::exception, public ExceptionData
{
public:
    ResourceException(int errorNumber, const std::string &message);
};

/**
 * \brief Same as bad_alloc but with a message
 */
class MemoryException: public std::bad_alloc, public ExceptionData
{
public:
    MemoryException(const std::string &message);
};


/**
 ** \brief Exception thrown when a file contains invalid data.
 **
 **/
class CorruptedFileException: public std::logic_error, public ExceptionData
{
public:
    CorruptedFileException(const std::string &type, const std::string &message);
};

/**
 ** \brief Exception thrown when the client supplied and unsupported version number.
 **
 **/
class UnsupportedVersionException: public std::logic_error, public ExceptionData
{
public:
    UnsupportedVersionException(const std::string &message);
};

/**
 ** \brief Exception thrown when the client supplied an invalid parameter.
 **
 **/
class InvalidParameterException: public std::logic_error, public ExceptionData
{
public:
    InvalidParameterException(const std::string &message);
};

/**
 ** \brief Exception thrown when an invalid command line option was detected.
 **
 **/
class InvalidOptionException: public std::logic_error, public ExceptionData
{
public:
    InvalidOptionException(const std::string &message);
};

/**
 ** \brief Exception thrown when a method invocation violates the pre-conditions.
 **
 **/
class PreConditionException: public std::logic_error, public ExceptionData
{
public:
    PreConditionException(const std::string &message);
};

/**
 ** \brief Exception thrown when a method invocation violates the post-conditions.
 **
 **/
class PostConditionException: public std::logic_error, public ExceptionData
{
public:
    PostConditionException(const std::string &message);
};

/**
 ** \brief Exception thrown when a method invocation generates an out-of-limits situation.
 **
 **/
class OutOfLimitsException: public std::logic_error, public ExceptionData
{
public:
    OutOfLimitsException(const std::string &message);
};


} // namespace common
} // namespace eagle

#endif // EAGLE_COMMON_EXCEPTIONS_HH
