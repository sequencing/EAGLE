/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Deadlock-free Semaphore class, able to reset itself after crashes
 **
 ** Deadlock-free semaphore class based on boost's named semaphore.
 ** Automatically increase its counter after a timed_wait to avoid deadlocks
 ** And automatically decrease it back to normal when possible
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_COMMON_SEMAPHORE_HH
#define EAGLE_COMMON_SEMAPHORE_HH

#include <string>
#include <boost/interprocess/sync/named_semaphore.hpp>


namespace eagle
{
namespace common
{

class Semaphore
{
public:
    Semaphore( const std::string& name, const unsigned int count );
    void wait();
    void post();

private:
    const unsigned int resourceCount_;
    boost::interprocess::named_semaphore semaphore_;
};


} // namespace common
} // namespace eagle

#endif // EAGLE_COMMON_SEMAPHORE_HH
