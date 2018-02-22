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

#include "common/Semaphore.hh"

#include <iostream>
#include <boost/date_time.hpp>

using namespace std;


namespace eagle
{
namespace common
{

Semaphore::Semaphore( const string& name, const unsigned int resourceCount )
    : resourceCount_( resourceCount )
    , semaphore_    ( boost::interprocess::open_or_create, name.c_str(), resourceCount )
{
}

void Semaphore::wait()
{
    try
    {
        // Ensure that the semaphore's count is not greater than expected
        unsigned int count = 0;
        while (semaphore_.try_wait())
        {
            count++;
            clog << "Semaphore level " << count << endl;
        }

        // Free up a maximum of 'count_' semaphore resources, minus the one we want to use
        if (count > resourceCount_)
        {
            clog << "Adjusting semaphore down" << endl;
            count = resourceCount_;
        }
        for (unsigned int i=1; i<count; ++i)
        {
            semaphore_.post();
        }

        // If there was no resource available, wait for one
        if (count == 0)
        {
            // We wait for X minutes, and then we decide to bypass the semaphore, in case it got locked somehow
            // This has the effect of increasing the total semaphore count, but it automatically gets readjusted by the initial count-check above
            clog << "Waiting for 5 minutes" << endl;
            boost::posix_time::ptime delay = boost::posix_time::microsec_clock::universal_time() + boost::posix_time::minutes(5);
            semaphore_.timed_wait( delay );
            clog << "Finished waiting for 5 minutes" << endl;
        }
    }
    catch (boost::interprocess::interprocess_exception& e)
    {
        clog << "Problem with Semaphore: " << e.what() << endl;
    }
}

void Semaphore::post()
{
    try
    {
        semaphore_.post();
    }
    catch (boost::interprocess::interprocess_exception& e)
    {
        clog << "Problem with Semaphore: " << e.what() << endl;
    }
}


} // namespace common
} // namespace eagle
