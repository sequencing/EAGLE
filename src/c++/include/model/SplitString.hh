/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Lazy-evaluated string split
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MODEL_SPLITSTRING_HH
#define EAGLE_MODEL_SPLITSTRING_HH

#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "common/Exceptions.hh"


namespace eagle
{
namespace model
{

class SplitString
{
public:
    SplitString( const std::string& line, const std::string& separators ) 
        : line_( line )
        , separators_( separators )
        , sizeLazilyEvaluated_( false )
    {
//        boost::split(tokens, line, boost::is_any_of( separators ));
    }

    unsigned int size() {
        if (!sizeLazilyEvaluated_)
        {
            size_ = 1;
            BOOST_FOREACH( const char c, separators_ )
            {
                size_ += std::count( line_.begin(), line_.end(), c );
            }
            sizeLazilyEvaluated_ = true;
            tokenPos.resize( size_+1, std::string::npos );
            tokenPos[0] = 0;
            tokenPos[size_] = line_.length()+1;
        }
        return size_;
    }

    std::string operator[]( const unsigned int index ) {
//        return tokens[index];
        size_t nextTokenPos = getTokenPos(index+1);
        size_t thisTokenPos = getTokenPos(index);
        return line_.substr( thisTokenPos, nextTokenPos - thisTokenPos - 1 );
    }

private:
    size_t getTokenPos( unsigned int index ) {
        size_t thisTokenPos = tokenPos[index];
        if (thisTokenPos != std::string::npos)
        {
            return thisTokenPos;
        }
        size_t previousTokenPos = getTokenPos( index-1 );
        thisTokenPos = line_.find_first_of( separators_, previousTokenPos );
        assert( thisTokenPos != std::string::npos );
        ++thisTokenPos;
        return thisTokenPos;
    }
//    std::vector<std::string> tokens;

    std::string line_;
    std::string separators_;
    bool sizeLazilyEvaluated_;
    unsigned int size_;
    std::vector<size_t> tokenPos;
};



} // namespace model
} // namespace eagle

#endif // EAGLE_MODEL_SPLITSTRING_HH
