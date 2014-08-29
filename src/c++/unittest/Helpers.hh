/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Some helper functions.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_UNIT_TEST_HELPERS_HH
#define EAGLE_UNIT_TEST_HELPERS_HH

#include <iostream>
#include <string>
#include <vector>
#include <boost/foreach.hpp>


inline std::vector<char> string2vector(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

inline std::string vector2string(const std::vector<char> &vec)
{
    return std::string(vec.begin(), vec.end());
}

inline std::string substr(const std::vector<char> &from, std::string::size_type __pos = 0, std::string::size_type __n = std::string::npos)
{
    return std::string(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::string &from, std::string::size_type __pos = 0, std::string::size_type __n = std::string::npos)
{
    return std::vector<char>(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos, std::string::size_type __n)
{
    return std::vector<char>(from.begin() + __pos, from.begin() + __pos + __n);
}

inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos)
{
    return std::vector<char>(from.begin() + __pos, from.end());
}

inline std::vector<char> operator+(const std::vector<char> &right, const std::vector<char> &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline std::vector<char> operator+(const std::vector<char> &right, const std::string &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline void show(const std::vector<char> &s)
{
    BOOST_FOREACH(char c, s) std::cerr << c;
}


#endif // EAGLE_UNIT_TEST_HELPERS_HH
