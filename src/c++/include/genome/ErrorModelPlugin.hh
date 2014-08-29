/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_ERROR_MODEL_PLUGIN_HH
#define EAGLE_GENOME_ERROR_MODEL_PLUGIN_HH

#include <fstream>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "model/Nucleotides.hh"
#include "common/Logger.hh"


namespace eagle
{
namespace genome
{

class ClusterErrorModelContext;

class ErrorModelPlugin
{
public:
    ErrorModelPlugin( const std::vector< std::string >& errorModelOptions, const std::string& errorModelId )
        : parsedUserOptions_( parseUserOptions( errorModelOptions, errorModelId ) )
    {
    }
    virtual ~ErrorModelPlugin()
    {
    }
    virtual void apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext ) = 0;

protected:
    template <class T>
    T getParsedValue( const std::string& key, const T& defaultValue )
    {
        std::map< std::string, std::pair< std::string, bool > >::iterator valueItr = parsedUserOptions_.find( key );
        if (valueItr != parsedUserOptions_.end())
        {
            T val;
            try
            {
                val = boost::lexical_cast<T>( valueItr->second.first );
            }
            catch (boost::bad_lexical_cast &)
            {
                EAGLE_ERROR( (boost::format("Invalid type for error model command line option \"%s\"") % key).str() );
            }
            valueItr->second.second = true; // Remember that this key+value has been used
            return val;
        }
        return defaultValue;
    }

    void reportUnusedCommandLineOptionsForThisPlugin();

private:
    std::map< std::string, std::pair< std::string, bool > > parseUserOptions( const std::vector< std::string >& errorModelOptions, const std::string& errorModelId );
    std::map< std::string, std::pair< std::string, bool > > parsedUserOptions_;
};


class LongreadBaseDuplicationModel: public ErrorModelPlugin
{
public:
    LongreadBaseDuplicationModel( const std::vector< std::string >& errorModelOptions );
    virtual void apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
    const double prob_;
};


class LongreadDeletionModel: public ErrorModelPlugin
{
public:
    LongreadDeletionModel( const std::vector< std::string >& errorModelOptions );
    void apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
    const double prob_;
    boost::random::discrete_distribution<> deletionLengthDist_;
    unsigned int basesLeftToDelete_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_ERROR_MODEL_PLUGIN_HH
