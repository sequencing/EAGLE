/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "io/Text.hh"
#include "genome/QualityModel.hh"
#include "genome/ErrorModelPlugin.hh"


using namespace std;


namespace eagle
{
namespace genome
{


map< string, pair< string, bool > > ErrorModelPlugin::parseUserOptions( const vector< string >& errorModelOptions, const string& errorModelId )
{
    map< string, pair< string, bool > > result;
    string prefix = errorModelId + ":";

    BOOST_FOREACH( const string& errorModelOption, errorModelOptions )
    {
        if (boost::starts_with( errorModelOption, prefix ))
        {
            vector<string> tokens;
            boost::split(tokens, errorModelOption, boost::is_any_of(":,"));
            for (unsigned int i=1; i<tokens.size(); ++i)
            {
                vector<string> tokens2;
                boost::split(tokens2, tokens[i], boost::is_any_of("="));
                if (tokens2.size() != 2)
                {
                    EAGLE_ERROR( (boost::format("Invalid error model option: %s should be of the form <key>=<value> in %s") % tokens[i] % errorModelOption).str() );
                }
                result[tokens2[0]] = make_pair( tokens2[1], false );
            }
        }
    }
    return result;
}

void ErrorModelPlugin::reportUnusedCommandLineOptionsForThisPlugin()
{
    typedef std::pair< std::string, std::pair< std::string, bool > > FullType;
    BOOST_FOREACH( const FullType & keyVal, parsedUserOptions_ )
    {
        if (keyVal.second.second != true)
        {
            EAGLE_WARNING( (boost::format("Unused plugin command line option: %s=%s") % keyVal.first % keyVal.second.first).str() );
        }
    }
}


LongreadBaseDuplicationModel::LongreadBaseDuplicationModel( const vector< string >& errorModelOptions )
    : ErrorModelPlugin( errorModelOptions, "LONGREAD-base-duplication" )
    , prob_( getParsedValue<double>( "prob", 0.0 ) )
{
    if (prob_)
    {
        clog << "LONGREAD Base Duplication Error Model initialised with:\n prob=" << prob_ << endl;
    }
    else
    {
        clog << "LONGREAD Base Duplication Error Model not in use" << endl;
    }
}

void LongreadBaseDuplicationModel::apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    if (prob_ == 0.0) { return; }

    if (randomErrorType == ErrorModel::NoError)
    {
        double randomValue = (double)randomGen() / (double)randomGen.max();
        if (randomValue < prob_)
        {
            randomErrorType = ErrorModel::BaseInsertion;
        }
    }
}


LongreadDeletionModel::LongreadDeletionModel( const vector< string >& errorModelOptions )
    : ErrorModelPlugin( errorModelOptions, "LONGREAD-deletion" )
    , prob_             ( getParsedValue<double>( "prob"     , 0.0 ) )
    , basesLeftToDelete_( 0 )
{
    string distTableFilename = getParsedValue<string>( "dist-file", ""  );
    if (prob_ > 0 && distTableFilename != "")
    {
        clog << "LONGREAD Deletion Error Model initialised with:\n";
        clog << " prob of deletion = " << prob_ << endl;

        // Parse tsv file
        io::DsvReader tsvReader( distTableFilename );
        vector< double > distValues;
        vector<string> tokens;
        while ( tsvReader.getNextLineFields<'\t'>(tokens) )
        {
            if ( tokens.size() == 0 ) { continue; }
            assert ( tokens.size() == 2 && "There should be 2 entries per line" );
            unsigned int index;
            vector<double> values;
            try
            {
                index = boost::lexical_cast<unsigned int>( tokens[0] );
                std::transform( tokens.begin()+1, tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
            }
            catch (const boost::bad_lexical_cast &e)
            {
                EAGLE_ERROR("Error while reading homopolymer indel table: a numerical field seems to contain non-numerical characters");
            }
            assert( values.size() == 1 );
            if (index >= distValues.size() )
            {
                distValues.resize( index+1 );
            }
            distValues[ index ] = values[0] ;
            clog << " sub-prob of deletion length " << index << " = " << values[0] << endl;
        }

        assert( !distValues.empty() );
        deletionLengthDist_ = boost::random::discrete_distribution<>( distValues );
    }
    else
    {
        clog << "LONGREAD Deletion Error Model not in use" << endl;
    }

    reportUnusedCommandLineOptionsForThisPlugin();
}

void LongreadDeletionModel::apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    if (prob_ == 0.0) { return; }

    if (basesLeftToDelete_ > 0)
    {
        randomErrorType = ErrorModel::BaseDeletion;
        --basesLeftToDelete_;
    }
    else
    {
        if (randomErrorType == ErrorModel::NoError)
        {
            // randomGen() can get values from 0 to max(), both included
            double randomValue = (double)randomGen() / ((double)randomGen.max() + 1.0);
            if (randomValue < prob_)
            {
                basesLeftToDelete_ = deletionLengthDist_(randomGen); //(int)( randomValue / prob_ * (double)(maxLength_ - minLength_ + 1) + minLength_ - 1 );
                if (basesLeftToDelete_ > 0)
                {
                    randomErrorType = ErrorModel::BaseDeletion;
                    --basesLeftToDelete_;
                }
            }
        }
    }
}


} // namespace genome
} // namespace eagle
