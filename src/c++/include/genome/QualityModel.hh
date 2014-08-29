/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_QUALITY_MODEL_HH
#define EAGLE_GENOME_QUALITY_MODEL_HH

#include <fstream>
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "model/Nucleotides.hh"
#include "genome/ErrorModelPlugin.hh"


namespace eagle
{
namespace genome
{

class MotifRepeatQualityDropInfo;

// ClusterErrorModelContext contains all the data that needs to be remembered to process consequent cycles of a read cluster
// It contains a substructure for each error plugin
struct ClusterErrorModelContext
{
    struct
    {
        unsigned int qualityLevel;
    } qualityModelContext;

    struct
    {
        int qualityDrop;
    } phasingContext;

    struct
    {
        char lastBase;
        signed char errorDirection; // insertion=+1, deletion=-1
        unsigned int homopolymerLength;
    } homopolymerModelContext;

    struct
    {
        uint64_t kmer;
        unsigned int kmerLength;
        double qualityDropLevel;
        MotifRepeatQualityDropInfo *shortTermEffect;
        float shortTermQualityDrop;
    } motifQualityDropModelContext;

    ClusterErrorModelContext()
    {
        initialiseForNewRead();
    }

    void initialiseForNewRead()
    {
        qualityModelContext.qualityLevel = 0;
        homopolymerModelContext.lastBase = 0;
        homopolymerModelContext.errorDirection = 0;
        homopolymerModelContext.homopolymerLength = 0;
        motifQualityDropModelContext.kmer = 0;
        motifQualityDropModelContext.kmerLength = 0;
        motifQualityDropModelContext.shortTermQualityDrop = 0;
        motifQualityDropModelContext.shortTermEffect = NULL;
        motifQualityDropModelContext.qualityDropLevel = 0;
        phasingContext.qualityDrop = 0;
    }
};

// QualityModel provides a way to generate the phred quality values
class QualityModel
{
public:
    QualityModel( const std::vector<boost::filesystem::path>& qualityTableFiles );

    unsigned int getQuality( boost::mt19937& randomGen, const unsigned int cycle, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext );
    unsigned int getQuality( boost::mt19937& randomGen, const unsigned int cycle, ClusterErrorModelContext& clusterErrorModelContext );
    double qualToProbError(unsigned int qual);

private:
    unsigned int parseQualityTableFile( const boost::filesystem::path& filename, const int cycleOffset = 0 );
    void createDiscreteDistributionPerCycle();

    std::vector< std::vector< boost::random::discrete_distribution<> > > qualityDistPerCyclePerLastQuality;

    // New stuff
    unsigned int parseBigQualityTableFile( const boost::filesystem::path& filename );
    bool useNewStuff_;
    std::vector< unsigned int > bigTable_;
};


class SequencingMismatchModel
{
public:
    SequencingMismatchModel( const boost::filesystem::path& mismatchTableFilename );
    void apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
    std::vector< boost::random::discrete_distribution<> > errorDistPerBase;
};


class HomopolymerIndelModel
{
public:
    HomopolymerIndelModel( const boost::filesystem::path& homopolymerIndelTableFilename );
    void apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
    std::vector< double > homoDeletionTable_;
    std::vector< double > homoInsertionTable_;
};


class MotifQualityDropModel
{
public:
    MotifQualityDropModel( const boost::filesystem::path& tableFilename );
    void applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext, const unsigned int cycle, boost::mt19937& randomGen );

private:
    MotifRepeatQualityDropInfo* getMotifRepeatQualityDrop( const uint64_t kmer1, const unsigned int repeatKmerLength, const unsigned int repeatCount );
    bool active_;
    std::vector< std::map< uint64_t, boost::shared_ptr< std::vector< MotifRepeatQualityDropInfo > > > > tableData_;
    eagle::model::IUPAC baseConverter_;
};


class RandomQualityDropModel
{
public:
    RandomQualityDropModel( /*const boost::filesystem::path& tableFilename*/ );
    void applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
};


class QualityGlitchModel
{
public:
    QualityGlitchModel( /*const boost::filesystem::path& tableFilename*/ );
    void applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
};


class HappyPhasingModel
{
public:
    HappyPhasingModel( /*const boost::filesystem::path& tableFilename*/ );
    void applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
};


class QQTable
{
public:
    QQTable( const boost::filesystem::path& qqTableFilename );
    double qualToErrorRate( const unsigned int quality );

private:
    std::vector< double > qualityToProbability_;
};


class ErrorModel
{
public:
    enum ErrorType { NoError, BaseSubstitution, BaseDeletion, BaseInsertion } ;

    ErrorModel( const std::vector<boost::filesystem::path>& qualityTableFiles, const boost::filesystem::path& mismatchTableFile, const boost::filesystem::path& homopolymerIndelTableFilename, const boost::filesystem::path& motifQualityDropTableFilename, const boost::filesystem::path& qqTableFilename, const std::vector< std::string >& errorModelOptions );
    void getQualityAndRandomError( boost::mt19937& randomGen, const unsigned int cycle, const char base, unsigned int& quality, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext );

private:
    QualityModel qualityModel_;
    SequencingMismatchModel sequencingMismatchModel_;
    HomopolymerIndelModel homopolymerIndelModel_;
    MotifQualityDropModel motifQualityDropModel_;
    RandomQualityDropModel randomQualityDropModel_;
    QualityGlitchModel qualityGlitchModel_;
    HappyPhasingModel happyPhasingModel_;
    LongreadBaseDuplicationModel longreadBaseDuplicationModel_;
    LongreadDeletionModel longreadDeletionModel_;
    QQTable qqTable_;
    eagle::model::IUPAC baseConverter_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_QUALITY_MODEL_HH
