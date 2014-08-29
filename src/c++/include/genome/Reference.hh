/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to deal with reference genomes.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_GENOME_REFERENCE_HH
#define EAGLE_GENOME_REFERENCE_HH


#include "common/FileSystem.hh"
#include "io/Fasta.hh"
#include "model/Contig.hh"

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>


namespace eagle
{
namespace genome
{


/*
 * \brief When the reference is a set of FASTA files
 */
class FastaReference: boost::noncopyable
{
public:
    FastaReference(   // read-only
        const eagle::io::FastaMetadata& metadata );
    FastaReference(   // write-only
        const boost::filesystem::path& outputDir,
        const bool overwrite = false );
    FastaReference(   // read-write
        const eagle::io::FastaMetadata& metadata,
        const boost::filesystem::path& outputDir,
        const bool overwrite);

    void load();
    void save();

    char get( const unsigned long globalPos,      const unsigned long offset, bool &overlapContigBoundary );
    char get( const eagle::model::Locus location, const unsigned long offset, bool &overlapContigBoundary );
    unsigned long read( eagle::model::Contig &contig, const std::string &contigName );

    unsigned long local2global(const eagle::model::Locus& position);
    eagle::model::Locus global2local(unsigned long globalPos);
    void convertFromGlobalPos( const unsigned long globalPos, int& refId, unsigned long& posInContig );

    std::string currentChromosome() const {return reader_.getContigName();}
    unsigned long estimatedLength() {return reader_.getContigSize();}
    std::vector<eagle::model::Contig>::const_iterator begin() const {return reference_.begin();}
    std::vector<eagle::model::Contig>::iterator begin() {return reference_.begin();}
    std::vector<eagle::model::Contig>::const_iterator end() const {return reference_.end();}
    std::vector<eagle::model::Contig>::iterator end() {return reference_.end();}
    eagle::model::Contig operator[](int i) const {return reference_[i];}
    eagle::model::Contig& operator[](int i) {return reference_[i];}
    //TODO:     {return *(reference_.find(contigName,sameName()));}
    eagle::model::Contig & getContig(const std::string &contigName);
    eagle::model::Contig const & getContig(const std::string &contigName) const;

    size_t contigCount() const {return reference_.size();}
    size_t fileCount() const {return ( (mode_ & std::ios_base::out) ? writer_.size() : reader_.size() );}
    void clear() {reference_.clear();}
    void resize(unsigned int n) {reference_.resize(n);}

    size_t length();
    size_t length(unsigned int i) const {return reference_[i].size();}
    std::vector<std::string> allContigNames() const;
    std::vector<unsigned long> allContigLengths() const;
    unsigned long getContigLength(const std::string &contigName) const;
    eagle::io::FastaMetadata & metadata()             {return metadata_;}
    eagle::io::FastaMetadata const & metadata() const {return metadata_;}

private:
    void inputStructure( const eagle::io::FastaMetadata& metadata );
    void outputStructure( const boost::filesystem::path& outputDir );
    void inputMode();
    void outputMode();

    eagle::io::MultiFastaReader reader_;
    eagle::io::MultiFastaWriter writer_;
    const std::ios_base::openmode mode_;

    eagle::io::FastaMetadata metadata_;
    std::vector< eagle::model::Contig > reference_;
    eagle::io::FastaInfo currentGetInfo_;

    unsigned int local2globalCache_roundRobin;
    std::vector<std::string> local2globalCache_chr;
    std::vector<unsigned long> local2globalCache_pos;

    eagle::io::FastaInfo global2localCache;
};

typedef std::vector<eagle::model::Contig>::const_iterator ReferenceIterator;
typedef std::pair< ReferenceIterator,ReferenceIterator > ReferenceBounds;

typedef std::vector<eagle::model::Contig>::iterator iterator;


/*
 * \brief When the reference is either a directory, or a set of FASTA files
 */
class MultiFastaReference : public FastaReference
{
public:
    MultiFastaReference(   // read-only (paths can be either dirs or files)
        const std::vector< boost::filesystem::path > &inputPaths )
    : FastaReference( initialize(inputPaths) )
    , overwrite_(false)    // N/A
    {}
    MultiFastaReference(   // read-only (path can be either dir or file)
        const boost::filesystem::path &inputPath )
    : FastaReference( initialize(inputPath) )
    , overwrite_(false)    // N/A
    {}
    MultiFastaReference(   // write-only
        const boost::filesystem::path outputDir,
        const bool overwrite)
    : FastaReference( outputDir,overwrite )
    , overwrite_(overwrite)
    {}
    void saveMetadata();

private:
    eagle::io::FastaMetadata initialize( const std::vector<boost::filesystem::path>& inputPaths );
    eagle::io::FastaMetadata initialize( const boost::filesystem::path& inputPath );

    const bool overwrite_;
};


/*****************/
/**  METADATA   **/
/*****************/


class GenomeSizeXml
{
public:
    GenomeSizeXml(const boost::filesystem::path& indexPath);
    GenomeSizeXml(const boost::filesystem::path& indexPath, bool overwrite);

    void load(eagle::io::FastaMetadata &index);
    void load(std::istream &is, eagle::io::FastaMetadata &index);
    void save(const eagle::io::FastaMetadata &index);
    void save(const boost::filesystem::path &file, const eagle::io::FastaMetadata &index);

    bool enabled() const {return !path_.empty();}
private:
    boost::filesystem::path path_;
    const std::ios_base::openmode mode_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_REFERENCE_HH
