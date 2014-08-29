/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef SHARED_FASTA_REFERENCE_HH
#define SHARED_FASTA_REFERENCE_HH

#include <vector>
#include <boost/filesystem.hpp>


#ifndef USE_TMP_FASTA_READER
  #include "genome/Reference.hh"
  typedef eagle::genome::MultiFastaReference PreferredFastaReader;
#else
  #include "genome/TmpFastaReader.hh"
  typedef eagle::genome::TmpFastaReader PreferredFastaReader;
  #warning using TmpFastaReader
#endif


namespace eagle
{
namespace genome
{

class SharedFastaReference
{
public:
    static void init( const boost::filesystem::path& sampleGenomeDir )
    {
        SharedFastaReference::fastaReference_ = new PreferredFastaReader( sampleGenomeDir );
        sampleGenomeDir_.clear();
        sampleGenomeDir_.push_back( sampleGenomeDir );
    }
    static void init( const std::vector<boost::filesystem::path>& sampleGenomeDir )
    {
        SharedFastaReference::fastaReference_ = new PreferredFastaReader( sampleGenomeDir );
        sampleGenomeDir_ = sampleGenomeDir;
    }
    static inline PreferredFastaReader *get()
    {
        assert (SharedFastaReference::fastaReference_);
        return SharedFastaReference::fastaReference_;
    }
    static inline void setActive( std::string id )
    {
        if (dictionary_.find(id) == dictionary_.end())
        {
            unsigned int newIdx = fastaReferenceArray_.size();
            PreferredFastaReader* newFastaReference = new PreferredFastaReader( sampleGenomeDir_ );
            fastaReferenceArray_.push_back( newFastaReference );
            dictionary_[id] = newIdx;
        }
        unsigned int idx = dictionary_[id];
        fastaReference_ = fastaReferenceArray_[ idx ];
    }
private:
    static PreferredFastaReader*  fastaReference_;
    static std::vector<boost::filesystem::path> sampleGenomeDir_;
    static std::vector< PreferredFastaReader* >  fastaReferenceArray_;
    static std::map< std::string, unsigned int > dictionary_;
};


} // namespace genome
} // namespace eagle

#endif // SHARED_FASTA_REFERENCE_HH
