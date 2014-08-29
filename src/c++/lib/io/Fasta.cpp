/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Component to read/write FASTA files.
 **
 ** \author Mauricio Varea
 **/

#include <algorithm>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "io/Fasta.hh"
#include "common/Exceptions.hh"
#include "common/Logger.hh"

namespace eagle
{
namespace io
{


std::ostream& operator<<( std::ostream& os, const FastaInfo& fi )
{
    os << fi.contigName << "\t"
       << fi.contigSize << "\t"
       << fi.position.second << "\t"
       << fi.contigWidth.first << "\t"
       << fi.contigWidth.second;
    return os;
}

std::ostream& operator<<( std::ostream& os, const FastaIndex& fi )
{
    for(std::vector< FastaInfo >::const_iterator it = fi.second.begin();
        it != fi.second.end();
        ++it)
    {
        if (fi.second.begin() != it) os << std::endl;
        os << " - " << fi.first << " " << *it << " (" << it->position.first << ")";
    }
    return os;
}

std::ostream& operator<<( std::ostream& os, const FastaMetadata& fm )
{
    bool first(true);
    BOOST_FOREACH(const FastaIndex &idx, fm)
    {
        if (first) {first = false;} else {os << std::endl;}
        os << idx;
    }
    return os;
}


void FastaMetadata::init( const std::vector< boost::filesystem::path > &P )
{
    using boost::lambda::_1;
    using boost::lambda::bind;
    this->clear();
    std::transform( P.begin(), P.end(), std::inserter(*this,this->begin()),
                    boost::bind( &std::make_pair< boost::filesystem::path, std::vector< eagle::io::FastaInfo > >,
                                 _1, std::vector< eagle::io::FastaInfo >() ));
}


void FastaMetadata::update( const boost::filesystem::path &p, const FastaInfo &fi )
{
    FastaInfo current;
    std::pair< FastaIterator, bool > status = this->insert(FastaIndex(
                                                           std::make_pair( p, std::vector<FastaInfo>(1,fi)) ));
    FastaIterator idx(status.first);
    unsigned int contigCount = 0;
    BOOST_FOREACH(const FastaIndex &idx2, *this)
    {
        contigCount += idx2.second.size();
    }
    
    if (status.second)
    {   // didn't find path => fi has been inserted as new element
        idx->second[0].contigNumber = contigCount;
        EAGLE_DEBUG(0, "[metadata] " << FastaIndex(std::make_pair( p, this->at(p) )) );
    } else {
        using boost::lambda::_1;
        using boost::lambda::bind;
        this->insert(FastaIndex( std::make_pair( p, std::vector<FastaInfo>(1,fi)) ));
        std::vector< FastaInfo >::iterator it = std::find_if( idx->second.begin(), idx->second.end(),
                                                              bind( &FastaInfo::sameName, _1, fi.contigName) );
        if (idx->second.end() != it)
        {   // some info about fi already exists => check that it's up-to-date
            if (fi.hasSize() && !it->sameSize(fi))
            {
                if (it->hasSize())
                {
                    EAGLE_WARNING_IF(     it->contigSize, "Metadata mismatch" );
                    EAGLE_WARNING_CONT_IF(it->contigSize, (boost::format( "    Chromosome '%s': Updating length to %lu (from %lu)" )
                                                        % it->contigName % fi.contigSize % it->contigSize).str() );
                }
                it->setSize( fi.contigSize );
            }
            if (fi.hasLogicalPosition() && !it->hasLogicalPosition())
            {
/*
                if (it->hasLogicalPosition())
                {
                    EAGLE_WARNING_IF(     fi.position.first, "Metadata mismatch" );
                    EAGLE_WARNING_CONT_IF(fi.position.first, (boost::format( "    Chromosome '%s': Updating global position (from %lu to %lu)" )
                                                           % it->contigName % it->position.first % fi.position.first).str() );
                }
*/
                it->setLogicalPosition( fi.position.first );
            }
            if (fi.hasPhysicalPosition() && !it->samePhysicalPosition(fi))
            {
                if (it->hasPhysicalPosition())
                {
                    EAGLE_WARNING_IF(     fi.position.second, "Metadata mismatch" );
                    EAGLE_WARNING_CONT_IF(fi.position.second, (boost::format( "                %s : Updating local indexing position (from %lu to %lu)" )
                                                            % std::string(it->contigName.length(),' ') % it->position.second % fi.position.second ).str() );
                }
                it->setPhysicalPosition( fi.position.second );
            }
            if (fi.hasWidth() && !it->sameWidth(fi))
            {
                if (it->hasWidth())
                {
                    EAGLE_WARNING_IF(     fi.hasWidth(), "Metadata mismatch" );
                    EAGLE_WARNING_CONT_IF(fi.contigWidth.first, (boost::format( "    Chromosome '%s': Updating number of bases (from %lu to %lu)" )
                                                              % it->contigName % it->contigWidth.first % fi.contigWidth.first).str() );
                    EAGLE_WARNING_CONT_IF(fi.contigWidth.second, (boost::format( "                %s : Updating number of bytes (from %lu to %lu)" )
                                                               % std::string(it->contigName.length(),' ') % it->contigWidth.second % fi.contigWidth.second ).str() );
                }
                it->setWidth( fi.contigWidth.first, fi.contigWidth.second );
            }
        } else {
            // if fi seems new to idx[p], then add it to the end...
            idx->second.push_back( fi );
            idx->second.back().contigNumber = contigCount+1;
        }
    }
}


/*************/
/**  READ   **/
/*************/


FastaReader &FastaReader::get(char &base, bool &newContig)
{
    char tmp;
    newContig = false;
    std::ifstream::get(tmp);
    while (*this && '>' == tmp)
    {
        char contigName[1000];
        std::ifstream::getline(contigName,1000);
        contigName_ = std::string(contigName);
        get(base, newContig);
        newContig = true;
        return *this;
    }
    while (*this && '\n' == tmp)
    {
        return get(base, newContig);
    }
    if (*this)
    {
        base = tmp;
    }
    return *this;
}


FastaReader& FastaReader::read( unsigned long pos, unsigned long length, unsigned long lineCount)
{
    cache_.resize( length+lineCount );
    if (! std::ifstream::seekg( pos, std::ios_base::beg ))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to access position %lu in FASTA file") % pos).str()));
    }
    if (! std::ifstream::read( &cache_[0], length+lineCount ))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read %lu bases from FASTA file") % length).str()));
    }
    return *this;
}

MultiFastaReader::MultiFastaReader(const FastaMetadata &index)
    : index_(index)
    , this_(index_.begin())
    , globalContigId_(-1)
    , localContigId_(-1)
{
    unsigned long globalPos(0L);
    BOOST_FOREACH(FastaIndex &idx, index_)
    {
        if( boost::filesystem::exists(idx.first.string() + ".fai") )
        {   // get indexing info from *.fai
            FastaInfo info;
            FaiReader fai( std::vector<boost::filesystem::path>(1, idx.first.string() + ".fai") );
            while (fai.getNextIndex(info))
            {
                info.position.first = globalPos;
                index_.update( idx.first, info );
                globalPos += info.contigSize;
            }
        } else { // check if genome_size.xml has written some info that may save us having to load the *.fa
            if (idx.second.end() != std::find_if( idx.second.begin(), idx.second.end(), FastaInfo::needsUpdating ))
            {   // otherwise, calculate indexing info from contents of *.fa
                boost::uintmax_t size = boost::filesystem::file_size( idx.first );
                this->open();
                FastaReader::read( 0, size, 0 );
                std::vector<char>::iterator it1 = cache().begin();
                char delim[] = {'>','\n',' '};
                while( it1 != cache().end() )
                {
                    std::vector<char>::iterator it2,it3,it4;
                    if( cache().end() == (it2 = std::find_first_of(it1+1,cache().end(),delim+1,delim+3)) )
                    {
                        BOOST_THROW_EXCEPTION(common::CorruptedFileException("FASTA","*** found open-ended header ***"));
                    }
                    if( cache().end() == (it3 = std::find(it2,cache().end(),delim[1])) )
                    {
                        BOOST_THROW_EXCEPTION(common::CorruptedFileException("FASTA","*** found open-ended header ***"));
                    }
                    if( cache().end() == (it4 = std::find(it3 + 1,cache().end(),delim[1])) )
                    {
                        BOOST_THROW_EXCEPTION(common::CorruptedFileException("FASTA","*** found open-ended body ***"));
                    }
                    FastaInfo info( std::string(it1 + 1, it2 ),
                                    0,
                                    globalPos,
                                    it3 + 1 - cache().begin(),
                                    it4 - it3 -1 );
                    unsigned int adj = static_cast<unsigned int>( cache().begin() == it1 );
                    if ( (it1 = std::find(it4,cache().end(),delim[0])) != cache().end() )
                    {
                        info.setSize( it1 - it3 - (it1-it3)/info.contigWidth.second - 1 );
                    } else {
                        info.setSize( cache().begin() + size - it3 - (cache().begin()+size-it3)/info.contigWidth.second - 1 - adj );
                    }
                    index_.update( idx.first, info );
                    globalPos += info.contigSize;
                }
            }
        }
        ++this_;
    }
    if (!index_.empty())
    {
        close();
        clear();
        cache().clear();
        this_ = index_.begin();
        this->open();
        if( !this_->second.empty() )
        {
            setContigName( this_->second[0].contigName );
        }
    }
}


void MultiFastaReader::open()
{
    if (boost::filesystem::is_directory( this_->first ))
    {
        EAGLE_ERROR( (boost::format("%s is a directory instead of being a fasta file. You may want to use --whole-genome instead of --reference-genome") % this_->first.string()).str() );
    }
    std::ifstream::open(this_->first.string().c_str());
    if (!is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open FASTA file %s for reading") % this_->first).str()));
    }
}

/*
 * (return value == false) => did not open any new file
 * (return value == true)  => opened a new file as a result of find() returning a different value
 */
bool MultiFastaReader::open(const boost::filesystem::path fastaPath)
{
    if (fastaPath.string().empty())
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException("FASTA filename cannot be an empty string"));
    }
    FastaIterator it = index_.find(fastaPath);
    if (it == index_.end())
    {
        std::stringstream message;
        message << (boost::format("*** Cannot open %s, as this name was not given at initialization time ***") % fastaPath).str() << std::endl;
        message << "    The following paths are available:" << std::endl;
        BOOST_FOREACH(const eagle::io::FastaIndex &index, index_)
        {
            message << (boost::format("       %s") % index.first).str() << std::endl;
        }
        BOOST_THROW_EXCEPTION(common::PreConditionException( message.str() ));
    }
    if (it != this_)
    {
        close();
        clear();
        this_ = it;
        this->open();
        return true;
    }
    return false;
}


boost::filesystem::path MultiFastaReader::find(const std::string& chr, FastaInfo& info)
{
    unsigned int j = 0;
    for (FastaMetadata::const_iterator idx = index_.begin();
         idx != index_.end(); ++idx)
    {
        for (unsigned int i = 0; i < idx->second.size(); i++)
        {
            if ( idx->second[i].sameName( chr ) )
            {
                info = idx->second[i];
                setContigName( info.contigName );
                setLocalContigId( i );
                setGlobalContigId( j );  // Needed for Bam ID
                return idx->first;
            }
            j++;
        }
    }
    EAGLE_WARNING( (boost::format("Contig '%s' not found!") % chr ).str() );
    EAGLE_WARNING( "Current contigs are:" );
    for (FastaMetadata::const_iterator idx = index_.begin();
         idx != index_.end(); ++idx)
    {
        EAGLE_WARNING( (boost::format("  - '%s':") % idx->first.string()).str() );
        for (unsigned int i = 0; i < idx->second.size(); i++)
        {
            EAGLE_WARNING( (boost::format("    - '%s' (%d)") % idx->second[i].contigName % i).str() );
            if ( idx->second[i].sameName( chr ) )
            {
                info = idx->second[i];
                setContigName( info.contigName );
                setLocalContigId( i );
                setGlobalContigId( j );  // Needed for Bam ID
                return idx->first;
            }
            j++;
        }
    }
    return boost::filesystem::path();
}


boost::filesystem::path MultiFastaReader::find(unsigned long pos, FastaInfo& info)
{
    unsigned int j = 0;
    for (FastaMetadata::const_iterator idx = index_.begin();
         idx != index_.end(); ++idx)
    {
        for (unsigned int i = 0; i < idx->second.size(); i++)
        {
            if ( idx->second[i].within( pos ) )
            {
                info = idx->second[i];
                setContigName( info.contigName );
                setLocalContigId( i );
                setGlobalContigId( j );  // Needed for Bam ID
                return idx->first;
            }
            j++;
        }
    }
    EAGLE_WARNING( (boost::format("Could not determine in which contig global pos %lu belongs to") % pos).str() );
    return boost::filesystem::path();
}


std::vector< FastaInfo > & MultiFastaReader::index(boost::filesystem::path fastaPath)
{
    FastaIterator idx = index_.find(fastaPath);
    if (index_.end() == idx) EAGLE_ERROR("Non-existent path!");
    return idx->second;
}
std::vector< FastaInfo > const & MultiFastaReader::index(boost::filesystem::path fastaPath) const
{
    FastaConstIterator idx = index_.find(fastaPath);
    if (index_.end() == idx) EAGLE_ERROR("Non-existent path!");
    return idx->second;
}


MultiFastaReader& MultiFastaReader::get(char &base, bool &newContig)
{
    char tmp;
    FastaReader::get(tmp, newContig);
    if (!*this)
    {
        close();
        //clear();
        ++this_;
        if (index_.end() != this_)
        {
            open();
            if (get(base, newContig))
            {
                newContig = true;
                ++globalContigId_;
                ++localContigId_;
            }
        }
    }
    else
    {
        base = tmp;
        if (newContig)
        {
            ++globalContigId_;
            ++localContigId_;
        }
    }
    return *this;
}


MultiFastaReader& MultiFastaReader::read( const FastaInfo &info, unsigned long skip, unsigned long size)
{
    if (0 == size) {size = info.contigSize;}
    if (skip > info.contigSize)
    {
        EAGLE_ERROR( (boost::format("cannot start reading from base number %lu in contig '%s', as it only has %lu bases")
                      % skip % info.contigName % info.contigSize).str() );
    }
    unsigned int eol = info.contigWidth.second - info.contigWidth.first;
    assert(1 == eol && "we do not support platforms that have EOL > 1 byte");
    if ((skip + size) <= info.contigSize )
    {
        FastaReader::read( info.position.second + skip + eol * (skip / info.contigWidth.first),
                           size,
                           eol * (size / info.contigWidth.first));
    } else {
        EAGLE_ERROR("Tried to read outside contig boundary");
    }
    return *this;
}


unsigned long MultiFastaReader::getContigSize()
{
    using boost::lambda::_1;
    using boost::lambda::bind;
    FastaInfo fi;
    boost::filesystem::path p = this->find(getContigName(),fi);
    std::vector< FastaInfo >::iterator it = std::find_if( index_[p].begin(), index_[p].end(),
                                                          bind( &FastaInfo::sameName, _1, fi.contigName) );
    if (index_[p].end() != it)
    {
        return it->contigSize;
    }
    return 0;
}


/**************/
/**  WRITE   **/
/**************/


void FastaWriter::open(const boost::filesystem::path fastaPath)
{
    fastaPath_ = fastaPath;
    if (fastaPath_.string().empty())
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException("FASTA filename cannot be an empty string"));
    }
    this->open();
}

void FastaWriter::open()
{
    std::ofstream::open(fastaPath_.string().c_str());
    if (!is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open FASTA file %s for writing") % fastaPath_).str()));
    }
}

FastaWriter& FastaWriter::write(const eagle::model::Contig& contig)
{
    std::ofstream::put('>');
    if (!std::ofstream::write(contig.name().c_str(), contig.name().length()))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to write contig name into %s") % fastaPath_).str()));
    }
    eagle::model::Contig::const_iterator base(contig.begin());
    unsigned long position(0);
    bool eol(false);
    while(contig.end() != base)
    {
        if (position % contigWidth_ || eol)
        {
            std::ofstream::put(*base);
            ++base;
            ++position;
            eol = false;
        } else {
            std::ofstream::put('\n');
            eol = true;
        }
    }
    std::ofstream::put('\n');
    return *this;
}

MultiFastaWriter& MultiFastaWriter::write(const eagle::model::Contig& contig, const FastaInfo& info)
{
    if (fastaDir_.empty())
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException("Path to FASTA output cannot be empty"));
    }
    boost::filesystem::path fastaPath = fastaDir_ / (boost::format("%s.fa") % contig.id()).str();
    boost::filesystem::path indexPath = fastaDir_ / (boost::format("%s.fa.fai") % contig.id()).str();
    // no need to check for indexPath, as base class DsvWriter already performs this check...
    if ( boost::filesystem::exists(fastaPath) )
    {
        if (overwrite_)
        {
            EAGLE_WARNING( "Overwriting " << fastaPath << " due to the --force switch." );
        } else {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Cannot write FASTA file %s: File already exists!") % fastaPath).str()));
        }
    }

    FastaWriter::open( fastaPath );
    FastaWriter::write( contig );
    close();
    FaiWriter fai( indexPath, overwrite_ );
    fai.open( 0 );  // always writing index for one contig at a time
    fai.write( info );
    index_.update( fastaPath, info );

    return *this;
}


/*****************/
/**  METADATA   **/
/*****************/


bool FaiReader::getNextIndex(FastaInfo &index)
{
    std::vector<std::string> tokens;
    while (getNextLineFields<'\t'>(tokens))
    {
        // Detect and warn about malformed lines
        if (5 == tokens.size())
        {
            index = FastaInfo(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4]);
        } else if (5 > tokens.size()) {
            EAGLE_WARNING_IF( pathList_.begin() <= thisPath_ && pathList_.end() > thisPath_,
                              (boost::format("Only %d tokens in %s:%lu") % tokens.size() % *thisPath_ % lineCount_ ).str() );
            EAGLE_WARNING_CONT( "*** LINE IGNORED ***" );
            continue;
        } else if (5 < tokens.size()) {
            EAGLE_WARNING_IF( pathList_.begin() <= thisPath_ && pathList_.end() > thisPath_,
                              (boost::format("More tokens (%d) than expected (5) at %s:%lu") % tokens.size() % *thisPath_ % lineCount_ ).str() );
            index = FastaInfo(tokens[0],tokens[1],tokens[2],tokens[3],tokens[4]);
        }
        return true;
    }
    return false;
}


FaiWriter& FaiWriter::write(const FastaInfo& index)
{
    if ( !(*this << index) )
    {
        std::stringstream message;
        message << "Failed to write line:  " << index << std::endl
                << "       into the FAI file" << std::endl;
//                << "       into file:  " << *thisPath_ << std::endl;
        BOOST_THROW_EXCEPTION(eagle::common::IoException(errno, message.str()));
    }
    std::ofstream::put('\n');
    return *this;
}



} // namespace io
} // namespace eagle
