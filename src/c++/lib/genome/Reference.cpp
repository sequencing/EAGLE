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

#include "common/Logger.hh"
#include "model/Struct.hh"
#include "genome/Reference.hh"

#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <numeric>
#include <algorithm>


namespace eagle
{
namespace genome
{

FastaReference::FastaReference (
    const eagle::io::FastaMetadata& metadata
    )
    : reader_(metadata)
    , mode_(std::ios_base::in)
    , local2globalCache_roundRobin( 0 )
    , local2globalCache_chr(2)
    , local2globalCache_pos(2)
{
    inputStructure(metadata);
#ifdef EAGLE_DEBUG_MODE
    std::cout << "FASTA index Metadata:" << std::endl
              << reader_.index() << std::endl;
#endif
}

FastaReference::FastaReference (
    const boost::filesystem::path& outputDir,
    const bool overwrite
    )
    : writer_(outputDir,overwrite)
    , mode_(std::ios_base::out)
    , local2globalCache_roundRobin( 0 )
    , local2globalCache_chr(2)
    , local2globalCache_pos(2)
{
    outputStructure(outputDir);
}

FastaReference::FastaReference (
    const eagle::io::FastaMetadata& metadata,
    const boost::filesystem::path& outputDir,
    const bool overwrite
    )
    : reader_(metadata)
    , writer_(outputDir,overwrite)
    , mode_(std::ios_base::in & std::ios_base::out)
    , local2globalCache_roundRobin( 0 )
    , local2globalCache_chr(2)
    , local2globalCache_pos(2)
{
    inputStructure(metadata);
    outputStructure(outputDir);
}

void FastaReference::inputStructure( const eagle::io::FastaMetadata& metadata )
{
    BOOST_FOREACH(const eagle::io::FastaIndex &idx, metadata)
    {
        std::clog << "+ Input reference genome: " << idx.first << std::endl;
        EAGLE_WARNING_IF(!boost::filesystem::exists( idx.first ), "The above file does not exist!")
    }
}

void FastaReference::outputStructure( const boost::filesystem::path& outputDir )
{
    if (!boost::filesystem::exists(outputDir) && !boost::filesystem::create_directory(outputDir) && !boost::filesystem::exists(outputDir))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to create directory %s for sample genome") % outputDir).str()));
    }
    std::clog << "+ Output path to sample genome: " << outputDir << std::endl;
}

void FastaReference::inputMode()
{
    if (! mode_ & std::ios_base::in)
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException(
                              (boost::format("Not allowed to read from %s" ) % reader_.file()).str() ));
    }
}
// TODO: print out the file name
void FastaReference::outputMode()
{
    if (! mode_ & std::ios_base::out)
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException(
                              (boost::format("Not allowed to write to FASTA file" )).str() ));
    }
}

unsigned long FastaReference::local2global( const eagle::model::Locus& location)
{
    assert( location.pos() > 0 );
    if (location.chr() == local2globalCache_chr[0])
    {
        return local2globalCache_pos[0] + location.pos() - 1;
    }
    else if (location.chr() == local2globalCache_chr[1])
    {
        return local2globalCache_pos[1] + location.pos() - 1;
    }
    else
    {
        std::clog << "local2global: updating cache pos " << local2globalCache_roundRobin << " from " <<  local2globalCache_chr[local2globalCache_roundRobin]<< " to " << location.chr() << std::endl;
        eagle::io::FastaInfo info;
        boost::filesystem::path file = reader_.find(location.chr(),info);
        if (file.empty())
        {
            std::stringstream message;
            message << "Could not convert local position " << location << " into global" << std::endl;
            EAGLE_ERROR( message.str() );
        }
        local2globalCache_chr[local2globalCache_roundRobin] = location.chr();
        local2globalCache_pos[local2globalCache_roundRobin] = info.position.first;
        local2globalCache_roundRobin = (local2globalCache_roundRobin+1)%2;
        return info.position.first + location.pos() - 1;
    }
}

eagle::model::Locus FastaReference::global2local(unsigned long globalPos)
{
    if ( !global2localCache.within(globalPos) )
    {
        std::clog << "global2local: updating cache from " << global2localCache.contigName;
        //    eagle::io::FastaInfo info;
        boost::filesystem::path file = reader_.find(globalPos,global2localCache);
        if (file.empty())
        {
            EAGLE_ERROR( (boost::format("Could not convert global location %lu into local") % globalPos).str() );
        }
        std::clog << " to " << global2localCache.contigName << std::endl;
    }
    return eagle::model::Locus( global2localCache.contigName, globalPos - global2localCache.position.first + 1);
}

void FastaReference::convertFromGlobalPos( const unsigned long globalPos, int& refId, unsigned long& posInContig )
{
    eagle::model::Locus location = global2local(globalPos);
    refId       = reader_.getGlobalContigId();
    posInContig = location.pos();
}

char FastaReference::get( const unsigned long globalPos, const unsigned long offset, bool &overlapContigBoundary )
{
    inputMode();
    if ( globalPos < currentGetInfo_.position.first || globalPos >= (currentGetInfo_.position.first + currentGetInfo_.contigSize) )
    {
        boost::filesystem::path file = reader_.find(globalPos,currentGetInfo_);
        if (reader_.open( file ) || reader_.cache().empty())
        {
            if (! reader_.read(currentGetInfo_) )
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read block of data from %s (cache has %lu bases)")
                                                                  % file % reader_.cache().size()).str()));
            }
        }
    }
    assert( globalPos >= currentGetInfo_.position.first && globalPos < (currentGetInfo_.position.first + currentGetInfo_.contigSize) && "Global position needs to be situated within contig's range" );
    overlapContigBoundary = !reader_.inCache( globalPos + offset );
    return reader_[globalPos + offset];
}


char FastaReference::get( const eagle::model::Locus location, const unsigned long offset, bool &overlapContigBoundary )
{
    inputMode();
    if (location.chr() != currentGetInfo_.contigName)
    {
        boost::filesystem::path file = reader_.find(location.chr(),currentGetInfo_);
        if (reader_.open( file ) || reader_.cache().empty())
        {
            if (! reader_.read(currentGetInfo_) )
            {
                BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read block of data from %s (cache has %lu bases)")
                                                                  % file % reader_.cache().size()).str()));
            }
        }
    }
    assert( location.pos() + currentGetInfo_.position.first + offset >= currentGetInfo_.position.first && "Global position needs to be situated within contig's range" );
    overlapContigBoundary = !reader_.inCache( location.pos() + offset );
    return reader_[location.pos() + offset];
}


// read entire contig
unsigned long FastaReference::read( eagle::model::Contig &contig, const std::string &contigName )
{
    inputMode();
    eagle::io::FastaInfo info;
    boost::filesystem::path file = reader_.find(contigName,info);
    if (reader_.open( file ))
    {
        if (! reader_.read(info) )
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to read block of data from %s") % file).str()));
        }
    }
    std::vector<char> segment;
    contig.reset();
    contig.name(contigName);
    segment.resize( std::remove_copy( reader_.cache().begin(), reader_.cache().end(), segment.begin(), '\n' ) - segment.begin() );
    contig.append(segment);
    return segment.size();
}

// TODO: use cache() and read() to make it go faster
void FastaReference::load()
{
    inputMode();
    char base;
    bool newContig = false;
    eagle::model::Contig contig;
    unsigned long position = 0L;
    unsigned long prevPos = position;
    if (! reader_.seekg( 0, std::ios_base::beg ))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to rewind FASTA file %s") % reader_.file()).str()));
    }
    boost::filesystem::path lastPath = reader_.file();
    while (reader_.get(base, newContig))
    {
        if (newContig)
        {
            if (contig.size())
            {
                reference_.push_back( contig );
                std::clog << "....loaded '" << contig.name() << "' in memory" << std::endl;
                reader_.index().update( lastPath,
                                        eagle::io::FastaInfo(contig.id(),contig.size(),prevPos) );
                lastPath = reader_.file();
                contig.reset();
            }
            std::clog << "..Contig #" << reader_.getLocalContigId() << " begins at position " << position << std::endl;
            prevPos = position;
            contig.name( reader_.getContigName() );
        }
        contig.put( base );
        ++position;
    }
    if (contig.size())
    {
        reference_.push_back(contig);
        std::clog << "....loaded '" << contig.name() << "' in memory" << std::endl;
        eagle::io::FastaInfo info;
        reader_.index().update( lastPath,
                                eagle::io::FastaInfo(contig.id(),contig.size(),prevPos) );
    }
#ifdef EAGLE_DEBUG_MODE
    std::cout << "Merged FASTA Metadata:" << std::endl
              << reader_.index() << std::endl;
#endif
}


void FastaReference::save()
{
    outputMode();
    ReferenceIterator contig(reference_.begin());
    unsigned int contigCount(0);
    std::pair<unsigned long, unsigned long> location(1,0);  // (logical,phisical)
    while(reference_.end() != contig)
    {
        std::clog << "..Writing contig (#" << contigCount << "): '" << contig->name() << "'" << std::endl;
        location.second += contig->name().size() + 2;
        writer_.write( *contig,
                       eagle::io::FastaInfo( contig->id(), contig->size(), location.first, location.second, FASTA_CONTIG_WIDTH ));
        // TODO: overload FastaInfo::operator+() to do this when chromosome name hasn't changed
        location.first  += contig->size();
        location.second = 0; // This was used when we wanted to dump multiple contigs in one file:
                             // location.second += 1 + contig->size() + contig->size() / FASTA_CONTIG_WIDTH;
        ++contig;
        ++contigCount;
    }
    metadata_ = writer_.index();  // transfer metadata, so that *.fai info can be replicated in genome_size.xml
}


size_t FastaReference::length()
{
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::bind;
    return std::accumulate(
        reference_.begin(), reference_.end(),
        size_t(0), bind<size_t>(std::plus<size_t>(), _1, bind(&eagle::model::Contig::size, _2)));
}


std::vector<std::string> FastaReference::allContigNames() const
{
    using boost::lambda::_1;
    using boost::lambda::bind;
    std::vector<std::string> contigNames;
    if (reader_.size())
    {
        for (eagle::io::FastaMetadata::const_iterator idx = reader_.index().begin();
             idx != reader_.index().end(); ++idx)
        {
            std::transform( idx->second.begin(), idx->second.end(), std::back_inserter(contigNames),
                            bind( &eagle::io::FastaInfo::contigName, _1 ));
        }
    } else {
        BOOST_FOREACH(const eagle::model::Contig& contig, reference_)
        {
            contigNames.push_back(contig.id());
        }
    }
    return contigNames;
}

std::vector<unsigned long> FastaReference::allContigLengths() const
{
    using boost::lambda::_1;
    using boost::lambda::bind;
    std::vector<unsigned long> contigLengths;
    if (reader_.size())
    {
        for (eagle::io::FastaMetadata::const_iterator idx = reader_.index().begin();
             idx != reader_.index().end(); ++idx)
        {
            std::transform( idx->second.begin(), idx->second.end(), std::back_inserter(contigLengths),
                            bind( &eagle::io::FastaInfo::contigSize, _1 ));
        }
    } else {
        BOOST_FOREACH(const eagle::model::Contig& contig, reference_)
        {
            contigLengths.push_back(contig.size());
        }
    }
    return contigLengths;
}

unsigned long FastaReference::getContigLength(const std::string &contigName) const
{
    const std::vector<std::string> contigNames = allContigNames();
    const std::vector<unsigned long> contigLengths = allContigLengths();
    for (unsigned int i = 0; i < contigNames.size(); ++i)
    {
        if (contigName == contigNames[i])
        {
            assert( i < contigLengths.size() );
            return contigLengths[i];
        }
    }
    BOOST_THROW_EXCEPTION(eagle::common::EagleException(0, (boost::format("Contig '%s' not found") % contigName).str()));
}

eagle::model::Contig & FastaReference::getContig(const std::string &contigName)
{
    iterator contig(reference_.begin());
    for (; reference_.end() != contig; ++contig)
    {
        if (contig->name() == contigName)
            return *contig;
    }
    BOOST_THROW_EXCEPTION(eagle::common::EagleException(0, (boost::format("Contig '%s' not found") % contigName).str()));
}

eagle::model::Contig const & FastaReference::getContig(const std::string &contigName) const
{
    ReferenceIterator contig(reference_.begin());
    for (; reference_.end() != contig; ++contig)
    {
        if (contig->name() == contigName)
            return *contig;
    }
    BOOST_THROW_EXCEPTION(eagle::common::EagleException(0, (boost::format("Contig '%s' not found") % contigName).str()));
}


/***************************/
/**  MultiFastaReference  **/
/***************************/

eagle::io::FastaMetadata MultiFastaReference::initialize( const std::vector<boost::filesystem::path>& inputPaths )
{
    eagle::io::FastaMetadata metadata;
    BOOST_FOREACH(const boost::filesystem::path &p, inputPaths)
    {
        if (boost::filesystem::is_directory(p))
        {
            GenomeSizeXml GS(p / GENOMESIZE_XML);
            if (GS.enabled())
            {
                GS.load( metadata );
            } else {
                eagle::common::Glob FS(".*\\.fa(sta)?$");
                metadata.init( FS.glob(p) );
            }
        } else {
            metadata.update( p, eagle::io::FastaInfo() );
        }
    }
#ifdef EAGLE_DEBUG_MODE
    std::cout << "Genome Metadata:" << std::endl
              << metadata << std::endl;
#endif
    return metadata;
}

eagle::io::FastaMetadata MultiFastaReference::initialize( const boost::filesystem::path& inputPath )
{
    return initialize( std::vector< boost::filesystem::path >( 1, inputPath ) );
}

void MultiFastaReference::saveMetadata()
{
    if (!metadata().empty())
    {   // Hack: we are writing genome_size.xml into the parent directory of the 1st chromosome
        GenomeSizeXml GS( metadata().begin()->first.parent_path() / GENOMESIZE_XML, overwrite_ );
        GS.save( metadata() );
    } else {
        EAGLE_WARNING("No metadata available!");
        EAGLE_WARNING_CONT("*** Will not write " GENOMESIZE_XML " ***");
    }
}



/*****************/
/**  METADATA   **/
/*****************/


GenomeSizeXml::GenomeSizeXml(
    const boost::filesystem::path& indexPath
    )
    : path_(indexPath)
    , mode_(std::ios_base::in)
{
    if (!path_.empty() && !boost::filesystem::exists(path_))
    {
        EAGLE_WARNING( (boost::format("Cannot read FASTA metadata from %s") % path_).str() );
        path_ = boost::filesystem::path();
    }
}
GenomeSizeXml::GenomeSizeXml(
    const boost::filesystem::path& indexPath,
    const bool overwrite
    )
    : path_(indexPath)
    , mode_(std::ios_base::out)
{
    if ( path_.empty() ) {
        BOOST_THROW_EXCEPTION(common::PreConditionException( "Path to '" GENOMESIZE_XML "' cannot be empty" ));
    } else if ( boost::filesystem::exists(path_) ) {
        if (overwrite)
        {
            EAGLE_WARNING( "Overwriting " << path_ << " due to the --force switch." );
        } else {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Cannot write FASTA metadata file %s: File already exists!") % path_).str()));
        }
    }
}

void GenomeSizeXml::load(eagle::io::FastaMetadata &metadata)
{
    if (! mode_ & std::ios_base::in)
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException( (boost::format("Not allowed to read from %s") % path_ ).str() ));
    }
    if (boost::filesystem::exists(path_))
    {
        std::ifstream is(path_.string().c_str());
        if (!is)
        {
            BOOST_THROW_EXCEPTION(common::IoException(errno, (boost::format("Failed to open %s for reading") % path_).str()));
        }
        this->load(is,metadata);
    } else {
        if (! path_.empty() )
        {
            EAGLE_WARNING( "Could not find " << path_ );
            EAGLE_WARNING_CONT( "*** Will not pre-load FASTA metadata ***" );
        }
    }
}

void GenomeSizeXml::load(std::istream &is, eagle::io::FastaMetadata &metadata)
{
    using namespace boost::property_tree;
    ptree tree;
    read_xml< ptree >(is, tree);
    unsigned long absolutePos(0);
    EAGLE_DEBUG(0, "Reading metadata from " << path_ );
    if( boost::optional< ptree & > sequenceSizes = tree.get_child_optional( "sequenceSizes" ) )
    {
        BOOST_FOREACH(const ptree::value_type& chromosome, sequenceSizes.get())
        {
            if (chromosome.first != "chromosome")
            {
                BOOST_THROW_EXCEPTION(common::CorruptedFileException("XML.GenomeSize",
                                     (boost::format("*** Expected <chromosome/> element in %s. Found <%s/> ***") % path_ % chromosome.first).str() ));
            }
            boost::filesystem::path fileName = chromosome.second.get<boost::filesystem::path>("<xmlattr>.fileName");
            std::string chrName = chromosome.second.get<std::string>("<xmlattr>.contigName");
            unsigned long chrSize = chromosome.second.get<unsigned long>("<xmlattr>.totalBases");
            metadata.update( path_.parent_path() / fileName,
                             eagle::io::FastaInfo( chrName, chrSize, absolutePos ) );
            absolutePos += chrSize;
        }
    } else if (boost::optional< ptree & > sequenceSizesLegacy = tree.get_child_optional( "SequenceSizes" )) {
        BOOST_FOREACH(const ptree::value_type& chromosome, sequenceSizesLegacy.get())
        {
            boost::filesystem::path fileName = chromosome.first;
            std::string chrName = fileName.stem().string();
            unsigned long chrSize = boost::lexical_cast<unsigned long>( chromosome.second.data() );
            metadata.update( path_.parent_path() / fileName,
                             eagle::io::FastaInfo( chrName, chrSize, absolutePos ) );
            absolutePos += chrSize;
        }
    } else {
        BOOST_THROW_EXCEPTION(common::CorruptedFileException("XML.GenomeSize","*** either <sequenceSizes/> or <SequenceSizes/> missing at the top level ***"));
    }
}

void GenomeSizeXml::save(const eagle::io::FastaMetadata &index)
{
    if (! mode_ & std::ios_base::out)
    {
        BOOST_THROW_EXCEPTION(common::PreConditionException( (boost::format("Not allowed to write to %s") % path_ ).str() ));
    }
    this->save(path_,index);
}

void GenomeSizeXml::save(const boost::filesystem::path &file, const eagle::io::FastaMetadata &index)
{
    using namespace boost::property_tree;
    ptree tree;
    tree.clear();
    tree.put("sequenceSizes", "");
    ptree &sequenceSizes = tree.get_child("sequenceSizes");
    EAGLE_DEBUG(0, "Writing metadata to " << path_ );
    for ( eagle::io::FastaMetadata::const_iterator chromosome = index.begin();
          chromosome != index.end();
          ++chromosome)
    {
        for (unsigned int i=0; i < chromosome->second.size(); i++)
        {
            ptree node;
            node.put( "<xmlattr>.fileName",   chromosome->first.filename().string() );
            node.put( "<xmlattr>.contigName", chromosome->second[i].contigName );
            node.put( "<xmlattr>.totalBases", chromosome->second[i].contigSize );
            sequenceSizes.add_child("chromosome",node);
        }
    }
    write_xml(file.string(), tree, std::locale(), xml_writer_settings<char>(' ', 4));
}



} // namespace genome
} // namespace eagle
