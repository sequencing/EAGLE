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

#ifndef EAGLE_IO_FASTA_HH
#define EAGLE_IO_FASTA_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>

#include "model/Contig.hh"
#include "io/Text.hh"


// NCBI suggests using width = 80
//      but I am using 70 due to the E_coli reference
//      used for unit/regression testing
#define FASTA_CONTIG_WIDTH 70

// TODO: improve how we handle these constants
#define GENOMESIZE_XML "genome_size.xml"


namespace eagle
{
namespace io
{


/*
 * \item when using numeric constructor, assume 1 byte EOL
 * \item not possible to initialize logical global location (set to 0xffffffffffffffffL) when using string constructor
 */
struct FastaInfo
{
    FastaInfo( std::string   name = "",
               unsigned long size = 0L,
               unsigned long global = 0xffffffffffffffffL,
               unsigned long absolute = 0xffffffffffffffffL,
               unsigned int  width = 0)
        : contigName(   name )
        , contigSize(   size )
        , position(     std::make_pair( global, absolute ))
        , contigWidth(  std::make_pair( width,  width+1 ))
        , contigNumber( 0 )
        {}
    FastaInfo( std::string name,
               std::string size,
               std::string absolute,
               std::string basesPerLine,
               std::string bytesPerLine)
        : contigName(   name )
        , contigSize(   boost::lexical_cast<unsigned long>(size) )
        , position(     std::make_pair( 0xffffffffffffffffL, boost::lexical_cast<unsigned long>(absolute) ))
        , contigWidth(  std::make_pair( boost::lexical_cast<unsigned int>(basesPerLine),
                                        boost::lexical_cast<unsigned int>(bytesPerLine) ))
        , contigNumber( 0 )
        {}  // This constructor facilitates the reading from *.fai
    void setName( const std::string& name )    {if ("" != name) contigName = name;}
    void setSize( unsigned long size )         {if (0L != size) contigSize = size;}
    void setLogicalPosition( unsigned long global )    {if (0xffffffffffffffffL != global)   position.first = global;}
    void setPhysicalPosition( unsigned long absolute ) {if (0xffffffffffffffffL != absolute) position.second = absolute;}
    void setPosition( unsigned long global,
                      unsigned long absolute ) {position = std::make_pair( (global   ? global   : this->position.first),
                                                                           (absolute ? absolute : this->position.second) );}
    void setWidth( unsigned int bases,
                   unsigned int bytes ) {contigWidth = std::make_pair( (bases ? bases : this->contigWidth.first),
                                                                       (bytes ? bytes : this->contigWidth.second) );}

    bool within(unsigned long pos) const {return (position.first <= pos) && (pos < (position.first + contigSize));}
    bool sameName( const FastaInfo& rhs ) const     {return this->contigName  == rhs.contigName;}
    bool sameSize( const FastaInfo& rhs ) const     {return this->contigSize  == rhs.contigSize;}
    bool sameLogicalPosition( const FastaInfo& rhs ) const {return this->position.first == rhs.position.first;}
    bool samePhysicalPosition( const FastaInfo& rhs ) const {return this->position.second == rhs.position.second;}
    bool samePosition( const FastaInfo& rhs ) const {return sameLogicalPosition(rhs) && samePhysicalPosition(rhs);}
    bool sameWidth( const FastaInfo& rhs ) const    {return this->contigWidth == rhs.contigWidth;}
    bool hasName() const             {return "" != this->contigName;}
    bool hasSize() const             {return 0L != this->contigSize;}
    bool hasLogicalPosition() const  {return 0xffffffffffffffffL != this->position.first;}
    bool hasPhysicalPosition() const {return 0xffffffffffffffffL != this->position.second;}
    bool hasPosition() const         {return hasLogicalPosition() || hasPhysicalPosition();}
    bool hasWidth() const            {return std::pair< unsigned int, unsigned int >(0,1) != this->contigWidth;}
    static bool needsUpdating( FastaInfo self ) {return !self.hasName() && !self.hasPosition();}

    std::string contigName;
    unsigned long contigSize;
    std::pair< unsigned long, unsigned long > position;   // logical(across files), physical (within file)
    std::pair< unsigned int, unsigned int > contigWidth;  // basesPerLine,          bytesPerLine
    unsigned int contigNumber;
};

typedef std::pair< boost::filesystem::path, std::vector< FastaInfo > > FastaIndexBase;
class FastaIndex : public FastaIndexBase
{
public:
    FastaIndex( const std::pair< boost::filesystem::path, std::vector< FastaInfo > > &p ) : FastaIndexBase(p) {}
    bool operator==( const FastaIndex &rhs) const { return first == rhs.first; }
    bool operator==( const boost::filesystem::path &p) const { return first == p; }
};

typedef std::vector< FastaIndex >::iterator FastaIterator;
typedef std::vector< FastaIndex >::const_iterator FastaConstIterator;

class FastaMetadata : public std::vector< FastaIndex >
{
public:
    void init( const std::vector< boost::filesystem::path > &P );
    void update( const boost::filesystem::path &p, const FastaInfo &fi );

    FastaIterator find( const boost::filesystem::path &key )
    {
        FastaIterator it = std::find( begin(), end(), key );
        return it;
    }

    FastaConstIterator find( const boost::filesystem::path &key ) const
    {
        FastaConstIterator it = std::find( begin(), end(), key );
        return it;
    }

    std::pair< FastaIterator, bool > insert( const FastaIndex &fi )
    {
        std::pair< FastaIterator, bool > result;
        FastaIterator it = find( fi.first );
        if (it == end())
        {
            push_back( fi );
            result.first = end()-1;
            result.second = true;
        }
        else
        {
            result.first = it;
            result.second = false;
        }
        return result;
    }

    FastaIterator insert( const FastaIterator& it, const FastaIndex &fi ) 
    {
        insert( fi );
        return it;
    }

    std::vector< FastaInfo >& at( const boost::filesystem::path &key )
    {
        FastaIterator it = find( key );
        assert (it != end());
        return it->second;
    }

    std::vector< FastaInfo >& operator[]( const boost::filesystem::path &key )
    {
        return at(key);
    }
/*
    size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }

  private:
    std::vector< FastaIndex > data_;
*/
private:
    mutable boost::filesystem::path cachedKey;
    mutable FastaIterator cachedIterator;
};


std::ostream& operator<<( std::ostream& os, const FastaInfo& fi );
//std::istream& operator>>( std::istream& is, FastaInfo& fi );
std::ostream& operator<<( std::ostream& os, const FastaIndex& fi );
std::ostream& operator<<( std::ostream& os, const FastaMetadata& fm );

//namespace fasta
//{
//    static unsigned int ContigWidth;
//}


class FastaReader : public std::ifstream, boost::noncopyable
{
public:
    FastaReader() {}
    /// read a base and set the 'newContig' flag when starting a new contig
    FastaReader& get(char &base, bool &newContig);
    FastaReader& read(unsigned long pos, unsigned long length, unsigned long lineCount);
    std::string getContigName() const {return contigName_;}
    void setContigName(std::string contigName) {contigName_ = contigName;}
    std::vector<char> & cache()             {return cache_;}
    std::vector<char> const & cache() const {return cache_;}

private:
    std::string contigName_;
    std::vector<char> cache_;
};


class MultiFastaReader: public FastaReader
{
public:
    MultiFastaReader(const FastaMetadata &index = FastaMetadata());
    MultiFastaReader& get(char &base, bool &newContig);
    int getGlobalContigId() const {return globalContigId_;}
    int getLocalContigId() const {return localContigId_;}
    void setGlobalContigId(int globalContigId) {globalContigId_ = globalContigId;}
    void setLocalContigId(int localContigId) {localContigId_ = localContigId;}
    unsigned long getContigSize();
    size_t size() const {return index_.size();}
    std::vector< FastaInfo > & index(boost::filesystem::path fastaPath);
    std::vector< FastaInfo > const & index(boost::filesystem::path fastaPath) const;
    FastaMetadata & index()             {return index_;}
    FastaMetadata const & index() const {return index_;}
    boost::filesystem::path file() const {assert(!index_.empty() && "Cannot perform MultiFastaReader::file() of an empty object!");
                                          return this_->first;}
    boost::filesystem::path find(const std::string& chr, FastaInfo& info);
    boost::filesystem::path find(unsigned long pos,      FastaInfo& info);
    bool open(const boost::filesystem::path fastaPath);
    // read block of data into cache
    MultiFastaReader& read( const FastaInfo &info, unsigned long skip = 0, unsigned long size = 0);

    //direct access to data in cache
    char operator[](unsigned long i)
    {
        unsigned long posInContig = i - this_->second[localContigId_].position.first;
        posInContig %= this_->second[localContigId_].contigSize;
        unsigned long fullLinesCount = posInContig / this_->second[localContigId_].contigWidth.first;
        unsigned long posInLine = posInContig % this_->second[localContigId_].contigWidth.first;
        return cache()[fullLinesCount * this_->second[localContigId_].contigWidth.second + posInLine];
    }
    bool inCache(unsigned long i)
    {
        return (i - this_->second[localContigId_].position.first) < this_->second[localContigId_].contigSize;
    }
private:
    void open();
    FastaMetadata index_;
    FastaIterator this_;
    int globalContigId_;
    int localContigId_;
};


class FastaWriter: public std::ofstream, boost::noncopyable
{
public:
    FastaWriter(const boost::filesystem::path fastaPath, const unsigned int width=FASTA_CONTIG_WIDTH)
    : fastaPath_(fastaPath), contigWidth_(width) {open();}
    FastaWriter(const unsigned int width=FASTA_CONTIG_WIDTH)
    : contigWidth_(width) {}
    void open(const boost::filesystem::path fastaPath);
    FastaWriter& write(const eagle::model::Contig& contig);
    size_t size() const {return 1;}
    boost::filesystem::path file() const {return fastaPath_;}
private:
    void open();
    boost::filesystem::path fastaPath_;
    unsigned int contigWidth_;
};

class MultiFastaWriter: public FastaWriter
{
public:
    MultiFastaWriter(const boost::filesystem::path &fastaDir = "", const bool overwrite = false)
    : FastaWriter(), fastaDir_(fastaDir), overwrite_(overwrite) {}
    MultiFastaWriter& write(const eagle::model::Contig& contig, const FastaInfo& info);
    FastaMetadata & index()             {return index_;}
    FastaMetadata const & index() const {return index_;}
private:
    const boost::filesystem::path fastaDir_;
    const bool overwrite_;
    FastaMetadata index_;
};

/*
 * to generate a *.fai file:
 *    > samtools faidx Homo_sapiens_assembly18.fasta
 *    108.446u 3.384s 2:44.61 67.9%   0+0k 0+0io 0pf+0w
 */
class FaiReader: public DsvReader
{
public:
    FaiReader(const std::vector<boost::filesystem::path> &faiPathList)
    : DsvReader(faiPathList) {}
    FaiReader() {}
    bool getNextIndex(FastaInfo &index);
};

class FaiWriter: public DsvWriter
{
public:
    FaiWriter(const std::vector<boost::filesystem::path> &faiPathList, const bool overwrite)
    : DsvWriter(faiPathList,overwrite) {}
    FaiWriter(const boost::filesystem::path faiPath, const bool overwrite)
    : DsvWriter(faiPath,overwrite)     {}
    FaiWriter(const bool overwrite)
    : DsvWriter(overwrite)             {}
    FaiWriter& write(const FastaInfo& info);
};


} // namespace io
} // namespace eagle

#endif // EAGLE_IO_FASTA_HH
