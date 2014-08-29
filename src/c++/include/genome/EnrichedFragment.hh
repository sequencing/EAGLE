/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Set of classes to represent DNA fragments enriched with sequencing primesrs, attachments, etc.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_GENOME_ENRICHED_FRAGMENT_HH
#define EAGLE_GENOME_ENRICHED_FRAGMENT_HH

#include <utility>
#include <iostream>
#include <boost/format.hpp>
#include "model/Fragment.hh"

using namespace std;


namespace eagle
{
namespace genome
{


struct FragmentComponent
{
    FragmentComponent() {}
    virtual ~FragmentComponent() {}//assert(false); } // should never deconstruct for the moment
    virtual char getBase( const unsigned int posInRead, const eagle::model::Fragment &fragment, const bool isForward ) const = 0;  // { assert(false); }

    unsigned int length_;
};

struct FragmentComponentHardcoded : public FragmentComponent
{
    FragmentComponentHardcoded( string bases);
    virtual ~FragmentComponentHardcoded() {}
    virtual char getBase( const unsigned int posInRead, const eagle::model::Fragment &fragment, const bool isForward ) const;

private:
    string bases_;
}; 

struct FragmentComponentRealDna : public FragmentComponent
{
    FragmentComponentRealDna();
    virtual ~FragmentComponentRealDna() {}
    virtual char getBase( const unsigned int posInRead, const eagle::model::Fragment &fragment, const bool isForward ) const;

private:

    //    bool directionIsForward_;
    //    signed int startPosInFragment_;
}; 

class FragmentStructure
{
public:
    vector< boost::shared_ptr< FragmentComponent > > components_;
    // TODO: destructor to free pointers
    vector< pair< unsigned int, bool > > reads_;

    char getBase( const unsigned int readNum, const unsigned int posInRead, const eagle::model::Fragment &fragment ) const;
    unsigned int getReadLength( const unsigned int readNum, const eagle::model::Fragment &fragment ) const;
    bool getReadInfo( const unsigned int readNum, /*bool &isIndex, unsigned long &startPos, */bool &directionIsForward ) const
        {
            if ( readNum < reads_.size() )
            {
//                isIndex = components_[reads_[readNum].first].isIndex_;
//                startPos = eFragment_.fragment_.startPos + components_[reads_[readNum].first];
                directionIsForward = reads_[readNum].second;
                return true;
            }
            return false;
        }
};

class FragmentStructureType1 : public FragmentStructure
{
public:
    FragmentStructureType1()
    {
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentRealDna( ) ) );

        reads_.push_back( make_pair(0,true) );
    }
};

class FragmentStructureType2Generic : public FragmentStructure
{
public:
    FragmentStructureType2Generic()
    {
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "AttachmentP5" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "SeqPrimer1" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentRealDna( ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "sEQpRIMER2" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "BarcadaCADA" ) ) ); // "Barcode" but 'o' and 'e' are invalid base values... and some extra values in case a deletion occurs
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "aTTACHMENTp7" ) ) );
    }

    void addRead( const unsigned int readNum ) { reads_.push_back( make_pair(2,readNum==1) ); }
    void addBarcode() { reads_.push_back( make_pair(4,true) ); }
};

class FragmentStructureType2GenericReverse : public FragmentStructure
{
public:
    FragmentStructureType2GenericReverse()
    {
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "AttachmentP5" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "SeqPrimer1" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentRealDna( ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "sEQpRIMER2" ) ) );
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "BarcadaCADA" ) ) ); // "Barcode" but 'o' and 'e' are invalid base values... and some extra values in case a deletion occurs
        components_.push_back( boost::shared_ptr<FragmentComponent>( new FragmentComponentHardcoded( "aTTACHMENTp7" ) ) );
    }

    void addRead( const unsigned int readNum ) { reads_.push_back( make_pair(2,readNum==2) ); }
    void addBarcode() { reads_.push_back( make_pair(4,true) ); }
};

class EnrichedFragment
{
public:
    EnrichedFragment( const eagle::model::Fragment &fragment, const vector<FragmentStructure> &multiplexedFragmentStructures, const unsigned int dnaFragmentDirection );
    char getBase( unsigned int read, unsigned int posInRead ) const;
    unsigned int getReadCount() const;
    unsigned int getReadLength( unsigned int r ) const;

    //unsigned long startPos_, fragmentLength_;
    const eagle::model::Fragment &fragment_;
    unsigned int dnaFragmentDirection_;
    const FragmentStructure &structure_;
};


} // namespace genome
} // namespace eagle

#endif // EAGLE_GENOME_ENRICHED_FRAGMENT_HH
