#ifndef ADD_TRANSLOCATIONS_H
#define ADD_TRANSLOCATIONS_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <seqan/file.h>
#include <seqan/random.h>

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/refinement.h> // includes interval tree and related tools

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "VcfVariant.h"

using namespace seqan;

typedef String<Dna5>                 TSequence;
typedef IntervalTree<int,CharString> TIntervalTree;
typedef std::map<std::string,TIntervalTree> TIntervalTreeMap;

//#define DEBUG

struct Options {

  bool showHelp;
  bool showVersion;
  unsigned verbosity;

  CharString gapFile;
  CharString refFile;
  CharString outputFile;
	CharString vcfInfile;

	CharString sampleId;

  unsigned numTranslocations;
  unsigned randomSeed;

  Options() {
    showHelp          = false;
    showVersion       = false;
    verbosity         = 1;
    numTranslocations = 500;
    randomSeed        = 0;
    sampleId          = "NA12878";
  }
};

#endif // ADD_TRANSLOCATIONS_H
