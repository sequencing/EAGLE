#include "AddTranslocations.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;

using boost::lexical_cast;

const unsigned MAX_CHR_NUM  = 25; // 22 chrs + X + Y + mito

size_t mapChr2Idx(const string& chr)
{
	if (chr == "chr1") {
		return 1;
	} else if (chr == "chr2") {
		return 2;
	} else if (chr == "chr3") {
		return 3;
	} else if (chr == "chr4") {
		return 4;
	} else if (chr == "chr5") {
		return 5;
	} else if (chr == "chr6") {
		return 6;
	} else if (chr == "chr7") {
		return 7;
	} else if (chr == "chr8") {
		return 8;
	} else if (chr == "chr9") {
		return 9;
	} else if (chr == "chr10") {
		return 10;
	} else if (chr == "chr11") {
		return 11;
	} else if (chr == "chr12") {
		return 12;
	} else if (chr == "chr13") {
		return 13;
	} else if (chr == "chr14") {
		return 14;
	} else if (chr == "chr15") {
		return 15;
	} else if (chr == "chr16") {
		return 16;
	} else if (chr == "chr17") {
		return 17;
	} else if (chr == "chr18") {
		return 18;
	} else if (chr == "chr19") {
		return 19;
	} else if (chr == "chr20") {
		return 20;
	} else if (chr == "chr21") {
		return 21;
	} else if (chr == "chr22") {
		return 22;
	} else if (chr == "chrX") {
		return 23;
	} else if (chr == "chrY") {
		return 24;
	} else if (chr == "chrM") {
		// don't think that we need this but nevertheless
		return 25;
	}	else {	
		cout <<"Unknown chromosome : " << chr << endl;
		return(EXIT_FAILURE);
	}
}

void setupCommandLineParser(CommandLineParser& parser, Options const & options) {

    addTitleLine(parser,"AddTranslocations - add translocations to an existing VCF file");
    addVersionLine(parser, "1.0");

    addSection(parser, "General Options");
    addOption(parser, CommandLineOption("v", "verbose", "Enable verbose mode (show steps).", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "Enable very verbose mode.", OptionType::Bool));
    addOption(parser, CommandLineOption("", "random-seed", "Random seed", OptionType::Int, options.randomSeed));

    addSection(parser, "Input Specification");
    addOption(parser, CommandLineOption("g", "gaps-file", "Path to UCSC gaps file", OptionType::String | OptionType::Mandatory));
    addOption(parser, CommandLineOption("r", "ref-file", "Path to reference sequence file", OptionType::String | OptionType::Mandatory));
    addOption(parser, CommandLineOption("c", "vcf-file", "Path to input vcf file", OptionType::String | OptionType::Mandatory));
    addOption(parser, CommandLineOption("o", "output-file", "Path to output file", OptionType::String | OptionType::Mandatory));

    addSection(parser, "Sampling options");
    addOption(parser, CommandLineOption("t", "num-translocations", "Number of translocations to be simulated", OptionType::Int, options.numTranslocations));
	
	// does not seem to work with this one switched on
    //requiredArguments(parser, 3);
}

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv) {
	bool stop = !parse(parser, argc, argv);
  if (stop && ! isSetLong(parser, "help")) {
		printHelp(parser,cerr);
  	return 1;
	}
  if (isSetLong(parser, "help")) {
  	options.showHelp = true;
    return 0;
  }
  if (isSetLong(parser, "version")) {
  	options.showVersion = true;
    return 0;
  }

	if (isSetLong(parser, "verbose"))
  	options.verbosity = 2;
  if (isSetLong(parser, "very-verbose"))
  	options.verbosity = 3;

  getOptionValueLong(parser, "gaps-file", options.gapFile);
  getOptionValueLong(parser, "ref-file", options.refFile);
  getOptionValueLong(parser, "vcf-file", options.vcfInfile);
  getOptionValueLong(parser, "output-file", options.outputFile);
  getOptionValueLong(parser, "num-translocations", options.numTranslocations);
  getOptionValueLong(parser, "random-seed", options.randomSeed);
  //options.gapFile = getArgumentValue(parser,0);
  //options.refFile = getArgumentValue(parser,1);

	return 0;
}

int loadRefSeqs(MultiSeqFile& multiSeqFile, 
								Options& options, 
								StringSet<int>& chrLengths, 
								StringSet<CharString>& chrIDs, 
								StringSet<TSequence>& chrSeqs
                                ) {

    if (!open(multiSeqFile.concat, toCString(options.refFile), OPEN_RDONLY)) {
  	    cout << "Cannot read from " << options.refFile << endl;
        return(-1);
    }

	if (options.verbosity >= 2) {
		cout << "Reading from " << options.refFile << endl;
	}
    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    unsigned seqCount = length(multiSeqFile);
	if (options.verbosity >= 2) {
  	    cout << "Found " << seqCount << " sequences. " << endl;
	}
  
	reserve(chrSeqs, seqCount, Exact());
    reserve(chrIDs, seqCount, Exact());
    TSequence seq;
    string id;

    cout << "Reading reference sequences " << endl;
    for (unsigned i=0;i<seqCount;++i) {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id
		if (options.verbosity >= 2) {
  		    cout << "id = " << id << endl;
        }
	    unsigned len = length(seq);
        appendValue(chrSeqs, seq, Generous());
        appendValue(chrIDs, id, Generous());
  	    appendValue(chrLengths,len,Generous());
    }
	return 0;
}

/*
TSequence generateRandomInsertionSequence(unsigned& eventLen) {
	TSequence seq;
	for(unsigned i=0;i<eventLen;++i) {
		//Rng<MersenneTwister> rng(SEED);
		Rng<MersenneTwister> rng; <- use options.randomSeed
  	Pdf<Uniform<int> > nucPdf(0,3);
		int rnd = pickRandomNumber(rng, nucPdf);
		switch (rnd) {
			case 0:
				appendValue(seq,'A');
				break;
			case 1:
				appendValue(seq,'C');
				break;
			case 2:
				appendValue(seq,'G');
				break;
			case 3:
				appendValue(seq,'T');
				break;
			default:
				cout << "Random number out of range: " << rnd << endl;
		}
	}
	assert(length(seq)==eventLen);
	return seq;
}
*/

//#CHROM POS    ID     REF ALT           QUAL FILT INFO
//2      321681 bnd_W  G   G[13:123457[  6    PASS SVTYPE=BND;MATEID=bnd_X;EVENT=RR0
//2      321682 bnd_V  T   ]13:123456]T  6    PASS SVTYPE=BND;MATEID=bnd_U;EVENT=RR0
//13     123456 bnd_U  C   C[2:321682[   6    PASS SVTYPE=BND;MATEID=bnd_V;EVENT=RR0
//13     123457 bnd_X  A   ]2:321681]A   6    PASS SVTYPE=BND;MATEID=bnd_W;EVENT=RR0

void simulateTranslocations(vcfStore& vcfLines, 
                            TIntervalTreeMap iForestGaps, 
                            Options& options, 
                            const StringSet<int>& chrLengths,
                            const StringSet<CharString>& chrIDs,
                            const StringSet<TSequence>& chrSeqs,
                            const int refCount)
{
	Rng<MersenneTwister> rng(options.randomSeed);
    Pdf<Uniform<int> > chrPdf(0, (refCount-1));

	int transLocLen(1);
	unsigned cnt(0);

	while (cnt<options.numTranslocations) {
        int chr1Idx   = pickRandomNumber(rng,chrPdf);
        int chr2Idx   = pickRandomNumber(rng,chrPdf);

        // we don't want translocations between the same chromosome
        if (chr1Idx==chr2Idx) continue;

        string chr1   = toCString(chrIDs[chr1Idx]);
        string chr2   = toCString(chrIDs[chr2Idx]);
            
        int max_pos1  = (transLocLen > chrLengths[chr1Idx] ) ? 0 : chrLengths[chr1Idx] - transLocLen;
        int max_pos2  = (transLocLen > chrLengths[chr2Idx] ) ? 0 : chrLengths[chr2Idx] - transLocLen;
        //cout << "max_pos = " << max_pos << " " << length(chrSeqs[chrIdx]) << endl;
        if (max_pos1==0 || max_pos2==0) continue;
        Pdf<Uniform<int> > pos1Pdf(0,max_pos1);
        Pdf<Uniform<int> > pos2Pdf(0,max_pos2);

        int eventStartPos1 = pickRandomNumber(rng,pos1Pdf);
        int eventStartPos2 = pickRandomNumber(rng,pos2Pdf);
        cout << "chrIdx=" << chr1Idx << " chr=" << chr1 << " pos=" << eventStartPos1 << " chrLen=" << length(chrSeqs[chr1Idx]) << endl;
        cout << "chrIdx=" << chr2Idx << " chr=" << chr2 << " pos=" << eventStartPos2 << " chrLen=" << length(chrSeqs[chr2Idx]) << endl;
        String<CharString> results;
        findIntervals(iForestGaps[chr1],eventStartPos1,(eventStartPos1+1),results); 
        if (length(results) != 0)
        {
            cerr << "First breakpoint overlaps with an existing variant, skipping it." << endl;
            continue;
        }
        findIntervals(iForestGaps[chr2],eventStartPos2,(eventStartPos2+1),results); 
        if (length(results) != 0)
        {
            cerr << "Second breakpoint overlaps with an existing variant, skipping it." << endl;
            continue;
        }
		
		// we need the reference nucleotides at each of the 4 breakpoint ends
		string bndWRefNuc(toCString(CharString(infix(chrSeqs[chr1Idx],eventStartPos1,(eventStartPos1+1)))));
		string bndVRefNuc(toCString(CharString(infix(chrSeqs[chr1Idx],(eventStartPos1+1),(eventStartPos1+2)))));
		string bndURefNuc(toCString(CharString(infix(chrSeqs[chr2Idx],eventStartPos2,(eventStartPos2+1)))));
		string bndXRefNuc(toCString(CharString(infix(chrSeqs[chr2Idx],(eventStartPos2+1),(eventStartPos2+2)))));

		string bndWAlt(bndWRefNuc + "[" + chr2 + ":" + lexical_cast<string>(eventStartPos2+1) + "[");
		string bndVAlt("]" + chr2 + ":" + lexical_cast<string>(eventStartPos2) + "]" + bndVRefNuc);
		string bndUAlt(bndURefNuc + "[" + chr1 + ":" + lexical_cast<string>(eventStartPos1+1) + "[");
		string bndXAlt("]" + chr1 + ":" + lexical_cast<string>(eventStartPos1) + "]" + bndXRefNuc);

		string bndWinfo(string("SVTYPE=BND;MATEID=bnd_X;EVENT=RR") + lexical_cast<string>(cnt));
		string bndVinfo(string("SVTYPE=BND;MATEID=bnd_U;EVENT=RR") + lexical_cast<string>(cnt));
		string bndUinfo(string("SVTYPE=BND;MATEID=bnd_V;EVENT=RR") + lexical_cast<string>(cnt));
		string bndXinfo(string("SVTYPE=BND;MATEID=bnd_W;EVENT=RR") + lexical_cast<string>(cnt));
	
		// create vcf record for each breakpoint
		VcfRecord bndWrec(chr1,eventStartPos1,transLocLen,("bnd_W_" + lexical_cast<string>(cnt)),bndWRefNuc,bndWAlt,30,"PASS",bndWinfo,"GT","1/0");
		VcfRecord bndVrec(chr1,(eventStartPos1+1),transLocLen,("bnd_V_" + lexical_cast<string>(cnt)),bndVRefNuc,bndVAlt,30,"PASS",bndVinfo,"GT","1/0");
		VcfRecord bndUrec(chr2,(eventStartPos2),transLocLen,("bnd_U_" + lexical_cast<string>(cnt)),bndURefNuc,bndUAlt,30,"PASS",bndUinfo,"GT","1/0");
		VcfRecord bndXrec(chr2,(eventStartPos2+1),transLocLen,("bnd_X_" + lexical_cast<string>(cnt)),bndXRefNuc,bndXAlt,30,"PASS",bndXinfo,"GT","1/0");

		vcfLines.push_back(bndWrec);
		vcfLines.push_back(bndVrec);
		vcfLines.push_back(bndUrec);
		vcfLines.push_back(bndXrec);

		// add translocation breakpoints to interval tree for this chromosome
		//addInterval(iForestVars[chrIdx],eventStartPos,eventEndPos,rec.id);
		++cnt; // count simulated translocations
	}
}

bool hasOverlap(VcfRecord& rec, TIntervalTreeMap iForestGaps, const bool lastBndPassed) {

	if (rec.info.find("bnd_V") != string::npos ||	
			rec.info.find("bnd_U") != string::npos ||
			rec.info.find("bnd_X") != string::npos ) {
		
		if (lastBndPassed) {
			return true;
		} else {
			cout << "Won't check " << rec.info << " for overlaps." << endl;
			return false;
		}
	}

	String<CharString> results;
	//size_t chrIdx = mapChr2Idx(rec.chr);
	findIntervals(iForestGaps[rec.chr],rec.pos,(rec.pos+rec.len),results); 
	if (length(results) != 0) {
		//cerr << "This variant overlaps with a UCSC gap or an existing variant, Removing it." << endl;
		return true;
	}
	return false;
}


int main(int argc, char const** argv) {

	// Setup command line parser.
    CommandLineParser clParser;
    Options options;
    setupCommandLineParser(clParser, options);
    VcfParser vcfParser;
    
    int res = parseCommandLineAndCheck(options, clParser, argc, argv);

  if (res != 0)
  	return (EXIT_FAILURE);
  if (options.showHelp || options.showVersion)
  	return (EXIT_SUCCESS);

	// allocating one interval tree per chromosome to check
	// for overlaps with UCSC gaps
	TIntervalTreeMap iForestGaps;
	if (options.verbosity >= 2) {
		cout << "Reading UCSC gaps from " << options.gapFile << endl;
	}

	//ifstream gapStream(toCString(options.gapFile),std::ios::in);
	ifstream gapStream(toCString(options.gapFile));
	if (!gapStream.is_open()) {
		cout << "Cannot read from " << options.gapFile << endl;
		return(EXIT_FAILURE);
	}

	string line;
  string chr;
  int start, end;
	unsigned ct = 0;	

    while(getline(gapStream,line)) {
        stringstream sstr(line);
        sstr >> chr >> start >> end;
        size_t chrIdx = mapChr2Idx(chr);
        cout << "Adding gap " << chr << " " << chrIdx << " " << start << " " << end << endl;
        addInterval(iForestGaps[chr],start,end);
        ++ct;
    }
    if (options.verbosity >= 2) {
        cout << "Read " << ct << " UCSC gaps. " << endl;
    }

	// loading reference genome
    MultiSeqFile multiSeqFile;
    StringSet<int> chrLengths;
    StringSet<CharString> chrIDs;
    StringSet<TSequence > chrSeqs;
	if (loadRefSeqs(multiSeqFile,options,chrLengths,chrIDs,chrSeqs) != 0 ) {
		cerr << "Cannot read from reference file " << options.refFile << endl;
		return(EXIT_FAILURE);
	}
    int refCount = length(chrSeqs);

	// allocating interval trees for simulated variants
	//TIntervalTree iForestVars[MAX_CHR_NUM];
	vcfStore vcfLines;	
	vector< string > vcfHeader;
	cout << "Reading existing variants from " << options.vcfInfile << endl;
	ifstream vcfInstream(toCString(options.vcfInfile));
	if (!vcfInstream.is_open()) {
		cout << "Cannot read from " << options.vcfInfile << endl;
	}
	unsigned varCount(0);
	unsigned overlapCount(0);
	bool lastBndPassed(false);

	while(getline(vcfInstream,line) && !line.empty()) {
		if (line[0] == '#') {
			//cout << line << endl;
			vcfHeader.push_back(line);
			continue;
		}
		VcfRecord r;
		vcfParser.buildVcfRecordFromString(line,r);
		++varCount;
		//cout << "read vcf=" << endl;
		//cout << r << endl;
		// only include into list if no overlap with UCSC gaps or existing variants
		bool ovl = hasOverlap(r,iForestGaps,lastBndPassed);

		if (!ovl) {
			vcfLines.push_back(r);
			//size_t chrIdx = mapChr2Idx(r.chr);

			if (r.info.find("bnd_W") != string::npos) {
				lastBndPassed=true;
			}
			if (r.info.find("bnd") == string::npos) {
				lastBndPassed=false;
			}

			if (r.len > 0)
			{
				addInterval(iForestGaps[r.chr],r.pos,(r.pos+r.len));
			} else {
				addInterval(iForestGaps[r.chr],(r.pos+r.len),r.pos);
			}
		}	else {
			cout << "overlap. Skipping this one. Overlap count : " << overlapCount << endl;
			++overlapCount;
		}
		cout << "varCount=" << varCount << " overlapCount=" << overlapCount << endl;
		cout << "Collected " << vcfLines.size() << " vcfEntries." << endl;
	}
	cout << "Done" << endl;
	cout << "varCount=" << varCount << " overlapCount=" << overlapCount << endl;
	cout << "Collected " << vcfLines.size() << " vcfEntries." << endl;
	cout << "Length of header " << vcfHeader.size() << endl;
	
	// simulate translocations
	if (options.numTranslocations>0) {
		cout << "Simulating " << options.numTranslocations << " translocations." << endl;
		simulateTranslocations(vcfLines,iForestGaps,options,chrLengths,chrIDs,chrSeqs,refCount);
	}

	// filter for UCSC gaps
	//filterForUcscGaps(vcfLines,iForestGaps);
	cout << "Writing " << vcfLines.size() << " variants to " << options.outputFile << endl;
	ofstream overAndOut(toCString(options.outputFile));
	for (vector<string>::const_iterator ct = vcfHeader.begin(); ct != vcfHeader.end();++ct) {
		overAndOut << *ct << endl;
	}
	overAndOut << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\tFORMAT\t" << options.sampleId << endl;
	for (vcfStore::const_iterator ct = vcfLines.begin(); ct != vcfLines.end(); ++ct) {
		overAndOut << *ct << endl;
	}
	overAndOut.close();
	return (EXIT_SUCCESS);
}

