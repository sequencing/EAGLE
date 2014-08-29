/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "genome/TmpFastaReader.hh"

namespace eagle
{
namespace genome
{

TmpFastaReader::TmpFastaReader( const boost::filesystem::path &refDir )
    : fastaRef_ ( refDir )
{
    unsigned long globalPos = 0;
    std::vector<std::string> contigNames = fastaRef_.allContigNames();
    for (unsigned int i=0; i<contigNames.size(); ++i)
    {
        struct TmpFastaFileInfo info;
        info.contigName = contigNames[i];
        info.baseCount = allContigLengths()[i];

        // file
        boost::filesystem::path fullFilename = refDir / (info.contigName + ".fa");
        info.file = new ifstream( fullFilename.string().c_str() );
        if( !info.file->good() )
        {
            fullFilename = refDir;
            info.file = new ifstream( fullFilename.string().c_str() );
            assert( info.file->good() );
        }

        // headerLength
        string line;
        getline( *info.file, line );
        info.headerLength = line.length() + 1;

        // basesPerLine
        getline( *info.file, line );
        info.basesPerLine = line.length();

        // sanity check: fileSize==headerLength+(baseCount/basesPerLine)*(basesPerLine+1) + basesOnLastLine
        info.fileSize = boost::filesystem::file_size( fullFilename );
        if (info.baseCount == 0)
        {
            clog << "No contig info filled by Reference.cpp" << endl;
            unsigned long fullLinesCount = (info.fileSize - info.headerLength) / (info.basesPerLine+1);
            unsigned long basesOnLastLine = info.fileSize - info.headerLength - fullLinesCount * (info.basesPerLine+1);
            info.baseCount = fullLinesCount * info.basesPerLine + basesOnLastLine - (basesOnLastLine?1:0);
        }

        unsigned long fullLinesCount = info.baseCount / info.basesPerLine;
        unsigned long basesOnLastLine = info.baseCount % info.basesPerLine;
        unsigned long expectedFileSize = info.headerLength + fullLinesCount * (info.basesPerLine+1) + basesOnLastLine + (basesOnLastLine?1:0);
        clog << "=============" << endl;
        clog << "Contig name: " << info.contigName << endl;
        clog << "baseCount: " << info.baseCount << endl;
        clog << "headerLength: " << info.headerLength << endl;
        clog << "basesPerLine: " << info.basesPerLine << endl;
        clog << "fileSize: " << info.fileSize << endl;
        clog << "fullLinesCount: " << fullLinesCount << endl;
        clog << "basesOnLastLine: " << basesOnLastLine << endl;
        clog << "expectedFileSize: " << expectedFileSize << endl;
        assert( info.fileSize == expectedFileSize );
        clog << "=============" << endl;

        info.globalPosMin = globalPos;
        info.globalPosMax = globalPos + info.baseCount - 1;
        globalPos += info.baseCount;

        fileInfos_.push_back(info);
    }
}

char TmpFastaReader::get( const unsigned long globalPos, const unsigned long offset, bool& overlapContigBoundary )
{
    static unsigned int lastContigNum = 0 ;
    static TmpFastaFileInfo *fileInfo = &fileInfos_[lastContigNum];
    static char *cachedFileData = 0;
    while (globalPos < fileInfo->globalPosMin)
    {
        fileInfo = &fileInfos_[--lastContigNum];

        if (cachedFileData)
        {
            free( cachedFileData );
            cachedFileData = 0;
        }
    }
    while (globalPos > fileInfo->globalPosMax)
    {
        fileInfo = &fileInfos_[++lastContigNum];

        if (cachedFileData)
        {
            free( cachedFileData );
            cachedFileData = 0;
        }
    }
    unsigned long posInContig = (globalPos-fileInfo->globalPosMin+offset)%fileInfo->baseCount;
    unsigned long fullLinesCount = posInContig / fileInfo->basesPerLine;
    unsigned long posInLine = posInContig % fileInfo->basesPerLine;
    unsigned long posInFile = fileInfo->headerLength + fullLinesCount * (fileInfo->basesPerLine+1) + posInLine;


    if (!cachedFileData)
    {
        clog << "TmpFastaReader: Reading file " << fileInfo->contigName << endl;
        cachedFileData = (char*)malloc( fileInfo->fileSize );
        fileInfo->file->seekg(0);
        fileInfo->file->read( cachedFileData, fileInfo->fileSize );
        clog << "TmpFastaReader: Done reading file" << endl;
    }

    /*
      fileInfo->file->seekg( posInFile );
      char result;
      fileInfo->file->read( &result, 1 );
    */
    char result = cachedFileData[posInFile];

    // Debugging
    /*
      static int debugCount = 500;
      if (debugCount>0)
      {
      clog << "TmpFastaReader: Reading from position " << globalPos << "+" << offset << " \t=> " << result << endl;
      debugCount--;
      }
    */
    return result;
}


void TmpFastaReader::convertFromGlobalPos( const unsigned long globalPos, int& refId, unsigned long& posInContig )
{
    static unsigned int lastContigNum = 0 ;
    static TmpFastaFileInfo *fileInfo = &fileInfos_[lastContigNum];
    while (globalPos < fileInfo->globalPosMin)
    {
        fileInfo = &fileInfos_[--lastContigNum];
    }
    while (globalPos > fileInfo->globalPosMax)
    {
        fileInfo = &fileInfos_[++lastContigNum];
    }
    refId = lastContigNum;
    posInContig = (globalPos - fileInfo->globalPosMin) % fileInfo->baseCount + 1;
}


} // namespace genome
} // namespace eagle
