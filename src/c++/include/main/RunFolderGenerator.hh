/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to induce variants in a reference.
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_RUN_FOLDER_GENERATOR_HH
#define EAGLE_MAIN_RUN_FOLDER_GENERATOR_HH

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "io/RunInfo.hh"
#include "main/RunFolderGeneratorOptions.hh"

namespace bfs = boost::filesystem;


namespace eagle
{
namespace main
{

class RunFolderGenerator
{
public:
    RunFolderGenerator( const RunFolderGeneratorOptions& options );
    void run();
private:
    void generateDirectoryStructure() const;
    void generateMetadata() const;
    void generateRunInfo() const;
    void generateConfig() const;
    void generateSampleSheet() const;
    void generateMatrix() const;
    void generatePhasing() const;

    const RunFolderGeneratorOptions &options_;
    eagle::io::RunInfo runInfo_;
    const bfs::path runFolderPath_;
    const bfs::path dataPath_;
    const bfs::path intensitiesPath_;
    const bfs::path baseCallsPath_;

    string runFolder_, runFolderDate_, runFolderId_, instrument_, flowcell_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_RUN_FOLDER_GENERATOR_HH
