/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Declaration of the skeleton of all c++ programs.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_COMMON_PROGRAM_HH
#define EAGLE_COMMON_PROGRAM_HH

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <typeinfo>
#include <boost/program_options.hpp>
#include <boost/noncopyable.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"

namespace eagle
{
namespace common
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;
typedef std::pair<bfs::path *, std::string> PathOption;

struct ProgramInfo
{
    ProgramInfo(std::string arg, unsigned int width)
    : tool(arg), toolName(tool.leaf().string()), project(EAGLE_NAME), version(EAGLE_VERSION), copyright(EAGLE_COPYRIGHT), optionsWidth(width) {}
    bfs::path tool;
    std::string toolName;
    std::string project;
    std::string version;
    std::string copyright;
    unsigned int optionsWidth;
    friend std::ostream& operator<< (std::ostream& os, ProgramInfo pi);
};

std::string naiveDemangling(const char *type);

/**
 ** Encapsulation of the processing of the command line options.
 **
 ** TODO: add config file and environment options
 **/
class Options: boost::noncopyable
{
public:
    enum Action
    {
        RUN, HELP, VERSION, ABORT
    };
    Options( bool anyOutput = true );
    virtual ~Options()
    {
    }
    Action parse(int argc, char *argv[]);
    std::string usage(bool full = false) const;
    unsigned int width() const {return bpo::options_description::m_default_line_length;}
protected:
    bpo::options_description parameters_;
    bpo::options_description namedOptions_;
    bpo::options_description unnamedOptions_;
    bpo::positional_options_description positionalOptions_;
    bpo::options_description general_;
private:
    virtual std::string usagePrefix() const = 0;
    virtual std::string usageSuffix() const
    {
        return "";
    }
    virtual void postProcess(bpo::variables_map &)
    {
    }
};

class OptionsHelper : bpo::variables_map
{
public:
    OptionsHelper(bpo::variables_map vm) : bpo::variables_map(vm) {}
    void requiredOptions(std::vector<std::string> requiredOptionList);
    std::string mutuallyExclusiveOptions(std::vector<std::string> mutuallyExclusiveList);

    void inputPathsExist();
    void outputFilesWriteable();
    void outputDirsWriteable();

    void addPathOptions(std::vector<bfs::path>& pathOptions, std::string label);
    void addPathOptions(bfs::path& pathOption, std::string label);
    void clearPathOptions() {pathOptions_.clear();}

    template<typename T>
    void inRange(std::pair<T,std::string> option, T minValue = std::numeric_limits<T>::min(), T maxValue = std::numeric_limits<T>::max())
    {
        if( minValue > option.first || option.first >= maxValue )
        {
            std::string t = naiveDemangling(typeid(T).name());
            const boost::format message = boost::format("\n   *** The '%s' option is out of range. Please specify a%s '%s' within [%s,%s)")
                                        % option.second
                                        % ((t[0] == 'a' || t[0] == 'e' || t[0] == 'i' || t[0] == 'o' || t[0] == 'u') ? "n" : "")
                                        % t
                                        % boost::lexical_cast<std::string>(minValue)
                                        % boost::lexical_cast<std::string>(maxValue);
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }

    }

private:
    std::vector<PathOption> pathOptions_;
};


/**
 ** Unified behavior of all programs.
 **/
template<class O>
void run(void(*callback)(const O &), int argc, char *argv[])
{
#ifdef EAGLE_DEBUG_MODE
    std::cout << "Command-line invocation:\n     ";
    for (int i = 0; i < argc; i++)
    {
        std::cout << std::string(argv[i]) << " ";
    }
    std::cout << std::endl;
#endif
    try
    {
        O options;
        ProgramInfo info( std::string(argv[0]), options.width() );
        const typename O::Action action = options.parse(argc, argv);
        if (O::RUN == action)
        {
            callback(options);
        }
        else
        {
            if (O::HELP == action || O::VERSION == action)
            {
                std::clog << info << std::endl;
                if (O::VERSION == action) exit(0);
            }
            std::clog << options.usage(O::HELP == action);
            exit(O::HELP != action);
        }
    }
    catch (const eagle::common::ExceptionData &exception)
    {
        std::clog << "Error: " << exception.getContext() << ": " << exception.getMessage() << std::endl;
        exit(1);
    }
    catch (const boost::exception &e)
    {
        std::clog << "Error: boost::exception: " << boost::diagnostic_information(e) << std::endl;
        exit(2);
    }
    catch (const std::bad_alloc &e)
    {
        std::clog << "memory allocation error: " << e.what() << std::endl;
        std::ifstream proc("/proc/self/status");
        std::string s;
        while(std::getline(proc, s), !proc.fail())
        {
            std::clog << "\t" << s << std::endl;
        }
        std::clog << "*** abandoned execution! ***" << std::endl;
        exit(3);
    }
    catch (const std::runtime_error &e)
    {
        std::clog << "runtime error: " << e.what() << std::endl;
        exit(4);
    }
    catch (const std::logic_error &e)
    {
        std::clog << "logic error: " << e.what() << std::endl;
        exit(5);
    }
}

} // namespace common
} // namespace eagle

#endif // EAGLE_COMMON_PROGRAM_HH
