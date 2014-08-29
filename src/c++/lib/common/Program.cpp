/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Implementation of the skeleton of all c++ programs.
 **
 ** \author Mauricio Varea
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace common
{


// TODO: improve!
std::string naiveDemangling(const char *type)
{
    if (0 == strcmp(type,"i"))
    {
        return std::string("int");
    } else if (0 == strcmp(type,"j")) {
        return std::string("unsigned int");
    }
    return std::string(type);
}


Options::Options( bool anyOutput )
    : parameters_("Command-line parameters")
    , namedOptions_("Command-line options")
    , general_()
{
    general_.add_options()("help,h",    "Produce help message and exit")
                          ("version,V", "Display version information");
    if (anyOutput)
    {
        general_.add_options()("force", "Overwrite output files");
    }
}

Options::Action Options::parse(int argc, char *argv[])
{
    try
    {
        bpo::options_description allOptions("Allowed options");
        allOptions.add(parameters_).add(namedOptions_).add(unnamedOptions_).add(general_);
        bpo::variables_map vm;
        bpo::store(
             bpo::command_line_parser(argc, argv)
                  .options(allOptions)
                  .positional(positionalOptions_)
                  .run(),
             vm);
        bpo::notify(vm);
        if (vm.count("help"))
        {
            return HELP;
        } else if (vm.count("version")) {
            return VERSION;
        } else {
            postProcess(vm);
            return RUN;
        }
    }
    catch (const std::exception &e)
    {
        std::clog << "Failed to parse the options: " << e.what() << std::endl;
        return ABORT;
    }
}

std::string Options::usage(bool full) const
{
    std::ostringstream os;
    os << this->usagePrefix() << std::endl << std::endl;
    if (parameters_.options().size())
    {
        os << parameters_ << std::endl;
    }
    os << namedOptions_ << std::endl;
    os << general_ << std::endl;
    if (full)
    {
        os << this->usageSuffix() << std::endl;
    }
    return os.str();
}

std::ostream& operator<< (std::ostream& os, ProgramInfo pi)
{
    os << std::setfill('=') << std::setw(  pi.optionsWidth ) << "=" << std::endl;
    os << std::setfill(' ') << std::setw( (pi.optionsWidth - pi.project.length())/2 ) << " ";
    os << pi.project << std::endl;
    os << std::setfill(' ') << std::setw( (pi.optionsWidth - pi.toolName.length() - pi.version.length() - 5)/2 ) << " ";
    os << pi.toolName << " -- v"<< pi.version << std::endl;
    os << std::setfill(' ') << std::setw( (pi.optionsWidth - pi.copyright.length() - 4)/2 ) << " ";
    os << "- " << pi.copyright << " -" << std::endl;
    os << std::setfill('=') << std::setw(  pi.optionsWidth ) << "="  << std::endl;
    return os;
}

/*
 * \brief check that no mandatory options are missing
 */
void OptionsHelper::requiredOptions(std::vector<std::string> requiredOptionList)
{
    BOOST_FOREACH(const std::string &required, requiredOptionList)
    {
        if(!count(required))
        {
            const boost::format message = boost::format("\n   *** Missing Option: The '%s' option is required ***\n") % required;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}

/*
 * \brief check that only one mutually exclusive option exist, and find out which one it is
 */
std::string OptionsHelper::mutuallyExclusiveOptions(std::vector<std::string> mutuallyExclusiveList)
{
    std::vector<std::string> onlyOneOption;
    BOOST_FOREACH(const std::string &exclusive, mutuallyExclusiveList)
    {
        if(count(exclusive))  { onlyOneOption.push_back(exclusive); }
    }
    if (1 != onlyOneOption.size())
    {
        std::stringstream message;
        message << std::endl<< "   *** One, and only one, of the following options is required ***" << std::endl;
        BOOST_FOREACH(const std::string &exclusive, mutuallyExclusiveList)
        {
            message << (boost::format("       '%s'") % exclusive).str() << std::endl;
        }
        message << (boost::format("   *** Found %d of them ***\n") % onlyOneOption.size() ).str() << std::endl;
        BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
    }
    return onlyOneOption[0];
}

void OptionsHelper::addPathOptions(std::vector<bfs::path>& pathOptions, std::string label)
{
    for (unsigned int i=0; i < pathOptions.size(); i++)
    {
        pathOptions_.push_back(PathOption(&pathOptions[i], label));
    }
}

void OptionsHelper::addPathOptions(bfs::path& pathOption, std::string label)
{
    if (! pathOption.empty())
    {
        pathOptions_.push_back(PathOption(&pathOption, label));
    }
}

/*
 * \brief check that all input paths exist
 */
void OptionsHelper::inputPathsExist()
{
    BOOST_FOREACH(const PathOption &pathOption, pathOptions_)
    {
        if(!bfs::exists(*pathOption.first))
        {
            const boost::format message = boost::format("\n   *** The '%s' path does not exist: %s ***\n") % pathOption.second % *pathOption.first;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}

/*
 * \brief check that all output files don't exist (or are forced to be re-written, if they are not dirs)
 */
void OptionsHelper::outputFilesWriteable()
{
    BOOST_FOREACH(const eagle::common::PathOption &pathOption, pathOptions_)
    {
        if (bfs::exists(*pathOption.first))
        {
            if( !count("force") )
            {
                const boost::format message = boost::format("\n   *** Option '%s' has an invalid value: ***"
                                                            "\n   ***        Cannot write into %s as it already exists! ***")
                                            % pathOption.second % *pathOption.first;
                BOOST_THROW_EXCEPTION(InvalidOptionException(message.str() +
                                                            "\n   ***        (you can use --force to overwrite this parameter, at your own risk) ***\n"));
            } else if ( bfs::is_directory(*pathOption.first) ) {
                const boost::format message = boost::format("\n   *** Option '%s' has an invalid value: ***"
                                                            "\n   ***     It should point to a file, but a directory already exists with name %s ***\n")
                                            % pathOption.second % *pathOption.first;
                BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
            }
        }
    }
}

/*
 * \brief check that all output dirs either are directories or don't exist at all
 */
void OptionsHelper::outputDirsWriteable()
{
    BOOST_FOREACH(const eagle::common::PathOption &pathOption, pathOptions_)
    {
        if(bfs::exists(*pathOption.first) && !bfs::is_directory(*pathOption.first))
        {
            const boost::format message = boost::format("\n   *** Option '%s' has an invalid value: ***"
                                                        "\n   ***     It should point to a directory, but a file already exists with name %s ***\n")
                                        % pathOption.second % *pathOption.first;
            BOOST_THROW_EXCEPTION(InvalidOptionException(message.str()));
        }
    }
}


} // namespace common
} // namespace eagle
