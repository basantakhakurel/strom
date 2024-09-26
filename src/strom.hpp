#pragma once

#include <iostream>
#include "data.hpp"
#include "tree_summary.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace strom
{
  class Strom
  {
  public:
    Strom();
    ~Strom();

    void clear();
    void processCommandLineOptions(int argc, const char *argv[]);
    void run();

  private:
    std::string _data_file_name;
    std::string _tree_file_name;

    Partition::SharedPtr _partition;
    Data::SharedPtr _data;

    TreeSummary::SharedPtr _tree_summary;

    static std::string _program_name;
    static unsigned _major_version;
    static unsigned _minor_version;
  };

  // member functions

  // constructor
  inline Strom::Strom()
  {
    // std::cout << "Constructing a Strom" << std::endl;
    clear();
  }

  // destructor
  inline Strom::~Strom()
  {
    // std::cout << "Destroying a Strom" << std::endl;
  }

  /**
   * @brief Sets each of the two non-static data members to an empty string
   *
   * @param name
   *
   * @returns void
   */
  inline void Strom::clear()
  {
    _data_file_name = "";
    _tree_file_name = "";
    _tree_summary = nullptr;
    _partition.reset(new Partition());
    _data = nullptr;
  }

  /**
   * Process command line options and configuration file
   *
   * @brief This function processes the command line options and configuration file
   * to sets the program's parameters. It supports options for help, version,
   * data file, and tree file.
   *
   * @param argc The number of command line arguments.
   * @param argv An array of command line arguments.
   *
   * @return None
   *
   * @throws boost::program_options::readinf_file If the configuration file cannot be read.
   */

  inline void Strom::processCommandLineOptions(int argc, const char *argv[])
  {
    std::vector<std::string> partition_subsets;
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")("version,v", "show program version")("datafile,d", boost::program_options::value(&_data_file_name)->required(), "name of a data file in NEXUS format")("treefile,t", boost::program_options::value(&_tree_file_name), "name of a tree file in NEXUS format")("subset", boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'");

    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

    try
    {
      const boost::program_options::parsed_options &parsed = boost::program_options::parse_config_file<char>("strom.conf", desc, false);
      boost::program_options::store(parsed, vm);
    }
    catch (boost::program_options::reading_file &x)
    {
      std::cout << "Note: configuration file (strom.conf) not found" << std::endl;
    }
    boost::program_options::notify(vm);

    // If the user specified "--help" on command line, output usage usage summary and quit
    if (vm.count("help") > 0)
    {
      std::cout << desc << "\n";
      std::exit(1);
    }

    // If the user specified "--version" on command line, output version and quit
    if (vm.count("version") > 0)
    {
      std::cout << boost::str(boost::format("This is %s version %d.%s") % _program_name % _major_version % _minor_version) << std::endl;
      std::exit(1);
    }

    // If user specified --subset on command line, break specified partition subset
    // definition into name and character set string and add to _partition
    if (vm.count("subset") > 0)
    {
      _partition.reset(new Partition());
      for (auto s : partition_subsets)
      {
        _partition->parseSubsetDefinition(s);
      }
    }
  }

  inline void Strom::run()
  {
    std::cout << "Starting......" << std::endl;
    std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;

    try
    {
      std::cout << "\n*** Reading and storing the data in the file" << _data_file_name << std::endl;
      _data = Data::SharedPtr(new Data());
      _data->setPartition(_partition);
      _data->getDataFromFile(_data_file_name);

      // Report information about data partition subsets
      unsigned nsubsets = _data->getNumSubsets();
      std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
      std::cout << "Number of partition susets: " << nsubsets << std::endl;
      for (unsigned subset = 0; subset < nsubsets; subset++)
      {
        DataType dt = _partition->getDataTypeForSubset(subset);
        std::cout << " Subset " << (subset + 1) << " (" << _data->getSubsetName(subset) << ")" << std::endl;
        std::cout << "    data type: " << dt.getDataTypeAsString() << std::endl;
        std::cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << std::endl;
        std::cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << std::endl;
      }
    }

    catch (XStrom &x)
    {
      std::cerr << "Strom encountered a problem:(\n " << x.what() << std::endl;
    }

    std::cout << "\nFinished!" << std::endl;
  }

}
