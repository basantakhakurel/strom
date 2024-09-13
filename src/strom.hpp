#pragma once

#include <iostream>
#include "tree_summary.hpp"
#include <boost/program_options.hpp>

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

    TreeSummary::SharedPtr _tree_summary;

    static std::string _program_name;
    static unsigned _major_version;
    static unsigned _minor_version;
  };

  inline Strom::Strom()
  {
    // std::cout << "Constructing a Strom" << std::endl;
  }

  inline Strom::~Strom()
  {
    // std::cout << "Destroying a Strom" << std::endl;
  }

  /**
   * @brief Sets each of the two non-static data members to an empty string
   *
   * @param None
   *
   * @returns void
   */
  inline void Strom::clear()
  {
    _data_file_name = "";
    _tree_file_name = "";
    _tree_summary = nullptr;
  }

  /**
   * Process command line options and configuration file.
   *
   * @brief This function processes the command line options and configuration file
   * to set the program's parameters. It supports options for help, version,
   * data file, and tree file.
   *
   * @param argc The number of command line arguments.
   * @param argv An array of command line arguments.
   *
   * @return None
   *
   * @throws boost::program_options::reading_file If the configuration file cannot be read.
   */
  inline void Strom::processCommandLineOptions(int argc, const char *argv[])
  {
    boost::program_options::variables_map vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()

        ("help,h", "produce help message")("version,v", "show program version")("datafile,d", boost::program_options::value(&_data_file_name), "name of a data file in NEXUS format")("treefile,t", boost::program_options::value(&_tree_file_name)->required(), "name of a tree file in NEXUS format");

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
  }

  inline void Strom::run()
  {
    std::cout << "Starting......." << std::endl;

    try
    {
      // create a new TreeSummary object and let _tree_summary point to it
      _tree_summary = TreeSummary::SharedPtr(new TreeSummary());

      // Read the tree file specified by the user
      _tree_summary->readTreefile(_tree_file_name, 0);
      Tree::SharedPtr tree = _tree_summary->getTree(0);

      // Summarize the trees read
      _tree_summary->showSummary();
    }

    catch (XStrom &x)
    {
      std::cerr << "Strom encountered a problem:(\n " << x.what() << std::endl;
    }

    std::cout << "\nFinished!" << std::endl;
  }

}
