#pragma once

#include <tuple>
#include <limits>
#include <cmath>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "genetic_code.hpp"
#include "datatype.hpp"
#include "xstrom.hpp"

namespace strom
{
  class Partition
  {
  public:
    typedef std::match_results<std::string::const_iterator>::const_reference regex_match_t;
    typedef std::tuple<unsigned, unsigned, unsigned, unsigned> subset_range_t;
    typedef std::vector<subset_range_t> partition_t;
    typedef std::vector<DataType> datatype_vect_t;
    typedef std::vector<unsigned> subset_sizes_vect_t;
    typedef std::vector<std::string> subset_names_vect_t;
    typedef std::shared_ptr<Partition> SharedPtr;

    Partition();
    ~Partition();

    unsigned getNumSites() const;
    unsigned getNumSubsets() const;
    std::string getSubsetName(unsigned subset) const;

    const partition_t &getSubsetRangeVect() const;

    unsigned findSubsetByName(const std::string &subset_name) const;
    unsigned findSubsetForSite(unsigned site_index) const;
    bool siteInSubset(unsigned site_index, unsigned subset_index) const;
    DataType getDataTypeForSubset(unsigned subset_index) const;
    const datatype_vector_t &getSubsetDataTypes() const;

    unsigned numSitesInSubset(unsigned subset_index) const;
    subset_sizes_vect_t calcSubsetSizes() const;

    void defaultPartition(unsigned nsites = std::numeric_limits<unsigned>::max());
    void parseSubsetDefinition(std::string &s);
    void finalize(unsigned nsites);

    void clear();

  private:
    int extractIntFromRegexMatch(regex_match_t s, unsigned min_value);
    void addSubsetRange(unsigned subset_index, std::string range_definition);
    void addSubset(unsigned subset_index, std::string subset_definition);

    unsigned _num_sites;
    unsigned _num_subsets;
    subset_names_vect_t _subset_names;
    partition_t _subset_ranges;
    datatype_vect_t _subset_data_types;

    const unsigned _infinity;
  };

  // constructor
  inline Partition::Partition() : _infinity(std::numeric_limits<unsigned> : max())
  {
    // std:;cout << "Constructing a Partition" << std::endl;
    clear();
  }

  // descructor
  inline Partition::~Partition()
  {
    // std::cout << "Destroying a Partition" << std::endl;
  }

  // accessor functions - provide access to the values stored in the private data members but do not allow you to change those variables
  inline unsigned Partition::getNumSites() const
  {
    return _num_sites;
  }

  inline unsigned Partition::getNumSubsets() const
  {
    return _num_subsets;
  }

  inline std::string Partition::getSubsetName(unsigned subset) const
  {
    assert(subset < _num_subsets);
    return _subset_names[subset];
  }

  inline const Partition::partition_t &Partition::getSubsetRangeVect() const
  {
    return _subset_ranges;
  }

  inline DataType Partition::getDataTypeForSubset(unsigned subset_index) const
  {
    assert(subset_index < _subset_data_types.size());
    return _subset_data_types[subset_index];
  }

} // namespace strom
