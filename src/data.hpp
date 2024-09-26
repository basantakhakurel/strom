#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <numeric>
#include <limits>
#include <map>
#include <boost/format.hpp>
#include "genetic_code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "xstrom.hpp"
#include "ncl/nxsmultiformat.h"
#include <boost/algorithm/string/join.hpp>

namespace strom
{

  class Data
  {
  public:
    typedef std::vector<std::string> taxon_names_t;
    typedef unsigned long long state_t;
    typedef std::vector<state_t> pattern_vect_t;
    typedef std::vector<state_t> monomorphic_vect_t;
    typedef std::vector<int> partition_key_t;
    typedef std::map<pattern_vect_t, unsigned> pattern_map_t;
    typedef std::vector<pattern_vect_t> data_matrix_t;
    typedef std::vector<pattern_map_t> pattern_map_vect_t;
    typedef std::vector<double> pattern_counts_t;
    typedef std::vector<unsigned> subset_end_t;
    typedef std::vector<unsigned> npatterns_vect_t;
    typedef std::pair<unsigned, unsigned> begin_end_pair_t;
    typedef std::shared_ptr<Data> SharedPtr;

    Data();
    ~Data();

    Partition::SharedPtr getPartition();
    void setPartition(Partition::SharedPtr partition);

    void getDataFromFile(const std::string filename);

    unsigned getNumSubsets() const;
    std::string getSubsetName(unsigned subset) const;

    unsigned getNumTaxa() const;
    const taxon_names_t &getTaxonNames() const;

    unsigned getNumPatterns() const;
    npatterns_vect_t calcNumPatternsVect() const;
    unsigned getNumPatternsInSubset(unsigned subset) const;
    unsigned getNumStatesForSubset(unsigned subset) const;
    unsigned calcSeqLen() const;
    unsigned calcSeqLenInSubset(unsigned subset) const;
    const data_matrix_t &getDataMatrix() const;
    begin_end_pair_t getSubsetBeginEnd(unsigned subset) const;
    const pattern_counts_t &getPatternCounts() const;
    const monomorphic_vect_t &getMonomorphic() const;
    const partition_key_t &getPartitionKey() const;

    void clear();

  private:
    unsigned storeTaxonNames(NxsTaxaBlock *taxaBlock, unsigned taxa_block_index);
    unsigned storeData(unsigned ntax, unsigned nchar, NxsCharactersBlock *charBlock, NxsCharactersBlock::DataTypesEnum datatype);
    unsigned buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets);
    void updatePatternMap(Data::pattern_vect_t &pattern, unsigned subset);
    void compressPatterns();

    Partition::SharedPtr _partition;
    pattern_counts_t _pattern_counts;
    monomorphic_vect_t _monomorphic;
    partition_key_t _partition_key;
    pattern_map_vect_t _pattern_map_vect;
    taxon_names_t _taxon_names;
    data_matrix_t _data_matrix;
    subset_end_t _subset_end;
  };

  // member functions

  // constructor and destructor
  inline Data::Data()
  {
    // std::cout << "Creating a Data Object" << std::endl;
    clear();
  }

  inline Data::~Data()
  {
    // std::cout << "Destryoing a Data Object" << std::endl;
  }

  // to access the partition information that the user supplies in the strom.conf file
  inline void Data::setPartition(Partition::SharedPtr partition)
  {
    _partition = partition;
  }

  // member functions to obtain various partition information

  // access to the Partition object through a shared pointer
  inline Partition::SharedPtr Data::getPartition()
  {
    return _partition;
  }

  // access to the number of data subsets in the partition
  inline unsigned Data::getNumSubsets() const
  {
    return (_partition ? _partition->getNumSubsets() : 1);
  }

  // access to the names assigned to the data subsets by the user
  inline std::string Data::getSubsetName(unsigned subset) const
  {
    return _partition ? _partition->getSubsetName(subset) : std::string("default");
  }

  // accesssor member functions

  inline const Data::partition_key_t &Data::getPartitionKey() const
  {
    return _partition_key;
  }

  inline const Data::pattern_counts_t &Data::getPatternCounts() const
  {
    return _pattern_counts;
  }

  inline const Data::monomorphic_vect_t &Data::getMonomorphic() const
  {
    return _monomorphic;
  }

  inline const Data::taxon_names_t &Data::getTaxonNames() const
  {
    return _taxon_names;
  }

  inline const Data::data_matrix_t &Data::getDataMatrix() const
  {
    return _data_matrix;
  }

  // function to obtain a pair with pattern index (beginning and end)
  inline Data::begin_end_pair_t Data::getSubsetBeginEnd(unsigned subset) const
  {
    assert(_subset_end.size() > subset);
    if (subset == 0)
      return std::make_pair(0, _subset_end[0]);
    else
      return std::make_pair(_subset_end[subset - 1], _subset_end[subset]);
  }

  // clear member function to empty data members that are vectors or maps
  inline void Data::clear()
  {
    _partition_key.clear();
    _pattern_counts.clear();
    _monomorphic.clear();
    _pattern_map_vect.clear();
    _taxon_names.clear();
    _data_matrix.clear();
    _subset_end.clear();
  }

  // function to get total number of patterns stored in all subsets (row length of _data_matrix vector)
  inline unsigned Data::getNumPatterns() const
  {
    if (_data_matrix.size() > 0)
      return (unsigned)_data_matrix[0].size();
    else
      return 0;
  }

  // function to obtain a vector containing the number of patterns in each subset
  inline Data::npatterns_vect_t Data::calcNumPatternsVect() const
  {
    unsigned nsubsets = (unsigned)_subset_end.size();
    std::vector<unsigned> num_patterns_vect(nsubsets, 0);
    for (unsigned s = 0; s < nsubsets; s++)
      num_patterns_vect[s] = getNumPatternsInSubset(s);
    return num_patterns_vect;
  }

  // function to get the number of states for each data subset
  inline unsigned Data::getNumStatesForSubset(unsigned subset) const
  {
    DataType data_type = _partition->getDataTypeForSubset(subset);
    return data_type.getNumStates();
  }

  // function to obtain the number of distinct patterns stored in a specified subset
  inline unsigned Data::getNumPatternsInSubset(unsigned subset) const
  {
    assert(_subset_end.size() > subset);
    return (unsigned)_subset_end[subset] - (subset == 0 ? 0 : _subset_end[subset - 1]);
  }

  // function to return the number of taxa
  inline unsigned Data::getNumTaxa() const
  {
    return (unsigned)_taxon_names.size();
  }

  // function that returns the total number of sites stored
  inline unsigned Data::calcSeqLen() const
  {
    return std::accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
  }

  // function that returns the number of sites in the subset whose index is provided
  inline unsigned Data::calcSeqLenInSubset(unsigned subset) const
  {
    begin_end_pair_t s = getSubsetBeginEnd(subset);
    return std::accumulate(_pattern_counts.begin() + s.first, _pattern_counts.begin() + s.second, 0);
  }

  // function that fills a map associating site patterns with counts
  inline unsigned Data::buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets)
  {
    pattern_vect_t pattern(ntaxa);

    _pattern_map_vect.clear();
    _pattern_map_vect.resize(nsubsets);

    const Partition::partition_t &tuples = _partition->getSubsetRangeVect();
    for (auto &t : tuples)
    {
      unsigned site_begin = std::get<0>(t);
      unsigned site_end = std::get<1>(t);
      unsigned site_skip = std::get<2>(t);
      unsigned site_subset = std::get<3>(t);
      for (unsigned site = site_begin; site <= site_end; site += site_skip)
      {
        // copy site into pattern
        for (unsigned taxon = 0; taxon < ntaxa; ++taxon)
        {
          pattern[taxon] = _data_matrix[taxon][site - 1];
        }

        // Add this pattern to _pattern_map_vect element corresponding to subset site_subset
        updatePatternMap(pattern, site_subset);
      }
    }

    // Tally total number of patterns across all subsets
    unsigned npatterns = 0;
    for (auto &map : _pattern_map_vect)
    {
      npatterns += (unsigned)map.size();
    }

    return npatterns;
  }

  // private member function to add 1 to the count stored for the supplied pattern in the supplied partition subset in the vector
  inline void Data::updatePatternMap(Data::pattern_vect_t &pattern, unsigned subset)
  {
    // If pattern is not already in the pattern_map, insert it and set the value to 1
    // If it does exist, increment its current value.
    // (see item 24 from Meyer's Efficient STL for more information on the technique used here)
    pattern_map_t::iterator lowb = _pattern_map_vect[subset].lower_bound(pattern);
    if (lowb != _pattern_map_vect[subset].end() && !(_pattern_map_vect[subset].key_comp()(pattern, lowb->first)))
    {
      // this pattern has already been seen
      lowb->second += 1;
    }
    else
    {
      // this pattern has not been seen yet
      _pattern_map_vect[subset].insert(lowb, pattern_map_t::value_type(pattern, 1));
    }
  }

  // function that takes all the pattern from the datamatrix and replaces the matrix with just unique
  // patterns and stores the number of sites exhibiting each pattern in the _pattern_counts data member
  inline void Data::compressPatterns()
  {
    // Perform sanity checks
    if (_data_matrix.empty())
      throw XStrom("Attempted to compress an empty data matrix");

    unsigned ntaxa = (unsigned)_data_matrix.size();
    unsigned seqlen = (unsigned)_data_matrix[0].size();

    // finalize partition
    unsigned nsubsets = getNumSubsets();
    _subset_end.resize(nsubsets);
    _partition->finalize(seqlen);

    // compact the data, storing it in _pattern_map_vect
    unsigned npatterns = buildSubsetSpecificMaps(ntaxa, seqlen, nsubsets);
    _pattern_counts.assign(npatterns, 0);
    _monomorphic.assign(npatterns, 0);
    _partition_key.assign(npatterns, -1);

    // rebuild the _data_matrix to hold compact data, storing counts in _pattern_counts
    _data_matrix.resize(ntaxa);
    for (auto &row : _data_matrix)
    {
      row.resize(npatterns);
    }

    unsigned p = 0;
    for (unsigned subset = 0; subset < nsubsets; subset++)
    {
      for (auto &pc : _pattern_map_vect[subset])
      {
        _pattern_counts[p] = pc.second; // record how many sites have pattern p
        _partition_key[p] = subset;     // record the subset to which pattern p belongs

        state_t constant_state = pc.first[0];
        unsigned t = 0;
        for (auto sc : pc.first)
        {
          assert(sc > 0);
          constant_state &= sc;
          _data_matrix[t][p] = sc;
          ++t;
        }
        // constant_state equals 0 if polymorphic or state code of state present if monomorphic
        _monomorphic[p] = constant_state;
        ++p;
      }

      _subset_end[subset] = p;

      // Everything for this susbet has been transferred to _data_matrix and _pattern_counts,
      // so we can now free this memory
      _pattern_map_vect[subset].clear();
    }
  }

  // function used by getDataFromFile to store the taxon labels found in a NEXUS taxa block
  inline unsigned Data::storeTaxonNames(NxsTaxaBlock *taxaBlock, unsigned taxa_block_index)
  {
    unsigned ntax = 0;
    if (taxa_block_index == 0)
    {
    }
  }

}
