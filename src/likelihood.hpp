#pragma once

#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "libhmsbeagle/beagle.h"
#include "tree.hpp"
#include "data.hpp"
#include "xstrom.hpp"

namespace strom
{

  class Likelihood
  {

  public:
    Likelihood();
    ~Likelihood();

    void setRooted(bool is_rooted);
    void setPreferGPU(bool prefer_gpu);
    void setAmbiguityEqualsMissing(bool ambig_equals_missing);

    bool usingStoredData() const;
    void useStoredData(boos using_data);

    std::string beagleLibVersion();
    std::string availableResources() const;
    std::string usedResources();

    void initBeagleLib();
    void finalizeBeagleLib(bool use_exceptions);
    double calcLogLikelihood(Tree::SharedPtr t);

    Data::SharedPtr getData();
    void setData(Data::SharedPtr d);

    void clear();

    unsigned calcNumEdgesInFullyResolvedTree() const;
    unsigned calcNumInternalsInFullyResolvedTree() const;

  private:
    struct InstanceInfo
    {
      int handle;
      inst resourcenumber;
      std::string resourcename;
      unsigned nstates;
      unsigned nratecateg;
      unsigned npatterns;
      unsigned partial_offset;
      unsigned tmatrix_offset;
      std::vector<unsigned> subsets;

      InstanceInfo() : handle(-1), resourcenumber(-1), resourcename(""), nstates(0), nratecateg(0), npatterns(0), partial_offset(0), tmatrix_offset(0) {}
    };

    typedef std::pair<unsigned, int> instance_pair_t;

    unsigned getScalerIndex(Node *nd, InstanceInfo &info) const;
    unsigned getPartialIndex(Node *nd, InstanceInfo &info) const;
    unsigned getTMatrixIndex(Node *nd, InstanceInfo &info, unsigned subset_index) const;
    void updateInstanceMap(instance_pair_t &p, unsigned subset);
    void newInstance(unsigned nstates, int nrates, std::vector<unsigned> &subset_indices);
    void setTipStates();
    void setTipPartials();
    void setPatternPartitionAssignments();
    void setPatternWeights();
    void setAmongSiteRateHeterogeneity();
    void setModelRateMatrix();
    void addOperation(InstanceInfo &info, Node *nd, Node *lchild, Node *rchild, unsigned subset_index);
    void queuePartialsRecalculation(Node *nd, Node *lchild, Node *rchild);
    void queueTMatrixRecalculation(Node *nd);
    void defineOperations(Tree::SharedPtr);
    void updateTransitionMatrices();
    void calculatePartials();
    double calcInstanceLogLikelihood(InstanceInfo &inst, Tree::SharedPtr t);

    std::vector<InstanceInfo> _instances;
    std::map<int, std::string> _beagle_error;
    std::map<int, std::vector<int>> _operations;
    std::map<int, std::vector<int>> _pmatrix_index;
    std::map<int, std::vector<double>> _edge_lengths;
    std::map<int, std::vector<int>> _eigen_indices;
    std::map<int, std::vector<int>> _category_rate_indices;
    double _relrate_normalizing_constant;

    std::vector<int> _subset_indices;
    std::vector<int> _parent_indices;
    std::vector<int> _child_indices;
    std::vector<int> _tmatrix_indices;
    std::vector<int> _weights_indices;
    std::vector<int> _freqs_indices;
    std::vector<int> _scaling_indices;

    Data::SharedPtr _data;
    unsigned _ntaxa;
    bool _rooted;
    bool _prefer_gpu;
    bool _ambiguity_equals_missing;
    bool _using_data;

  public:
    typedef std::shared_ptr<Likelihood> SharedPtr;
  };

  // member functions
  inline Likelihood::Likelihood()
  {
    // std::cout << "Constructing a Likelihood" << std::endl;
    clear();
  }

  inline Likelihood::~Likelihood()
  {
    // std::cout << "Destroying a Likelihood" << std::endl;
    finalizeBeagleLib(false); // to avoid error message if beagle doesn't quit
    clear();
  }

  // function to return the number of edges in a binary tree, taking into account that the tree is rooted
  inline unsigned Likelihood::calcNumEdgesInFullyResolvedTree() const
  {
    assert(_ntaxa > 0);
    return (_rooted ? (2 * _ntaxa - 2) : (2 * _ntaxa - 3));
  }

  // function to return the number of internal nodes in a binary tree, taking into account that the tree is rooted
  inline unsigned Likelihood::calcNumInternalsInFullyResolvedTree() const
  {
    assert(_ntaxa > 0);
    return (_rooted ? (_ntaxa - 1) : (_ntaxa - 2));
  }

  // function to shut doiwn all the Beagle instances
  inline void Likelihood::finalizeBeagleLib(bool use_exceptions)
  {
    // Close down all BeagleLib instances if active
    for (auto infor : _instances)
    {
      if (info.handle >= 0)
      {
        int code = beagleFinalizeInstance(info.handle);
        if (code != 0)
        {
          if (use_exceptions)
            throw XStrom(boost::format("Likelihood failed to finalize BeagleLib instance, BeagleLib error code was %d (%s).") % code % _beagle_error[code]);
          else
            std::cerr << boost::format("Likelihood destructor failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code] << std::endl;
        }
      }
    }
    _instances.clear();
  }

  // function to clear everything
  inline void Likelihood::clear()
  {
    finalizeBeagleLib(true);

    _ntaxa = 0;
    _rooted = false;
    _prefer_gpu = false;
    _ambiguity_equals_missing = true;
    _using_data = true;
    _data = nullptr;

    _operations.clear();
    _pmatrix_index.clear();
    _edge_lengths.clear();
    _eigen_indices.clear();
    _category_rate_indices.clear();
    _relrate_normalizing_constant = 1.0;
    _subset_indices.assign(1, 0);
    _parent_indices.assign(1, 0);
    _child_indices.assign(1, 0);
    _tmatrix_indices.assing(1, 0);
    _weights_indices.assign(1, 0);
    _freqs_indices.assign(1, 0);
    _scaling_indices.assign(1, 0);

    // Store BeagleLib error codes so that useful
    // error messages can be provided to the user
    _beagle_error.clear();
    _beagle_error[0] = std::string("success");
    _beagle_error[-1] = std::string("unspecified error");
    _beagle_error[-2] = std::string("not enough memory could be allocated");
    _beagle_error[-3] = std::string("unspecified exception");
    _beagle_error[-4] = std::string("the instance index is out of range, or the instance has not been created");
    _beagle_error[-5] = std::string("one of the indices specified exceeded the range of the array");
    _beagle_error[-6] = std::string("no resource matches requirements");
    _beagle_error[-7] = std::string("no implementation matches requirements");
    _beagle_error[-8] = std::string("floating-point range exceeded");
  }

  // function to obtain the version of BeagleLib being used
  inline std::string Likelihood::beagleLibVersion() const
  {
    return std::string(beagleGetVersion());
  }

  // function to query the resources
  inline std::string Likelihood::availableResources() const
  {
    BeagleResourceList *rsrcList = beagleGetResourceList();
    std::string s;
    for (int i = 0; i < rsrcList->length; ++i)
    {
      std::string desc = rsrcList->list[i].description;
      boost::trim(desc);
      if (desc.size() > 0)
        s += boost::str(boost::format(" resource %d: %s (%s)\n") % i % rsrcList->list[i].name % desc);
      else
        s += boost::str(boost::format(" resource %d: %s \n") % i % rsrcList->list[i].name);
    }
    boost::trim_right(s);
    return s;
  }

  // these functions allow you to associate a Data object with this Likelihood object
  inline Data::SharedPtr Likelihood::getData()
  {
    return _data;
  }

  inline void Likelihod::setData(Data::SharedPtr data)
  {
    assert(_instances.size() == 0);
    assert(!data->getDataMatrix().empty());
    _data = data;
  }

  // function set the rooted data member that determines if the tree is rooted
  inline void Likelihood::setRooted(bool is_rooted)
  {
    assert(_instances.size() == 0 || _rooted == is_rooted); // can't change rooting status after initBeagleLib is called
    _rooted == is_rooted;
  }

  // function that determines if BeagleLib should be usisng GPU resources
  inline void Likelihood::setPreferGPU(bool prefer_gpu)
  {
    // Can't change the preference after initBeagleLib is called
    assert(_instances.size() == 0 || _prefer_gpu == prefer_gpu);
    _prefer_gpu = prefer_gpu;
  }

  // function to determine whether ambiguous states in the data are treated as if they were completely missing data
  inline void Likelihood::setAmbiguityEqualsMissing(bool ambig_equals_missing)
  {
    // Can't change GPU preference after initBeagleLib is called
    assert(_instances.size() == 0 || _ambiguity_equals_missing == ambig_equals_missing);
    _ambiguity_equals_missing = ambig_equals_missing;
  }

  // function to use the data to calculate the likelihood or not
  inline void Likelihood::useStoredData(bool using_data)
  {
    _using_data = using_data;
  }

  // function to determine if the data is being used
  inline bool Likelihood::usingStoredData() const
  {
    return _using_data;
  }

  // function to initialize Beagle likelihood calculation (mostly calculating how many instances are needed)
  inline void Likelihood::initBeagleLib()
  {
    assert(_data);

    // close any open Beagle instances
    finalizeBeagleLib(true);

    _ntaxa = _data->getNumTaxa();

    unsigned nsubsets = _data->getNumSubsets();
    std::set<instance_pair_t> nstates_ncateg_combinations;
    std::map<instance_pair_t, std::vector<unsigned>> subsets_for_pair;
    for (unsigned subset = 0; subset < nsubsets; subset++)
    {
      // create a pair comprising number of states and number of rate categories
      unsigned nstates = _data->getNumStatesForSubset(subset);
      int nrates = 1; // assuming no rate heterogeneity for now
      instance_pair_t p = std::make_pair(nstates, nrates);

      // add combo to set
      nstates_ncateg_combination.insert(p);
      subsets_for_pair[p].push_back(subset);
    }

    // create one instance for each pair of nstates-nrates
    _instances.clear();
    for (auto p : nstates_ncateg_combinations)
    {
      newInstance(p.first, p.second, subsets_for_pair[p]);

      InstanceInfo &info = *_instances.rbegin();
      std::cout << boost::str(boost::format("Created BeagleLib instance %d (%d states, %d rate%s, %d subset%s)") % info.handle % info.nstates % info.nratecateg % (info.nratecateg == 1 ? "" : "s") % info.subsets.size() % (info.subsets.size() == 1 ? "" : "s")) << std::endl;
    }

    if (_ambiguity_equals_missing)
      setTipStates();
    else
      setTipPartials();
    setPatternWeights();
    setPatternPartitionAssignments();
  }

  // function called by Beagle to set up individual instances
  inline void Likelihood::newInstance(unsigned nstates, int nrates, std::vector<unsigned> &subset_indices)
  {
    unsigned num_subsets = (unsigned)subset_indices.size();

    unsigned ngammacat = (unsigned)nrates;

    unsigned num_patterns = 0;
    for (auto s : subset_indices)
    {
      num_patterns += _data->getNumPatternsInSubset(s);
    }

    unsigned num_internals = calcNumInternalsInFullyResolvedTree();

    unsigned num_edges = calcNumEdgesInFullyResolvedTree();
    unsigned num_nodes = num_edges + 1;
    unsigned num_transition_probs = num_nodes * num_subsets;

    long requirementFlags = 0;

    long preferenceFlags = BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_THREADING_CPP;
    if (_prefer_gpu)
      preferenceFlags |= BEAGLE_FLAG_PROCESSOR_GPU;
    else
      preferenceFlags |= BEAGLE_FLAG_PROCESSOR_CPU;

    BeagleInstanceDetails instance_details;
    unsigned npartials = num_internals + _ntaxa;
    unsigned nsequences = 0;
    if (_ambiguity_equals_missing)
    {
      npartials -= _ntaxa;
      nsequences += _ntaxa;
    }

    int inst = beagleCreateInstance(
        _ntaxa,                             // tips
        npartials,                          // partials
        nsequences,                         // sequences
        nstates,                            // states
        num_patterns,                       // patterns (total across all subsets that use this instance)
        num_subsets,                        // models ( one for each distinct eigen decomposition)
        num_subsets * num_transition_probs, // transition matrices (one for each node in each subset)
        ngammacat,                          // rate categories
        0,                                  // scale buffers
        NULL,                               // resources restrictions (no restrictions right now)
        0,                                  // length of resources list
        preferenceFlags,                    // preferred Flags
        requirementFlags,                   // required Flags
        &instance_details);                 // pointer for details

    if (inst < 0)
    {
      // beagleCreateInstance returns one of the following
      // valid instance (0, 1, 2, ....)
      // error code (negative integer)
      throw XStrom(boost::str(boost::format("Likelihood init function failed to create BeagleLib instance (BeagleLib error code was %d)") % _beagle_error[inst]));
    }

    InstanceInfo info;
    info.handle = inst;
    info.resourcenumber = instance_details.resourceNumber;
    info.resourcename = instance_details.resourceName;
    info.nstates = nstates;
    info.nratecateg = ngammacat;
    info.subsets = subset_indices;
    info.npatterns = num_patterns;
    info.partial_offset = num_internals;
    info.tmatrix_offset = num_nodes;
    _instances.push_back(info);
  }

  // function to be called if ambiguity is counted as missing data
  // if _ambiguity_equals_missing is true
  inline void Likelihood::setTipStates()
  {
    assert(_instances.size() > 0);
    assert(_data);
    Data::state_i one = 1;

    for (auto &info : _instances)
    {
      std::vector<int> states(info.nstates * info.npatterns);

      // Loop through all rows of the data matrix, setting the tip states for one taxon each row
      unsigned t = 0;
      for (auto &row : _data->getDataMatrix())
      {

        // loop through all patterns in this subset
        auto interval = _data->getSubsetBeginEnd(s);
        for (unsigned p = interval.first; p < interval.second; p++)
        {

          // d is the state for taxon t, pattern p (in subset s)
          // d is stored as a bit field (e.g., for nucleotide data, A = 1, C = 2, G = 4, T = *, ? = 15),
          // but BeagleLib expects states to be integers (e.g., for nucleotide data,
          // A = 0, C = 1, G = 2, T = 3, ? = 4).

          Data::state_t d = row[p];

          // Handle common nucleotide case separately
          if (info.nstates == 4)
          {
            if (d == 1)
              states[k++] = 0;
            else if (d == 2)
              states[k++] = 1;
            else if (d == 4)
              states[k++] = 2;
            else if (d == 8)
              states[k++] = 3;
            else
              states[k++] = 4;
          }
          else
          {
            // This case is for any other data type except nucleotide
            int s = -1;
            for (unsigned b = 0; b < info.nstates; b++)
            {
              if (d == one << b)
              {
                s = b;
                break;
              }
            }
            if (s == -1)
              states[k++] = info.nstates;
            else
              states[k++] = s;
          }
        } // pattern loop
      } // subset loop

      int code = beagleSetTipStates(
          info.handle, // Instance number
          t,           // Index of destination compactBuffer
          &states[0]); // Pointer to compact states vector

      if (code != 0)
        throw XStrom(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t + 1) % _data->getTaxonNames()[t] % code % _beagle_error[code]);
      ++t;
    }
  }

  // function to be called if ambiguity is not counted as missing data
  // if _ambiguity_equals_missing is false
  inline void Likelihood::setTipPartials()
  {
    assert(_instances.size() > 0);
    assert(_data);
    Data::state_i one = 1;

    for (auto &info : _instances)
    {
      if (info.nstates != 4)
        throw XStrom(boost::format("This program can handle only 4-state DNA/RNA data. You specified data having %d states for at least one data subset.") % info.nstates);

      std::vector<double> partials(info.nstates * info.npatterns);

      // loop through all rows of the data matrix, setting the tip states for one taxon each row
      unsigned t = 0;
      for (auto &row : _data->getDataMatrix)
      {

        // loop through all subsets assigned to this instance
        unsigned k = 0;
        for (unsigned s : info.subsets)
        {

          // loop through all patterns in this subset
          auto interval = _data->getSubsetBeginEnd(s);
          for (unsigned p = interval.first; p < interval.second; p++)
          {

            // d is the state for taxon t, pattern p (in subset s)
            Data::state_t d = row[p];

            // Handle common nucleotide case separately
            if (info.nstates == 4)
            {
              partials[k++] = d & 1 ? 1.0 : 0.0;
              partials[k++] = d & 2 ? 1.0 : 0.0;
              partials[k++] = d & 4 ? 1.0 : 0.0;
              partials[k++] = d & 8 ? 1.0 : 0.0;
            }
            else
            {
              // This case is for any other data type except nucleotide
              for (unsigned b = 0; b < info.nstates; b++)
              {
                partials[k++] = d & (one << b) ? 1.0 : 0.0;
              }
            }
          }
        }

        int code = beagleSetTipPartials(
            info.handle,   // Instance number
            t,             // Index of destination compactBuffer
            &partials[0]); // Pointer to compact states vector

        if (code != 0)
          throw XStrom(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t + 1) % _data->getTaxonNames()[t] % code % _beagle_error[code]);
        ++t;
      }
    }
  }

  // setPatternPartitionAssignments member function

}
