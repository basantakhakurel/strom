#pragma once

#include <memory>
#include <boost/format.hpp>
#include "lot.hpp"
#include "data.hpp"
#include "tree.hpp"
#include "model.hpp"
#include "likelihood.hpp"
#include "tree_manip.hpp"
#include "updater.hpp"
#include "gamma_ratevar_updater.hpp"
#include "omega_updater.hpp"
#include "pinvar_updater.hpp"
#include "statefreq_updater.hpp"
#include "exchangeability_updater.hpp"
#include "subset_relrate_updater.hpp"
#include "tree_updater.hpp"
#include "tree_length_updater.hpp"

namespace strom
{

  class Chain
  {

    friend class Likelihood;

  public:
    typedef std::vector<Updater::SharedPtr> updater_vect_t;
    typedef std::shared_ptr<Chain> SharedPtr;

    Chain();
    ~Chain();

    void clear();

    void startTuning();
    void stopTuning();

    void setTreeFromNewick(std::string &newick);
    unsigned createUpdaters(Model::SharedPtr model, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood);

    TreeManip::SharedPtr getTreeManip();
    Model::SharedPtr getModel();
    double getLogLikelihood() const;

    void setHeatingPower(double p);
    double getHeatingPower() const;

    void setChainIndex(unsigned idx);
    double getChainIndex() const;

    Updater::SharedPtr findUpdaterByName(std::string name);
    std::vector<std::string> getUpdaterNames() const;
    std::vector<double> getAcceptPercentages() const;
    std::vector<unsigned> getNumUpdates() const;
    std::vector<double> getLambdas() const;
    void setLambdas(std::vector<double> &v);

    double calcLogLikelihood() const;
    double calcLogJointPrior() const;

    void start();
    void stop();
    void nextStep(int iteration);

  private:
    Model::SharedPtr _model;
    Lot::SharedPtr _lot;
    TreeManip::SharedPtr _tree_manipulator;

    updater_vect_t _updaters;

    unsigned _chain_index;
    double _heating_power;
    double _log_likelihood;
  };

  // constructo and destructor
  inline Chain::Chain()
  {
    // std::cout << "Chain being created" << std::endl;
    clear();
  }

  inline Chain::~Chain()
  {
    // std::cout << "Chain being destroyed" << std::endl;
  }

  // clear member function
  inline void Chain::clear()
  {
    _log_likelihood = 0.0;
    _updaters.clear();
    _chain_index = 0;
    setHeatingPower(1.0);
    startTuning();
  }

  // functions to start and stop tuning
  inline void Chain::startTuning()
  {
    for (auto u : _updaters)
      u->setTuning(true);
  }

  inline void Chain::stopTuning()
  {
    for (auto u : _updaters)
      u->setTuning(false);
  }

  // function to set the tree from a newick string
  inline void Chain::setTreeFromNewick(std::string &newick)
  {
    assert(_updaters.size() > 0);
    if (!_tree_manipulator)
      _tree_manipulator.reset(new TreeManip);
    _tree_manipulator->buildFromNewick(newick, false, false);
    for (auto u : _updaters)
      u->setTreeManip(_tree_manipulator);
  }

  /**
   * @brief Create updaters for various model parameters and add them to the chain.
   *
   * This function initializes and configures updaters for different model
   * parameters such as state frequencies, exchangeabilities, rate variances,
   * pinvar, omega, subset relative rates, and tree parameters. Each updater is
   * assigned a likelihood, a random number generator (lot), lambda value, target
   * acceptance rate, prior parameters, and a weight. The sum of the weights is
   * calculated and used to determine the probability for each updater.
   *
   * @param model Shared pointer to the Model object containing parameters.
   * @param lot Shared pointer to the Lot object for random number generation.
   * @param likelihood Shared pointer to the Likelihood object.
   * @return The total number of updaters created and added to the chain.
   */
  inline unsigned Chain::createUpdaters(Model::SharedPtr model, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood)
  {
    _model = model;
    _lot = lot;
    _updaters.clear();

    double wstd = 1.0;
    double wtreelength = 1.0;
    double wtreetopology = 19.0;
    double sum_weights = 0.0;

    // Add state frequency parameter updaters to _updaters
    Model::state_freq_params_t &statefreq_shptr_vect = _model->getStateFreqParams();
    for (auto statefreq_shptr : statefreq_shptr_vect)
    {
      Updater::SharedPtr u = StateFreqUpdater::SharedPtr(new StateFreqUpdater(statefreq_shptr));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(1.0);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters(std::vector<double>(statefreq_shptr->getStateFreqsSharedPtr()->size(), 1.0));
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add exchangeability parameter updaters to _updaters
    Model::exchangeability_params_t &exchangeability_shptr_vect = _model->getExchangeabilityParams();
    for (auto exchangeability_shptr : exchangeability_shptr_vect)
    {
      Updater::SharedPtr u = ExchangeabilityUpdater::SharedPtr(new ExchangeabilityUpdater(exchangeability_shptr));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(1.0);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add rate variance parameter updater to _updaters
    Model::ratevar_params_t &ratevar_shptr_vect = _model->getRateVarParams();
    for (auto ratevar_shptr : ratevar_shptr_vect)
    {
      Updater::SharedPtr u = GammaRateVarUpdater::SharedPtr(new GammaRateVarUpdater(ratevar_shptr));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(1.0);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({1.0, 1.0});
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add pinvar parameter updaters to _updaters
    Model::pinvar_params_t &pinvar_shptr_vect = _model->getPinvarParams();
    for (auto pinvar_shptr : pinvar_shptr_vect)
    {
      Updater::SharedPtr u = PinvarUpdater::SharedPtr(new PinvarUpdater(pinvar_shptr));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(0.5);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({1.0, 1.0});
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add omega parameter updaters to _updaters
    Model::omega_params_t &omega_shptr_vect = _model->getOmegaParams();
    for (auto omega_shptr : omega_shptr_vect)
    {
      Updater::SharedPtr u = OmegaUpdater::SharedPtr(new OmegaUpdater(omega_shptr));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(1.0);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({1.0, 1.0});
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add subset relative rate parameter updater to _updaters
    if (!_model->isFixedSubsetRelRates())
    {
      Updater::SharedPtr u = SubsetRelRateUpdater::SharedPtr(new SubsetRelRateUpdater(_model));
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(1.0);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters(std::vector<double>(_model->getNumSubsets(), 1.0));
      u->setWeight(wstd);
      sum_weights += wstd;
      _updaters.push_back(u);
    }

    // Add tree updater and tree length updater to _updaters
    if (!_model->isFixedTree())
    {
      double tree_length_shape = 1.0;
      double tree_length_scale = 10.0;
      double dirichlet_param = 1.0;

      Updater::SharedPtr u = TreeUpdater::SharedPtr(new TreeUpdater());
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(0.5);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
      u->setWeight(wtreetopology);
      sum_weights += wtreetopology;
      _updaters.push_back(u);

      u = TreeLengthUpdater::SharedPtr(new TreeLengthUpdater());
      u->setLikelihood(likelihood);
      u->setLot(lot);
      u->setLambda(0.2);
      u->setTargetAcceptanceRate(0.3);
      u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
      u->setWeight(wtreelength);
      sum_weights += wtreelength;
      _updaters.push_back(u);
    }

    for (auto u : _updaters)
    {
      u->calcProb(sum_weights);
    }

    return (unsigned)_updaters.size();
  }

  // accessor for the tree manipulator data member
  inline TreeManip::SharedPtr Chain::getTreeManip()
  {
    return _tree_manipulator;
  }

  // accessor function for the model data member
  inline Model::SharedPtr Chain::getModel()
  {
    return _model;
  }

  // accessor function for the log likelihood
  inline double Chain::getLogLikelihood() const
  {
    return _log_likelihood;
  }

  // functions to set and get the heating power
  inline double Chain::getHeatingPower() const
  {
    return _heating_power;
  }

  inline void Chain::setHeatingPower(double p)
  {
    _heating_power = p;
    for (auto u : _updaters)
      u->setHeatingPower(p);
  }

  // functions to set and get the chain index
  inline double Chain::getChainIndex() const
  {
    return _chain_index;
  }

  inline void Chain::setChainIndex(unsigned idx)
  {
    _chain_index = idx;
  }

  // function to loop through all the updaters and return the ones with the name
  inline Updater::SharedPtr Chain::findUpdaterByName(std::string name)
  {
    Updater::SharedPtr retval = nullptr;
    for (auto u : _updaters)
    {
      if (u->getUpdaterName() == name)
      {
        retval = u;
        break;
      }
    }
    assert(retval != nullptr);
    return retval;
  }

  // function to return the names of the updaters
  inline std::vector<std::string> Chain::getUpdaterNames() const
  {
    std::vector<std::string> v;
    for (auto u : _updaters)
      v.push_back(u->getUpdaterName());
    return v;
  }

  // function to provide the performance of each updater
  inline std::vector<double> Chain::getAcceptPercentages() const
  {
    std::vector<double> v;
    for (auto u : _updaters)
      v.push_back(u->getAcceptPct());
    return v;
  }

  // function to report the number of attempts to update each parameter
  inline std::vector<unsigned> Chain::getNumUpdates() const
  {
    std::vector<unsigned> v;
    for (auto u : _updaters)
      v.push_back(u->getNumUpdates());
    return v;
  }

  // functions to get and set the lambdas
  inline std::vector<double> Chain::getLambdas() const
  {
    std::vector<double> v;
    for (auto u : _updaters)
      v.push_back(u->getLambda());
    return v;
  }

  /**
   * Sets the lambda values for all of the updaters.
   *
   * Each element of the input vector v is used to set the lambda value for the
   * corresponding updater (as returned by getUpdaterNames()). The size of v
   * must match the size of the vector returned by getUpdaterNames(). The
   * lambda values are used to scale the proposals for each updater.
   */
  inline void Chain::setLambdas(std::vector<double> &v)
  {
    assert(v.size() == _updaters.size());
    unsigned index = 0;
    for (auto u : _updaters)
    {
      u->setLambda(v[index++]);
    }
  }

  // function to calculate the log likelihood
  inline double Chain::calcLogLikelihood() const
  {
    return _updaters[0]->calcLogLikelihood();
  }

  // function to calculate the log prior for the model parameter that it is responsible for updating
  inline double Chain::calcLogJointPrior() const
  {
    double lnP = 0.0;
    for (auto u : _updaters)
    {
      if (u->_name != "Tree Length")
        lnP += u->calcLogPrior();
    }
    return lnP;
  }

  // functions to start and stop the chain
  inline void Chain::start()
  {
    _tree_manipulator->selectAllPartials();
    _tree_manipulator->selectAllTMatrices();
    _log_likelihood = calcLogLikelihood();
  }

  inline void Chain::stop()
  {
  }

  // function to move to the next step in MCMC chain
  inline void Chain::nextStep(int iteration)
  {
    assert(_lot);
    double u = _lot->uniform();
    double cumprob = 0.0;
    unsigned i = 0;
    for (auto updater : _updaters)
    {
      cumprob += updater->_prob;
      if (u <= cumprob)
        break;
      i++;
    }
    assert(i < _updaters.size());
    _log_likelihood = _updaters[i]->update(_log_likelihood);
  }
}
