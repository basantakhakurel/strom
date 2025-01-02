#pragma once

#include "updater.hpp"

namespace strom
{

  class TreeLengthUpdater : public Updater
  {

  public:
    typedef std::shared_ptr<TreeLengthUpdater> SharedPtr;

    TreeLengthUpdater();
    ~TreeLengthUpdater();

    virtual void clear();
    virtual void proposeNewState();
    virtual void revert();

    virtual double calcLogPrior();

    void pullFromModel();
    void pushToModel() const;

  private:
    double _prev_point;
    double _curr_point;
  };

  /**
   * @brief Constructor for TreeLengthUpdater.
   *
   * This constructor for TreeLengthUpdater sets the name of the
   * updater to "Tree Length", and clears the Updater by calling
   * the clear() function inherited from the Updater class.
   */
  inline TreeLengthUpdater::TreeLengthUpdater()
  {
    // std::cout << "Creating a TreeLengthUpdater..." << std::endl;
    clear();
    _name = "Tree Length";
  }

  /**
   * @brief Destructor for TreeLengthUpdater.
   *
   * This destructor for TreeLengthUpdater simply outputs a message
   * to the console indicating that a TreeLengthUpdater is being
   * destroyed.
   */
  inline TreeLengthUpdater::~TreeLengthUpdater()
  {

    // std::cout << "Destroying a TreeLengthUpdater..." << std::endl;
  }

  /**
   * @brief Resets the state of the TreeLengthUpdater.
   *
   * This function clears the state of the TreeLengthUpdater by setting
   * the previous and current points to 0.0 and calling the reset function.
   * It also invokes the clear function of the base class Updater to reset
   * its own state.
   *
   */
  inline void TreeLengthUpdater::clear()
  {
    Updater::clear();
    _prev_point = 0.0;
    _curr_point = 0.0;
    reset();
  }

  /**
   * @brief Calculates the prior probability of the current tree length.
   *
   * @return The prior probability of the current tree length as a double.
   *
   * This function returns the prior probability of the current tree length
   * by calling the calcLogEdgeLengthPrior() function of the base class Updater.
   */
  inline double TreeLengthUpdater::calcLogPrior()
  {
    return Updater::calcLogEdgeLengthPrior();
  }

  /**
   * @brief Pulls the current tree length from the model.
   *
   * This function pulls the current tree length from the model by calling
   * the calcTreeLength() function of the TreeManipulator class and storing
   * the result in the _curr_point member variable.
   */
  inline void TreeLengthUpdater::pullFromModel()
  {
    _curr_point = _tree_manipulator->calcTreeLength();
  }

  /**
   * @brief Pushes the current tree length to the model.
   *
   * @param[in] scaler The scaling factor to apply to the tree lengths.
   *
   * This function pushes the current tree length to the model by scaling
   * all edge lengths in the tree by the given scaler.
   */
  inline void TreeLengthUpdater::pushToModel() const
  {
    double scaler = _curr_point / _prev_point;
    _tree_manipulator->scaleAllEdgeLengths(scaler);
  }

  /**
   * @brief Proposes a new state for the tree length.
   *
   * This function proposes a new tree length by multiplying the current tree length
   * by a random multiplier derived from a uniform distribution. It calculates the
   * log of the Hastings ratio under the GammaDir parameterization and updates the
   * current tree length in the model. The proposal invalidates all transition matrices
   * and partials.
   */
  inline void TreeLengthUpdater::proposeNewState()
  {
    // Save copy of _curr_point in case revert is necessary.
    pullFromModel();
    _prev_point = _curr_point;

    // Let _curr_point be proposed value
    double m = exp(_lambda * (_lot->uniform() - 0.5));
    _curr_point = m * _prev_point;
    pushToModel();

    // calculate log of Hastings ratio under GammaDir parameterization
    _log_hastings_ratio = log(m);

    // This proposal invalidates all transition matrices and partials
    _tree_manipulator->selectAllPartials();
    _tree_manipulator->selectAllTMatrices();
  }

  /**
   * @brief Reverts the current state of the tree length to the previous state.
   *
   * This function swaps the current and previous tree lengths and pushes the
   * previous tree length to the model. It is called when a proposed tree length
   * is rejected.
   */
  inline void TreeLengthUpdater::revert()
  {
    // swap _curr_point and _prev_point so that edge length scaler
    // in pushCurrentStateToModel will be correctly calculated
    double tmp = _curr_point;
    _curr_point = _prev_point;
    _prev_point = tmp;
    pushToModel();
  }
}
