#pragma once

#include <deque>
#include <filesystem>
#include <vector>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Core/Stresses.hpp"
#include "Core/box.hpp"
#include "Misc/Tags/TagMap.hpp"
#include "Misc/Tags/TagStorage.hpp"
#include "Misc/Math/Tensor2.hpp"
#include "Misc/Math/Vector.hpp"

namespace networkV4
{

using breakInfo = std::tuple<bonded::BondInfo, bonded::bondTypes, bonded::breakTypes, Utils::Tags::tagFlags>;
using bondQueue = std::deque<breakInfo>;

class network
{
public:
  network() = delete;
  network(const box& _box, const size_t _N, const size_t _B);

public:
  auto getNodes()
      -> nodes&;  // TODO: make so Input and intergrators can access nodes
  auto getNodes() const -> const nodes&;

  auto getBonds() -> bonded::bonds&;  // TODO: make so Input and intergrators
                                      // can access nodes
  auto getBonds() const -> const bonded::bonds&;

  auto getEnergy() const -> const double;

  auto getRestBox() const -> const box;
  auto getBox() const -> const box;

  auto getTags() -> Utils::Tags::tagMap&;
  auto getTags() const -> const Utils::Tags::tagMap&;

public:
  auto getStresses()
      -> stresses&;  // TODO: make so Input and intergrators can access nodes
  auto getStresses() const -> const stresses&;

  auto getBreakQueue()
      -> bondQueue&;  // TODO: make so Input and intergrators can access nodes
  auto getBreakQueue() const -> const bondQueue&;

public:
  double getShearStrain() const;
  auto getElongationStrain() const -> Utils::Math::vec2d;

  void shear(double _step);
  void setBox(const box& _box);

  void wrapNodes();

public:
  template <bool _evalBreak = false, bool _evalStress = false>
  void computeForces();

  auto computeEnergy() -> double;
  void computeBreaks();

private:
  void evalBreak(const Utils::Math::vec2d& _dist,
                 const bonded::BondInfo& _binfo,
                 bonded::bondTypes& _type,
                 bonded::breakTypes& _break,
                 Utils::Tags::tagFlags& _tags);

  template <bool _evalBreak = false>
  void applyforce(const bonded::BondInfo& _binfo,
                  const Utils::Math::vec2d& _dist,
                  const Utils::Math::vec2d& _force,
                    const Utils::Tags::tagFlags& _tags);

#if defined(_OPENMP)
private:
  template <bool _evalBreak = false, bool _evalStress = false>
  void computePass(auto _parts);
#endif

private:
  box m_box;
  box m_restbox;

  double m_energy;
  stresses m_stresses;

  nodes m_nodes;
  bonded::bonds m_bonds;

  bondQueue m_breakQueue;

  Utils::Tags::tagMap m_tags;
};

inline void merge(bondQueue& _b1, const bondQueue& _b2)
{
  _b1.insert(_b1.end(), _b2.begin(), _b2.end());
}

#pragma omp declare reduction(+ : bondQueue : merge(omp_out, omp_in))

}  // namespace networkV4
