#pragma once

#include <cstdint>

#include "Core/Network.hpp"
#include "IO/DataOut.hpp"
#include "IO/NetworkOut.hpp"
#include "Misc/Enums.hpp"
#include "Protocols/Protocol.hpp"

namespace networkV4
{

class propogator : public protocol
{
public:
  propogator();
  propogator(std::vector<double> _strains, StrainType _strainType);
  propogator(const std::vector<double>& _strains,
             StrainType _strainType,
             bondType _breakType,
             double _esp,
             double _tol,
             double _maxStep,
             bool _respectLambda);
  ~propogator();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  void runLambda(network& _network);
  void runStrain(network& _network);

private:
  void evalStrain(network& _network, double _strain);
  void relax(network& _network);

  auto getMostStrained(network& _network, bondType _type) -> size_t;
  void breakMostStrained(network& _network, bondType _type);

  auto findSingleBreak(network& _network) -> bool;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   const bond& _bond) -> std::vector<writeableTypes>;

private:
  std::vector<double> m_strains = {};
  bondType m_breakType = bondType::any;
  bool m_respectLambda = false;

  double m_esp = config::integrators::adaptiveIntegrator::esp;
  double m_tol = config::integrators::miminizer::tol;
  double m_maxStep = 0.1;

  size_t m_strainCount = 0;
};

}  // namespace networkV4
