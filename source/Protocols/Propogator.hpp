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
             double _maxStep);
  ~propogator();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  void evalStrain(network& _network, double _strain);
  void breakMostStrained(network& _network);
  void breakMostStrained(network& _network, bondType _type);
  void relax(network& _network);

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   size_t _bondIndex) -> std::vector<writeableTypes>;

private:
  std::vector<double> m_strains;
  bondType m_breakType;

  double m_esp;
  double m_tol;
  double m_maxStep;

  size_t m_strainCount;
};

}  // namespace networkV4
