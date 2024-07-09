#pragma once

#include <cstdint>

#include "DataOut.hpp"
#include "Network.hpp"
#include "Enums.hpp"

namespace networkV4
{

class protocol
{
public:
  protocol();
  protocol(StrainType _strainType, BreakType _breakType);
  virtual ~protocol();

public:
  double getStrain(const network& _network) const;
  void strain(network& _network, double _step);

public:
  void breakData(const network& _network,
                 BreakType _type,
                 double& _maxVal,
                 std::size_t& _count);
  auto breakBonds(network& _network,
                  BreakType _type) const -> std::vector<std::size_t>;

public:
  virtual void run(network& _network);

  virtual void initIO(const network& _network,
                      std::unique_ptr<dataOut>& _dataOut);

protected:
  StrainType m_strainType;
  BreakType m_breakType;
  dataOut* m_dataOut;
};

class quasiStaticStrain : public protocol
{
public:
  quasiStaticStrain();
  quasiStaticStrain(double _maxStrain,
                    StrainType _strainType,
                    BreakType _breakType);
  quasiStaticStrain(double _maxStrain,
                    StrainType _strainType,
                    BreakType _breakType,
                    double _esp,
                    double _tol);
  ~quasiStaticStrain();

public:
  void run(network& _network) override;
  void initIO(const network& _network, std::unique_ptr<dataOut>& _dataOut);

private:
  void evalStrain(network& _network,
                  double _step,
                  double& _maxVal,
                  std::size_t& _count);
  auto findSingleBreak(network& _network) -> bool;
  auto relaxBreak(network& _network) -> std::size_t;

  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network, size_t _bondIndex) -> std::vector<writeableTypes>;

private:
  double m_maxStrain;

  double m_esp;
  double m_tol;

  double m_t;

  std::size_t m_strainCount;
};

}  // namespace networkV4
