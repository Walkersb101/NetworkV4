#pragma once

#include <cstdint>

#include "DataOut.hpp"
#include "Enums.hpp"
#include "Network.hpp"
#include "NetworkOut.hpp"

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
    void strainBreakData(const network& _network,
                    double& _maxVal,
                    std::size_t& _count) const;
    auto strainBreak(network& _network) const -> std::vector<std::size_t>;

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
                      std::unique_ptr<dataOut>& _dataOut,
                      std::unique_ptr<networkOut>& _networkOut);

protected:
  StrainType m_strainType;
  BreakType m_breakType;
  dataOut* m_dataOut;
  networkOut* m_networkOut;
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
                    double _tol,
                    bool _errorOnNotSingleBreak);
  ~quasiStaticStrain();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  void evalStrain(network& _network,
                  double _step,
                  double& _maxVal,
                  std::size_t& _count);
  auto converge(network& _baseNetwork,
                network& _networkB,
                double& _a,
                double& _b,
                double& _fa,
                double& _fb,
                size_t& _breakCountB,
                double _tol) -> bool;
  auto findSingleBreak(network& _network) -> size_t;

  auto relaxBreak(network& _network) -> std::size_t;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   size_t _bondIndex) -> std::vector<writeableTypes>;

private:
  double m_maxStrain;

  double m_esp;
  double m_tol;
  bool m_errorOnNotSingleBreak;

  double m_t;

  std::size_t m_strainCount;
};

}  // namespace networkV4
