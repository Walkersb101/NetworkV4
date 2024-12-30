#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Protocol.hpp"

#include "Integration/Integrator.hpp"

#include "Misc/Config.hpp"
#include "Misc/Tools.hpp"
#include "Misc/Roots.hpp"

networkV4::protocol::protocol() {}

networkV4::protocol::protocol(StrainType _strainType)
    : m_strainType(_strainType)
    , m_dataOut(nullptr)
    , m_networkOut(nullptr)
{
}

networkV4::protocol::~protocol() {}

double networkV4::protocol::getStrain(const network& _network) const
{
  switch (m_strainType) {
    case StrainType::Shear:
      return _network.getShearStrain();
    case StrainType::Elongation:
      return 0.0;//_network.getElongationStrain(); //TODO: Fix this
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::protocol::strain(network& _network, double _step)
{
  switch (m_strainType) {
    case StrainType::Shear: {
        _network.shear(_step);
      break;
    }
    case StrainType::Elongation: {
      //_network.elongate(_step); //TODO: Fix this
      break;
    }
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::protocol::run(network& _network) {}

void networkV4::protocol::initIO(const network& _network,
                                 std::unique_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
                                 std::unique_ptr<IO::networkDumps::networkDump>& _networkOut)
{
}