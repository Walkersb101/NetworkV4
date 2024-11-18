#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Protocol.hpp"

#include "Core/Intergrator.hpp"
#include "Core/Network.hpp"
#include "Core/Roots.hpp"
#include "IO/DataOut.hpp"
#include "Misc/Config.hpp"
#include "Misc/Tools.hpp"

networkV4::protocol::protocol() {}

networkV4::protocol::protocol(StrainType _strainType,
                              std::unique_ptr<BreakTypes>& _breakType)
    : m_strainType(_strainType)
    , m_breakProtocol(_breakType.get())
    , m_dataOut(nullptr)
    , m_networkOut(nullptr)
{
}

networkV4::protocol::protocol(StrainType _strainType)
    : m_strainType(_strainType)
    , m_breakProtocol(nullptr)
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
      return _network.getElongationStrain();
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
      _network.elongate(_step);
      break;
    }
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::protocol::run(network& _network) {}

void networkV4::protocol::initIO(const network& _network,
                                 std::unique_ptr<dataOut>& _dataOut,
                                 std::unique_ptr<networkOut>& _networkOut)
{
}