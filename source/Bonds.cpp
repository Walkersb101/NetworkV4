#include <cstddef>
#include <cmath>

#include "Bonds.hpp"

network::bond::bond()
    : m_src(0)
    , m_dst(0)
    , m_connected(false)
    , m_type(network::bondType::single)
{
}
network::bond::bond(std::size_t _src,  // NOLINT
                    std::size_t _dst,  // NOLINT
                    bool _connected)
    : m_src(_src)
    , m_dst(_dst)
    , m_connected(_connected)
    , m_type(network::bondType::single)
{
}
network::bond::bond(std::size_t _src, // NOLINT
                    std::size_t _dst, // NOLINT
                    bool _connected,
                    network::bondType _type)
    : m_src(_src)
    , m_dst(_dst)
    , m_connected(_connected)
    , m_type(_type)
{
}

auto network::bond::src() const -> size_t
{
  return m_src;
}
auto network::bond::dst() const -> size_t
{
  return m_dst;
}

auto network::bond::connected() const -> bool
{
  return m_connected;
}
auto network::bond::connected() -> bool&
{
  return m_connected;
}

auto network::bond::type() const -> network::bondType
{
  return m_type;
}

/*
void bond::serialize(std::ostream& outStream) const {
    outStream.write((char*)&m_connected, sizeof(bool));
    outStream.write((char*)&m_src, sizeof(size_t));
    outStream.write((char*)&m_dst, sizeof(size_t));
}

void bond::deserialize(std::istream& inStream) {
    inStream.read((char*)&m_connected, sizeof(bool));
    inStream.read((char*)&m_src, sizeof(size_t));
    inStream.read((char*)&m_dst, sizeof(size_t));
}
*/

network::WLCBond::WLCBond()
    : m_l0(0.0)
    , m_mu(0.0)
{
}

network::WLCBond::WLCBond(std::size_t _src,
                          std::size_t _dst,
                          bool _connected,
                          network::bondType _type,
                          double _l0,
                          double _mu)
    : network::bond::bond(_src, _dst, _connected, _type)
    , m_l0(_l0)
    , m_mu(_mu)
{
}

auto network::WLCBond::force(double _r) const -> double
{
    return m_mu * ((_r / m_l0) - 1.0);
}
auto network::WLCBond::energy(double _r) const -> double
{
    return 0.5 * m_mu * std::pow(_r - m_l0, 2) / m_l0;
}
