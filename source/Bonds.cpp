#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "Bonds.hpp"

networkV4::bond::bond()
    : bond(0, 0, false, 0.0, 0.0, 0.0)
{
}
networkV4::bond::bond(std::size_t _src,
                      std::size_t _dst,
                      bool _connected,
                      double _l0,
                      double _mu,
                      double _lambda)
    : bond(_src,
           _dst,
           _connected,
           networkV4::bondType::single,
           _l0,
           _mu,
           _lambda)
{
}

networkV4::bond::bond(std::size_t _src,
                      std::size_t _dst,
                      bool _connected,
                      networkV4::bondType _type,
                      double _l0,
                      double _mu,
                      double _lambda)
    : m_src(_src)
    , m_dst(_dst)
    , m_connected(_connected)
    , m_type(_type)
    , m_l0(_l0)
    , m_mu(_mu)
    , m_lambda(_lambda)
{
}

auto networkV4::bond::src() const -> std::size_t
{
  return m_src;
}
auto networkV4::bond::dst() const -> std::size_t
{
  return m_dst;
}

auto networkV4::bond::connected() const -> bool
{
  return m_connected;
}
auto networkV4::bond::connected() -> bool&
{
  return m_connected;
}

auto networkV4::bond::type() const -> networkV4::bondType
{
  return m_type;
}

auto networkV4::bond::naturalLength() const -> double
{
  return m_l0;
}

auto networkV4::bond::naturalLength() -> double&
{
  return m_l0;
}

auto networkV4::bond::mu() const -> double
{
  return m_mu;
}

auto networkV4::bond::mu() -> double&
{
  return m_mu;
}

auto networkV4::bond::lambda() const -> double
{
  return m_lambda;
}

auto networkV4::bond::lambda() -> double&
{
  return m_lambda;
}

auto networkV4::bond::force(double _r) const -> double
{
  return m_mu * ((_r / m_l0) - 1.0);
}
auto networkV4::bond::energy(double _r) const -> double
{
  return 0.5 * m_mu * std::pow(_r - m_l0, 2) / m_l0;
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

// ----------------------------------------------------------------------------

networkV4::bonds::bonds()
    : m_bonds()
{
}

void networkV4::bonds::clear()
{
  m_bonds.clear();
}

void networkV4::bonds::resize(std::size_t _size)
{
  m_bonds.resize(_size);
}

void networkV4::bonds::reserve(std::size_t _size)
{
  m_bonds.reserve(_size);
}

void networkV4::bonds::shrink_to_fit()
{
  m_bonds.shrink_to_fit();
}

auto networkV4::bonds::size() const -> std::size_t
{
  return m_bonds.size();
}

auto networkV4::bonds::empty() const -> bool
{
  return m_bonds.empty();
}

void networkV4::bonds::add(const bond& _bond)
{
  m_bonds.emplace_back(_bond);
}

void networkV4::bonds::set(std::size_t _index, const bond& _bond)
{
  m_bonds[_index] = _bond;
}

void networkV4::bonds::remove(std::size_t _index)
{
  m_bonds.erase(m_bonds.begin() + _index);
}

auto networkV4::bonds::get(std::size_t _index) -> networkV4::bond&
{
  return m_bonds[_index];
}

auto networkV4::bonds::get(std::size_t _index) const -> const networkV4::bond&
{
  return m_bonds[_index];
}

auto networkV4::bonds::operator[](std::size_t _index) -> networkV4::bond&
{
  return m_bonds[_index];
}

auto networkV4::bonds::operator[](std::size_t _index) const
    -> const networkV4::bond&
{
  return m_bonds[_index];
}

auto networkV4::bonds::connectedCount() const -> std::size_t
{
  return std::count_if(m_bonds.begin(),
                       m_bonds.end(),
                       [](const auto& b) { return b.connected(); });
}

auto networkV4::bonds::connectedCount(bondType _type) const -> std::size_t
{
  return std::count_if(m_bonds.begin(),
                       m_bonds.end(),
                       [_type](const auto& b) { return b.connected() && b.type() == _type; });
}

void networkV4::bonds::sort()
{
  std::sort(m_bonds.begin(),
            m_bonds.end(),
            [this](const auto& _lhs, const auto& _rhs)
            { return srcDestOrder(_lhs, _rhs); });
}

auto networkV4::bonds::begin() -> std::vector<networkV4::bond>::iterator
{
  return m_bonds.begin();
}

auto networkV4::bonds::begin() const
    -> std::vector<networkV4::bond>::const_iterator
{
  return m_bonds.begin();
}

auto networkV4::bonds::cbegin() const -> std::vector<bond>::const_iterator
{
  return m_bonds.cbegin();
}

auto networkV4::bonds::end() -> std::vector<networkV4::bond>::iterator
{
  return m_bonds.end();
}

auto networkV4::bonds::end() const
    -> std::vector<networkV4::bond>::const_iterator
{
  return m_bonds.end();
}

auto networkV4::bonds::cend() const
    -> std::vector<networkV4::bond>::const_iterator
{
  return m_bonds.cend();
}

inline bool networkV4::bonds::srcDestOrder(const networkV4::bond& _lhs,
                                           const networkV4::bond& _rhs) const
{
  return _lhs.src() < _rhs.src()
      || (_lhs.src() == _rhs.src() && _lhs.dst() < _rhs.dst());
}
