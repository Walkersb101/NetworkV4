#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "Misc/Enums.hpp"

namespace networkV4
{
class bond
{
public:
  bond();
  bond(std::size_t _src,
       std::size_t _dst,
       bool _connected,
       double _l0,
       double _mu,
       double _lambda);
  bond(std::size_t _src,
       std::size_t _dst,
       bool _connected,
       bondType _type,
       double _l0,
       double _mu,
       double _lambda);
  // bond(const bond&) = default;
  // auto operator=(const bond&) -> bond& = default;
  // bond(bond&&) = default;
  // auto operator=(bond&&) -> bond& = default;

  auto src() const -> std::size_t;
  auto dst() const -> std::size_t;

  auto connected() const -> bool;
  auto connected() -> bool&;

  auto type() const -> bondType;

  auto naturalLength() const -> double;
  auto naturalLength() -> double&;

  auto mu() const -> double;
  auto mu() -> double&;

  auto lambda() const -> double;
  auto lambda() -> double&;

  auto force(double _r) const -> double;
  auto energy(double _r) const -> double;

  // virtual void serialize(std::ostream& outStream) const = 0;
  // virtual void deserialize(std::istream& inStream) = 0;

private:
  std::size_t m_src;
  std::size_t m_dst;
  bool m_connected;
  bondType m_type;
  double m_l0;
  double m_mu;
  double m_lambda;
};

// ------------------------------------------------------------

class bonds
{
public:
  bonds();
  ~bonds() = default;
  bonds(const bonds&) = default;
  auto operator=(const bonds&) -> bonds& = default;
  bonds(bonds&&) = default;

public:
  void clear();
  void resize(std::size_t _size);
  void reserve(std::size_t _size);
  void shrink_to_fit();

  auto size() const -> std::size_t;
  auto empty() const -> bool;

public:
  void add(const bond& _bond);

  void set(std::size_t _index, const bond& _bond);

  void remove(std::size_t _index);

  auto get(std::size_t _index) -> bond&;
  auto get(std::size_t _index) const -> const bond&;

  auto operator[](std::size_t _index) -> bond&;
  auto operator[](std::size_t _index) const -> const bond&;

public:
  auto connectedCount() const -> std::size_t;
  auto connectedCount(bondType _type) const -> std::size_t;

public:
  void sort();

public:
    auto getConnectedVector() const -> std::vector<bool>;
    auto getNaturalLengthVector() const -> std::vector<double>;
    auto getMuVector() const -> std::vector<double>;
    auto getLambdaVector() const -> std::vector<double>;

public:
    void setConnectedVector(const std::vector<bool>& _connected);
    void setNaturalLengthVector(const std::vector<double>& _l0);
    void setMuVector(const std::vector<double>& _mu);
    void setLambdaVector(const std::vector<double>& _lambda);

public:
  auto begin() -> std::vector<bond>::iterator;
  auto begin() const -> std::vector<bond>::const_iterator;
  auto cbegin() const -> std::vector<bond>::const_iterator;

  auto end() -> std::vector<bond>::iterator;
  auto end() const -> std::vector<bond>::const_iterator;
  auto cend() const -> std::vector<bond>::const_iterator;

private:
  inline bool srcDestOrder(const bond& _lhs, const bond& _rhs) const;

private:
  std::vector<bond> m_bonds;
};
}  // namespace networkV4
