#pragma once

#include <cstdint>
#include <memory>
#include <vector>

namespace network
{
enum class bondType : std::uint8_t
{
  single,
  sacrificial,
  matrix
};

class bond
{
public:
  bond();
  bond(std::size_t _src, std::size_t _dst, bool _connected);
  bond(std::size_t _src, std::size_t _dst, bool _connected, bondType _type);
  virtual ~bond() = default;
  bond(const bond&) = default;
  auto operator=(const bond&) -> bond& = default;
  bond(bond&&) = default;
  auto operator=(bond&&) -> bond& = default;

  auto src() const -> std::size_t;
  auto dst() const -> std::size_t;

  auto connected() const -> bool;
  auto connected() -> bool&;

  auto type() const -> bondType;

  virtual auto force(double _r) const -> double = 0;
  virtual auto energy(double _r) const -> double = 0;

  // virtual void serialize(std::ostream& outStream) const = 0;
  // virtual void deserialize(std::istream& inStream) = 0;

private:
  std::size_t m_src;
  std::size_t m_dst;
  bool m_connected;
  bondType m_type;
};

class WLCBond : public bond  // NOLINT
{
public:
  WLCBond();
  WLCBond(std::size_t _src,
          std::size_t _dst,
          bool _connected,
          bondType _type,
          double _l0,
          double _mu);
  ~WLCBond() override = default;
  WLCBond(const WLCBond&) = default;
  auto operator=(const WLCBond&) -> WLCBond& = default;
  WLCBond(WLCBond&&) = default;
  auto operator=(WLCBond&&) -> WLCBond& = default;

  auto force(double _r) const -> double override;
  auto energy(double _r) const -> double override;

  // void serialize(std::ostream& outStream) const override;
  // void deserialize(std::istream& inStream) override;

private:
  double m_l0;
  double m_mu;
};

// ------------------------------------------------------------

class bonds
{
public:
  bonds();
  ~bonds() = default;

public:
  void clear();
  void resize(std::size_t _size);
  void reserve(std::size_t _size);
  void shrink_to_fit();

  auto size() const -> std::size_t;
  auto empty() const -> bool;

public:
  void add(const bond& _bond);
  void add(bond&& _bond);

  void set(std::size_t _index, const bond& _bond);
  void set(std::size_t _index, bond&& _bond);

  void remove(std::size_t _index);

  auto get(std::size_t _index) -> bond&;
  auto get(std::size_t _index) const -> const bond&;

  auto operator[](std::size_t _index) -> bond&;
  auto operator[](std::size_t _index) const -> const bond&;

public:
  auto connectedCount() const -> std::size_t;

public:
  void sort();

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
  std::vector<std::unique_ptr<bond>> m_bonds;
};
}  // namespace network
