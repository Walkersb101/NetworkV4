#pragma once

#include <cstdint>

namespace network
{
    enum class bondType : std::uint8_t {
        single,
        sacrificial,
        matrix
    };

class bond
{
public:
  bond();
  bond(std::size_t _src, std::size_t _dst,
       bool _connected);
  bond(std::size_t _src,
       std::size_t _dst,
       bool _connected,
       bondType _type);
  virtual ~bond() = default;
  bond(const bond&) = default;
  auto operator=(const bond&) -> bond& = default;
  bond(bond&&) = default;
  auto operator=(bond&&) -> bond& = default;

  auto src() const -> size_t;
  auto dst() const -> size_t;

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

class WLCBond : public bond //NOLINT
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

}  // namespace network
