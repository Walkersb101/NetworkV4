#pragma once

#include "Misc/Vec2.hpp"

namespace networkV4
{

class box
{
  box(const double _Lx, const double _Ly, const double _xy)
      : m_Lx(_Lx)
      , m_Ly(_Ly)
      , m_xy(_xy)
      , m_invLx(1.0 / _Lx)
      , m_invLy(1.0 / _Ly)
      , m_invxy(-_xy / (_Lx * _Ly))
  {
  }

public:
  auto getLx() const -> double { return m_Lx; }
  auto getLy() const -> double { return m_Ly; }
  auto getxy() const -> double { return m_xy; }

  void setLx(const double _Lx)
  {
    m_Lx = _Lx;
    m_invLx = 1.0 / _Lx;
    m_invxy = -m_xy / (_Lx * m_Ly);
    m_halfLx = 0.5 * _Lx;
  }
  void setLy(const double _Ly)
  {
    m_Ly = _Ly;
    m_invLy = 1.0 / _Ly;
    m_invxy = -m_xy / (m_Lx * _Ly);
    m_halfLy = 0.5 * _Ly;
  }

  void setxy(const double _xy)
  {
    m_xy = _xy;
    m_invxy = -_xy / (m_Lx * m_Ly);
  }

public:
  inline auto lambda2x(const Utils::vec2d& _pos) const -> Utils::vec2d
  {
    return Utils::vec2d(m_Lx * _pos.x + m_xy * _pos.y, m_Ly * _pos.y);
  };

  inline auto x2Lambda(const Utils::vec2d& _pos) const -> Utils::vec2d
  {
    return Utils::vec2d(m_invLx * _pos.x + m_invxy * _pos.y, m_invLy * _pos.y);
  };

public:
  inline auto wrapPosition(const Utils::vec2d& _pos) const -> Utils::vec2d
  {
    Utils::vec2d lambda = x2Lambda(_pos);
    while (lambda.x > 1.0)
      lambda.x -= 1.0;
    while (lambda.x < 0.0)
      lambda.x += 1.0;
    while (lambda.y > 1.0)
      lambda.y -= 1.0;
    while (lambda.y < 0.0)
      lambda.y += 1.0;
    return lambda2x(lambda);
  }

  inline auto minDist(const Utils::vec2d& _pos1,
                      const Utils::vec2d& _pos2) const -> Utils::vec2d
  {
    Utils::vec2d dist = _pos2 - _pos1;
    while (std::abs(dist.y) > m_halfLy) {
      if (dist.y > 0.0) {
        dist.y -= m_Ly;
        dist.x -= m_xy;
      } else {
        dist.y += m_Ly;
        dist.x += m_xy;
      }
    }
    while (std::abs(dist.x) > m_halfLx) {
      if (dist.x > 0.0)
        dist.x -= m_Lx;
      else
        dist.x += m_Lx;
    }
    return dist;
  }

private:
  double m_Lx;
  double m_Ly;
  double m_xy;

  double m_invLx;
  double m_invLy;
  double m_invxy;

  double m_halfLx;
  double m_halfLy;
};

}  // namespace networkV4