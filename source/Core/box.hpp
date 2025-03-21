#pragma once

#include <stdexcept>
#include <vector>

#include "Misc/Math/Vector.hpp"

namespace networkV4
{

class box
{
public:
  box(const double _Lx, const double _Ly, const double _xy = 0.0)
      : m_Lx(_Lx)
      , m_Ly(_Ly)
      , m_xy(_xy)
      , m_halfLx(0.5 * _Lx)
      , m_halfLy(0.5 * _Ly)
      , m_area(_Lx * _Ly)
  {
    if (_Lx <= 0.0 || _Ly <= 0.0) {
      throw std::runtime_error("box::box: invalid box size");
    }
    updateDependent();
  }
  box(const Utils::Math::vec2d& _domain, const double _xy = 0.0)
      : box(_domain[0], _domain[1], _xy)
  {
  }

public:
  auto getLx() const -> double { return m_Lx; }
  auto getLy() const -> double { return m_Ly; }
  auto getxy() const -> double { return m_xy; }

  void setLx(const double _Lx)
  {
    if (_Lx <= 0.0) {
      throw std::runtime_error("box::setLx: invalid box size");
    }

    m_Lx = _Lx;
    updateDependent();
  }
  void setLy(const double _Ly)
  {
    if (_Ly <= 0.0) {
      throw std::runtime_error("box::setLy: invalid box size");
    }

    m_Ly = _Ly;
    updateDependent();
  }

  void setxy(const double _xy)
  {
    m_xy = _xy;
    updateDependent();
  }

public:
  auto shearStrain() const -> const double { return m_xy / m_Ly; }
  auto area() const -> const double { return m_area; }
  auto invArea() const -> const double { return m_invArea; }
  auto getDomain() const -> const Utils::Math::vec2d
  {
    return {m_Lx, m_Ly};
  }
  auto getBox() const -> std::vector<std::vector<double>> const
  {
    return {{m_Lx, 0.0}, {m_xy, m_Ly}};
  }

public:
  inline auto lambda2x(const Utils::Math::vec2d& _pos) const -> Utils::Math::vec2d
  {
    return {m_Lx * _pos[0] + m_xy * _pos[1], m_Ly * _pos[1]};
  };

  inline auto x2Lambda(const Utils::Math::vec2d& _pos) const -> Utils::Math::vec2d
  {
    return {m_invLx * _pos[0] + m_invxy * _pos[1], m_invLy * _pos[1]};
  };

public:
  inline auto wrapPosition(const Utils::Math::vec2d& _pos) const -> Utils::Math::vec2d
  {
    Utils::Math::vec2d lambda = x2Lambda(_pos);
    while (lambda[0] > 1.0)
      lambda[0] -= 1.0;
    while (lambda[0] < 0.0)
      lambda[0] += 1.0;
    while (lambda[1] > 1.0)
      lambda[1] -= 1.0;
    while (lambda[1] < 0.0)
      lambda[1] += 1.0;
    return lambda2x(lambda);
  }

  inline auto minDist(const Utils::Math::vec2d& _pos1,
                      const Utils::Math::vec2d& _pos2) const -> Utils::Math::vec2d
  {
    Utils::Math::vec2d dist = _pos1 - _pos2;
    while (std::abs(dist[1]) > m_halfLy) {
      if (dist[1] > 0.0) {
        dist[1] -= m_Ly;
        dist[0] -= m_xy;
      } else {
        dist[1] += m_Ly;
        dist[0] += m_xy;
      }
    }
    while (std::abs(dist[0]) > m_halfLx) {
      if (dist[0] > 0.0)
        dist[0] -= m_Lx;
      else
        dist[0] += m_Lx;
    }
    return dist;
  }

private:
    void updateDependent(){
        m_invLy = 1.0 / m_Ly;
        m_invLx = 1.0 / m_Lx;
        m_invxy = -m_xy / (m_Lx * m_Ly);
        m_halfLy = 0.5 * m_Ly;
        m_halfLx = 0.5 * m_Lx;
        m_area = m_Lx * m_Ly;
        m_invArea = 1.0 / m_area;
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

  double m_area;
  double m_invArea;
};

}  // namespace networkV4