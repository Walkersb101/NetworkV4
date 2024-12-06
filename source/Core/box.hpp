#pragma once

#include "Misc/Vec2.hpp"

namespace networkV4
{

class box
{
  box(const Utils::vec2d& _a, const Utils::vec2d& _b)
      : m_a(_a)
      , m_b(_b)
  {
    calculateInverseDeterminant();
  }

public:
  auto getA() const -> const Utils::vec2d& { return m_a; }
  auto getB() const -> const Utils::vec2d& { return m_b; }

  void setA(const Utils::vec2d& _a)
  {
    m_a = _a;
    calculateInverseDeterminant();
  }
  void setB(const Utils::vec2d& _b)
  {
    m_b = _b;
    calculateInverseDeterminant();
  }

public:
  auto toFractional(const Utils::vec2d& _pos) const -> Utils::vec2d
  {
    return Utils::vec2d(m_invDetertminant * (_pos.x * m_b.y - _pos.y * m_b.x),
                        m_invDetertminant * (m_a.x * _pos.y - m_a.y * _pos.x));
  };

  auto toCartesian(const Utils::vec2d& _pos) const -> Utils::vec2d
  {
    return Utils::vec2d(m_a.x * _pos.x + m_b.x * _pos.y,
                        m_a.y * _pos.x + m_b.y * _pos.y);
  };

private:
void calculateInverseDeterminant()
{
  const double determinant = m_a.x * m_b.y - m_a.y * m_b.x;
  if (determinant == 0.0) {
    throw("box::calculateInverseDeterminant: determinant is zero");
  }
  m_invDetertminant = 1.0 / determinant;
}

private:
Utils::vec2d m_a;
Utils::vec2d m_b;

double m_invDetertminant;
}

}  // namespace networkV4