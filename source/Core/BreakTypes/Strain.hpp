#pragma once

#include <optional>

#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace BreakTypes
{
class StrainBreak
{
public:
    StrainBreak(double _lambda);

public:
    bool checkBreak(const Utils::vec2d& _r) const;
    std::optional<double> thresholdData(const Utils::vec2d& _r) const;

private:
    const double m_lambda;
};

StrainBreak::StrainBreak(double _lambda)
    : m_lambda(_lambda)
{}

inline bool StrainBreak::checkBreak(const Utils::vec2d& _r) const
{
    return _r.norm() > m_lambda;
}

inline std::optional<double> StrainBreak::thresholdData(const Utils::vec2d& _r) const
{
    return _r.norm() - m_lambda;
}

}  // namespace BreakTypes
}  // namespace networkV4