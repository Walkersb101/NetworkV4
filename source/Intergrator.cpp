#include <algorithm>

#include "Intergrator.hpp"

#include "Config.hpp"

network::intergrator::intergrator() {}

network::intergrator::~intergrator() {}

void network::intergrator::intergrate(network& _network) {}

network::overdampedEuler::overdampedEuler()
    : m_dt(defaultConfig.default_dt)
{
}

network::overdampedEuler::overdampedEuler(double _dt)
    : m_dt(_dt)
{
}

network::overdampedEuler::~overdampedEuler() {}

void network::overdampedEuler::intergrate(network& _network)
{
  nodes& networkNodes = _network.getNodes();
  networkNodes.clearForces();
  _network.computeForces();
  std::transform(PAR networkNodes.positions_begin(),
                 networkNodes.positions_end(),
                 networkNodes.forces_begin(),
                 networkNodes.positions_begin(),
                 [this](const vec2d& _pos, const vec2d& _force)
                 { return _pos + _force * m_dt; });
}

network::OverdampedEulerHeun::OverdampedEulerHeun()
    : m_dt(defaultConfig.default_dt)
{
}

network::OverdampedEulerHeun::OverdampedEulerHeun(double _dt)
    : m_dt(_dt)
{
}

network::OverdampedEulerHeun::~OverdampedEulerHeun() {}

void network::OverdampedEulerHeun::intergrate(network& _network)
{
  nodes& networkNodes = _network.getNodes();
  networkNodes.clearForces();
  _network.clearStresses();
  _network.computeForces();
  std::copy(networkNodes.forces_begin(),
            networkNodes.forces_end(),
            m_tempForces.begin());
  std::transform(PAR networkNodes.positions_begin(),
                 networkNodes.positions_end(),
                 networkNodes.forces_begin(),
                 networkNodes.positions_begin(),
                 [this](const vec2d& _pos, const vec2d& _force)
                 { return _pos + _force * m_dt; });
  _network.clearStresses();
  m_tempStresses = _network.getStresses();
  networkNodes.clearForces();
  _network.computeForces();

  double halfdt = m_dt * 0.5;

    // Averages the forces and positions
  std::transform(PAR networkNodes.forces_begin(),
                 networkNodes.forces_end(),
                 m_tempForces.begin(),
                 networkNodes.forces_begin(),
                 [this](const vec2d& _force, const vec2d& _tempForce)
                 { return (_force + _tempForce) * 0.5; });
    //Steps based on the average forces x1 = x0 + [f(x0) + f(x')] * dt/2
    // x' = x0 + f(x0) * dt -> x1 = x' + [f(x') - f(x0)] * dt/2
  std::transform(PAR networkNodes.positions_begin(),
                 networkNodes.positions_end(),
                 networkNodes.forces_begin(),
                 networkNodes.positions_begin(),
                 [halfdt](const vec2d& _pos, const vec2d& _force)
                 { return _pos + _force * halfdt; });
  std::transform(PAR networkNodes.positions_begin(),
                 networkNodes.positions_end(),
                 m_tempForces.begin(),
                 networkNodes.positions_begin(),
                 [halfdt](const vec2d& _pos, const vec2d& _force)
                 { return _pos - _force * halfdt; });
}
