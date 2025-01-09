#pragma once 

namespace networkV4
{
namespace protocols
{

class quasiStaticStrainDouble : public protocolBase
{
public:
  quasiStaticStrainDouble() = delete;
  quasiStaticStrainDouble(
      std::shared_ptr<deform::deformBase>& _deform,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
      std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
      const network& _network)

    