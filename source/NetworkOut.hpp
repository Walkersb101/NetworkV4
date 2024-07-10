#pragma once

#include <filesystem>

#include "Network.hpp"

class networkOut
{
public:
  networkOut();
  networkOut(const std::filesystem::path& _path);
  virtual ~networkOut();

public:
  virtual void save(const networkV4::network& _net, const std::string& _name) const = 0;

protected:
  std::filesystem::path m_path;
};

class noNetworkOut : public networkOut
{
public:
  noNetworkOut();
  ~noNetworkOut();

public:
  void save(const networkV4::network& _net, const std::string& _name) const override;
};

class networkOutBinV2 : public networkOut
{
public:
  networkOutBinV2(const std::filesystem::path& _path);
  ~networkOutBinV2();

public:
  void save(const networkV4::network& _net, const std::string& _name) const override;
};