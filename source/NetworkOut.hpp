#pragma once

#include <filesystem>

#include <bxzstr.hpp>

#include "Network.hpp"

template<typename T>
bool write(const std::stringstream& _ss,
           const std::filesystem::path& _path,
           T _file)
{
  if (!_file.is_open()) {
    return false;
  }
  _file << _ss.str() << std::flush;
  _file.close();
  return true;
}

template<>
bool write<bxz::Compression>(const std::stringstream& _ss,
                             const std::filesystem::path& _path,
                             bxz::Compression _compression);

class networkOut
{
public:
  networkOut();
  networkOut(const std::filesystem::path& _path);
  virtual ~networkOut();

public:
  virtual void save(const networkV4::network& _net,
                    const std::size_t _step,
                    const double _time,
                    const std::string& _type) const;

protected:
  std::filesystem::path m_path;
};

class noNetworkOut : public networkOut
{
public:
  noNetworkOut();
  ~noNetworkOut();

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override;
};

class networkOutBinV2 : public networkOut
{
public:
  networkOutBinV2(const std::filesystem::path& _path,
                  bxz::Compression _compression);
  ~networkOutBinV2();

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override;

private:
  bxz::Compression m_compression;
};

class networkOutHDF5 : public networkOut
{
public:
  networkOutHDF5(const std::filesystem::path& _path);
  ~networkOutHDF5();

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override;

private:
    HighFive::File m_file;
};