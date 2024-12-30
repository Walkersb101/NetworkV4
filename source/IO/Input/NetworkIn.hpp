#pragma once

#include <filesystem>

#include "Core/Network.hpp"

namespace IO
{
namespace NetworkIn
{

class networkIn
{
public:
  networkIn();
  virtual ~networkIn() = 0;

public:
  virtual networkV4::network load() = 0;
};

}  // namespace NetworkIn
}  // namespace IO
