#include <fstream>

#include "Core/Simulation.hpp"

auto main(int argc, char* argv[]) -> int
{
  std::filesystem::path path =
      argc == 2 ? argv[1] : "/home/sam/Documents/Code/NetworkV4/test/test.toml";

  networkV4::Simulation input(path);
  
  input.run();

  return 0;
}
