#include <fstream>

#include "DataOut.hpp"
#include "Intergrator.hpp"
#include "Network.hpp"
#include "Protocol.hpp"
#include "TomlLoad.hpp"

auto main(int argc, char* argv[]) -> int
{
  std::filesystem::path path =
      argc == 2 ? argv[1] : "/home/sam/Documents/Code/NetworkV4/test/test.toml";

  tomlIn input(path);
  auto net = input.readNetwork();
  auto problem = input.readProblem();
  auto dataOut = input.readDataOut();
  auto netOut = input.readNetworkOut();

  problem->initIO(net, dataOut, netOut);
  problem->run(net);

  return 0;
}
