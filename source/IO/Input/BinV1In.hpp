void networkV4::network::loadFromBinV1(std::ifstream& _file, double _lambda)
{
  std::size_t N;
  std::size_t B;
  _file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  _file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  double domainX;
  _file.read(reinterpret_cast<char*>(&domainX), sizeof(double));
  setDomain(vec2d(domainX, domainX));
  _file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));
  m_restSize = m_domain;

  m_nodes.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    vec2d pos;
    _file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
    m_nodes.addNode(pos);
  }

  m_bonds.resize(B);
  for (auto& b : m_bonds) {
    bool connected;
    std::size_t index1;
    std::size_t index2;
    double naturalLength;
    double constant;

    _file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    _file.read(reinterpret_cast<char*>(&constant), sizeof(double));

    b = bond(index1,
             index2,
             connected,
             bondType::single,
             naturalLength,
             constant,
             _lambda);
  }
  initStresses();
}