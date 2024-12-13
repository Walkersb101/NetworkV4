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

void networkV4::network::loadFromBinV2(std::ifstream& _file)
{
  std::size_t N;
  std::size_t B;
  _file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  _file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  vec2d domain;  
  _file.read(reinterpret_cast<char*>(&domain), sizeof(vec2d));
  _file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));
  setDomain(domain);
  m_restSize = m_domain;

  m_nodes.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    vec2d pos;
    _file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
    m_nodes.addNode(pos);
  }

  m_bonds.reserve(B);
  for (size_t i = 0; i < B; ++i) {
    std::size_t index1;
    std::size_t index2;
    bool connected;
    bool matrix;
    double naturalLength;
    double constant;
    double lambda;

    _file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&matrix), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    _file.read(reinterpret_cast<char*>(&constant), sizeof(double));
    _file.read(reinterpret_cast<char*>(&lambda), sizeof(double));

    if (index1 >= N || index2 >= N) {
      throw std::runtime_error("Invalid bond indices: " + std::to_string(index1)
                               + ", " + std::to_string(index2));
    }

    const bondType type = matrix ? bondType::matrix : bondType::sacrificial;
    m_bonds.add(
        bond(index1, index2, connected, type, naturalLength, constant, lambda));
  }
  initStresses();
}