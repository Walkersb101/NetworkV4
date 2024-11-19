cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    NetworkV4
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    NetworkV4_lib OBJECT
    source/Core/Bonds.cpp
    source/Core/BreakTypes.cpp
    source/Core/Network.cpp
    source/Core/Simulation.cpp
    source/Core/Nodes.cpp

    source/Integration/Integrator.cpp
    source/Integration/Minimization.cpp
    source/Integration/OverDamped.cpp

    source/IO/DataOut.cpp
    source/IO/hdf5Out.cpp
    source/IO/NetworkOut.cpp
    
    source/Protocols/Protocol.cpp
    source/Protocols/StepStrain.cpp
    source/Protocols/Quasistatic.cpp
    source/Protocols/Propogator.cpp
    
    source/Misc/Config.cpp
    source/Misc/EnumString.cpp
    source/Misc/Roots.cpp
)

target_include_directories(
    NetworkV4_lib ${warning_guard}
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(NetworkV4_lib PUBLIC cxx_std_20)

# ---- FetchContent ----

include(FetchContent)

FetchContent_Declare(
    toml
    GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
    GIT_TAG v4.0.2
)

FetchContent_MakeAvailable(toml)

target_link_libraries(NetworkV4_lib PUBLIC toml11::toml11)

FetchContent_Declare(
    bxzstr
    GIT_REPOSITORY https://github.com/tmaklin/bxzstr.git
    GIT_TAG v1.2.2
)

FetchContent_MakeAvailable(bxzstr)

target_link_libraries(NetworkV4_lib PUBLIC bxzstr)

FetchContent_Declare(
  trng4
  GIT_REPOSITORY https://github.com/rabauke/trng4.git
  GIT_TAG 1ef58a7
  )
FetchContent_MakeAvailable(trng4)

target_link_libraries(NetworkV4_lib PUBLIC trng4)

FetchContent_Declare(
  HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
  GIT_SHALLOW ON
  GIT_TAG v3.0.0-beta1
  )
FetchContent_GetProperties(highfive)
FetchContent_MakeAvailable(highfive)

include_directories(BEFORE SYSTEM ${highfive_SOURCE_DIR} ${highfive_BINARY_DIR}/include)
target_link_libraries(NetworkV4_lib PUBLIC HighFive)

# ---- Declare executable ----

add_executable(NetworkV4_exe source/main.cpp)
add_executable(NetworkV4::exe ALIAS NetworkV4_exe)

set_property(TARGET NetworkV4_exe PROPERTY OUTPUT_NAME NetworkV4)

target_compile_features(NetworkV4_exe PRIVATE cxx_std_20)

target_link_libraries(NetworkV4_exe PRIVATE NetworkV4_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT NetworkV4_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of NetworkV4"
  )
endif()

include(cmake/dev-mode.cmake)
