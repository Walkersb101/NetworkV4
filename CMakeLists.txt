cmake_minimum_required(VERSION 3.14)

project(
    networkV4
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)
# ---- Declare library ----

add_library(
    NetworkV4_lib OBJECT
    source/Core/Bonds.cpp
    source/Core/Network.cpp
    source/Core/Nodes.cpp
    source/Core/Simulation.cpp

    source/Integration/Integrator.cpp
    #source/Integration/Minimization.cpp
    #source/Integration/OverDamped.cpp

    source/IO/BaseIO.cpp
    source/IO/TimeSeries/DataOut.cpp
    #source/IO/hdf5Out.cpp
    source/IO/NetworkDump/NetworkOut.cpp
    source/IO/Input/NetworkIn.cpp
    
    source/Protocols/Protocol.cpp
    #source/Protocols/StepStrain.cpp
    #source/Protocols/Quasistatic.cpp
    #source/Protocols/Propogator.cpp
    
    #source/Misc/Roots.cpp
)

target_include_directories(
    NetworkV4_lib
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
  trng4
  GIT_REPOSITORY https://github.com/rabauke/trng4.git
  GIT_TAG 1ef58a7
  )
FetchContent_MakeAvailable(trng4)

target_link_libraries(NetworkV4_lib PUBLIC trng4)

FetchContent_Declare(
    zstd
    URL "https://github.com/facebook/zstd/releases/download/v1.5.5/zstd-1.5.5.tar.gz"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR build/cmake
)

FetchContent_MakeAvailable(zstd)

target_link_libraries(NetworkV4_lib PRIVATE libzstd_static)

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

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.15.0
)
FetchContent_MakeAvailable(spdlog)
target_link_libraries(NetworkV4_lib PUBLIC spdlog::spdlog_header_only)

FetchContent_Declare(
    libmorton
    GIT_REPOSITORY https://github.com/Forceflow/libmorton.git
    GIT_TAG v0.2.12
)
FetchContent_MakeAvailable(libmorton)
target_link_libraries(NetworkV4_lib PUBLIC libmorton)

FetchContent_Declare(
    range-v3
    GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
    GIT_TAG 0.12.0
)
FetchContent_MakeAvailable(range-v3)
target_link_libraries(NetworkV4_lib PUBLIC range-v3)

# ---- Declare executable ----

add_executable(NetworkV4_exe source/main.cpp)
add_executable(networkV4::exe ALIAS NetworkV4_exe)

set_property(TARGET NetworkV4_exe PROPERTY OUTPUT_NAME networkV4)

target_compile_features(NetworkV4_exe PRIVATE cxx_std_20)

target_link_libraries(NetworkV4_exe PRIVATE NetworkV4_lib)

# ---- Developer mode ----

if(NOT NetworkV4_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of NetworkV4"
  )
endif()

#include(cmake/dev-mode.cmake)
