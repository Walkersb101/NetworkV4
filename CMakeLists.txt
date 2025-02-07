cmake_minimum_required(VERSION 3.14)

project(
    networkV4
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

# ---- Options ----
option(USE_OPENMP "Enable OpenMP support" OFF)

# ---- Declare library ----

add_library(
    NetworkV4_lib OBJECT
    source/Core/Bonds.cpp
    source/Core/Network.cpp
    source/Core/Nodes.cpp
    source/Core/Simulation.cpp

    source/Core/Threading/NetworkOMP.cpp

    source/IO/BaseIO.cpp
    source/IO/TimeSeries/DataOut.cpp
    source/IO/NetworkDump/NetworkOut.cpp
    source/IO/Input/NetworkIn.cpp
    
    source/Protocols/DoubleNetworks/Quasistatic.cpp
    source/Protocols/DoubleNetworks/Propogator.cpp

    #source/Protocols/StepStrain.cpp
    #source/Protocols/Propogator.cpp
    
    #source/Misc/Roots.cpp
)

target_include_directories(
    NetworkV4_lib
    PUBLIC
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(NetworkV4_lib PUBLIC cxx_std_20)
target_compile_options(NetworkV4_lib PRIVATE -std=c++20)

# ---- FetchContent ----

include(FetchContent)

FetchContent_Declare(
    toml
    GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
    GIT_TAG v4.0.2
)

FetchContent_MakeAvailable(toml)

target_link_libraries(NetworkV4_lib PUBLIC toml11::toml11)

#FetchContent_Declare(
#  trng4
#  GIT_REPOSITORY https://github.com/rabauke/trng4.git
#  GIT_TAG 1ef58a7
#  )
#FetchContent_MakeAvailable(trng4)
#
#target_link_libraries(NetworkV4_lib PUBLIC trng4)

set(ZSTD_BUILD_STATIC ON)
set(ZSTD_BUILD_SHARED OFF)

FetchContent_Declare(
    zstd
    URL "https://github.com/facebook/zstd/releases/download/v1.5.5/zstd-1.5.5.tar.gz"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR build/cmake
    #GIT_REPOSITORY    https://github.com/facebook/zstd.git
    #GIT_TAG           v1.5.6
)

FetchContent_MakeAvailable(zstd)
target_link_libraries(NetworkV4_lib PUBLIC libzstd_static)

target_include_directories(
    NetworkV4_lib
    PUBLIC
    ${zstd_SOURCE_DIR}/lib
)

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

#FetchContent_Declare(
#    spdlog
#    GIT_REPOSITORY https://github.com/gabime/spdlog.git
#    GIT_TAG v1.15.0
#)
#FetchContent_MakeAvailable(spdlog)
#target_link_libraries(NetworkV4_lib PUBLIC spdlog::spdlog_header_only)

set(BUILD_TESTING OFF)
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

FetchContent_Declare(
    expected
    GIT_REPOSITORY https://github.com/TartanLlama/expected.git
    GIT_TAG v1.1.0
)
FetchContent_MakeAvailable(expected)
target_link_libraries(NetworkV4_lib PUBLIC expected)

# ---- OpenMP ----

if(USE_OPENMP)
    message(STATUS "OpenMP requested")
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(NetworkV4_lib PUBLIC OpenMP::OpenMP_CXX)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
        message(WARNING "OpenMP requested but not found")
    endif()
endif()

# ---- Declare executable ----

add_executable(NetworkV4_exe source/main.cpp)
add_executable(networkV4::exe ALIAS NetworkV4_exe)

if (USE_OPENMP)
    set_property(TARGET NetworkV4_exe PROPERTY OUTPUT_NAME networkV4_omp)
else()
    set_property(TARGET NetworkV4_exe PROPERTY OUTPUT_NAME networkV4)
endif()

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
