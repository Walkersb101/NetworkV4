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
    source/Bonds.cpp
    source/Config.cpp
    source/Intergrator.cpp
    source/DataOut.cpp
    source/Network.cpp
    source/NetworkOut.cpp
    source/Nodes.cpp
    source/Protocol.cpp
    source/Roots.cpp
    source/TomlLoad.cpp
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
