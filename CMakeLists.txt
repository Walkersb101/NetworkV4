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
    source/Vec2.cpp
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
    toml11
    GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
)

FetchContent_MakeAvailable(toml11)

# ---- Declare executable ----

add_executable(NetworkV4_exe source/main.cpp)
add_executable(NetworkV4::exe ALIAS NetworkV4_exe)

set_property(TARGET NetworkV4_exe PROPERTY OUTPUT_NAME NetworkV4)

target_compile_features(NetworkV4_exe PRIVATE cxx_std_20)

target_link_libraries(NetworkV4_exe PRIVATE NetworkV4_lib)

target_link_libraries(NetworkV4_exe PRIVATE toml11::toml11)

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
