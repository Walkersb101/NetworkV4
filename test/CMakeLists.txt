# Parent project does not export its library target, so this CML implicitly
# depends on being added from it, i.e. the testing is done only from the build
# tree and is not feasible from an install location

project(NetworkV4Tests LANGUAGES CXX)

# ---- Catch2 ----

include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.6.0
)

FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

# ---- Tests ----

add_executable(NetworkV4_tests source/BondsTest.cpp source/Vec2Test.cpp)
target_link_libraries(NetworkV4_tests PRIVATE NetworkV4_lib)
target_compile_features(NetworkV4_tests PRIVATE cxx_std_20)

target_link_libraries(NetworkV4_tests PRIVATE Catch2::Catch2WithMain)

#add_test(NAME NetworkV4_tests COMMAND NetworkV4_tests)

catch_discover_tests(NetworkV4_tests)

# ---- End-of-file commands ----

add_folders(Test)
