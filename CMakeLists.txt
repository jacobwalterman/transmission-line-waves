cmake_minimum_required(VERSION 3.15...4.4)
project(transmission_line_waves)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_VARIABLE NB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

if (CMAKE_VERSION VERSION_LESS 3.18)
    set(DEV_MODULE Development)
else()
    set(DEV_MODULE Development.Module)
endif()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.4 QUIET)
if(NOT Eigen3_FOUND)
  include(FetchContent)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(eigen)
endif()

add_executable(transmission_line_waves src/main.cpp)
target_link_libraries(transmission_line_waves PRIVATE Eigen3::Eigen)
nanobind_add_module(_transmission_line_waves src/bindings.cpp src/solver.cpp src/transmission_line_solver.cpp)
target_link_libraries(_transmission_line_waves PRIVATE Eigen3::Eigen)

install(TARGETS _transmission_line_waves
        LIBRARY DESTINATION .
        RUNTIME DESTINATION .)
