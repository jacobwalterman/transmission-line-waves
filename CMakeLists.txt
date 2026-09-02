cmake_minimum_required(VERSION 3.16)
project(coupled_piezo_heat)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

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

find_package(Eigen3 REQUIRED)

add_executable(coupled_piezo_heat src/main.cpp)
target_link_libraries(coupled_piezo_heat PRIVATE Eigen3::Eigen)
nanobind_add_module(_coupled_piezo_heat src/bindings.cpp src/solver.cpp src/transmission_line_solver.cpp)
target_link_libraries(_coupled_piezo_heat PRIVATE Eigen3::Eigen)
