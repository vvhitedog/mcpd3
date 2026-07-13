# mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
# Copyright (C) 2021 Matt Gara
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.
cmake_minimum_required(VERSION 2.8)

project(mcpd3)

set(PROJECT_NAME "mcpd3 - minimum cut primal dual dual decomposition")
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH} )
set(PROJECT_AUTHOR Matt Gara )
set(PROJECT_AUTHOR_EMAIL gara.matt@gmail.com )

# Add custom cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CTest)

# Set standard to c++11
if (${CMAKE_VERSION} VERSION_GREATER 3.0)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
else ()
  ADD_DEFINITIONS( -std=c++17 )
endif ()

# Required packages
find_package (Threads REQUIRED)

# Options
option(GPERF_PROFILER_BUILD "Build with gperftools CPU profiling enabled." OFF )
option(TRAP_SIGNED_INTEGER_OVERFLOW "Trap signed and unsigned integer overflow using the `ftrapv` flag." OFF )
option(METIS_ENABLED "Use the metis library for graph partitioning (highly recommended.)" OFF )
set(MCPD_CAPACITY_MODE "32" CACHE STRING
  "Capacity precision: 32, 64, 128, or gmp")
set_property(CACHE MCPD_CAPACITY_MODE PROPERTY STRINGS 32 64 128 gmp)

add_library(mcpd3_capacity_config INTERFACE)
target_include_directories(mcpd3_capacity_config INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR})
find_path(MCPD_BOOST_INCLUDE_DIR NAMES boost/multiprecision/cpp_int.hpp
  HINTS ENV BOOST_ROOT PATH_SUFFIXES include)
if (NOT MCPD_BOOST_INCLUDE_DIR)
  message(FATAL_ERROR
    "mcpd3 capacity support requires Boost.Multiprecision headers")
endif ()
target_include_directories(mcpd3_capacity_config INTERFACE
  ${MCPD_BOOST_INCLUDE_DIR})
if (MCPD_CAPACITY_MODE STREQUAL "32")
  target_compile_definitions(mcpd3_capacity_config INTERFACE
    MCPD_CAPACITY_MODE_32=1)
elseif (MCPD_CAPACITY_MODE STREQUAL "64")
  target_compile_definitions(mcpd3_capacity_config INTERFACE
    MCPD_CAPACITY_MODE_64=1)
elseif (MCPD_CAPACITY_MODE STREQUAL "128")
  target_compile_definitions(mcpd3_capacity_config INTERFACE
    MCPD_CAPACITY_MODE_128=1)
elseif (MCPD_CAPACITY_MODE STREQUAL "gmp")
  find_path(MCPD_GMPXX_INCLUDE_DIR NAMES gmpxx.h
    HINTS ENV GMP_ROOT PATH_SUFFIXES include)
  find_path(MCPD_GMP_INCLUDE_DIR NAMES gmp.h
    HINTS ENV GMP_ROOT PATH_SUFFIXES include include/x86_64-linux-gnu)
  find_library(MCPD_GMPXX_LIBRARY NAMES gmpxx
    HINTS ENV GMP_ROOT PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)
  find_library(MCPD_GMP_LIBRARY NAMES gmp
    HINTS ENV GMP_ROOT PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)
  if (NOT MCPD_GMPXX_INCLUDE_DIR OR NOT MCPD_GMP_INCLUDE_DIR OR
      NOT MCPD_GMPXX_LIBRARY OR NOT MCPD_GMP_LIBRARY)
    message(FATAL_ERROR
      "MCPD_CAPACITY_MODE=gmp requires the GMP C and C++ development files")
  endif ()
  target_compile_definitions(mcpd3_capacity_config INTERFACE
    MCPD_CAPACITY_MODE_GMP=1)
  target_include_directories(mcpd3_capacity_config INTERFACE
    ${MCPD_GMPXX_INCLUDE_DIR} ${MCPD_GMP_INCLUDE_DIR})
  target_link_libraries(mcpd3_capacity_config INTERFACE
    ${MCPD_GMPXX_LIBRARY} ${MCPD_GMP_LIBRARY})
else ()
  message(FATAL_ERROR
    "MCPD_CAPACITY_MODE must be one of: 32, 64, 128, gmp")
endif ()

# Set flags based on options
if (GPERF_PROFILER_BUILD) 
  set(PROFILER_LINK_FLAG profiler)
else ()
  set(PROFILER_LINK_FLAG)
endif (GPERF_PROFILER_BUILD)
if (TRAP_SIGNED_INTEGER_OVERFLOW)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -ftrapv")
endif (TRAP_SIGNED_INTEGER_OVERFLOW)
if (METIS_ENABLED)
  find_package(METIS REQUIRED)
  include_directories(${METIS_INCLUDE_DIRS})
  find_library(GKLIB_LIBRARY
    NAMES GKlib
    HINTS ${METIS_LIB_DIR} ${METIS_DIR}/lib ${METIS_LIBRARY_DIRS})
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAVE_METIS")
endif (METIS_ENABLED)

include_directories(.)
set(MAXFLOW_SRC maxflow/graph.cpp  maxflow/maxflow.cpp)

## Simple example
add_executable(simple_example example/simple.cpp ${MAXFLOW_SRC})
target_link_libraries(simple_example PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT} )

## DIMACS maxflow format example
add_executable(dimacs_example example/dimacs.cpp ${MAXFLOW_SRC})
target_link_libraries(dimacs_example PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT} )

## DIMACS maxflow format example for dual decomposition
add_executable(dimacs_dual_decomp_example example/dimacs_dual_decomp.cpp ${MAXFLOW_SRC})
target_link_libraries(dimacs_dual_decomp_example PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${METIS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
if (METIS_ENABLED AND GKLIB_LIBRARY)
  target_link_libraries(dimacs_dual_decomp_example PUBLIC ${GKLIB_LIBRARY} m)
endif ()

## Native monolith benchmark. This intentionally uses DualDecomposition
## directly and does not go through the worker/distributed package path.
add_executable(mcpd3_native_monolith_benchmark benchmark/native_monolith.cpp ${MAXFLOW_SRC})
target_link_libraries(mcpd3_native_monolith_benchmark PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${METIS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
if (METIS_ENABLED AND GKLIB_LIBRARY)
  target_link_libraries(mcpd3_native_monolith_benchmark PUBLIC ${GKLIB_LIBRARY} m)
endif ()

if (BUILD_TESTING)
  add_executable(capacity_precision_test tests/capacity_precision_test.cpp ${MAXFLOW_SRC})
  target_link_libraries(capacity_precision_test PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})
  add_test(NAME capacity_precision_test COMMAND capacity_precision_test)

  add_executable(partition_worker_test tests/partition_worker_test.cpp ${MAXFLOW_SRC})
  target_link_libraries(partition_worker_test PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})
  add_test(NAME partition_worker_test COMMAND partition_worker_test)
endif ()
