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
option(MCPD_LEGACY_32BIT_DD_REPLAY
  "Replay historical narrow 32-bit DD arithmetic with defined wrapping. Benchmark compatibility only."
  OFF)
include(cmake/McpdCapacity.cmake)
mcpd_configure_capacity(mcpd3_capacity_config ${CMAKE_CURRENT_SOURCE_DIR})

if (MCPD_LEGACY_32BIT_DD_REPLAY)
  if (NOT MCPD_CAPACITY_MODE STREQUAL "32")
    message(FATAL_ERROR
      "MCPD_LEGACY_32BIT_DD_REPLAY requires MCPD_CAPACITY_MODE=32")
  endif ()
  message(WARNING
    "MCPD_LEGACY_32BIT_DD_REPLAY deterministically reproduces historical "
    "wrapped DD arithmetic. It is for benchmark replay only and does not "
    "provide production-valid solver arithmetic.")
  target_compile_definitions(mcpd3_capacity_config INTERFACE
    MCPD_LEGACY_32BIT_DD_REPLAY=1)
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

add_executable(mcpd3_pdhg_benchmark benchmark/pdhg_undirected.cpp
  primaldual/pdhg_undirected.cpp ${MAXFLOW_SRC})
target_link_libraries(mcpd3_pdhg_benchmark PUBLIC mcpd3_capacity_config
  ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})

if (BUILD_TESTING)
  add_executable(capacity_precision_test tests/capacity_precision_test.cpp ${MAXFLOW_SRC})
  target_link_libraries(capacity_precision_test PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})
  add_test(NAME capacity_precision_test COMMAND capacity_precision_test)

  add_executable(partition_worker_test tests/partition_worker_test.cpp ${MAXFLOW_SRC})
  target_link_libraries(partition_worker_test PUBLIC mcpd3_capacity_config ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})
  add_test(NAME partition_worker_test COMMAND partition_worker_test)

  add_executable(pdhg_undirected_test tests/pdhg_undirected_test.cpp
    primaldual/pdhg_undirected.cpp ${MAXFLOW_SRC})
  target_link_libraries(pdhg_undirected_test PUBLIC mcpd3_capacity_config
    ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT})
  add_test(NAME pdhg_undirected_test COMMAND pdhg_undirected_test)
  add_test(NAME pdhg_benchmark_existing_smoke COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --solver existing)
  add_test(NAME pdhg_benchmark_pdhg_smoke COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --solver pdhg --pdhg-check-interval 1 --pdhg-max-iterations 10000)
  add_test(NAME pdhg_benchmark_both_smoke COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --solver both --pdhg-check-interval 1 --pdhg-max-iterations 10000)
  add_test(NAME pdhg_benchmark_options_smoke COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --solver pdhg
    --pdhg-max-iterations 10 --pdhg-check-interval 1
    --pdhg-tau 0.1 --pdhg-sigma 0.1 --pdhg-theta 0.5
    --pdhg-step-size-scale 0.9 --pdhg-step-balance 2
    --capacity-quantum 0 --pdhg-eps-abs 1000 --pdhg-eps-rel 0
    --pdhg-time-limit-seconds 1
    --pdhg-stagnation-checks 2 --pdhg-stagnation-tolerance 0
    --pdhg-use-ergodic-primal true --pdhg-use-ergodic-dual true
    --pdhg-lower-bound-safety-factor 128
    --print-history true --verbose true)
  add_test(NAME pdhg_benchmark_invalid_solver COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --solver invalid)
  set_tests_properties(pdhg_benchmark_invalid_solver PROPERTIES WILL_FAIL TRUE)
  add_test(NAME pdhg_benchmark_negative_size COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --pdhg-max-iterations -1)
  set_tests_properties(pdhg_benchmark_negative_size PROPERTIES WILL_FAIL TRUE)
  add_test(NAME pdhg_benchmark_invalid_boolean COMMAND mcpd3_pdhg_benchmark
    --dimacs ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/pdhg_path.max
    --print-history maybe)
  set_tests_properties(pdhg_benchmark_invalid_boolean PROPERTIES WILL_FAIL TRUE)
endif ()
