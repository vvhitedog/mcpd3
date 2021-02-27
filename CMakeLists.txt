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

# Set standard to c++11
if (${CMAKE_VERSION} VERSION_GREATER 3.0)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
else ()
  ADD_DEFINITIONS( -std=c++17 )
endif ()

# Required packages
find_package (Threads REQUIRED)
find_package(Boost 1.69.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

# Options
option(GPERF_PROFILER_BUILD "Build with gperftools CPU profiling enabled." OFF )
option(TRAP_SIGNED_INTEGER_OVERFLOW "Trap signed and unsigned integer overflow using the `ftrapv` flag." OFF )
option(METIS_ENABLED "Use the metis library for graph partitioning (highly recommended.)" OFF )

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
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAVE_METIS")
endif (METIS_ENABLED)

include_directories(.)
set(MAXFLOW_SRC maxflow/graph.cpp  maxflow/maxflow.cpp)

## Simple example
add_executable(simple_example example/simple.cpp ${MAXFLOW_SRC})
target_link_libraries(simple_example PUBLIC ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT} )

## DIMACS maxflow format example
add_executable(dimacs_example example/dimacs.cpp ${MAXFLOW_SRC})
target_link_libraries(dimacs_example PUBLIC ${PROFILER_LINK_FLAG} ${CMAKE_THREAD_LIBS_INIT} )

## DIMACS maxflow format example for dual decomposition
add_executable(dimacs_dual_decomp_example example/dimacs_dual_decomp.cpp ${MAXFLOW_SRC})
target_link_libraries(dimacs_dual_decomp_example PUBLIC ${PROFILER_LINK_FLAG} ${METIS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
