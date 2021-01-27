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

# Dependencies
find_package(METIS REQUIRED)

include_directories(.)

set(MAXFLOW_SRC maxflow/graph.cpp  maxflow/maxflow.cpp)

## Simple example
add_executable(simple_example example/simple.cpp ${MAXFLOW_SRC})

## DIMACS maxflow format example
add_executable(dimacs_example example/dimacs.cpp ${MAXFLOW_SRC})

#set(DDPD_SRC src)
#set(SHARED ${DDPD_SRC}/util/timer.cpp ${DDPD_SRC}/maxflow/graph.cpp  ${DDPD_SRC}/maxflow/maxflow.cpp)
#
## Compile the dual decomp init binary
#set(DUAL_DECOMP_INIT ${DDPD_SRC}/dual_decomp_init.cpp )
#add_executable(dual_decomp_init ${DUAL_DECOMP_INIT} ${SHARED})
#target_include_directories( dual_decomp_init PRIVATE ${DDPD_SRC} ${METIS_INCLUDE_DIR} )
#target_link_libraries(dual_decomp_init ${METIS_LIBRARIES} )
#
## Compile the dual decomp init binary
#set(DUAL_DECOMP_SOLVE ${DDPD_SRC}/dual_decomp_solve.cpp )
#add_executable(dual_decomp_solve ${DUAL_DECOMP_SOLVE} ${SHARED})
#target_include_directories( dual_decomp_solve PRIVATE ${DDPD_SRC} )
#target_link_libraries(dual_decomp_solve pthread profiler)
#
#set(DUAL_DECOMP_flow_test ${DDPD_SRC}/dual_decomp_flow_test.cpp )
#add_executable(dual_decomp_flow_test ${DUAL_DECOMP_flow_test} ${SHARED})
#target_include_directories( dual_decomp_flow_test PRIVATE ${DDPD_SRC} )
#target_link_libraries(dual_decomp_flow_test pthread profiler)
