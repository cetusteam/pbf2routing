#   This file is part of the osmpbf library.
#
#   Copyright(c) 2014 Oliver Groß.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 3 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, see
#   <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1)

project(pbf2routing)

set(DEBUG_FLAGS " ${DEBUG_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${DEBUG_FLAGS}")

find_package(Protobuf REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

include_directories(
	${PROTOBUF_INCLUDE_DIRS}
	${ZLIB_INCLUDE_DIRS}
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_BINARY_DIR}
)

set(OSMPBF_LIBRARIES osmpbf protos-ch)

add_executable(${PROJECT_NAME} pbf2routing.cpp)
add_dependencies(${PROJECT_NAME} osmpbf)
target_link_libraries(${PROJECT_NAME} ${OSMPBF_LIBRARIES} ${ZLIB_LIBRARIES} ${PROTOBUF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
