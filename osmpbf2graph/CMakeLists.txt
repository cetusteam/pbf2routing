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

project(osmpbf2route)

MACRO(ADD_TARGET_PROPERTIES _target _name)
	SET(_properties)
	FOREACH(_prop ${ARGN})
		SET(_properties "${_properties} ${_prop}")
	ENDFOREACH(_prop)
	GET_TARGET_PROPERTY(_old_properties ${_target} ${_name})
	#MESSAGE(STATUS "adding property to ${_target} ${_name}:" ${_properties})
	IF(NOT _old_properties)
		# in case it's NOTFOUND
		SET(_old_properties)
	ENDIF(NOT _old_properties)
	SET_TARGET_PROPERTIES(${_target} PROPERTIES ${_name} "${_old_properties} ${_properties}")
ENDMACRO(ADD_TARGET_PROPERTIES)

set(DEBUG_FLAGS " ${DEBUG_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
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
	${CMAKE_SOURCE_DIR}
)

set(OSMPBF_LIBRARIES ${osmpbf_BINARY_DIR}/libosmpbf.a)

add_executable(pbf2route pbf2route.cpp)
add_dependencies(pbf2route osmpbf)
target_link_libraries(pbf2route ${OSMPBF_LIBRARIES} ${ZLIB_LIBRARIES} ${PROTOBUF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(pbf2route_openmp pbf2route.cpp)
add_dependencies(pbf2route_openmp osmpbf)
target_link_libraries(pbf2route_openmp ${OSMPBF_LIBRARIES} ${ZLIB_LIBRARIES} ${PROTOBUF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_target_properties(pbf2route_openmp COMPILE_FLAGS "-fopenmp")
add_target_properties(pbf2route_openmp LINK_FLAGS "-fopenmp")


add_executable(pbf2route_parsehelpers pbf2route_parsehelpers.cpp)
add_dependencies(pbf2route_parsehelpers osmpbf)
target_link_libraries(pbf2route_parsehelpers ${OSMPBF_LIBRARIES} ${ZLIB_LIBRARIES} ${PROTOBUF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_target_properties(pbf2route_parsehelpers COMPILE_FLAGS "-fopenmp")
add_target_properties(pbf2route_parsehelpers LINK_FLAGS "-fopenmp")

add_executable(osmpbf_route osmpbf_route.cpp)
add_dependencies(osmpbf_route osmpbf)
target_link_libraries(osmpbf_route ${OSMPBF_LIBRARIES} ${ZLIB_LIBRARIES} ${PROTOBUF_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
