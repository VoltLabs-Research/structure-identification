#
# Copyright (C) Rodolfo Herrera Hernandez. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(structure-identification VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

# Find CoreToolkit dependency
if(NOT TARGET coretoolkit::coretoolkit)
	find_package(coretoolkit REQUIRED)
endif()
if(NOT EMSCRIPTEN)
	if(NOT TARGET TBB::tbb)
		find_package(TBB REQUIRED)
	endif()
	if(NOT TARGET duckdb::duckdb)
		find_package(duckdb REQUIRED)
	endif()
endif()
if(NOT TARGET Boost::headers AND NOT TARGET Boost::boost AND NOT TARGET boost::headers AND NOT TARGET boost::boost)
	find_package(Boost REQUIRED)
endif()
if(NOT TARGET spdlog::spdlog)
	find_package(spdlog REQUIRED)
endif()
if(NOT TARGET nlohmann_json::nlohmann_json)
	find_package(nlohmann_json REQUIRED)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp)
	find_package(yaml-cpp REQUIRED)
endif()

set(STRUCTURE_IDENTIFICATION_BOOST_TARGET "")
if(TARGET Boost::headers)
	set(STRUCTURE_IDENTIFICATION_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
	set(STRUCTURE_IDENTIFICATION_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
	set(STRUCTURE_IDENTIFICATION_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
	set(STRUCTURE_IDENTIFICATION_BOOST_TARGET boost::boost)
else()
	message(FATAL_ERROR "StructureIdentification requires a Boost headers target")
endif()

set(LIBRARY_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/structure_analysis_context.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/structure_identification_export.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/reconstructed_dump_utils.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/cluster_graph_io.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/cluster_graph_builder.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/orientation_cluster_rule_provider.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/crystal_topology_library.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/reconstructed_state_canonicalizer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/structure_analysis.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/core/reconstructed_structure.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/structures/crystal_topology_registry.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/structures/cluster_graph.cpp
)

add_library(${PROJECT_NAME} STATIC ${LIBRARY_SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../OpenDXA/lattices")
	target_compile_definitions(${PROJECT_NAME} PUBLIC
		VOLT_OPENDXA_LATTICE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../OpenDXA/lattices"
	)
endif()

set(_SI_LINK
	${STRUCTURE_IDENTIFICATION_BOOST_TARGET}
	coretoolkit::coretoolkit
	spdlog::spdlog
	nlohmann_json::nlohmann_json
	yaml-cpp::yaml-cpp
)
if(NOT EMSCRIPTEN)
	list(APPEND _SI_LINK TBB::tbb duckdb::duckdb)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC ${_SI_LINK})

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)
