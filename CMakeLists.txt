cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

#Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)


SET(Prob3plusplus_VERSION 3.10.2)
project(Prob3plusplus VERSION ${Prob3plusplus_VERSION} LANGUAGES CXX C)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
elseif(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

SET(HEADERS EarthDensity.h 
	BargerPropagator.h 
	mosc.h
	mosc3.h)

SET(SOURCE EarthDensity.cc 
	BargerPropagator.cc 
	mosc.c
	mosc3.c)

add_library(Prob3plusplus SHARED ${SOURCE})

set_target_properties(Prob3plusplus PROPERTIES 
	PUBLIC_HEADER "${HEADERS}"
	EXPORT_NAME Prob3plusplus)

target_include_directories(
  Prob3plusplus PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                       $<INSTALL_INTERFACE:include>
)

install(TARGETS Prob3plusplus 
		EXPORT Prob3plusplus-target
		LIBRARY DESTINATION lib/
		PUBLIC_HEADER DESTINATION include/)

install(EXPORT Prob3plusplus-target
  FILE Prob3plusplusTargets.cmake
  NAMESPACE Prob3plusplus::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  Prob3plusplusConfig.cmake.in ${CMAKE_BINARY_DIR}/Prob3plusplusConfig.cmake
  INSTALL_DESTINATION 
  	/this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

install(FILES 
		${CMAKE_BINARY_DIR}/Prob3plusplusConfig.cmake 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)

write_basic_package_version_file(${CMAKE_INSTALL_PREFIX}/cmake/Prob3plusplusConfigVersion.cmake
  VERSION ${Prob3plusplus_VERSION}
  COMPATIBILITY AnyNewerVersion)