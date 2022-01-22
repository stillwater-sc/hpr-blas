########################################################################
# CMakeLists.txt
#
# Authors: Matthias Moller, Nauman Ahmed, Theodore Omtzigt
########################################################################
include(cmake/Banners.cmake)
print_header()

# Set minimum version of CMake. Since we are using the VERSION option of the
# project command, we need at least 3.0. To honor the amazing work that
# the folks at KitWare have done for the open-source community, we are
# going to specify a recent version. 
# As of UNIVERSAL v3.0 December 2019
# Ubuntu 16.04 LTS runs cmake 3.5.1
# Ubuntu 18.04 LTS runs cmake 3.10.2
# container runs cmake 3.7.1
cmake_minimum_required(VERSION 3.16)

####
## Enable project() command to manage VERSION variables
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif(POLICY CMP0048)

####
# Set project variables
if(NOT DEFINED HPR_BLAS_VERSION_MAJOR)
  set(HPR_BLAS_VERSION_MAJOR 1)
endif()
if(NOT DEFINED UNIVERSAL_VERSION_MINOR)
  set(HPR_BLAS_VERSION_MINOR 11)
endif()
if(NOT DEFINED UNIVERSAL_VERSION_PATCH)
  set(HPR_BLAS_VERSION_PATCH 1)
endif()

########################################################################
# Project name
########################################################################
# HPR-BLAS: high-performance reproducible linear algebra
project (hpr-blas 
	DESCRIPTION "A header only C++ library for High-Performance Reproducible BLAS"
	VERSION "${HPR_BLAS_VERSION_MAJOR}.${HPR_BLAS_VERSION_MINOR}.${HPR_BLAS_VERSION_PATCH}" 
	LANGUAGES C CXX
	HOMEPAGE_URL "https://github.com/stillwater-sc/hpr-blas")

####
# Change default build type to Release
#
# The CACHE STRING logic here and elsewhere is needed to force CMake
# to pay attention to the value of these variables.
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No default build type specified: setting CMAKE_BUILD_TYPE=Release")
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the build type: options are: Debug Release RelWithDebInfo MinSizeRel"
		FORCE)
else(NOT CMAKE_BUILD_TYPE)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("====================================================================================")
		message(STATUS "Build type is set to Debug: Performance will be negatively impacted")
		message(STATUS "Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build")
		message("====================================================================================")
	endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
endif(NOT CMAKE_BUILD_TYPE)

########################################################################
# We do not support in-source build
########################################################################
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not permitted.\nPlease make a separate folder for building, otherwise type \nmake\nthat will create a ./build folder and will compile in that folder. Cmake has created garbage files/dirs (you may manually remove them):\nCMakeCache.txt CMakeFiles")
endif()

####
# Create the library target
set(project_library_target_name ${PROJECT_NAME})
set(PACKAGE_NAME hprblas)

####
# Set environmental options for tracing, testing, and verbosity
option(HPR_BLAS_CMAKE_TRACE "Tracing CMake results, i.e. printing variable settings." OFF)
option(HPR_BLAS_ENABLE_TESTS "Enable the build and run of tests." ON)
option(HPR_BLAS_VERBOSE_TESTS "Always print test output, otherwise only errors. Only relevant when tests enabled." OFF)

########################################################################
# Set RPATH on MacOSX
########################################################################
if (APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  SET(CMAKE_SKIP_BUILD_RPATH FALSE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF ("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF ()
endif()

########################################################################
# Use solution folders for Visual Studio
########################################################################
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

macro(trace_variable variable)
    if (HPR_BLAS_CMAKE_TRACE)
        message(STATUS "${variable} = ${${variable}}")
    endif()
endmacro()
 
########################################################################
# Augment the CMAKE_MODULE_PATH to pick up custom cmake modules
########################################################################       
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

########################################################################
# Configure the compiler optimization flags
########################################################################
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  include( OptimizeForArchitecture )
  OptimizeForArchitecture()
  foreach (flag ${OFA_ARCHITECTURE_FLAGS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  endforeach()
endif()

# Global configuration
include(Config)

# Must be located in root dir, doesn't work in tests
if (HPR_BLAS_ENABLE_TESTS)
    enable_testing()
    # include(Dart)
endif()

####
## Enable package_ROOT variables
if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)
find_package(Boost)
message(STATUS "Boost_FOUND " ${Boost_FOUND})
if (Boost_FOUND)
	message(STATUS "Including Boost_INCLUDE_DIR " ${Boost_INCLUDE_DIR})
	include_directories(${Boost_INCLUDE_DIR})
endif (Boost_FOUND)

cmake_policy(SET CMP0054 NEW)
find_package(MTL REQUIRED)
include_directories(${MTL_INCLUDE_DIRS})
add_definitions(${MTL_CXX_DEFINITIONS})
message(STATUS "MTL INCLUDE DIR " ${MTL_INCLUDE_DIRS})

# Set posit include directory
find_package(UNIVERSAL REQUIRED)
include_directories(${UNIVERSAL_INCLUDE_DIRS})
message(STATUS "UNIVERSAL INCLUDE DIR " ${UNIVERSAL_INCLUDE_DIRS})

# Set hpr-blas include directory
set(HPR_BLAS_INCLUDE_DIR "./include")
include_directories(${HPR_BLAS_INCLUDE_DIR})
message(STATUS "HPR-BLAS INCLUDE DIR " ${HPR_BLAS_INCLUDE_DIR})

####
# Configure the compiler options
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "C++20 support has been enabled by default")

macro (compile_all testing prefix folder)
    # cycle through the sources
    # For the according directories, we assume that each cpp file is a separate test
    # and create a executable target and an associated test target for that file
    foreach (source ${ARGN})
        get_filename_component (test ${source} NAME_WE)
        string(REPLACE " " ";" new_source ${source})
        set(test_name ${prefix}_${test})
        #message(STATUS "Add test ${test_name} from source ${new_source}.")
        add_executable (${test_name} ${new_source})
	#message(STATUS "args: ${testing} - ${prefix} - ${folder}")
	set_target_properties(${test_name} PROPERTIES FOLDER ${folder})
        if (${testing} STREQUAL "true")
			if (HPR_BLAS_CMAKE_TRACE)
                message(STATUS "testing: ${test_name} ${RUNTIME_OUTPUT_DIRECTORY}/${test_name}")
            endif()
            add_test(${test_name} ${RUNTIME_OUTPUT_DIRECTORY}/${test_name})
        endif()
    endforeach (source)
endmacro (compile_all)

# Set up install directories. INCLUDE_INSTALL_DIR and
# CMAKECONFIG_INSTALL_DIR must not be absolute paths.
if(WIN32)
    set(include_install_dir Include)
    set(include_install_dir_full Include)
    set(config_install_dir CMake)
elseif(UNIX)
    set(include_install_dir include)
    set(include_install_dir_postfix "${project_library_target_name}")
    set(include_install_dir_full    "${include_install_dir}/${include_install_dir_postfix}")

    set(config_install_dir share/${PACKAGE_NAME})
else()
    message(FATAL_ERROR "Not supported system type. Options: UNIX or WIN32.")
endif()

add_subdirectory("tools/benchmark")
#add_subdirectory("tools/characterization")
add_subdirectory("tools/matrix_generation")
add_subdirectory("blas")
add_subdirectory("applications")
add_subdirectory("benchmarks")
add_subdirectory("solvers")
add_subdirectory("c_api/lib")
add_subdirectory("c_api/test")

