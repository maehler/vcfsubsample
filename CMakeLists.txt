cmake_minimum_required(VERSION 3.6.1)
project(vcfsubsample)
set(vcfsubsample_MAJOR 0)
set(vcfsubsample_MINOR 1)
set(vcfsubsample_PATCH 0)
set(vcfsubsample_VERSION_STRING "\"${vcfsubsample_MAJOR}.${vcfsubsample_MINOR}.${vcfsubsample_PATCH}\"")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/cmake-argp/cmake" ${CMAKE_MODULE_PATH})

find_package(argp REQUIRED)
find_path(HTSLIB_INCLUDE_DIR htslib/vcf.h)
find_library(HTSLIB_LIBRARIES hts)

if (NOT HTSLIB_INCLUDE_DIR)
  message(FATAL_ERROR "could not find htslib headers")
else()
  message("htslib include dir = ${HTSLIB_INCLUDE_DIR}")
endif()

if (NOT HTSLIB_LIBRARIES)
  message(FATAL_ERROR "could not find libhts")
else()
  message("htslib lib = ${HTSLIB_LIBRARIES}")
endif()

set(LIBS ${LIBS} ${HTSLIB_LIBRARIES} ${ARGP_LIBRARIES})
include_directories("${CMAKE_SOURCE_DIR}/include")

configure_file("${CMAKE_SOURCE_DIR}/include/config.h.in"
               "${CMAKE_SOURCE_DIR}/include/config.h")

file(GLOB SOURCES "src/*.c")

add_executable(vcfsubsample ${SOURCES})
target_link_libraries(vcfsubsample ${LIBS})
