cmake_minimum_required(VERSION 2.6)
project(isONclust2)

set(Tutorial_VERSION_MAJOR 0)
set(Tutorial_VERSION_MINOR 1)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

configure_file(
        "${PROJECT_SOURCE_DIR}/src/Config.h.in"
        "${PROJECT_BINARY_DIR}/Config.h"
    )

include_directories("${PROJECT_BINARY_DIR}")

set(SOURCE 
    src/main.cpp
    )

add_executable(isONclust2 ${SOURCE})
