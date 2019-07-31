cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_FLAGS "-std=c++11")

project (OpenLBM)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

file(GLOB SRC src/*.cpp)
include_directories(head)
add_executable(lbm ${SRC})
