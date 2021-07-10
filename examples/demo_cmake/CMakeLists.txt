# This is demo on how to link libpointmatcher to an
# external project.

cmake_minimum_required (VERSION 3.10)
project(convert)

set(CMAKE_CXX_STANDARD 11)

find_package(libpointmatcher REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${libpointmatcher_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})

add_executable(${PROJECT_NAME} convert.cpp)

target_link_libraries(${PROJECT_NAME} ${libpointmatcher_LIBRARIES} ${Boost_LIBRARIES})
