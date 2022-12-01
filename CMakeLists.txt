cmake_minimum_required(VERSION 3.20)
project(Blinn)

set(CMAKE_CXX_STANDARD 20)

include_directories(${PROJECT_SOURCE_DIR}/)
add_subdirectory(excerpt)
link_directories(${PROJECT_SOURCE_DIR}/excerpt)
add_executable(${CMAKE_PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC excerpt)
