cmake_minimum_required(VERSION 3.20)
project(Blinn)

add_compile_options(--std=c++20)
add_compile_options(-mfma4)
add_compile_options(-Wall)
add_compile_options(-Wextra)


set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/)
add_subdirectory(excerpt)
link_directories(${PROJECT_SOURCE_DIR}/excerpt)
add_executable(${CMAKE_PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC excerpt Eigen3::Eigen)

