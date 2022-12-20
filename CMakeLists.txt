cmake_minimum_required(VERSION 3.20)
project(Blinn)

#Flag
set(CMAKE_CXX_STANDARD 20)
add_compile_options(--std=c++20)
add_compile_options(-mfma4)
add_compile_options(-Wall)
add_compile_options(-Wextra)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)

#Submodule
include_directories(${PROJECT_SOURCE_DIR}/)
add_subdirectory(excerpt)
link_directories(${PROJECT_SOURCE_DIR}/excerpt)

add_executable(${PROJECT_NAME} app.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC excerpt BlinnSolver)



