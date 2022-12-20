set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/BlinnSolver/solver.hpp")
add_library(BlinnSolver solver.cpp ${HEADER_LIST})
target_include_directories(BlinnSolver PUBLIC ../include)
include_directories(${PROJECT_SOURCE_DIR}/)
target_link_libraries(BlinnSolver PUBLIC excerpt Eigen3::Eigen)
target_compile_features(BlinnSolver PUBLIC cxx_std_20)

