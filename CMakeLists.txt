cmake_minimum_required(VERSION 3.0)
project(optimization_framework)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_executable(optmF optFrameWork.cpp)
#target_link_libraries(optmF  optmMethod)




