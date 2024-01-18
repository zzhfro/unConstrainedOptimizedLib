# 设置项目名称
project(MyProject)

# 设置 CMake 的最低版本
cmake_minimum_required(VERSION 3.10)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
# 添加头文件目录
include_directories(include)
include_directories("/usr/include/eigen3")
# 查找所有源文件
file(GLOB SOURCES "src/*.cpp")

# 生成可执行文件
add_executable(testOpt ${SOURCES} testOpt.cpp)
add_executable(testExample  testExample.cpp)
# 如果有其他的库文件，可以使用以下方式链接
# target_link_libraries(my_executable my_library)




