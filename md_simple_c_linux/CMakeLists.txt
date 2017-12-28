cmake_minimum_required(VERSION 3.9)
project(md_simple_c_linux C)

set(CMAKE_C_STANDARD 99)

add_executable(md_simple_c_linux main.c)

target_link_libraries(md_simple_c_linux m)