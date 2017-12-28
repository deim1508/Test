cmake_minimum_required(VERSION 3.9)
project(md_simple_c_linux_thread C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(md_simple_c_linux_thread main.c thpool.c)

target_link_libraries(md_simple_c_linux_thread m)