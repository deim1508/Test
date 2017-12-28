cmake_minimum_required(VERSION 3.9)
project(c_thread C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(c_thread main.c thpool.c)

target_link_libraries(c_thread pthread)