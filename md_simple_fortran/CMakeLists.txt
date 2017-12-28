cmake_minimum_required(VERSION 3.9)
project(md_simple_fortran)

enable_language (Fortran)

add_executable(md_simple_fortran main.f90 init_random_seed.f90)

