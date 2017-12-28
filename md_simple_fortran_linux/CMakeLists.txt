cmake_minimum_required(VERSION 3.9)
project(md_simple_fortran)

enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_FLAGS_RELEASE "-fbounds-check -fbacktrace")
set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace")

add_executable(md_simple_fortran main.f90 init_random_seed.f90)

