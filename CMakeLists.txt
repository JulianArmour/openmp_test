cmake_minimum_required(VERSION 3.10)
project(openmp_test) # you can change the project name

find_package(OpenMP)

#add_executable(openmptest openmptest.c)
add_executable(openmptest q2.c)

if(OpenMP_CXX_FOUND)
    message("found openMP")
    target_link_libraries(openmptest PUBLIC OpenMP::OpenMP_CXX)
endif()