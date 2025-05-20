
cmake_minimum_required(VERSION 3.10)

project(TEST1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Debug build settings
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-g -O0")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

# Release build settings
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -flto -ffast-math -funroll-loops -fno-math-errno")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

 set(SOURCE_DIR ${CMAKE_SOURCE_DIR})



add_executable(Pruning
    ${SOURCE_DIR}/PruningToolMain.cpp

)



target_compile_options(Pruning PUBLIC -march=native -mtune=native -ffast-math)

# Include header files
include_directories(${CMAKE_SOURCE_DIR}/include)


