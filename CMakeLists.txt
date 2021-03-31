project(wf_solver CXX)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "--std=c++11")
add_executable(wf_solver libint_adapter.cpp HartreeFock.cpp InputParser.cpp main.cpp )

target_link_libraries(wf_solver int2)

# Copy the example files to output directory
#set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_BINARY_DIR})

#add_subdirectory(example)
#add_subdirectory(test)

