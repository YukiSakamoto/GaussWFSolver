project(wf_solver CXX)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "--std=c++11 -O3")
add_executable(wf_solver libint_adapter.cpp HartreeFock.cpp InputParser.cpp main.cpp )

target_link_libraries(wf_solver int2)

# Copy the example files to output directory
#set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_BINARY_DIR})

#add_subdirectory(example)
#add_subdirectory(test)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
file(GLOB EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/example/*)
file(COPY ${EXAMPLES} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
