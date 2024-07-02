# project/CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(PET-Project)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add subdirectories
add_subdirectory(sim)
add_subdirectory(digi)
