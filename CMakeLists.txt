cmake_minimum_required(VERSION 3.15)
project(Directional LANGUAGES CXX)

add_library(directional INTERFACE)

target_include_directories(directional INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Add polyscope dependency
add_subdirectory(external/polyscope)
target_link_libraries(directional INTERFACE polyscope)
