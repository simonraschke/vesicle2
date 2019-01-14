cmake_minimum_required(VERSION 3.1)
project(vesicle VERSION 2.0.0)

message(WARNING "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(WARNING "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(WARNING "BUILD TYPE: ${CMAKE_BUILD_TYPE}")


set(CMAKE_VERBOSE_MAKEFILE OFF)

add_subdirectory(lib)
add_subdirectory(src)

set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

### vesicle
set(BASIC_SOURCE src/vesicle.cpp)
add_executable(vesicle ${BASIC_SOURCE})

set_target_properties(vesicle PROPERTIES RUNTIME_OUTPUT_DIRECTORY "~/bin/")
set_target_properties(vesicle PROPERTIES CXX_STANDARD 17)
target_compile_options(vesicle PRIVATE -Wall -Wextra -Wpedantic -ferror-limit=2 -march=nocona -std=c++17)

# Link
target_link_libraries(vesicle vesicle_lib)

install( TARGETS vesicle RUNTIME DESTINATION bin )