cmake_minimum_required(VERSION 3.1)
project(vesicle2 VERSION 2.0.0)

message(WARNING "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(WARNING "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(WARNING "BUILD TYPE: ${CMAKE_BUILD_TYPE}")


set(CMAKE_VERBOSE_MAKEFILE OFF)

add_subdirectory(lib)
add_subdirectory(src)

set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math -mtune=generic")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math -mtune=generic")

### vesicle
set(BASIC_SOURCE src/vesicle.cpp)
add_executable(vesicle2 ${BASIC_SOURCE})

set_target_properties(vesicle2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "~/bin/")
set_target_properties(vesicle2 PROPERTIES CXX_STANDARD 17)
target_compile_options(vesicle2 PRIVATE 
    -Wall 
    -Wextra 
    -Wpedantic 
    -Weffc++
    -march=nocona
    -mtune=nocona
    -std=c++17
)

# Link
target_link_libraries(vesicle2 vesicle_lib)

install( TARGETS vesicle2 RUNTIME DESTINATION bin )