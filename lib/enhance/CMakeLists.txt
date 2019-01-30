message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB enhance_source_files
    *.cpp
    *.hpp
)

add_library(enhance STATIC ${enhance_source_files})

target_link_libraries(enhance dl)
target_compile_options(enhance PRIVATE -std=c++14 -O2 -g0 -Wall -Wextra -Wpedantic -march=nocona )
target_compile_definitions(enhance PRIVATE NDEBUG EIGEN_NO_DEBUG)