message(${CMAKE_CURRENT_SOURCE_DIR})
# project(vesicle_lib)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

file(GLOB vesicle_source_files
    */*.cpp
    */*/*.cpp
)

add_library(vesicle_lib STATIC ${vesicle_source_files})

target_compile_definitions(vesicle_lib PRIVATE 
    USE_MATH_DEFINES 
    EIGEN_FAST_MATH=1 
    EIGEN_MAX_CPP_VER=17
# #     EIGEN_DONT_PARALLELIZE
)

set(CMAKE_CXX_FLAGS_DEBUG "-Og -ggdb -D_DEBUG -DDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

target_compile_options(vesicle_lib PRIVATE 
    -Wall 
    -Wextra 
    -Wpedantic 
    -Weffc++
    -march=nocona 
    -std=c++17
)

target_link_libraries(vesicle_lib enhance tbb pthread boost_program_options stdc++fs)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/lib)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/lib/h5xx)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/lib/observer-ptr/include/nonstd)