# set the cmake version (>= 3.17 needed for CUDAToolkit)
cmake_minimum_required(VERSION 3.17...3.24)

# set project name
project(fasterjet VERSION 1.0 DESCRIPTION "faster than fastjet" LANGUAGES CXX)

# write binaries to the bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# enable CUDA language support
enable_language(CUDA)

# set the CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 61)

# set project header files
set(HEADER_FILES "include/fasterjet/DetailedException.hpp"
                 "include/fasterjet/jet.hpp"
                 "include/fasterjet/jetalgorithmsCPU.hpp"
                 "include/fasterjet/jetalgorithmsGPU.hpp"
                 "include/fasterjet/jetPrivate.hpp"
                 "include/fasterjet/precompiledheader.hpp"
                 "include/fasterjet/threadpool.hpp"
                 "include/fasterjet/timer.hpp"
                 "include/fasterjet/utils_avx.hpp"
                 "include/fasterjet/utils.hpp" )

# set project source files
set(SOURCE_FILES "src/jetalgorithmsCPU.cpp"
                 "src/jetPrivate.cpp"
                 "src/main.cpp" 
                 "src/utils.cpp" )

# set cuda source files
set(SOURCE_FILES_CUDA "src/jetalgorithmsGPU.cu" )

# setup fasterjet cuda target
add_library(fjcu ${SOURCE_FILES_CUDA})
find_package(CUDAToolkit REQUIRED)
target_link_libraries(fjcu PUBLIC CUDA::cudart)
set_target_properties(fjcu PROPERTIES CUDA_STANDARD 20)
target_include_directories(fjcu PUBLIC "include/" "${PROJECT_BINARY_DIR}/include/")

# setup fasterjet cpp target
add_executable(fj ${SOURCE_FILES} ${HEADER_FILES})

# doesnt seem to work with gcc
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
target_precompile_headers(fj PUBLIC "include/fasterjet/precompiledheader.hpp")
endif()

target_compile_features(fj PUBLIC cxx_std_20)
set_target_properties(fj PROPERTIES CXX_EXTENSIONS OFF)
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(fj PUBLIC Eigen3::Eigen)
target_link_libraries(fj PUBLIC fjcu)
target_include_directories(fj PUBLIC "include/" "${PROJECT_BINARY_DIR}/include/" "contrib/" "contrib/fastjet/include/")

add_library(fastjettools STATIC IMPORTED)
add_library(fastjet STATIC IMPORTED)
add_library(fastjetplugins STATIC IMPORTED)
add_library(siscone_spherical STATIC IMPORTED)
add_library(siscone STATIC IMPORTED)

set(tmpA "${CMAKE_CURRENT_SOURCE_DIR}/contrib/fastjet/buildd/lib/")
set(tmpB "${CMAKE_CURRENT_SOURCE_DIR}/contrib/fastjet/include/")

set_target_properties(fastjettools PROPERTIES IMPORTED_LOCATION "${tmpA}/libfastjettools.a" INTERFACE_INCLUDE_DIRECTORIES "${tmpB}")
set_target_properties(fastjet PROPERTIES IMPORTED_LOCATION "${tmpA}/libfastjet.a" INTERFACE_INCLUDE_DIRECTORIES "${tmpB}")
set_target_properties(fastjetplugins PROPERTIES IMPORTED_LOCATION "${tmpA}/libfastjetplugins.a" INTERFACE_INCLUDE_DIRECTORIES "${tmpB}")
set_target_properties(siscone_spherical PROPERTIES IMPORTED_LOCATION "${tmpA}/libsiscone_spherical.a" INTERFACE_INCLUDE_DIRECTORIES "${tmpB}")
set_target_properties(siscone PROPERTIES IMPORTED_LOCATION "${tmpA}/libsiscone.a" INTERFACE_INCLUDE_DIRECTORIES "${tmpB}")



# link_directories("contrib/fastjet/buildd/lib/")
# link_directories("${CMAKE_CURRENT_SOURCE_DIR}/contrib/fastjet/buildd/lib/")
target_link_libraries(fj PUBLIC fastjettools)
target_link_libraries(fj PUBLIC fastjet)
target_link_libraries(fj PUBLIC fastjetplugins)
target_link_libraries(fj PUBLIC siscone_spherical)
target_link_libraries(fj PUBLIC siscone)


if (UNIX)
  target_link_libraries(fj PUBLIC -lpthread -ldl)
endif (UNIX)

# clang
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -mavx2")
  message("using Clang compiler")

# gcc
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -mavx2 -march=native")
  message("using gcc compiler")
  
# icc
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND UNIX)
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -mavx2")
  message("using icc (Intel-gcc) compiler")
  
# icl
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND NOT UNIX)
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /arch:AVX2")
  message("using icl (Intel-MSVC) compiler")

# cl
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /arch:AVX2")
  message("using cl (MSVC) compiler")

endif()
