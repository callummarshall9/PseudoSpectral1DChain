cmake_minimum_required(VERSION 3.14)
project(single_chain)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

add_executable(single_chain main.cpp Propagator.cpp Propagator.h SCFT.cpp SCFT.h qPropagator.cpp qPropagator.h qDaggerPropagator.cpp qDaggerPropagator.h FieldMethod.h mathutils.h)
find_package(PkgConfig) # With pkg-config
pkg_check_modules(FFTW3 REQUIRED fftw3)
target_include_directories(single_chain PUBLIC ${FFTW3_INCLUDE_DIRS})
target_link_libraries(single_chain ${FFTW3_LINK_LIBRARIES})
