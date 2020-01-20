cmake_minimum_required(VERSION 3.14)
project(single_chain)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

add_executable(single_chain main.cpp Propagator.cpp Propagator.h SCFT.cpp SCFT.h qAPropagator.cpp qAPropagator.h qBPropagator.cpp qBPropagator.h)
target_link_libraries(single_chain fftw3)