cmake_minimum_required(VERSION 3.1.0)
project(xgtopw VERSION 1.1.0 LANGUAGES C CXX)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
    message("debug mode")
endif()


add_executable(${PROJECT_NAME} ${sources} ${headers})

set_target_properties(xgtopw PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/headers)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -fext-numeric-literals)

target_compile_definitions(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:GTOPW_DEBUG>")

find_package(OpenMP REQUIRED)
find_library(Quadmath NAMES libquadmath.so libquadmath.so.0)
find_library(Math NAMES libm.so)

target_link_libraries(xgtopw PUBLIC ${Math} ${Quadmath} OpenMP::OpenMP_CXX)
