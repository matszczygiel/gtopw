cmake_minimum_required(VERSION 3.1.0)
project(xgtopw VERSION 1.1.0 LANGUAGES C CXX)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h)

add_executable(${PROJECT_NAME} ${sources} ${headers})

set_target_properties(xgtopw PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/headers)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -pipe -fext-numeric-literals)

find_package(OpenMP REQUIRED)
target_link_libraries(xgtopw -lm -lquadmath OpenMP::OpenMP_CXX)
