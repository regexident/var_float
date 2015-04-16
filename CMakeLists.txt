cmake_minimum_required(VERSION 3.1)
project(var_float)

include_directories(
    ${PROJECT_SOURCE_DIR}/tests/catch/include
    ${PROJECT_SOURCE_DIR}/include
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(TESTS_SOURCE_FILES
	src/var_float.c
	tests/var_float_tests.cpp
)

add_executable(tests
	${TESTS_SOURCE_FILES}
)
