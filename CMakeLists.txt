cmake_minimum_required(VERSION 3.5)

project(nelder-mead-test LANGUAGES CXX)

#set(CMAKE_VERBOSE_MAKEFILE ON)

set(example_programs
    beale
    booth
    cosine
    distance
    example
    rosenbrock
    squared
    )

file(GLOB HEADERS
    ${CMAKE_SOURCE_DIR}/*.h
    )

foreach(example_program ${example_programs})
    set(source ${example_program}.cpp)
    add_executable(${example_program} ${source} ${HEADERS})

target_compile_options(${example_program} PUBLIC
    -Wall -Wextra
    $<$<CONFIG:Release>:-O3 -Wall -Werror>
    $<$<CONFIG:Debug>:-g3 -O0 -fno-inline -Wall -Werror>
    $<$<CONFIG:Profile>:-O3 -pg -Wall -Werror>
    )
set_property(TARGET ${example_program} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${example_program} PROPERTY CXX_EXTENSIONS OFF)
set_property(TARGET ${example_program} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(${example_program} PUBLIC
    ${CMAKE_SOURCE_DIR}
    )

target_link_libraries(${example_program})
endforeach()
