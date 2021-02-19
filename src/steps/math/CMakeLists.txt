add_library(stepsmath STATIC
    tetrahedron.cpp
    tools.cpp
    linsolve.cpp
    triangle.cpp
    ghk.cpp
)

set_property(TARGET stepsmath PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsmath
    PUBLIC .
    PUBLIC ../../)
