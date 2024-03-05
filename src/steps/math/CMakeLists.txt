add_library(stepsmath STATIC tetrahedron.cpp tools.cpp linsolve.cpp triangle.cpp ghk.cpp sphere.cpp)

set_property(TARGET stepsmath PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsmath PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsmath stepsutil)
