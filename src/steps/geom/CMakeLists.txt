add_library(stepsgeo STATIC
    tetmesh.cpp
    comp.cpp
    geom.cpp
    patch.cpp
    tmcomp.cpp
    tmpatch.cpp
    sdiffboundary.cpp
    memb.cpp
    diffboundary.cpp
)

set_property(TARGET stepsgeo PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsgeo
    PUBLIC .
    PUBLIC ../../)

target_link_libraries(stepsgeo stepsmath stepsutil)