add_library(stepsgeom STATIC
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

set_property(TARGET stepsgeom PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsgeom PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsgeom stepsmath stepsutil stepsmodel)

if(STEPS_USE_DIST_MESH)
    add_subdirectory(dist)
endif()
