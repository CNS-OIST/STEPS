add_library(stepsmodel STATIC
    model.cpp
    diff.cpp
    chan.cpp
    reac.cpp
    spec.cpp
    sreac.cpp
    surfsys.cpp
    volsys.cpp
    chanstate.cpp
    ohmiccurr.cpp
    ghkcurr.cpp
    vdeptrans.cpp
    vdepsreac.cpp
)

set_property(TARGET stepsmodel PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsmodel PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsmodel stepsutil)