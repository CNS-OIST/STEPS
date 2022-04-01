add_library(stepswmrssa OBJECT
    comp.cpp
    kproc.cpp
    patch.cpp
    reac.cpp
    sreac.cpp
    wmrssa.cpp
)

set_property(TARGET stepswmrssa PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepswmrssa PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepswmrssa stepsutil)
