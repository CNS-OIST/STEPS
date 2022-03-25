add_library(stepswmrk4 STATIC
    wmrk4.cpp
)

set_property(TARGET stepswmrk4 PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepswmrk4 PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepswmrk4 stepsutil stepssolver)