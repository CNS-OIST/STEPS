add_library(stepswmdirect STATIC
    comp.cpp
    kproc.cpp
    patch.cpp
    reac.cpp
    sreac.cpp
    wmdirect.cpp
)

set_property(TARGET stepswmdirect PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepswmdirect PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepswmdirect stepsutil stepssolver)