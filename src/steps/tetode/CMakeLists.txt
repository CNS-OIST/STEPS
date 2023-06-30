add_library(stepstetode STATIC comp.cpp patch.cpp tet.cpp tri.cpp tetode.cpp)

set_property(TARGET stepstetode PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepstetode PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepstetode stepsutil stepssolver SUNDIALS::SUNDIALS)
