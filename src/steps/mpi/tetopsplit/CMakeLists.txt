add_library(
  stepstetopsplit STATIC
  comp.cpp
  diff.cpp
  sdiff.cpp
  kproc.cpp
  patch.cpp
  reac.cpp
  sreac.cpp
  tet.cpp
  tetopsplit.cpp
  tri.cpp
  ghkcurr.cpp
  vdepsreac.cpp
  diffboundary.cpp
  wmvol.cpp
  sdiffboundary.cpp)

set_property(TARGET stepstetopsplit PROPERTY POSITION_INDEPENDENT_CODE ON)

set_property(TARGET stepstetopsplit PROPERTY COMPILE_FLAGS "-Wno-old-style-cast")

target_include_directories(stepstetopsplit PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepstetopsplit PUBLIC stepsutil stepsmpi stepssolver)
