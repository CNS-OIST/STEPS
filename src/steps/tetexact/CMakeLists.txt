add_library(
  stepstetexact STATIC
  comp.cpp
  diff.cpp
  sdiff.cpp
  kproc.cpp
  patch.cpp
  reac.cpp
  sreac.cpp
  tet.cpp
  tetexact.cpp
  tri.cpp
  ghkcurr.cpp
  vdepsreac.cpp
  diffboundary.cpp
  wmvol.cpp
  sdiffboundary.cpp)

set_property(TARGET stepstetexact PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsmodel PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepstetexact stepsutil stepssolver)
