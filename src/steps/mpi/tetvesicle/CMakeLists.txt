add_library(
  stepstetvesicle STATIC
  comp_rdef.cpp
  comp_vesraft.cpp
  diffboundary.cpp
  diff.cpp
  endocytosis.cpp
  exocytosis.cpp
  ghkcurr.cpp
  kproc.cpp
  linkspec.cpp
  linkspecpair.cpp
  patch_rdef.cpp
  patch_vesraft.cpp
  path.cpp
  pointspec.cpp
  qtable.cpp
  raft.cpp
  raftdis.cpp
  raftendocytosis.cpp
  raftgen.cpp
  raftproxy.cpp
  raftsreac.cpp
  reac.cpp
  sdiffboundary.cpp
  sdiff.cpp
  sreac.cpp
  tet_rdef.cpp
  tetvesicle_rdef.cpp
  tetvesicle_vesraft.cpp
  tet_vesraft.cpp
  tri_rdef.cpp
  tri_vesraft.cpp
  vdepsreac.cpp
  vesbind.cpp
  vesicle.cpp
  vesproxy.cpp
  vesreac.cpp)

set_property(TARGET stepstetvesicle PROPERTY POSITION_INDEPENDENT_CODE ON)

set_property(TARGET stepstetvesicle PROPERTY COMPILE_FLAGS
                                             "-Wno-old-style-cast -Wno-unused-parameter")

target_include_directories(stepstetvesicle PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepstetvesicle PUBLIC stepsutil stepsmpi stepssolver GSL::gsl GSL::gslcblas)
