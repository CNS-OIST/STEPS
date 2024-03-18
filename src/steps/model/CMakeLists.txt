add_library(
  stepsmodel STATIC
  chan.cpp
  chanstate.cpp
  complex.cpp
  complexevents.cpp
  complexreac.cpp
  complexsreac.cpp
  diff.cpp
  endocytosis.cpp
  exocytosis.cpp
  ghkcurr.cpp
  linkspec.cpp
  model.cpp
  ohmiccurr.cpp
  raft.cpp
  raftdis.cpp
  raftendocytosis.cpp
  raftgen.cpp
  raftsreac.cpp
  raftsys.cpp
  reac.cpp
  spec.cpp
  sreac.cpp
  surfsys.cpp
  vdepsreac.cpp
  vesbind.cpp
  vesicle.cpp
  vessdiff.cpp
  vessreac.cpp
  vessurfsys.cpp
  vesunbind.cpp
  volsys.cpp)

set_property(TARGET stepsmodel PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepsmodel PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

target_link_libraries(stepsmodel stepsutil)
