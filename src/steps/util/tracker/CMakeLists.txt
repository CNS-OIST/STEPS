add_library(stepstracker STATIC memory_tracker.cpp peak_rss.cpp region_tracker.cpp time_tracker.cpp)

set_property(TARGET stepstracker PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(stepstracker PUBLIC "${PROJECT_SOURCE_DIR}/src/steps")

if(MPI_FOUND)
  target_compile_definitions(stepstracker PUBLIC STEPS_USE_MPI)
  target_link_libraries(stepstracker MPI::MPI_CXX)
endif()
