add_library(instrumentor INTERFACE)
target_include_directories(instrumentor INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if(MPI_FOUND)
    target_compile_definitions(instrumentor INTERFACE STEPS_USE_MPI)
    target_link_libraries(instrumentor INTERFACE MPI::MPI_CXX)
endif()

# =============================================================================
# Profiling/Instrumentation Options with Caliper
# =============================================================================
if(STEPS_USE_CALIPER_PROFILING)
  find_package(caliper REQUIRED)
  target_include_directories(instrumentor INTERFACE ${caliper_INCLUDE_DIR})
  target_compile_definitions(instrumentor INTERFACE STEPS_CALIPER)
  target_link_libraries(instrumentor INTERFACE caliper)
endif()
message(STATUS "Caliper: ${STEPS_USE_CALIPER_PROFILING}")

# =============================================================================
# Profiling/Instrumentation Options with Likwid
# =============================================================================
if(STEPS_USE_LIKWID_PROFILING)
  find_package(likwid 5 REQUIRED)
  target_include_directories(instrumentor INTERFACE ${likwid_INCLUDE_DIR})
  target_compile_definitions(instrumentor INTERFACE LIKWID_PERFMON)
  target_link_libraries(instrumentor INTERFACE ${likwid_LIBRARIES})
endif()
message(STATUS "Likwid: ${STEPS_USE_LIKWID_PROFILING}")

# =============================================================================
# Profiling/Instrumentation with STEPS region tracker
# =============================================================================
if(STEPS_USE_NATIVE_PROFILING)
  target_compile_definitions(instrumentor INTERFACE STEPS_REGION_TRACKER)
endif()
message(STATUS "Region tracker: ${STEPS_USE_NATIVE_PROFILING}")
