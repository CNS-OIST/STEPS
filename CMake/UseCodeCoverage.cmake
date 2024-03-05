# Enable Code Coverage
# ~~~
# Variables you may define are:
#  CODECOV_HTMLOUTPUTDIR - the name of the directory where HTML results are placed. Defaults to "coverage_html"
#  CODECOV_XMLOUTPUTFILE - the name of the directory where HTML results are placed. Defaults to "coverage.xml"
#  CODECOV_GCOVR_OPTIONS - additional options given to gcovr commands.
# ~~~

option(ENABLE_CODECOVERAGE "Enable code coverage testing support")

if(ENABLE_CODECOVERAGE)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Code coverage results with an optimised (non-Debug) build may be misleading")
  endif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")

  if(NOT DEFINED CODECOV_XMLOUTPUTFILE)
    set(CODECOV_XMLOUTPUTFILE coverage.xml)
  endif(NOT DEFINED CODECOV_XMLOUTPUTFILE)

  if(NOT DEFINED CODECOV_HTMLOUTPUTDIR)
    set(CODECOV_HTMLOUTPUTDIR coverage_html)
  endif(NOT DEFINED CODECOV_HTMLOUTPUTDIR)

  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCXX)
    foreach(exe gcov lcov gcovr)
      string(TOUPPER ${exe} exe_var)
      set(exe_var "CODECOV_${exe_var}")
      find_program(${exe_var} ${exe})
      if(NOT ${exe_var})
        message(SEND_ERROR "Could not find ${exe} executable")
      endif()
    endforeach()
    add_definitions(-fprofile-arcs -ftest-coverage)
    link_libraries(gcov)
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} --coverage)
    add_custom_target(
      coverage_init
      ${CODECOV_LCOV}
      --gcov-tool
      ${CODECOV_GCOV}
      --base-directory
      .
      --directory
      ${CMAKE_BINARY_DIR}
      --capture
      --initial)
    string(REPLACE " " ";" codecov_gcovr_opts "${CODECOV_GCOVR_OPTIONS}")
    set(CODECOV_GCOVR_BASECMD ${CODECOV_GCOVR} ${codecov_gcovr_opts} --gcov-executable
                              ${CODECOV_GCOV} --root ${PROJECT_SOURCE_DIR})
    add_custom_target(
      coverage
      ${CODECOV_GCOVR_BASECMD} --xml --output ${CODECOV_XMLOUTPUTFILE}
      COMMAND mkdir -p ${CODECOV_HTMLOUTPUTDIR}
      COMMAND ${CODECOV_GCOVR_BASECMD} --html --html-details --output
              ${CODECOV_HTMLOUTPUTDIR}/index.html)
  endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCXX)
endif(ENABLE_CODECOVERAGE)
