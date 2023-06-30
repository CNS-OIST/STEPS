# Utility function to mark as advanced all variables matching a regular expression pattern. This is
# useful to do when importing projects.

function(mark_as_advanced_all)
  get_cmake_property(_variableNames VARIABLES)
  foreach(_variableName ${_variableNames})
    if(ARGV0)
      unset(MATCHED)
      string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
      if(NOT MATCHED)
        continue()
      endif()
      mark_as_advanced(${_variableName})
    endif()
  endforeach()
endfunction()
