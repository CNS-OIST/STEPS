# /FindCLANG_TIDY.cmake
#
# This CMake script will search for clang-tidy and set the following
# variables
#
# CLANG_TIDY_FOUND : Whether or not clang-tidy is available on the target system
# CLANG_TIDY_EXECUTABLE : Fully qualified path to the clang-tidy executable
#
# The following variables will affect the operation of this script
# CLANG_TIDY_SEARCH_PATHS : List of directories to search for clang-tidy in,
#                           before searching any system paths. This should be
#                           the prefix to which clang-tidy was installed, and
#                           not the path that contains the clang-tidy binary.
#                           Eg /opt/ not /opt/bin/

# Search for the canonical executable, then search for ones with
# a version from newest to oldest.
find_program(CLANG_TIDY_EXECUTABLE
  NAMES clang-tidy clang-tidy-6.0 clang-tidy-5.0
        clang-tidy-4.0 clang-tidy-3.9 clang-tidy-3.8
        clang-tidy-3.7 clang-tidy-3.6 clang-tidy-3.5
  HINTS ${CLANG_TIDY_SEARCH_PATHS}
)

if(CLANG_TIDY_EXECUTABLE)
  set (CLANG_TIDY_FOUND TRUE)
else()
  set (CLANG_TIDY_FOUND FALSE)
endif()
