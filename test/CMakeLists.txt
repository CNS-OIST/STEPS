include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0)
FetchContent_MakeAvailable(Catch2)

include_directories("${CMAKE_SOURCE_DIR}/src")

add_subdirectory(unit)
add_subdirectory(validation)
