
# gather all steps related portability configuration


# C++ support in g++ pre-version 4.8 is insufficient
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8)
	add_definitions( -DGNU_FORCE_INLINE=[[gnu::always_inline]] ) 
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    message(FATAL_ERROR "C++ compiler is g++ ${CMAKE_CXX_COMPILER_VERSION}; require version 4.7 or higher.")
endif()



# add some OSX helpers for BLAS / Lapack
## include OpenBLAS path by default for libpath
if(APPLE)
	list(APPEND CMAKE_SYSTEM_LIBRARY_PATH "/usr/local/opt/openblas/lib")
	list(APPEND CMAKE_SYSTEM_INCLUDE_PATH "/usr/local/opt/openblas/include")
endif()
