# Set variables representing corresponding command-line options for
# various compiler functionality

# Defaults follow GNU conventions

set(C_DIALECT_OPT_C89    "-std=c89")
set(C_DIALECT_OPT_C89EXT "-std=gnu89")
set(C_DIALECT_OPT_C99    "-std=c99")
set(C_DIALECT_OPT_C99EXT "-std=gnu99")
set(C_DIALECT_OPT_C11    "-std=c11")
set(C_DIALECT_OPT_C11EXT "-std=gnu11")

set(CXX_DIALECT_OPT_CXX03    "-std=c++03")
set(CXX_DIALECT_OPT_CXX03EXT "-std=gnu++03")
set(CXX_DIALECT_OPT_CXX11    "-std=c++11")
set(CXX_DIALECT_OPT_CXX11EXT "-std=gnu++11")
set(CXX_DIALECT_OPT_CXX14    "-std=c++14")
set(CXX_DIALECT_OPT_CXX14EXT "-std=gnu++14")


if( (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^ppc" ) OR ( CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^power" ) )
     ## ppc arch do not support -march= syntax 
     set(COMPILER_OPT_ARCH_NATIVE "-mcpu=native")
     add_definitions(-D_BLUEGENE)
else()
     set(COMPILER_OPT_ARCH_NATIVE "-march=native")
endif()



if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        set(CXX_DIALECT_OPT_CXX14    "-std=c++1y")
        set(CXX_DIALECT_OPT_CXX14EXT "-std=gnu++1y")
    endif()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
        set(CXX_DIALECT_OPT_CXX11    "-std=c++0x")
        set(CXX_DIALECT_OPT_CXX11EXT "-std=gnu++0x")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        set(CXX_DIALECT_OPT_CXX14    "-std=c++1y")
        set(CXX_DIALECT_OPT_CXX14EXT "-std=gnu++1y")
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES "XL")
    set(C_DIALECT_OPT_C89    "-qlanglvl=stdc89")
    set(C_DIALECT_OPT_C89EXT "-qlanglvl=extc89")
    set(C_DIALECT_OPT_C99    "-qlanglvl=stdc99")
    set(C_DIALECT_OPT_C99EXT "-qlanglvl=extc99")
    set(C_DIALECT_OPT_C11    "-qlanglvl=extc1x")
    set(C_DIALECT_OPT_C11EXT "-qlanglvl=extc1x")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "XL")
    set(CXX_DIALECT_OPT_CXX03    "-qlanglvl=extended")
    set(CXX_DIALECT_OPT_CXX03EXT "-qlanglvl=extended")
    set(CXX_DIALECT_OPT_CXX11    "-qlanglvl=extended0x")
    set(CXX_DIALECT_OPT_CXX11EXT "-qlanglvl=extended0x")

    set(COMPILER_OPT_ARCH_NATIVE "-qarch=auto")
endif()
