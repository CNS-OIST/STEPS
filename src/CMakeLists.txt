####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

set(lib_soversion "${VERSION_MAJOR}.${VERSION_MINOR}")
set(lib_version "${lib_soversion}.${VERSION_PATCH}")

set(cvode   "third_party/cvode-2.6.0/src")
set(nr      "third_party/nr")
set(superlu "third_party/superlu_dist-4.1/src")


set(superlu_sources
    "${superlu}/GetDiagU.c"
    "${superlu}/comm.c"
#    "${superlu}/cublas_utils.c"
    "${superlu}/ddistribute.c"    "${superlu}/dgsequ.c"
    "${superlu}/dlamch.c"         "${superlu}/dlangs.c"
    "${superlu}/dlaqgs.c"         "${superlu}/dldperm_dist.c"
    "${superlu}/dmemory.c"        "${superlu}/dmyblas2.c"
    "${superlu}/dsp_blas2.c"      "${superlu}/dsp_blas3.c"
    "${superlu}/dutil.c"          "${superlu}/etree.c"
    "${superlu}/get_perm_c.c"     "${superlu}/get_perm_c_parmetis.c"
    "${superlu}/lsame.c"          "${superlu}/mc64ad_dist.c"
    "${superlu}/memory.c"         "${superlu}/mmd.c"
    "${superlu}/pddistribute.c"   "${superlu}/pdgsequ.c"
    "${superlu}/pdgsmv.c"         "${superlu}/pdgsmv_AXglobal.c"
    "${superlu}/pdgsrfs.c"        "${superlu}/pdgsrfs_ABXglobal.c"
    "${superlu}/pdgssvx.c"        "${superlu}/pdgssvx_ABglobal.c"
    "${superlu}/pdgstrf.c"        "${superlu}/pdgstrf2.c"
    "${superlu}/pdgstrs.c"        "${superlu}/pdgstrs1.c"
    "${superlu}/pdgstrs_Bglobal.c"
    "${superlu}/pdgstrs_lsum.c"   "${superlu}/pdlangs.c"
    "${superlu}/pdlaqgs.c"        "${superlu}/pdsymbfact_distdata.c"
    "${superlu}/pdutil.c"         "${superlu}/psymbfact.c"
    "${superlu}/psymbfact_util.c" "${superlu}/pxerbla.c"
    "${superlu}/slamch.c"         "${superlu}/sp_colorder.c"
    "${superlu}/sp_ienv.c"        "${superlu}/static_schedule.c"
    "${superlu}/superlu_grid.c"   "${superlu}/superlu_timer.c"
    "${superlu}/symbfact.c"       "${superlu}/util.c"
    "${superlu}/xerbla.c"
)

set(lib_sources
    "steps/init.cpp"                           "steps/error.cpp"
    "steps/finish.cpp"
    "steps/geom/tetmesh.cpp"                   "steps/geom/comp.cpp"
    "steps/geom/geom.cpp"                      "steps/geom/patch.cpp"
                                               "steps/geom/tmcomp.cpp"
    "steps/geom/tmpatch.cpp"                   "steps/geom/sdiffboundary.cpp"
    "steps/geom/memb.cpp"                      "steps/geom/diffboundary.cpp"
    "steps/model/model.cpp"                    "steps/model/diff.cpp"
    "steps/model/chan.cpp"                     "steps/model/reac.cpp"
    "steps/model/spec.cpp"                     "steps/model/sreac.cpp"
    "steps/model/surfsys.cpp"                  "steps/model/volsys.cpp"
    "steps/model/chanstate.cpp"                "steps/model/ohmiccurr.cpp"
    "steps/model/ghkcurr.cpp"                  "steps/model/vdeptrans.cpp"
    "steps/model/vdepsreac.cpp"                "steps/math/tetrahedron.cpp"
    "steps/math/tools.cpp"                     "steps/math/linsolve.cpp"
    "steps/math/triangle.cpp"                  "steps/math/ghk.cpp"
    "steps/tetode/comp.cpp"                    "steps/tetode/patch.cpp"
    "steps/tetode/tet.cpp"                     "steps/tetode/tri.cpp"
    "steps/tetode/tetode.cpp"                  "steps/solver/api_comp.cpp"
    "steps/solver/api_main.cpp"                "steps/solver/api_memb.cpp"
    "steps/solver/api_patch.cpp"               "steps/solver/api_tet.cpp"
    "steps/solver/api_vert.cpp"                "steps/solver/api_tri.cpp"
    "steps/solver/api_diffboundary.cpp"        "steps/solver/api_recording.cpp"
    "steps/solver/api_batchdata.cpp"           "steps/solver/api_roidata.cpp"
    "steps/solver/compdef.cpp"                 "steps/solver/diffdef.cpp"
    "steps/solver/patchdef.cpp"                "steps/solver/api_sdiffboundary.cpp"
    "steps/solver/reacdef.cpp"                 "steps/solver/specdef.cpp"
    "steps/solver/sreacdef.cpp"                "steps/solver/statedef.cpp"
    "steps/solver/chandef.cpp"                 "steps/solver/ghkcurrdef.cpp"
    "steps/solver/diffboundarydef.cpp"         "steps/solver/ohmiccurrdef.cpp"
    "steps/solver/vdeptransdef.cpp"            "steps/solver/vdepsreacdef.cpp"
    "steps/solver/sdiffboundarydef.cpp"
    "steps/solver/efield/dVsolver.cpp"
    "steps/solver/efield/bdsystem.cpp"
    "steps/solver/efield/dVsolver.cpp"
    "steps/solver/efield/efield.cpp"           "steps/solver/efield/matrix.cpp"
    "steps/solver/efield/tetcoupler.cpp"       "steps/solver/efield/tetmesh.cpp"
    "steps/solver/efield/vertexconnection.cpp" "steps/solver/efield/vertexelement.cpp"
    "steps/tetexact/comp.cpp"                  "steps/tetexact/diff.cpp"
    "steps/tetexact/sdiff.cpp"                 "steps/tetexact/kproc.cpp"
    "steps/tetexact/patch.cpp"                 "steps/tetexact/reac.cpp"
    "steps/tetexact/sreac.cpp"                 "steps/tetexact/tet.cpp"
    "steps/tetexact/tetexact.cpp"              "steps/tetexact/tri.cpp"
    "steps/tetexact/ghkcurr.cpp"               "steps/tetexact/vdeptrans.cpp"
    "steps/tetexact/vdepsreac.cpp"             "steps/tetexact/diffboundary.cpp"
    "steps/tetexact/wmvol.cpp"                 "steps/tetexact/sdiffboundary.cpp"
    "steps/wmdirect/comp.cpp"
    "steps/wmdirect/kproc.cpp"                 "steps/wmdirect/patch.cpp"
    "steps/wmdirect/reac.cpp"                  "steps/wmdirect/sreac.cpp"
    "steps/wmdirect/wmdirect.cpp"              "steps/wmrk4/wmrk4.cpp"
    "steps/rng/rng.cpp"                        "steps/rng/mt19937.cpp"
    "steps/rng/r123.cpp"
    "steps/rng/create.cpp"
    "steps/util/checkid.cpp"
    #
    "${cvode}/cvode/cvode_band.cpp"            "${cvode}/cvode/cvode_bandpre.cpp"
    "${cvode}/cvode/cvode_bbdpre.cpp"          "${cvode}/cvode/cvode_dense.cpp"
    "${cvode}/cvode/cvode_diag.cpp"            "${cvode}/cvode/cvode_direct.cpp"
    "${cvode}/cvode/cvode_io.cpp"              "${cvode}/cvode/cvode_spbcgs.cpp"
    "${cvode}/cvode/cvode_spgmr.cpp"           "${cvode}/cvode/cvode_spils.cpp"
    "${cvode}/cvode/cvode_sptfqmr.cpp"         "${cvode}/cvode/cvode.cpp"
    "${cvode}/nvec_ser/fnvector_serial.cpp"    "${cvode}/nvec_ser/nvector_serial.cpp"
    "${cvode}/sundials/sundials_band.cpp"      "${cvode}/sundials/sundials_dense.cpp"
    "${cvode}/sundials/sundials_direct.cpp"    "${cvode}/sundials/sundials_iterative.cpp"
    "${cvode}/sundials/sundials_math.cpp"      "${cvode}/sundials/sundials_nvector.cpp"
    "${cvode}/sundials/sundials_spbcgs.cpp"    "${cvode}/sundials/sundials_spgmr.cpp"
    "${cvode}/sundials/sundials_sptfqmr.cpp"
    #
    "${nr}/gamma.cpp"                           "${nr}/incgammabeta.cpp"
)

set_source_files_properties(${superlu_sources} PROPERTIES COMPILE_FLAGS "-std=c99 -DSUPERLU_NO_METIS")
set_source_files_properties("${superlu}/slamch.c" "${superlu}/dlamch.c" PROPERTIES COMPILE_FLAGS "-O0")

set(lib_public_headers
    "steps/common.h"
    "steps/error.hpp"                          "steps/init.hpp"
    "steps/finish.hpp"
    #
    "steps/geom/comp.hpp"                      "steps/geom/diffboundary.hpp"
    "steps/geom/geom.hpp"                      "steps/geom/memb.hpp"
    "steps/geom/patch.hpp"                     "steps/geom/sdiffboundary.hpp"
    "steps/geom/tetmesh.hpp"                   "steps/geom/tetmesh_rw.hpp"
    "steps/geom/tmcomp.hpp"                    "steps/geom/tmpatch.hpp"
    #
    "steps/util/collections.hpp"               "steps/util/fnv_hash.hpp"
    "steps/util/type_traits.hpp"               "steps/util/checkid.hpp"
    #
    "steps/math/constants.hpp"                 "steps/math/ghk.hpp"
    "steps/math/linsolve.hpp"                  "steps/math/tetrahedron.hpp"
    "steps/math/tools.hpp"                     "steps/math/triangle.hpp"
    "steps/math/point.hpp"                     "steps/math/bbox.hpp"
    #
    "steps/model/chan.hpp"                     "steps/model/chanstate.hpp"
    "steps/model/diff.hpp"                     "steps/model/ghkcurr.hpp"
    "steps/model/model.hpp"                    "steps/model/ohmiccurr.hpp"
    "steps/model/reac.hpp"                     "steps/model/spec.hpp"
    "steps/model/sreac.hpp"                    "steps/model/surfsys.hpp"
    "steps/model/vdepsreac.hpp"                "steps/model/vdeptrans.hpp"
    "steps/model/volsys.hpp"
    #
    "steps/rng/mt19937.hpp"                    "steps/rng/rng.hpp"
    "steps/rng/r123.hpp"
    "steps/rng/create.hpp"
    #
    "steps/solver/api.hpp"                     "steps/solver/chandef.hpp"
    "steps/solver/compdef.hpp"
    "steps/solver/diffboundarydef.hpp"         "steps/solver/diffdef.hpp"
    "steps/solver/sdiffboundarydef.hpp"
    "steps/solver/efield/bdsystem_lapack.hpp"  "steps/solver/efield/bdsystem.hpp"
    "steps/solver/efield/dVsolver.hpp"         "steps/solver/efield/efield.hpp"
    "steps/solver/efield/dVsolver_slu.hpp"
    "steps/solver/efield/efieldsolver.hpp"     "steps/solver/efield/linsystem.hpp"
    "steps/solver/efield/matrix.hpp"           "steps/solver/efield/tetcoupler.hpp"
    "steps/solver/efield/tetmesh.hpp"          "steps/solver/efield/vertexconnection.hpp"
    "steps/solver/efield/vertexelement.hpp"
    "steps/solver/ghkcurrdef.hpp"              "steps/solver/ohmiccurrdef.hpp"
    "steps/solver/patchdef.hpp"                "steps/solver/reacdef.hpp"
    "steps/solver/specdef.hpp"                 "steps/solver/sreacdef.hpp"
    "steps/solver/statedef.hpp"
    "steps/solver/types.hpp"                   "steps/solver/vdepsreacdef.hpp"
    "steps/solver/vdeptransdef.hpp"
    #
    "steps/tetexact/comp.hpp"                  "steps/tetexact/crstruct.hpp"
    "steps/tetexact/diff.hpp"                  "steps/tetexact/diffboundary.hpp"
    "steps/tetexact/ghkcurr.hpp"               "steps/tetexact/kproc.hpp"
    "steps/tetexact/patch.hpp"                 "steps/tetexact/reac.hpp"
    "steps/tetexact/sdiff.hpp"                 "steps/tetexact/sreac.hpp"
    "steps/tetexact/tet.hpp"                   "steps/tetexact/tetexact.hpp"
    "steps/tetexact/tri.hpp"                   "steps/tetexact/vdepsreac.hpp"
    "steps/tetexact/vdeptrans.hpp"             "steps/tetexact/wmvol.hpp"
    "steps/tetexact/sdiffboundary.hpp"
    #
    "steps/tetode/comp.hpp"                    "steps/tetode/patch.hpp"
    "steps/tetode/tet.hpp"                     "steps/tetode/tetode.hpp"
    "steps/tetode/tri.hpp"
    #
    "steps/wmdirect/comp.hpp"                  "steps/wmdirect/kproc.hpp"
    "steps/wmdirect/patch.hpp"                 "steps/wmdirect/reac.hpp"
    "steps/wmdirect/sreac.hpp"                 "steps/wmdirect/wmdirect.hpp"
    #
    "steps/wmrk4/wmrk4.hpp"
)

if(LAPACK_FOUND)
    list(APPEND lib_sources "steps/solver/efield/bdsystem_lapack.cpp")
    list(APPEND lib_public_headers "steps/solver/efield/bdsystem_lapack.hpp")
endif()

# Come back to lapack version later ...
#if(USE_BDSYSTEM_LAPACK AND LAPACK_FOUND)
#    add_definitions(-DUSE_BDSYSTEM_LAPACK)
#    message(STATUS "E-field solver: using BDSystemLapack implementation")
#else()
#    add_definitions(-DUSE_BDSYSTEM)
#    message(STATUS "E-field solver: using BDSystem implementation")
#endif()

include_directories(".")

if(MPI_FOUND)
    list(APPEND lib_sources
    "steps/solver/efield/slusystem.cpp"
    "steps/mpi/mpi_init.cpp"                    "steps/mpi/mpi_finish.cpp"
    "steps/mpi/tetopsplit/comp.cpp"             "steps/mpi/tetopsplit/diff.cpp"
    "steps/mpi/tetopsplit/sdiff.cpp"            "steps/mpi/tetopsplit/kproc.cpp"
    "steps/mpi/tetopsplit/patch.cpp"            "steps/mpi/tetopsplit/reac.cpp"
    "steps/mpi/tetopsplit/sreac.cpp"            "steps/mpi/tetopsplit/tet.cpp"
    "steps/mpi/tetopsplit/tetopsplit.cpp"       "steps/mpi/tetopsplit/tri.cpp"
    "steps/mpi/tetopsplit/ghkcurr.cpp"          "steps/mpi/tetopsplit/vdeptrans.cpp"
    "steps/mpi/tetopsplit/vdepsreac.cpp"        "steps/mpi/tetopsplit/diffboundary.cpp"
    "steps/mpi/tetopsplit/wmvol.cpp"			"steps/mpi/tetopsplit/sdiffboundary.cpp"
    #
    ${superlu_sources}
    )

    list(APPEND lib_public_headers
    "steps/solver/efield/slusystem.hpp"
    "steps/solver/efield/dVsolver_slu.hpp"
    "steps/mpi/mpi_common.hpp"
    "steps/mpi/mpi_init.hpp"                    "steps/mpi/mpi_finish.hpp"
    "steps/mpi/tetopsplit/comp.hpp"             "steps/mpi/tetopsplit/crstruct.hpp"
    "steps/mpi/tetopsplit/diff.hpp"             "steps/mpi/tetopsplit/diffboundary.hpp"
    "steps/mpi/tetopsplit/ghkcurr.hpp"          "steps/mpi/tetopsplit/kproc.hpp"
    "steps/mpi/tetopsplit/patch.hpp"            "steps/mpi/tetopsplit/reac.hpp"
    "steps/mpi/tetopsplit/sdiff.hpp"            "steps/mpi/tetopsplit/sreac.hpp"
    "steps/mpi/tetopsplit/tet.hpp"              "steps/mpi/tetopsplit/tetopsplit.hpp"
    "steps/mpi/tetopsplit/tri.hpp"              "steps/mpi/tetopsplit/vdepsreac.hpp"
    "steps/mpi/tetopsplit/vdeptrans.hpp"        "steps/mpi/tetopsplit/wmvol.hpp"
    "steps/mpi/tetopsplit/sdiffboundary.hpp"
    )

    include_directories("." ${MPI_C_INCLUDE_PATH})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif()


if (PETSC_FOUND)
    include_directories(${PETSC_INCLUDE_DIRS})
    list(APPEND lib_sources "steps/solver/efield/dVsolver_petsc.cpp")
    list(APPEND lib_public_headers "steps/solver/efield/dVsolver_petsc.hpp")
endif()

add_definitions(-D_ELPP_NO_DEFAULT_LOG_FILE=1 -D_ELPP_STL_LOGGING=1 -D_ELPP_DISABLE_DEFAULT_CRASH_HANDLING=1)

# enable below to turn on MPI profiling
add_definitions(-DMPI_PROFILING=1)

# enable debug file
#add_definitions(-UNDEBUG -D_DEBUG=1)

# enable surface diffusion boundary debug
#add_definitions(-DSDIFFBOUNDARY_DEBUG=1)

# only enable if you know what you are doing :)
#add_definitions(-DMPI_DEBUG=1 -UNDEBUG)
#add_definitions(-DSERIAL_EFIELD_DEBUG=1)
#add_definitions(-D_ELPP_STL_LOGGING=1 -D_ELPP_STACKTRACE_ON_CRASH=1)

add_library(libsteps-obj OBJECT ${lib_sources})
set_target_properties(libsteps-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(libsteps_static STATIC $<TARGET_OBJECTS:libsteps-obj>)
set_target_properties(libsteps_static PROPERTIES OUTPUT_NAME steps)

add_library(libsteps SHARED $<TARGET_OBJECTS:libsteps-obj>)
target_link_libraries(libsteps ${libsteps_link_libraries})
set_target_properties(libsteps PROPERTIES VERSION "${lib_version}" SOVERSION "${lib_soversion}")
set_target_properties(libsteps PROPERTIES OUTPUT_NAME steps)


# ==========================================================================
# Installation of the headers and shared lib
#
option(INSTALL_DEV_LIBS "Install steps lib and headers which are required for development only" OFF)

if(INSTALL_DEV_LIBS)
install(TARGETS libsteps
    RUNTIME DESTINATION bin COMPONENT lib
    ARCHIVE DESTINATION lib COMPONENT lib
    LIBRARY DESTINATION lib COMPONENT lib)

foreach(header ${lib_public_headers})
    get_filename_component(inc_subdir ${header} DIRECTORY)
    install(FILES ${header} DESTINATION include/${inc_subdir})
endforeach()
endif()
