set(GRIDTYPE  ALUGRID_SIMPLEX)
set(GRIDDIM  3)
set(WORLDDIM  3)
set(POLORDER  3)
set(WANT_ISTL  1)
set(USE_OLD_GRIDPART 0)

set(CLEANFILES  *.log)

set(check_PROGRAMS  test-surface test-coupled )

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL} -DGRIDDIM=${GRIDDIM} -DWORLDDIM=${WORLDDIM} -DUSE_OLD_GRIDPART=${USE_OLD_GRIDPART}")

add_executable(test-surface test-surface.cc)
add_dune_alberta_flags(test-surface GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-surface )

add_executable(test-coupled test-coupled.cc)
add_dune_alberta_flags(test-coupled GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-coupled )

dune_symlink_to_source_files(FILES parameter sphere2.dgf ball.dgf)
