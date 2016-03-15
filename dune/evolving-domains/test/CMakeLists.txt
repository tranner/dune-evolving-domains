set(GRIDTYPE  ALUGRID_SIMPLEX)
set(GRIDDIM  3)
set(WORLDDIM  3)
set(POLORDER  3)
set(WANT_ISTL  1)

set(CLEANFILES  *.log)

set(check_PROGRAMS  test-surface )

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL} -DGRIDDIM=${GRIDDIM} -DWORLDDIM=${WORLDDIM}")

add_executable(test-surface test-surface.cc)
add_dune_alberta_flags(test-surface GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-surface )

# add_executable(main-coupled main-coupled.cc)
# add_dune_alberta_flags(main-coupled GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
# target_link_dune_default_libraries( main-coupled )

dune_symlink_to_source_files(FILES parameter sphere2.dgf ball.dgf)
