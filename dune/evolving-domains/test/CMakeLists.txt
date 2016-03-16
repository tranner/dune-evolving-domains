set(GRIDTYPE  ALUGRID_SIMPLEX)
set(WORLDDIM  3)
set(POLORDER  1)
set(WANT_ISTL  1)
set(USE_OLD_GRIDPART 1)

set(CLEANFILES  *.log)

set(check_PROGRAMS  test-surface test-coupled )

## test surface ##

math(EXPR SURFACE_GRIDDIM '${WORLDDIM}-1')
add_executable(test-surface test-surface.cc)
add_dune_alberta_flags(test-surface GRIDDIM ${SURFACE_GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-surface )
target_compile_definitions( test-surface PRIVATE ${GRIDTYPE};GRIDDIM=${SURFACE_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART} )

## test coupled ##

add_executable(test-coupled test-coupled.cc)
set(BULK_GRIDDIM ${WORLDDIM})
add_dune_alberta_flags(test-coupled GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-coupled )
target_compile_definitions( test-coupled PRIVATE ${GRIDTYPE};GRIDDIM=${BULK_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART} )

dune_symlink_to_source_files(FILES parameter sphere2.dgf ball.dgf)
