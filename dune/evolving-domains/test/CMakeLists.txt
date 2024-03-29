set(GRIDTYPE  ALBERTAGRID)
set(WORLDDIM  3)
set(POLORDER  3)
set(WANT_ISTL  1)
set(USE_OLD_GRIDPART 0)

set(CLEANFILES  *.log)

set(check_PROGRAMS  test-surface test-coupled )

## test surface ##

math(EXPR SURFACE_GRIDDIM '${WORLDDIM}-1')
add_executable(test-surface test-surface.cc)
add_dune_alberta_flags(test-surface GRIDDIM ${SURFACE_GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-surface )
target_compile_definitions( test-surface PRIVATE ${GRIDTYPE};GRIDDIM=${SURFACE_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=0 )

add_executable(test-surface-old test-surface.cc)
add_dune_alberta_flags(test-surface-old GRIDDIM ${SURFACE_GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-surface-old )
target_compile_definitions( test-surface-old PRIVATE ${GRIDTYPE};GRIDDIM=${SURFACE_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=1;USE_OLD_GRIDPART=1 )


## test coupled ##

add_executable(test-coupled test-coupled.cc)
set(BULK_GRIDDIM ${WORLDDIM})
add_dune_alberta_flags(test-coupled GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( test-coupled )
target_compile_definitions( test-coupled PRIVATE ${GRIDTYPE};GRIDDIM=${BULK_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART} )

dune_symlink_to_source_files(FILES parameter sphere2.dgf ball.dgf)

## test refinement projection ##

add_executable( test-refinement-projection test.cc)
dune_target_enable_all_packages( test-refinement-projection )
add_dune_alberta_flags( test-refinement-projection GRIDDIM 3 WORLDDIM 3 )
