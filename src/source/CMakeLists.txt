set(GRIDTYPE  ALBERTAGRID)
set(GRIDDIM  3)
set(WORLDDIM  3)
set(POLORDER  2)
set(WANT_ISTL  1)
set(USE_OLD_GRIDPART 0)

set(CLEANFILES  *.log)

set(check_PROGRAMS  main-surface-${POLORDER} main-coupled-${POLORDER})

add_executable(main-surface-${POLORDER} main-surface.cc)
math(EXPR SURFACE_GRIDDIM '${WORLDDIM}-1')
add_dune_alberta_flags(main-surface-${POLORDER} GRIDDIM ${SURFACE_GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( main-surface-${POLORDER} )
target_compile_definitions( main-surface-${POLORDER} PRIVATE ${GRIDTYPE};GRIDDIM=${SURFACE_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART} )


add_executable(main-coupled-${POLORDER} main-coupled.cc)
set(BULK_GRIDDIM ${GRIDDIM})
add_dune_alberta_flags(main-coupled-${POLORDER} GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( main-coupled-${POLORDER} )
target_compile_definitions( main-coupled-${POLORDER} PRIVATE ${GRIDTYPE};GRIDDIM=${BULK_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART} )

add_executable(poisson-surface poisson-surface.cc)
math(EXPR SURFACE_GRIDDIM '${WORLDDIM}-1')
add_dune_alberta_flags(poisson-surface GRIDDIM ${SURFACE_GRIDDIM} WORLDDIM ${WORLDDIM})
target_link_dune_default_libraries( poisson-surface )
target_compile_definitions( poisson-surface PRIVATE ${GRIDTYPE};GRIDDIM=${SURFACE_GRIDDIM};WORLDDIM=${WORLDDIM};WANT_ISTL=${WANT_ISTL};POLORDER=${POLORDER};USE_OLD_GRIDPART=${USE_OLD_GRIDPART};FEMTIMER=1 )


set(EXTRA_DIST   main.cc deformation.hh dirichletconstraints.hh elliptic.hh femscheme.hh heat.hh heatmodel.hh heatscheme.hh model.hh probleminterface.hh rhs.hh temporalprobleminterface.hh heat2.dox)

#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

dune_symlink_to_source_files(FILES run.sh)
